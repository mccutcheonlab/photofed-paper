"""Build the photometry pickles used by the figure notebooks from raw TDT tanks.

Writes only the 24 pickles listed in scripts/photometry_files.py: pellet onset
and offset for each of the 12 FR1 and progressive-ratio conditions.

Run from the repository root:

    pixi run python src/raw_to_pickles.py --tanks-dir <folder with tanks>

Tanks are read from --tanks-dir, which defaults to tanks/ in the repository.
The downloader will populate that folder once raw data are published.
"""

import argparse
import csv
import math
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import tdt
import trompy as tp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from photometry_files import PICKLE_NAMES  # noqa: E402

METADATA_PATH = ROOT / "metadata" / "fotofeda_photometry_metadata_selected.xlsx"
DEFAULT_TANKS_DIR = ROOT / "tanks"
DEFAULT_OUT_DIR = ROOT / "data"

EVENT_KEY = "ttl_pellet"
# Snipper window, matching the original notebook call: trompy defaults
# (10 s baseline, 30 s total, i.e. -10 to +20 s around the event).
SNIP_BASELINE = 10
SNIP_TRIAL = 30
# The original notebook binned the 30 s snip into int(15 * fs) bins, giving the
# deposit's 15258 columns. Each bin averages about two samples.
BIN_SECONDS = 15
N_BOOTSTRAP = 1000
CI_PERCENT = 95
BOOTSTRAP_SEED = 0


def wanted_outputs() -> dict[tuple[str, str], str]:
    """Map (condition prefix, alignment) to the pickle name, from PICKLE_NAMES."""
    wanted = {}
    for name in PICKLE_NAMES:
        stem = name[: -len(".pickle")]
        prefix, _, align = stem.partition("_ttl_pellet_")
        wanted[(prefix, align)] = name
    return wanted


def bootstrap_mean_ci(data: np.ndarray, n: int, ci: float, seed: int):
    """Bootstrap the mean trace over events. Returns mean, lower and upper bounds."""
    rng = np.random.default_rng(seed)
    n_events, n_bins = data.shape
    boot_means = np.empty((n, n_bins), dtype=np.float32)
    for i in range(n):
        indices = rng.integers(0, n_events, n_events)
        boot_means[i] = data[indices].mean(axis=0)
    mean = data.mean(axis=0)
    lower = np.percentile(boot_means, (100 - ci) / 2, axis=0)
    upper = np.percentile(boot_means, 100 - (100 - ci) / 2, axis=0)
    return mean, lower, upper


def kept_event_times(times, fs: float, n_samples: int) -> np.ndarray:
    """Return the event times that trompy.snipper keeps.

    Mirrors the filter inside trompy.snipper (finite values, and a full snip
    must fit inside the recording) so that event_times stays row-aligned with
    the snips.
    """
    times = np.asarray(times, dtype=float)
    times = times[np.isfinite(times)]
    samples = np.array([math.ceil(t * fs) for t in times], dtype=np.int64)
    baseline = int(SNIP_BASELINE * fs)
    trial = int(SNIP_TRIAL * fs)
    keep = (samples - baseline > 0) & (samples - baseline + trial < n_samples)
    return times[keep]


def snip_events(corrected: np.ndarray, times, fs: float) -> tuple[np.ndarray, np.ndarray]:
    """Snip a corrected trace around each event, as the original notebook did.

    The window runs from -SNIP_BASELINE to +(SNIP_TRIAL - SNIP_BASELINE) seconds,
    binned into int(BIN_SECONDS * fs) columns. Baseline subtraction is done by
    trompy.snipper (adjust_baseline=True).
    """
    kept = kept_event_times(times, fs, len(corrected))
    snips, _ = tp.snipper(
        corrected,
        kept,
        fs=fs,
        baseline_length=SNIP_BASELINE,
        trial_length=SNIP_TRIAL,
        bins=int(BIN_SECONDS * fs),
    )
    if snips.shape[0] != len(kept):
        raise RuntimeError(f"snipper returned {snips.shape[0]} snips for {len(kept)} events")
    return snips.astype(np.float32), kept


def read_tank(tank_path: Path):
    """Read one TDT tank and return the block."""
    return tdt.read_block(str(tank_path))


def corrected_signal(block, row: pd.Series):
    """Return the isosbestic-corrected blue signal and its sampling rate."""
    blue = getattr(block.streams, f"_{row['blue']}").data
    uv_stream = getattr(block.streams, f"_{row['uv']}")
    fs = uv_stream.fs
    return tp.processdata(blue, uv_stream.data, fs=fs), fs


def build_pickle(label: str, parts: list[dict], fs: float) -> dict:
    """Assemble one pickle payload from per-tank snips."""
    all_snips = np.vstack([p["snips"] for p in parts])
    event_times = np.concatenate([p["event_times"] for p in parts]).astype(np.float32)
    mouse_ids = np.concatenate([p["mouse_ids"] for p in parts])
    session_numbers = np.concatenate([p["session_numbers"] for p in parts])
    event_idx = np.concatenate([p["event_idx"] for p in parts])

    mean, lower, upper = bootstrap_mean_ci(
        all_snips, n=N_BOOTSTRAP, ci=CI_PERCENT, seed=BOOTSTRAP_SEED
    )
    n_bins = all_snips.shape[1]
    # Labels the true window. The deposit labels its axis -5 to +10 s instead.
    time_vector = np.linspace(-SNIP_BASELINE, SNIP_TRIAL - SNIP_BASELINE, n_bins)

    return {
        "label": label,
        "fs": fs,
        "pre_time": SNIP_BASELINE,
        "post_time": SNIP_TRIAL - SNIP_BASELINE,
        "time_vector": time_vector.astype(np.float32),
        "snips": all_snips,
        "event_times": event_times,
        "mouse_ids": mouse_ids,
        "session_numbers": session_numbers,
        "event_index_within_session": event_idx,
        "n_events": all_snips.shape[0],
        "summary": {
            "mean": mean.astype(np.float32),
            "lower": lower.astype(np.float32),
            "upper": upper.astype(np.float32),
            "n_bootstrap": N_BOOTSTRAP,
            "ci_percent": CI_PERCENT,
        },
    }


def assemble(metadata: pd.DataFrame, tanks_dir: Path, out_dir: Path) -> list[dict]:
    """Build every wanted pickle and return a list of skipped-tank records."""
    wanted = wanted_outputs()
    skipped: list[dict] = []
    out_dir.mkdir(parents=True, exist_ok=True)

    metadata = metadata.assign(
        condition=metadata["session_number"].astype(str)
        + "_" + metadata["diet"].astype(str)
        + "_" + metadata["session_pellet"].astype(str)
        + "_" + metadata["session_mode"].astype(str)
    )

    prefixes = sorted({prefix for prefix, _ in wanted})
    for prefix in prefixes:
        cond_df = metadata[metadata["condition"] == prefix]
        if cond_df.empty:
            skipped.append({"condition": prefix, "mouseid": "", "align": "", "reason": "no metadata rows"})
            continue

        parts = {align: [] for align in ("onset", "offset") if (prefix, align) in wanted}
        fs_used = None

        for _, row in cond_df.iterrows():
            try:
                block = read_tank(tanks_dir / row["tdt_tank"])
                corrected, fs = corrected_signal(block, row)
            except Exception as error:  # record and keep going; reported at the end
                skipped.append({"condition": prefix, "mouseid": row["mouseid"], "align": "all", "reason": str(error)})
                continue

            fs_used = fs
            ttl_label = f"_{row[EVENT_KEY]}_"
            for align in parts:
                try:
                    if not hasattr(block.epocs, ttl_label):
                        raise ValueError("missing TTL")
                    times = np.asarray(getattr(getattr(block.epocs, ttl_label), align))
                    if times.size == 0:
                        raise ValueError("no TTL events")
                    snips, kept = snip_events(corrected, times, fs)
                except Exception as error:
                    skipped.append({"condition": prefix, "mouseid": row["mouseid"], "align": align, "reason": str(error)})
                    continue

                n = snips.shape[0]
                parts[align].append({
                    "snips": snips,
                    "event_times": kept,
                    "mouse_ids": np.array([row["mouseid"]] * n, dtype=object),
                    "session_numbers": np.array([row["session_number"]] * n, dtype=object),
                    "event_idx": np.arange(n, dtype=np.int32),
                })

        for align, align_parts in parts.items():
            name = wanted[(prefix, align)]
            if not align_parts:
                continue
            label = name[: -len(".pickle")]
            payload = build_pickle(label, align_parts, fs_used)
            with open(out_dir / name, "wb") as handle:
                pickle.dump(payload, handle, protocol=pickle.HIGHEST_PROTOCOL)
            print(f"Wrote {name}: {payload['n_events']} events")

    return skipped


def main() -> None:
    """Assemble the figure pickles and check that every expected file was written."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--tanks-dir", type=Path, default=DEFAULT_TANKS_DIR, help="Folder containing the TDT tanks")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR, help="Folder for the output pickles")
    parser.add_argument("--allow-skips", action="store_true",
                        help="Do not fail when a tank or event type is skipped (see skipped_log.csv)")
    args = parser.parse_args()

    existing = [name for name in PICKLE_NAMES if (args.out_dir / name).exists()]
    if existing:
        raise SystemExit("Refusing to overwrite existing files:\n" + "\n".join(existing))

    metadata = pd.read_excel(METADATA_PATH)
    skipped = assemble(metadata, args.tanks_dir, args.out_dir)

    log_path = args.out_dir / "skipped_log.csv"
    with open(log_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["condition", "mouseid", "align", "reason"])
        writer.writeheader()
        writer.writerows(skipped)

    missing = [name for name in PICKLE_NAMES if not (args.out_dir / name).exists()]
    if missing:
        raise SystemExit("Missing outputs:\n" + "\n".join(missing))
    if skipped and not args.allow_skips:
        raise SystemExit(f"{len(skipped)} skipped item(s); see {log_path} (or pass --allow-skips)")
    print(f"Wrote {len(PICKLE_NAMES)} pickles to {args.out_dir}")


if __name__ == "__main__":
    main()
