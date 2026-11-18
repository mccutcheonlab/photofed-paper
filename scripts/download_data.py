"""Download and verify the figure inputs from Zenodo."""

import argparse
import hashlib
import json
import shutil
import tempfile
import urllib.error
import urllib.request
import zipfile
from pathlib import Path
from urllib.parse import urlparse

from photometry_files import ARCHIVE_NAME, PICKLE_NAMES

RECORD = "24315309"
ROOT = Path(__file__).resolve().parents[1]
SEPARATE_FILES = (
    "food_and_body_weight.xlsx",
    "operant_behaviour_summary.csv",
    "pellet_retrieval_latency.csv",
)


def download_file(record: str, target: Path, name: str = ARCHIVE_NAME) -> None:
    """Download a named record file and verify its Zenodo checksum.

    Args:
        record: Numeric, published Zenodo record ID.
        target: Temporary path for the downloaded file.
        name: Exact filename in the record.
    """
    if not record.isascii() or not record.isdigit():
        raise ValueError("Zenodo record ID must contain only ASCII digits")
    try:
        with urllib.request.urlopen(f"https://zenodo.org/api/records/{record}") as response:
            metadata = json.load(response)
    except urllib.error.HTTPError as error:
        if error.code == 404:
            raise ValueError(f"Zenodo record {record} was not found or is not published") from error
        raise

    entries = [entry for entry in metadata["files"] if entry["key"] == name]
    if len(entries) != 1:
        raise ValueError(f"Expected exactly one {name} in Zenodo record {record}")
    entry = entries[0]
    algorithm, separator, checksum = entry["checksum"].partition(":")
    if algorithm != "md5" or not separator or len(checksum) != 32:
        raise ValueError("Zenodo did not provide a valid MD5 checksum")
    url = entry["links"]["self"]
    if urlparse(url).scheme != "https":
        raise ValueError("Zenodo download URL must use HTTPS")

    with urllib.request.urlopen(url) as response, target.open("xb") as destination:
        shutil.copyfileobj(response, destination)
    with target.open("rb") as stream:
        actual = hashlib.file_digest(stream, lambda: hashlib.md5(usedforsecurity=False)).hexdigest()
    if actual.lower() != checksum.lower():
        raise ValueError(f"Downloaded {name} does not match the Zenodo checksum")


def extract_archive(archive_path: Path, destination: Path) -> None:
    """Extract only the expected files, without overwriting existing files.

    Args:
        archive_path: Verified ZIP archive.
        destination: Directory to receive the pickle files.
    """
    expected = {f"data/{name}" for name in PICKLE_NAMES}
    with zipfile.ZipFile(archive_path) as archive:
        entries = archive.namelist()
        if len(entries) != len(expected) or set(entries) != expected:
            raise ValueError("Archive contents do not match the expected photometry inputs")
        existing = [name for name in PICKLE_NAMES if (destination / name).exists()]
        if existing:
            raise FileExistsError("Refusing to overwrite files:\n" + "\n".join(existing))

        destination.mkdir(parents=True, exist_ok=True)
        created: list[Path] = []
        try:
            for name in PICKLE_NAMES:
                path = destination / name
                with archive.open(f"data/{name}") as source, path.open("xb") as output:
                    created.append(path)
                    shutil.copyfileobj(source, output)
        except BaseException:
            for path in created:
                path.unlink()
            raise


def download_inputs(record: str, destination: Path) -> None:
    """Download and install all figure inputs without overwriting files.

    Args:
        record: Numeric, published Zenodo record ID.
        destination: Directory to receive all figure inputs.
    """
    names = (*PICKLE_NAMES, *SEPARATE_FILES)
    existing = [name for name in names if (destination / name).exists() or (destination / name).is_symlink()]
    if existing:
        raise FileExistsError("Refusing to overwrite files:\n" + "\n".join(existing))

    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        archive = root / ARCHIVE_NAME
        download_file(record, archive)
        for name in SEPARATE_FILES:
            download_file(record, root / name, name)
        payload = root / "payload"
        extract_archive(archive, payload)

        destination.mkdir(parents=True, exist_ok=True)
        created: list[Path] = []
        try:
            for name in names:
                source_path = payload / name if name in PICKLE_NAMES else root / name
                path = destination / name
                with source_path.open("rb") as source, path.open("xb") as output:
                    created.append(path)
                    shutil.copyfileobj(source, output)
        except BaseException:
            for path in created:
                path.unlink()
            raise


def main() -> None:
    """Download and extract the published figure inputs."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record", default=RECORD, help="Published Zenodo record ID")
    parser.add_argument("--dest", type=Path, default=ROOT / "data", help="Directory for the figure inputs")
    args = parser.parse_args()
    download_inputs(args.record, args.dest)
    print(f"Downloaded {len(PICKLE_NAMES)} pickles and {len(SEPARATE_FILES)} other files to {args.dest}")


if __name__ == "__main__":
    main()
