# photofed_paper

This repository is a slimmed down version of the codebase developed by Hamid Taghipourbibalan and James McCutcheon to analyse a fibre photometry experiment examining dopamine responses under protein restriction. The full paper is available on bioRxiv (https://www.biorxiv.org/content/10.64898/2026.07.31.741998v1). Data are available on Zenodo. Citation and download details follow.

The original repository that includes alternative analysis avenues, temporary binary files etc is available at https://github.com/Htbibalan/photofed_paper.

## Funding

This work was supported by a Tromsø Research Foundation Starting Grant to J.E. McCutcheon (19-SG-JMcC).

## Getting started

1. Install [Pixi](https://pixi.sh) and run `pixi install`.
2. Download the processed inputs from the public
   [Zenodo record](https://doi.org/10.5281/zenodo.21024650), as described
   below.
3. Run `pixi run jupyter lab`, open a figure notebook and run its cells.

Run the download and helper tests with
`pixi run python -m unittest discover -s tests -v`.

| Notebook | Analysis |
| --- | --- |
| `fig2_food_and_bw.ipynb` | Food intake and body weight |
| `fig3_pellet_behaviour.ipynb` | Operant behaviour and latency panels |
| `fig4_fr1_photo.ipynb` | FR1 photometry at pellet delivery |
| `fig5_progratio_photo.ipynb` | Progressive-ratio photometry at pellet delivery |
| `figSx_fr1_photo_pellet_offset.ipynb` | FR1 photometry at pellet retrieval |
| `figSx_progratio_photo_offset.ipynb` | Progressive-ratio photometry at pellet retrieval |
| `figSx_retrieval_latencies.ipynb` | Retrieval latency analyses |

## Citation

If you use this code or data, please cite the preprint in
[CITATION.cff](CITATION.cff) and the
[Zenodo dataset](https://doi.org/10.5281/zenodo.21024650). The preprint
reference will be updated when the paper is published.

## License

The code is released under the GNU General Public License v3.0 or later (GPL-3.0-or-later).
The data are available under CC BY 4.0 license as published on Zenodo.

## Photometry data

The figure notebooks require processed photometry pickles. The public Zenodo
record [10.5281/zenodo.21024650](https://doi.org/10.5281/zenodo.24315309)
contains `photometry_pickles.zip` and the three supplementary inputs. Restore
all inputs with:

```shell
pixi run python scripts/download_data.py
```

The downloader verifies Zenodo's checksums and places the 24 pickles and three
separate inputs into `data/`. It refuses to overwrite existing files. Use
`--dest` to choose another directory or `--record` for a different published
record ID. The notebooks expect the files in the repository's `data/`
directory. Only load pickles from a trusted source: unpickling can execute
arbitrary code.

Additional inputs currently referenced by the notebooks:

| Input | Current expected location |
| --- | --- |
| `food_and_body_weight.xlsx` | `data/` |
| `operant_behaviour_summary.csv` | `data/` |
| `pellet_retrieval_latency.csv` | `data/` |

## Notebook Path Setup

Canonical helper location is `src/repo_paths.py`.

Use this once at the top of a notebook to get repository-absolute paths:

```python
from pathlib import Path
import runpy

ROOT = next(p for p in [Path.cwd(), *Path.cwd().parents] if (p / "pixi.toml").exists())
globals().update(runpy.run_path(str(ROOT / "src" / "repo_paths.py")))
```

After running that cell, these variables are available:

- `ROOT`
- `DATA_DIR`
- `FIGS_DIR`
- `RESULTS_DIR`
- `SRC_DIR`

Examples:

```python
df_path = DATA_DIR / "operant_behaviour_summary.csv"
out_path = FIGS_DIR / "panels" / "my_output.png"
```

## Rebuilding the photometry pickles from raw tanks (optional)

The figure notebooks do not need raw data. To rebuild the 24 photometry
pickles from the raw TDT tanks, place the tanks in `tanks/` (or pass
`--tanks-dir`) and run:

```shell
pixi run python src/raw_to_pickles.py --tanks-dir <folder with tanks>
```

`src/raw_to_pickles.py` applies the isosbestic correction (trompy
`processdata`) and snips each TTL pellet onset and offset, using the session
metadata in `metadata/`. It writes to `data/` and refuses to overwrite
existing pickles. Events from the first 10 s or last 20 s of a recording are
dropped, and any skipped tank or event type is written to `skipped_log.csv`.
The run fails on skips unless `--allow-skips` is given. Its output matches the
published pickles' event counts and array shapes; snip values agree to within
about 3e-3 but are not bit-identical.

`src/pickles_to_snips.py` loads the pickles and builds the snip arrays used
by the figure notebooks.

