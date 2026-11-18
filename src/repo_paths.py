"""Repository-relative path helpers shared across scripts and notebooks."""

from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"
FIGS_DIR = ROOT / "figs"
RESULTS_DIR = ROOT / "results"
SRC_DIR = ROOT / "src"
STATS_DIR = ROOT / "stats"


def repo_path(*parts: str) -> Path:
    """Return a repository-absolute path by joining to ROOT."""
    return ROOT.joinpath(*parts)


# Ensure imports from the repository and src/ are available.
for p in (ROOT, SRC_DIR):
    p_str = str(p)
    if p_str not in sys.path:
        sys.path.insert(0, p_str)
