"""Tests for importing the publication plotting helpers."""

import runpy
import sys
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch


class FigurePlottingImportTests(unittest.TestCase):
    """Keep publication plotting independent of exploratory scripts."""

    def test_import_without_legacy_timing_script(self) -> None:
        """The plotting helpers can load without the legacy timing module."""
        plotting = Path(__file__).resolve().parents[1] / "src" / "figure_plotting.py"
        dependencies: dict[str, object] = {
            name: MagicMock()
            for name in (
                "matplotlib",
                "matplotlib.pyplot",
                "matplotlib.colors",
                "matplotlib.transforms",
                "numpy",
                "seaborn",
            )
        }
        with patch.dict(sys.modules, dependencies):
            module = runpy.run_path(str(plotting))
        self.assertTrue(callable(module["make_auc_plot"]))


if __name__ == "__main__":
    unittest.main()
