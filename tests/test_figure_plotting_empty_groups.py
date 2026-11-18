"""Tests for photometry plots with empty phase/pellet groups."""

import warnings
import unittest

import matplotlib
matplotlib.use("Agg")

import numpy as np
import pandas as pd

from src.figure_plotting import make_heatplot_and_lineplot_per_session


class EmptyPlotGroupTests(unittest.TestCase):
    """Keep empty exploratory plot groups explicit and warning-free."""

    def test_empty_group_is_labeled_without_mean_warnings(self) -> None:
        """A missing pellet group is labeled rather than plotted as an empty heatmap."""
        metadata = pd.DataFrame(
            {
                "pellettype": ["grain"],
                "dietphase": ["NR2"],
                "max_pellets": [7],
            }
        )
        snips = np.ones((1, 150))

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            figure, axes, _ = make_heatplot_and_lineplot_per_session(
                snips, metadata, session="NR2", max_pellets_filter=6
            )

        self.addCleanup(matplotlib.pyplot.close, figure)
        self.assertEqual(axes[1].texts[0].get_text(), "No data")
        self.assertFalse(caught)


if __name__ == "__main__":
    unittest.main()
