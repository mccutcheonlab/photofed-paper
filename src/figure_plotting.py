import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.transforms import blended_transform_factory

FIG4_DIVERGING_CMAP = LinearSegmentedColormap.from_list(
	"custom_diverging",
	["#016895", "white", "#C74632"],
	N=256,
)

FIG4_PALETTE_OPTIONS = {
	"option_1_moss_ochre": {
		"NR": {"grain": "#4F772D", "sucrose": "#90A955"},
		"PR": {"grain": "#A06A00", "sucrose": "#E0B44C"},
	},
	"option_2_pine_sand": {
		"NR": {"grain": "#2D6A4F", "sucrose": "#95D5B2"},
		"PR": {"grain": "#8C613C", "sucrose": "#DDB892"},
	},
	"option_3_olive_clay": {
		"NR": {"grain": "#556B2F", "sucrose": "#A3B18A"},
		"PR": {"grain": "#9A6324", "sucrose": "#DDBEA9"},
	},
	"option_4_grey_reddish": {
		"NR": {"grain": "#8C8C8C", "sucrose": "#C2C2C2"},
		"PR": {"grain": "#C74632", "sucrose": "#DE9387"},
	},
}


def get_palette(active_palette_name="option_1_moss_ochre"):
	active_palette = FIG4_PALETTE_OPTIONS[active_palette_name]
	colors_nr = [active_palette["NR"]["grain"], active_palette["NR"]["sucrose"]]
	colors_pr = [active_palette["PR"]["grain"], active_palette["PR"]["sucrose"]]
	phase_pellet_colors = {
		("NR", "grain"): active_palette["NR"]["grain"],
		("NR", "sucrose"): active_palette["NR"]["sucrose"],
		("PR", "grain"): active_palette["PR"]["grain"],
		("PR", "sucrose"): active_palette["PR"]["sucrose"],
	}
	return active_palette, colors_nr, colors_pr, phase_pellet_colors


def sort_by_signal(snips, window=(50, 100)):
	signal = snips[:, window[0] : window[1]].mean(axis=1)
	sorted_indices = np.argsort(signal)
	return snips[sorted_indices]


def make_heatplot_and_lineplot_per_session(
	snips,
	x_array,
	session,
	vmin=-1,
	vmax=1,
	cmap=FIG4_DIVERGING_CMAP,
	max_pellets_filter=999,
):
	figure, axes = plt.subplots(
		nrows=3,
		figsize=(1.8, 4),
		sharex=True,
		gridspec_kw={"height_ratios": [1, 1, 0.5]},
	)
	cbar_figure, cbar_axis = plt.subplots(figsize=(0.2, 0.5))

	grain_snips = snips[
		x_array.query(
			"pellettype == 'grain' and dietphase == @session and max_pellets > @max_pellets_filter"
		).index,
		:,
	]
	if len(grain_snips):
		sns.heatmap(sort_by_signal(grain_snips), cmap=cmap, vmin=vmin, vmax=vmax, ax=axes[0], cbar_ax=cbar_axis)
		axes[2].plot(grain_snips.mean(axis=0), color="k")
	else:
		axes[0].text(0.5, 0.5, "No data", ha="center", va="center", transform=axes[0].transAxes)

	sucrose_snips = snips[
		x_array.query(
			"pellettype == 'sucrose' and dietphase == @session and max_pellets > @max_pellets_filter"
		).index,
		:,
	]
	if len(sucrose_snips):
		sns.heatmap(sort_by_signal(sucrose_snips), cmap=cmap, vmin=vmin, vmax=vmax, ax=axes[1], cbar_ax=cbar_axis)
		axes[2].plot(sucrose_snips.mean(axis=0), color="r")
	else:
		axes[1].text(0.5, 0.5, "No data", ha="center", va="center", transform=axes[1].transAxes)

	for axis in axes[:2]:
		axis.invert_yaxis()

	for axis in axes:
		axis.set_xticks([])
		axis.set_yticks([])

	sns.despine(ax=axes[2], left=True, bottom=True)
	return figure, axes, cbar_axis


def make_heatplot_and_lineplot_per_session_mean_of_rats(
	grain_snips,
	sucrose_snips,
	vmin=-1,
	vmax=1,
	cmap=FIG4_DIVERGING_CMAP,
	colors=("k", "r"),
	sort_snips=True,
):
	figure, axes = plt.subplots(
		nrows=3,
		figsize=(1.8, 4),
		sharex=True,
		gridspec_kw={"height_ratios": [1, 1, 0.8]},
	)
	cbar_figure, cbar_axis = plt.subplots(figsize=(0.2, 1.5))

	if sort_snips:
		grain_snips = sort_by_signal(grain_snips)
		sucrose_snips = sort_by_signal(sucrose_snips)

	sns.heatmap(grain_snips, cmap=cmap, vmin=vmin, vmax=vmax, ax=axes[0], cbar_ax=cbar_axis)
	axes[2].plot(grain_snips.mean(axis=0), color=colors[0])
	axes[2].fill_between(
		np.arange(grain_snips.shape[1]),
		grain_snips.mean(axis=0) - (grain_snips.std(axis=0) / np.sqrt(grain_snips.shape[0])),
		grain_snips.mean(axis=0) + grain_snips.std(axis=0) / np.sqrt(grain_snips.shape[0]),
		color=colors[0],
		alpha=0.3,
	)

	sns.heatmap(sucrose_snips, cmap=cmap, vmin=vmin, vmax=vmax, ax=axes[1], cbar_ax=cbar_axis)
	axes[2].plot(sucrose_snips.mean(axis=0), color=colors[1])
	axes[2].fill_between(
		np.arange(sucrose_snips.shape[1]),
		sucrose_snips.mean(axis=0) - (sucrose_snips.std(axis=0) / np.sqrt(sucrose_snips.shape[0])),
		sucrose_snips.mean(axis=0) + sucrose_snips.std(axis=0) / np.sqrt(sucrose_snips.shape[0]),
		color=colors[1],
		alpha=0.3,
	)

	axes[2].axvspan(50, 100, color="grey", alpha=0.1)

	for axis in axes[:2]:
		axis.invert_yaxis()

	for axis in axes:
		axis.set_xticks([])
		axis.set_yticks([])

	axes[2].plot(
		[125, 150],
		[0, 0],
		color="k",
		linestyle="-",
		linewidth=2,
		alpha=1,
		transform=blended_transform_factory(axes[2].transData, axes[2].transAxes),
	)
 
	axes[2].plot(
		[0, 0],
		[0.5, 1.5],
		color="k",
		linestyle="-",
		linewidth=2,
		alpha=1,
	)
	
	sns.despine(ax=axes[2], left=True, bottom=True)
	cbar_axis.set_yticks([])
	return figure, axes, cbar_figure

def make_auc_plot(auc_df, colors_nr=None, colors_pr=None, phases=["NR1", "PR", "NR2"]):
	fig, ax = plt.subplots(
		figsize=(3, 1.8), gridspec_kw={"left": 0.3, "bottom": 0.25}
	)

	phase_spacing = 2.4
	phases = auc_df["dietphase"].unique()
	pellets = auc_df["pellettype"].unique()
	jitter = 0.08
	rng = np.random.default_rng(42)

	for i, phase in enumerate(phases):
		if phase.startswith("NR"):
			colors = colors_nr
		else:
			colors = colors_pr
		x_jittered_both, yvals_both = [], []
		for j, pellet in enumerate(pellets):
			mask = (auc_df["dietphase"] == phase) & (auc_df["pellettype"] == pellet)
			x_center = i * phase_spacing + j
			y_vals = auc_df.loc[mask, "auc"].to_numpy()

			ax.bar(x_center, y_vals.mean(), color=colors[j], alpha=1)

			x_jittered = x_center + rng.uniform(-jitter, jitter, size=len(y_vals))
			ax.scatter(x_jittered, y_vals, s=14, color="white", edgecolor=colors[j], linewidth=0.7, zorder=3, alpha=0.5)

			yvals_both.append(y_vals)
			x_jittered_both.append(x_jittered)

		ax.plot(x_jittered_both, yvals_both, color="0.5", alpha=0.3, linewidth=0.75, zorder=1)

	ax.set_xticks([i * phase_spacing + 0.5 for i in range(len(phases))])
	ax.set_xticklabels(phases)
	ax.set_ylabel("Dopamine (AUC)")
	sns.despine(ax=ax, offset=5)

	return fig, ax


def make_difference_plot(diff_stats, colors=None, phases=["NR1", "PR", "NR2"]):
	phase_diff = diff_stats["phase_diff"].reindex(phases)
	sem = diff_stats["diff_sem"].reindex(phases).fillna(0)

	fig, ax = plt.subplots(
		figsize=(1.8, 1.8), gridspec_kw={"left": 0.35, "bottom": 0.25}
	)
	xpos = np.arange(len(phase_diff))

	# Error bars extend in the direction of the mean (as non-negative magnitudes)
	lower_err = np.where(phase_diff.values > 0, 0, sem.to_numpy())
	upper_err = np.where(phase_diff.values > 0, sem.to_numpy(), 0)
	ratio_yerr = np.vstack([lower_err, upper_err])

	ax.bar(
		xpos,
		phase_diff.values,
		yerr=ratio_yerr,
		color=colors,
		ecolor="k",
		capsize=0,
	)

	ax.axhline(0, color="0.4", linestyle="--", linewidth=1, zorder=-20)
	ax.set_xticks(xpos)
	ax.set_xticklabels(phase_diff.index)
	ax.set_ylabel("Chow - Sucrose")

	sns.despine(ax=ax, offset=5)
 
	print(phase_diff)
 
	return fig, ax

 
 