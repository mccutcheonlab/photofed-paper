import numpy as np
import pandas as pd


def add_pellet_number(df, id_col="id", out_col="pellet_number"):
	return df.assign(**{out_col: df.groupby(id_col).cumcount() + 1})


def add_group_max(df, group_cols, value_col, out_col):
	return df.assign(**{out_col: df.groupby(group_cols)[value_col].transform("max")})


def add_group_min(df, group_cols, value_col, out_col):
	return df.assign(**{out_col: df.groupby(group_cols)[value_col].transform("min")})


def add_auc(df, snips, bins=(50, 100), out_col="auc", scale=10):
	aucs = np.trapezoid(snips[:, bins[0] : bins[1]], axis=1) / scale
	return df.assign(**{out_col: aucs})


def get_snips_and_xarray(filename, pellettype, dietphase, nbins=150, zscore=False):
	import dill

	data = dill.load(open(filename, "rb"))
	snips = data["snips"]

	groups = np.array_split(np.arange(snips.shape[1]), nbins)
	snips = np.stack([snips[:, group].mean(axis=1) for group in groups], axis=1)

	if zscore:
		baseline = snips[:, :30]
		snips = (snips - baseline.mean(axis=1, keepdims=True)) / baseline.std(axis=1, keepdims=True)

	x_array = (
		pd.DataFrame(
			{
				"id": data["mouse_ids"],
				"pellettype": pellettype,
				"dietphase": dietphase,
				"pellet_number": np.arange(snips.shape[0]),
			},
			index=np.arange(snips.shape[0]),
		)
		.pipe(add_pellet_number)
		.pipe(add_group_max, group_cols=["id"], value_col="pellet_number", out_col="max_pellets")
		.pipe(add_auc, snips=snips, bins=(50, 100), out_col="auc")
	)

	return snips, x_array


def assemble_snips_and_metadata(session_specs, nbins=150, zscore=False):
	all_snips = []
	all_x_arrays = []

	for spec in session_specs:
		snips, x_array = get_snips_and_xarray(
			spec["filename"],
			spec["pellettype"],
			spec["dietphase"],
			nbins=nbins,
			zscore=zscore,
		)
		all_snips.append(snips)
		all_x_arrays.append(x_array)

	snips = np.concatenate(all_snips, axis=0)
	x_array = (
		pd.concat(all_x_arrays, axis=0)
		.reset_index(drop=True)
		.pipe(
			add_group_min,
			group_cols=["id", "dietphase"],
			value_col="max_pellets",
			out_col="min_pellets_by_session",
		)
		.pipe(
			add_group_min,
			group_cols=["id", "pellettype"],
			value_col="max_pellets",
			out_col="min_pellets_by_pellettype",
		)
		.pipe(
			add_group_min,
			group_cols=["id"],
			value_col="max_pellets",
			out_col="min_pellets_all_sessions",
		)
	)

	return snips, x_array


def mean_snips_by_group(
	x_array,
	snips,
	group_cols,
	dietphase=None,
	pellettype=None,
	extra_query=None,
	equalize_within=None,
):
	if isinstance(group_cols, str):
		group_cols = [group_cols]

	if snips.ndim != 2:
		raise ValueError("snips must be a 2D array")

	if len(x_array) != snips.shape[0]:
		raise ValueError(
			f"Row mismatch: len(x_array)={len(x_array)} vs snips.shape[0]={snips.shape[0]}"
		)

	work = x_array.copy()
	work["_row"] = np.arange(len(work))

	clauses = []
	if dietphase is not None:
		clauses.append("dietphase == @dietphase")
	if pellettype is not None and equalize_within != "pellettype":
		clauses.append("pellettype == @pellettype")
	if extra_query is not None and str(extra_query).strip():
		clauses.append(f"({extra_query})")

	if clauses:
		work = work.query(" and ".join(clauses))

	if work.empty:
		return pd.DataFrame(columns=group_cols), np.empty((0, snips.shape[1]))

	grouped_rows = work.groupby(group_cols, sort=False)["_row"].apply(np.array)

	if len(group_cols) == 1:
		group_df = pd.DataFrame({group_cols[0]: grouped_rows.index.to_numpy()})
	else:
		group_df = grouped_rows.index.to_frame(index=False)

	if len(grouped_rows) == 0:
		return pd.DataFrame(columns=group_cols), np.empty((0, snips.shape[1]))

	mean_traces = np.vstack([snips[idx].mean(axis=0) for idx in grouped_rows.values])

	# Filter to groups with both pellet types if equalizing within pellettype
	if equalize_within == "pellettype":
		if "id" in group_df.columns and "pellettype" in group_df.columns:
			balance_cols = [c for c in ["id", "dietphase"] if c in group_df.columns]

			# Check which balance combos have both grain and sucrose
			combo_has_both = group_df.groupby(balance_cols)["pellettype"].apply(
				lambda x: {"grain", "sucrose"}.issubset(set(x))
			)

			# For each row in group_df, check if its balance combo has both
			mask = group_df.set_index(balance_cols).index.map(combo_has_both).values

			group_df = group_df[mask].reset_index(drop=True)
			mean_traces = mean_traces[mask]

	# Apply pellettype filter if specified (after equalize_within check)
	if pellettype is not None and equalize_within == "pellettype":
		mask = group_df["pellettype"] == pellettype
		group_df = group_df[mask].reset_index(drop=True)
		mean_traces = mean_traces[mask]

	return group_df, mean_traces


def get_paired_animals(x_array, dietphase, extra_query=None):
	"""
	Identify animals with data in both grain and sucrose pellet types.

	Parameters
	----------
	x_array : pd.DataFrame
		Metadata dataframe with id, pellettype, dietphase columns
	dietphase : str
		Diet phase to analyze
	extra_query : str, optional
		Additional filter criteria applied before pairing check

	Returns
	-------
	set
		Animal IDs present in both grain and sucrose for the given diet phase
	"""
	work = x_array[x_array["dietphase"] == dietphase].copy()

	if extra_query is not None and str(extra_query).strip():
		work = work.query(extra_query)

	if work.empty:
		return set()

	grain_ids = set(work[work["pellettype"] == "grain"]["id"].unique())
	sucrose_ids = set(work[work["pellettype"] == "sucrose"]["id"].unique())

	return grain_ids & sucrose_ids


def mean_snips_by_group_paired(
	x_array,
	snips,
	group_cols,
	dietphase=None,
	extra_query=None,
):
	"""
	Get mean snips for paired data across grain and sucrose pellet types.

	Only includes animals that have data in both pellet types for each
	specified diet phase. Returns separate results for each pellet type.

	Parameters
	----------
	x_array : pd.DataFrame
		Metadata dataframe
	snips : np.ndarray
		2D array of snip traces (rows=observations, cols=time bins)
	group_cols : str or list
		Column(s) to group by
	dietphase : str, optional
		Diet phase to analyze. If None, returns results for all diet phases.
	extra_query : str, optional
		Additional query string to filter x_array before pairing check

	Returns
	-------
	dict
		If dietphase specified: {"grain": (df, traces), "sucrose": (df, traces)}
		If dietphase is None: {phase: {"grain": (df, traces), "sucrose": (df, traces)}, ...}
	"""
	if isinstance(group_cols, str):
		group_cols = [group_cols]

	if dietphase is None:
		diet_phases = x_array["dietphase"].unique()
		return {
			phase: mean_snips_by_group_paired(
				x_array,
				snips,
				group_cols,
				dietphase=phase,
				extra_query=extra_query,
			)
			for phase in diet_phases
		}

	paired_animals = get_paired_animals(x_array, dietphase, extra_query)

	if not paired_animals:
		empty_df = pd.DataFrame(columns=group_cols)
		empty_traces = np.empty((0, snips.shape[1]))
		return {
			"grain": (empty_df, empty_traces),
			"sucrose": (empty_df, empty_traces),
		}

	paired_mask = x_array["id"].isin(paired_animals)
	x_filtered = x_array[paired_mask].reset_index(drop=True)
	snips_filtered = snips[paired_mask]

	grain_df, grain_traces = mean_snips_by_group(
		x_filtered,
		snips_filtered,
		group_cols,
		dietphase=dietphase,
		pellettype="grain",
		extra_query=extra_query,
	)

	sucrose_df, sucrose_traces = mean_snips_by_group(
		x_filtered,
		snips_filtered,
		group_cols,
		dietphase=dietphase,
		pellettype="sucrose",
		extra_query=extra_query,
	)

	return {
		"grain": (grain_df, grain_traces),
		"sucrose": (sucrose_df, sucrose_traces),
	}

def reduce_to_matched_mice(
    x_array,
    snips,
    extra_query,
    phase_col="dietphase",
    pellet_col="pellettype",
    id_col="id",
    require_both_pellet_types=True,
):
    """Apply extra_query and optionally enforce IDs present in both pellet types per phase."""
    extra_mask = x_array.eval(extra_query)
    x_extra = x_array.loc[extra_mask].copy()
    snips_extra = snips[extra_mask.to_numpy(), :]

    before_counts = (
        x_extra.groupby([phase_col, pellet_col], sort=False)[id_col]
        .nunique()
        .rename("n_ids_before")
    )

    if require_both_pellet_types:
        common_ids_by_phase = {}
        for phase in x_extra[phase_col].unique():
            grain_ids = set(
                x_extra.query(f"{phase_col} == @phase and {pellet_col} == 'grain'")[id_col].unique()
            )
            sucrose_ids = set(
                x_extra.query(f"{phase_col} == @phase and {pellet_col} == 'sucrose'")[id_col].unique()
            )
            common_ids_by_phase[phase] = grain_ids.intersection(sucrose_ids)

        matched_mask = x_extra.apply(
            lambda r: r[id_col] in common_ids_by_phase.get(r[phase_col], set()), axis=1
        )
    else:
        matched_mask = pd.Series(True, index=x_extra.index)

    x_reduced = x_extra.loc[matched_mask].copy()
    snips_red = snips_extra[matched_mask.to_numpy(), :]

    after_counts = (
        x_reduced.groupby([phase_col, pellet_col], sort=False)[id_col]
        .nunique()
        .rename("n_ids_after")
    )

    match_summary = (
        pd.concat([before_counts, after_counts], axis=1)
        .fillna(0)
        .astype(int)
        .reset_index()
    )

    return x_reduced, snips_red, match_summary


def calculate_grain_sucrose_difference(
    auc_df,
    id_col: str = "id",
    phase_col: str = "dietphase",
    pellet_col: str = "pellettype",
    value_col: str = "auc",
) -> dict:
	"""
	Calculate grain/sucrose difference per rat, then aggregate by diet phase.

	Each rat's difference is computed independently, treating each rat as one
	observation. The phase-level statistics are then computed from these
	per-rat differences.

	Parameters
	----------
	auc_df : pd.DataFrame
		DataFrame with columns for id, dietphase, pellettype, and value
	id_col : str
		Column name for animal ID (default: "id")
	phase_col : str
		Column name for diet phase (default: "dietphase")
	pellet_col : str
		Column name for pellet type (default: "pellettype")
	value_col : str
		Column name for the value to ratio (default: "auc")

	Returns
	-------
	dict
		'diff_points': DataFrame with per-rat differences indexed by id and phase
		'phase_diff': Series of mean differences by phase
		'diff_sem': Series of SEM of differences by phase
	"""


	diff_points = (
		auc_df.pivot_table(
			index=[id_col, phase_col],
			columns=pellet_col,
			values=value_col,
			aggfunc="mean",
		)
		.dropna(subset=["grain", "sucrose"])
		.assign(diff=lambda d: d["grain"] - d["sucrose"])
	)


	phase_diff = diff_points.groupby(phase_col, sort=False)["diff"].mean()
	diff_sem = diff_points.groupby(phase_col, sort=False)["diff"].sem()

	return {
		"diff_points": diff_points,
		"phase_diff": phase_diff,
		"diff_sem": diff_sem,
	}
