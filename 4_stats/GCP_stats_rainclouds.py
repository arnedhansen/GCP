# %% GCP Stats Rainclouds — trial-level GED, gaze, and behavior
# Raincloud figures with half-kernel densities, boxplots, jittered trials.
# Significance brackets are optional (commented out by default).
#
# Run:
#   python GCP_stats_rainclouds.py

# %% Imports
import os
import sys
import warnings
from typing import Optional

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.io
from scipy.stats import gaussian_kde
from statsmodels.tools.sm_exceptions import ConvergenceWarning

sys.path.insert(0, os.path.expanduser("/Users/Arne/Documents/GitHub"))

from functions.rainclouds_plotting_helpers import add_stat_brackets
from functions.stats_helpers import iqr_outlier_filter, mixedlm_pairwise_contrasts, p_to_signif

warnings.filterwarnings("ignore", category=ConvergenceWarning, module="statsmodels")

# %% GCP colours (color_def.m)
PAL = ["#FFE680", "#E69966", "#E66666", "#000000"]
CONDITION_ORDER = ["25%", "50%", "75%", "100%"]
COND_NUM_TO_LABEL = {1: "25%", 2: "50%", 3: "75%", 4: "100%"}
COMPARISONS = [
    ("25%", "50%"),
    ("25%", "75%"),
    ("25%", "100%"),
    ("50%", "75%"),
    ("50%", "100%"),
    ("75%", "100%"),
]

VARIABLES = [
    ("GammaFrequency", "Peak Gamma Frequency [Hz]", "gamma_freq"),
    ("GammaPower", "Peak Gamma Power [dB]", "gamma_power"),
    ("dBMSRate", "Microsaccade Rate [dB]", "ms"),
    ("dBPupilSize", "Pupil Size [dB]", "pupil"),
    ("dBVelV", "Gaze Velocity Y [dB]", "vely"),
    ("ReactionTime", "Reaction Time [s]", "rt"),
]

FIGURE_SAVE_DPI = 600
YLABEL_GRID_X = -0.15
YLIM_PAD_FRAC = 0.06
SHOW_SIGNIFICANCE_BRACKETS = False  # set True to restore MixedLM brackets

mpl.rcParams.update({
    "figure.dpi": 160,
    "savefig.dpi": FIGURE_SAVE_DPI,
    "savefig.transparent": False,
    "savefig.facecolor": "white",
    "savefig.bbox": None,
    "ps.fonttype": 42,
    "font.size": 18,
    "axes.titlesize": 18,
    "axes.labelsize": 15,
    "legend.fontsize": 15,
    "xtick.labelsize": 15,
    "ytick.labelsize": 15,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "mathtext.default": "regular",
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "figure.edgecolor": "white",
    "axes.edgecolor": "white",
})


def _as_str_list(values) -> list[str]:
    if isinstance(values, str):
        return [values]
    return [str(v) for v in np.atleast_1d(values).tolist()]


def reconstruct_trial_peak_power(
    peak_freq: np.ndarray,
    powratio: np.ndarray,
    scan_freqs: np.ndarray,
    halfwidth_hz: float = 5.0,
) -> np.ndarray:
    n_trl = peak_freq.shape[0]
    power = np.full(n_trl, np.nan, dtype=float)
    scan_freqs = np.asarray(scan_freqs, dtype=float).ravel()
    for t in range(n_trl):
        if not np.isfinite(peak_freq[t]):
            continue
        band = np.abs(scan_freqs - peak_freq[t]) <= halfwidth_hz
        if not np.any(band):
            continue
        power[t] = np.nanmean(powratio[t, band])
    return power


def load_ged_trial_metrics(mat_path: str, peak_power_halfwidth_hz: float = 5.0) -> pd.DataFrame:
    ged = scipy.io.loadmat(mat_path, squeeze_me=True, struct_as_record=False)
    trials_peaks = ged["trials_peaks"]
    trials_powratio = ged["trials_powratio_fullscan"]
    scan_freqs = np.asarray(ged["scan_freqs"], dtype=float).ravel()
    subjects = _as_str_list(ged["subjects"])

    power_masks = ged.get("trials_outlier_mask_power_full")
    n_cond, n_subj = trials_peaks.shape
    rows: list[dict] = []

    for c in range(n_cond):
        for s in range(n_subj):
            pf = trials_peaks[c, s]
            pr = trials_powratio[c, s]
            if pf is None or pr is None:
                continue

            pf = np.asarray(pf, dtype=float).ravel()
            pr = np.asarray(pr, dtype=float)
            if pr.ndim == 1:
                pr = pr.reshape(1, -1)
            if pf.size == 0 or pr.shape[0] != pf.size:
                continue

            pp = reconstruct_trial_peak_power(pf, pr, scan_freqs, peak_power_halfwidth_hz)
            if power_masks is not None:
                mask = power_masks[c, s]
                if mask is not None:
                    mask = np.asarray(mask, dtype=bool).ravel()
                    if mask.size == pf.size:
                        pp[mask] = np.nan

            cond_label = CONDITION_ORDER[c] if c < len(CONDITION_ORDER) else str(c + 1)
            sid = subjects[s] if s < len(subjects) else str(s + 1)
            for trial_idx in range(pf.size):
                if not np.isfinite(pf[trial_idx]) and not np.isfinite(pp[trial_idx]):
                    continue
                rows.append(
                    {
                        "ID": sid,
                        "Condition": cond_label,
                        "Trial": trial_idx + 1,
                        "GammaFrequency": float(pf[trial_idx]) if np.isfinite(pf[trial_idx]) else np.nan,
                        "GammaPower": float(pp[trial_idx]) if np.isfinite(pp[trial_idx]) else np.nan,
                    }
                )

    return pd.DataFrame(rows)


def load_merged_trial_metrics(csv_path: str) -> pd.DataFrame:
    dat = pd.read_csv(csv_path)
    dat["ID"] = dat["ID"].astype(str)
    dat["Condition"] = dat["Condition"].map(COND_NUM_TO_LABEL)
    dat = dat.rename(
        columns={
            "Gaze_dBMSRate": "dBMSRate",
            "Gaze_dBPupilSize": "dBPupilSize",
            "Gaze_dBVelV": "dBVelV",
            "Behavior_ReactionTime": "ReactionTime",
        }
    )
    return dat


def filter_gcp_analysis_cohort(dat: pd.DataFrame, controls_dir: str) -> pd.DataFrame:
    if "Include" in dat.columns:
        out = dat.loc[dat["Include"].astype(bool)].copy()
        print(f"GED cohort filter (Include column): {out['ID'].nunique()} subjects kept.")
        return out

    inclusion_path = os.path.join(controls_dir, "GCP_subject_inclusion.mat")
    if not os.path.isfile(inclusion_path):
        print(f"WARNING: inclusion file not found: {inclusion_path}")
        return dat

    inc = scipy.io.loadmat(inclusion_path, squeeze_me=True, struct_as_record=False)
    tbl = inc["subject_inclusion"]
    subj_ids = np.atleast_1d(tbl.SubjID).astype(int)
    include = np.atleast_1d(tbl.Include).astype(bool)
    included = {str(sid) for sid, flag in zip(subj_ids, include) if flag}
    out = dat.loc[dat["ID"].isin(included)].copy()
    print(f"GED cohort filter (controls mat): {out['ID'].nunique()} subjects kept.")
    return out


def label_condition(dat: pd.DataFrame) -> pd.DataFrame:
    out = dat.copy()
    out["Condition"] = pd.Categorical(out["Condition"], categories=CONDITION_ORDER, ordered=True)
    return out


def data_ylim(yvals: np.ndarray, pad_frac: float = YLIM_PAD_FRAC) -> tuple[float, float]:
    ymin = float(np.nanmin(yvals))
    ymax = float(np.nanmax(yvals))
    yr = ymax - ymin
    if not np.isfinite(yr) or yr == 0:
        yr = max(abs(ymax), abs(ymin), 1.0) * 0.1
    pad = yr * pad_frac
    return ymin - pad, ymax + pad


def bracket_labels_from_mixedlm(dvar: pd.DataFrame, var: str, comparisons: list) -> list[str]:
    df = dvar.rename(columns={var: "value"}).copy()
    try:
        pw = mixedlm_pairwise_contrasts(
            df, value_col="value", group_col="Condition", id_col="ID", p_adjust="fdr_bh"
        )
    except Exception as exc:
        print(f"WARNING: MixedLM brackets failed for {var}: {exc}")
        return ["n.s."] * len(comparisons)

    labels = []
    for g1, g2 in comparisons:
        row = pw.loc[(pw["group1"] == g1) & (pw["group2"] == g2)]
        labels.append("n.s." if row.empty else p_to_signif(float(row["p_adj"].iloc[0])))
    return labels


def plot_raincloud(
    dvar: pd.DataFrame,
    var: str,
    ylab: str,
    sname: str,
    pal_dict: dict,
    output_dir: str,
) -> None:
    fig, ax = plt.subplots(figsize=(8, 6), facecolor="white")
    ax.set_facecolor("white")

    viol_alpha = 0.60
    dot_alpha = 0.18
    dot_size = 24
    box_width = 0.20
    cloud_offset = -0.20
    max_violsw = 0.40
    bw_method = 0.15

    xpos = {c: i for i, c in enumerate(CONDITION_ORDER)}
    rng = np.random.default_rng(12345)

    yvals_all = dvar[var].dropna().to_numpy()
    ymin_plot, ymax_plot = data_ylim(yvals_all)
    kde_ymin, kde_ymax = ymin_plot, ymax_plot

    for cond_lab in CONDITION_ORDER:
        yvals = dvar.loc[dvar["Condition"] == cond_lab, var].dropna().to_numpy()
        if yvals.size == 0:
            continue

        kde = gaussian_kde(yvals, bw_method=bw_method)
        y_grid = np.linspace(kde_ymin, kde_ymax, 400)
        dens = kde(y_grid)
        scale = (max_violsw / np.nanmax(dens)) if np.nanmax(dens) > 0 else 0.0

        x_left = xpos[cond_lab] + cloud_offset - dens * scale
        x_right = np.full_like(y_grid, xpos[cond_lab] + cloud_offset)
        poly_x = np.concatenate([x_right, x_left[::-1]])
        poly_y = np.concatenate([y_grid, y_grid[::-1]])
        ax.fill(
            poly_x,
            poly_y,
            facecolor=pal_dict[cond_lab],
            edgecolor="none",
            alpha=viol_alpha,
            clip_on=True,
        )

        x_jit = xpos[cond_lab] + rng.uniform(-box_width / 2, box_width / 2, size=yvals.size)
        ax.scatter(
            x_jit,
            yvals,
            s=dot_size,
            alpha=dot_alpha,
            color=pal_dict[cond_lab],
            linewidths=0,
            zorder=3,
        )

        bp = ax.boxplot(
            [yvals],
            positions=[xpos[cond_lab]],
            widths=box_width,
            vert=True,
            patch_artist=True,
            showfliers=False,
            whis=(5, 95),
            medianprops=dict(color="black", linewidth=1.5),
            boxprops=dict(linewidth=1.0, edgecolor="black"),
            whiskerprops=dict(linewidth=1.0, color="black"),
            capprops=dict(linewidth=1.0, color="black"),
            meanline=False,
            showmeans=False,
        )
        for patch in bp["boxes"]:
            patch.set_facecolor(mpl.colors.to_rgba(pal_dict[cond_lab], 0.05))
            patch.set_edgecolor("black")

    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.yaxis.grid(True, linewidth=1, alpha=0.35)
    ax.xaxis.grid(False)

    ax.set_title("")
    ax.set_xticks(range(len(CONDITION_ORDER)))
    ax.set_xticklabels(CONDITION_ORDER)
    ax.set_xlabel("")
    ax.annotate(
        "Contrast",
        xy=(xpos[CONDITION_ORDER[1]], 0),
        xycoords=("data", "axes fraction"),
        xytext=(0, -28),
        textcoords="offset points",
        ha="center",
        va="top",
    )

    ymid = 0.5 * (ymin_plot + float(np.nanmax(yvals_all)))
    ax.set_ylabel("")
    ax.yaxis.get_label().set_visible(False)
    ax.text(
        YLABEL_GRID_X,
        ymid,
        ylab,
        transform=ax.get_yaxis_transform(which="grid"),
        rotation=90,
        ha="center",
        va="center",
    )

    range_y = ymax_plot - ymin_plot
    step = 0.10 * range_y
    y_positions = [float(np.nanmax(yvals_all)) + range_y * YLIM_PAD_FRAC + i * step for i in range(len(COMPARISONS))]

    # Significance brackets (re-enable with SHOW_SIGNIFICANCE_BRACKETS = True)
    if SHOW_SIGNIFICANCE_BRACKETS:
        bracket_extra = range_y * 0.12
        ymax_plot = ymax_plot + bracket_extra + (len(COMPARISONS) - 1) * (0.10 * range_y)
        labels = bracket_labels_from_mixedlm(dvar, var, COMPARISONS)
        add_stat_brackets(
            ax=ax,
            xcats=CONDITION_ORDER,
            comparisons=COMPARISONS,
            y_positions=y_positions,
            labels=labels,
            xmap=xpos,
        )

    ax.set_ylim(ymin_plot, ymax_plot)
    ax.set_xlim(-0.6, len(CONDITION_ORDER) - 0.4)

    fig.tight_layout()
    fig.subplots_adjust(left=0.17)
    out_path = os.path.join(output_dir, f"GCP_stats_rainclouds_{sname}.png")
    fig.savefig(
        out_path,
        dpi=FIGURE_SAVE_DPI,
        transparent=False,
        facecolor=fig.get_facecolor(),
        edgecolor="white",
    )
    plt.close(fig)
    print(f"Saved raincloud fig. -> {os.path.basename(out_path)}")


def main() -> None:
    base_dir = "/Volumes/g_psyplafor_methlab$/Students/Arne/GCP"
    features_dir = os.path.join(base_dir, "data", "features")
    controls_dir = os.path.join(base_dir, "data", "controls")
    output_dir = os.path.join(base_dir, "figures", "stats", "rainclouds")
    os.makedirs(output_dir, exist_ok=True)

    merged_csv = os.path.join(features_dir, "GCP_merged_data_trials.csv")
    ged_mat = os.path.join(features_dir, "GCP_eeg_GED.mat")

    if not os.path.isfile(merged_csv):
        raise FileNotFoundError(f"Merged trial CSV not found: {merged_csv}")
    if not os.path.isfile(ged_mat):
        raise FileNotFoundError(f"GED MAT file not found: {ged_mat}")

    merged = filter_gcp_analysis_cohort(
        label_condition(load_merged_trial_metrics(merged_csv)),
        controls_dir,
    )
    ged = filter_gcp_analysis_cohort(
        label_condition(load_ged_trial_metrics(ged_mat)),
        controls_dir,
    )

    data_by_var = {
        "GammaFrequency": ged,
        "GammaPower": ged,
        "dBMSRate": merged,
        "dBPupilSize": merged,
        "dBVelV": merged,
        "ReactionTime": merged,
    }

    pal_dict = dict(zip(CONDITION_ORDER, PAL))

    for var, ylab, sname in VARIABLES:
        src = data_by_var[var]
        if var not in src.columns:
            print(f"WARNING: {var} missing; skipping.")
            continue

        dvar = src.loc[src[var].notna(), ["ID", "Condition", var]].copy()
        if dvar.empty:
            print(f"WARNING: no data for {var}; skipping.")
            continue

        dvar = iqr_outlier_filter(dvar, [var], by="Condition")
        dvar = label_condition(dvar)
        dvar = dvar.loc[dvar[var].notna()].copy()
        if dvar.empty:
            print(f"WARNING: no data for {var} after outlier filtering; skipping.")
            continue

        print(f"{var}: n={len(dvar)}, subjects={dvar['ID'].nunique()}")
        plot_raincloud(dvar, var, ylab, sname, pal_dict, output_dir)

    print("\nRaincloud figures saved to:", output_dir)


if __name__ == "__main__":
    main()
