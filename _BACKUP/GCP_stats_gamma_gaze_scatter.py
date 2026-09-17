# %% GCP Gamma–Gaze Contrast Scatter
# Subject-level scatter of change (100% − 25% contrast) in gamma metrics vs.
# oculomotor metrics. Grey shading marks the hypothesis-consistent quadrant:
#   gamma frequency/power increase with contrast (Δ > 0 on x)
#   microsaccade rate decreases with contrast (Δ < 0 on y)
#   eye velocity increases with contrast (Δ > 0 on y)
#
# Gamma values are peaks from condition-averaged GED spectra (not trial-peak
# means). Gaze values are subject-condition means of trial metrics.
#
# Run:
#   python GCP_stats_gamma_gaze_scatter.py

from __future__ import annotations

import os
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

sys.path.insert(0, os.path.expanduser("/Users/Arne/Documents/GitHub"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from functions.stats_helpers import iqr_outlier_filter

from GCP_stats_rainclouds import (
    CONDITION_ORDER,
    filter_gcp_analysis_cohort,
    label_condition,
    load_ged_condition_metrics,
    load_merged_trial_metrics,
)

LOW_CONTRAST = "25%"
HIGH_CONTRAST = "100%"

FIGURE_SAVE_DPI = 600
HYPOTHESIS_FACE_COLOR = (0.5, 0.5, 0.5)
HYPOTHESIS_ALPHA = 0.30

PANELS = [
    {
        "x": "GammaFrequency",
        "y": "MSRate_bl",
        "xlab": r"$\Delta$ Peak Gamma Frequency [Hz] (100% $-$ 25%)",
        "ylab": r"$\Delta$ Microsaccade Rate [%] (100% $-$ 25%)",
        "hypothesis": "x_pos_y_neg",
        "title": "Gamma Frequency vs. Microsaccades",
    },
    {
        "x": "GammaPower",
        "y": "MSRate_bl",
        "xlab": r"$\Delta$ Peak Gamma Power [dB] (100% $-$ 25%)",
        "ylab": r"$\Delta$ Microsaccade Rate [%] (100% $-$ 25%)",
        "hypothesis": "x_pos_y_neg",
        "title": "Gamma Power vs. Microsaccades",
    },
    {
        "x": "GammaFrequency",
        "y": "BCEA_bl",
        "xlab": r"$\Delta$ Peak Gamma Frequency [Hz] (100% $-$ 25%)",
        "ylab": r"$\Delta$ BCEA [%] (100% $-$ 25%)",
        "hypothesis": "x_pos_y_pos",
        "title": "Gamma Frequency vs. BCEA",
    },
    {
        "x": "GammaPower",
        "y": "BCEA_bl",
        "xlab": r"$\Delta$ Peak Gamma Power [dB] (100% $-$ 25%)",
        "ylab": r"$\Delta$ BCEA [%] (100% $-$ 25%)",
        "hypothesis": "x_pos_y_pos",
        "title": "Gamma Power vs. BCEA",
    },
    {
        "x": "GammaFrequency",
        "y": "Vel2D_bl",
        "xlab": r"$\Delta$ Peak Gamma Frequency [Hz] (100% $-$ 25%)",
        "ylab": r"$\Delta$ Eye Velocity [%] (100% $-$ 25%)",
        "hypothesis": "x_pos_y_pos",
        "title": "Gamma Frequency vs. Eye Velocity",
    },
    {
        "x": "GammaPower",
        "y": "Vel2D_bl",
        "xlab": r"$\Delta$ Peak Gamma Power [dB] (100% $-$ 25%)",
        "ylab": r"$\Delta$ Eye Velocity [%] (100% $-$ 25%)",
        "hypothesis": "x_pos_y_pos",
        "title": "Gamma Power vs. Eye Velocity",
    },
]

mpl.rcParams.update({
    "figure.dpi": 160,
    "savefig.dpi": FIGURE_SAVE_DPI,
    "savefig.transparent": False,
    "savefig.facecolor": "white",
    "ps.fonttype": 42,
    "font.size": 14,
    "axes.titlesize": 15,
    "axes.labelsize": 13,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
})


def _rename_gaze_columns(dat: pd.DataFrame) -> pd.DataFrame:
    out = dat.copy()
    rename_map = {
        "Gaze_MSRate_bl": "MSRate_bl",
        "Gaze_BCEA_bl": "BCEA_bl",
        "Gaze_Vel2D_bl": "Vel2D_bl",
        "Gaze_VelV_bl": "VelV_bl",
        "GED_GammaFrequency": "GammaFrequency",
        "GED_GammaPower": "GammaPower",
    }
    for old, new in rename_map.items():
        if old in out.columns and new not in out.columns:
            out = out.rename(columns={old: new})
    return out


def _trial_table_with_outliers(
    dat: pd.DataFrame,
    variables: list[str],
) -> pd.DataFrame:
    keep = ["ID", "Condition", *variables]
    missing = [v for v in variables if v not in dat.columns]
    if missing:
        raise KeyError(f"Missing columns in trial table: {missing}")

    out = dat.loc[:, keep].copy()
    out["ID"] = out["ID"].astype(str)
    out = label_condition(out)

    for var in variables:
        dvar = out.loc[out[var].notna(), ["ID", "Condition", var]].copy()
        if dvar.empty:
            continue
        dvar = iqr_outlier_filter(dvar, [var], by="Condition")
        keep_mask = dvar[var].notna()
        drop_idx = dvar.index[~keep_mask]
        out.loc[drop_idx, var] = np.nan

    return out


def _subject_condition_means(dat: pd.DataFrame, variables: list[str]) -> pd.DataFrame:
    rows: list[dict] = []
    for sid, sub_df in dat.groupby("ID", sort=False):
        for cond in CONDITION_ORDER:
            cond_df = sub_df.loc[sub_df["Condition"] == cond]
            if cond_df.empty:
                continue
            row = {"ID": sid, "Condition": cond}
            for var in variables:
                row[var] = float(np.nanmean(cond_df[var].to_numpy(dtype=float)))
            rows.append(row)
    return pd.DataFrame(rows)


def _contrast_delta(means: pd.DataFrame, variables: list[str]) -> pd.DataFrame:
    rows: list[dict] = []
    for sid, sub_df in means.groupby("ID", sort=False):
        low = sub_df.loc[sub_df["Condition"] == LOW_CONTRAST]
        high = sub_df.loc[sub_df["Condition"] == HIGH_CONTRAST]
        if low.empty or high.empty:
            continue
        row: dict[str, object] = {"ID": sid}
        valid = True
        for var in variables:
            lo = float(low[var].iloc[0])
            hi = float(high[var].iloc[0])
            if not np.isfinite(lo) or not np.isfinite(hi):
                valid = False
                break
            row[var] = hi - lo
        if valid:
            rows.append(row)
    return pd.DataFrame(rows)


def _symmetric_limits(values: np.ndarray, pad_frac: float = 0.08) -> tuple[float, float]:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return -1.0, 1.0
    lim = float(np.max(np.abs(finite)))
    if lim == 0:
        lim = 1.0
    pad = lim * pad_frac
    return -lim - pad, lim + pad


def _shade_hypothesis_quadrant(ax: plt.Axes, hypothesis: str, xlim: tuple[float, float], ylim: tuple[float, float]) -> None:
    xmin, xmax = xlim
    ymin, ymax = ylim

    if hypothesis == "x_pos_y_neg":
        patch_x = [0.0, xmax, xmax, 0.0]
        patch_y = [0.0, 0.0, ymin, ymin]
    elif hypothesis == "x_pos_y_pos":
        patch_x = [0.0, xmax, xmax, 0.0]
        patch_y = [0.0, 0.0, ymax, ymax]
    else:
        raise ValueError(f"Unknown hypothesis code: {hypothesis}")

    ax.add_patch(
        plt.Polygon(
            np.column_stack([patch_x, patch_y]),
            closed=True,
            facecolor=HYPOTHESIS_FACE_COLOR,
            edgecolor="none",
            alpha=HYPOTHESIS_ALPHA,
            zorder=0,
        )
    )


def _correlation_text(x: np.ndarray, y: np.ndarray) -> str:
    mask = np.isfinite(x) & np.isfinite(y)
    if np.count_nonzero(mask) < 3:
        return "N < 3"
    r_p, p_p = pearsonr(x[mask], y[mask])
    r_s, p_s = spearmanr(x[mask], y[mask])
    return (
        f"Pearson r = {r_p:.3f}, p = {p_p:.3g}\n"
        f"Spearman rho = {r_s:.3f}, p = {p_s:.3g}"
    )


def _hypothesis_count(x: np.ndarray, y: np.ndarray, hypothesis: str) -> int:
    mask = np.isfinite(x) & np.isfinite(y)
    if hypothesis == "x_pos_y_neg":
        return int(np.count_nonzero(mask & (x > 0) & (y < 0)))
    if hypothesis == "x_pos_y_pos":
        return int(np.count_nonzero(mask & (x > 0) & (y > 0)))
    return 0


def plot_gamma_gaze_scatter(delta: pd.DataFrame, output_dir: str) -> str:
    fig, axes = plt.subplots(3, 2, figsize=(15.12, 12.5), facecolor="white")
    axes = axes.ravel()

    for ax, panel in zip(axes, PANELS):
        x = delta[panel["x"]].to_numpy(dtype=float)
        y = delta[panel["y"]].to_numpy(dtype=float)

        xlim = _symmetric_limits(x)
        ylim = _symmetric_limits(y)
        _shade_hypothesis_quadrant(ax, panel["hypothesis"], xlim, ylim)

        ax.axhline(0.0, color="0.35", linewidth=0.8, zorder=1)
        ax.axvline(0.0, color="0.35", linewidth=0.8, zorder=1)
        ax.scatter(x, y, s=70, c="k", alpha=0.85, edgecolors="none", zorder=3)

        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_xlabel(panel["xlab"])
        ax.set_ylabel(panel["ylab"])
        ax.set_title(panel["title"])

    fig.suptitle(
        "Contrast Related Change in Gamma and Oculomotor Metrics (100% - 25%)",
        fontsize=16,
        y=0.98,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.96])

    out_path = os.path.join(output_dir, "GCP_stats_gamma_gaze_scatter.png")
    fig.savefig(
        out_path,
        dpi=FIGURE_SAVE_DPI,
        transparent=False,
        facecolor=fig.get_facecolor(),
        edgecolor="white",
    )
    plt.close(fig)
    return out_path


def main() -> None:
    base_dir = "/Volumes/g_psyplafor_methlab$/Students/Arne/GCP"
    features_dir = os.path.join(base_dir, "data", "features")
    controls_dir = os.path.join(base_dir, "data", "controls")
    output_dir = os.path.join(base_dir, "figures", "stats", "gamma_gaze_scatter")
    os.makedirs(output_dir, exist_ok=True)

    merged_csv = os.path.join(features_dir, "GCP_merged_data_trials.csv")
    ged_mat = os.path.join(features_dir, "GCP_eeg_GED.mat")

    if not os.path.isfile(merged_csv):
        raise FileNotFoundError(f"Merged trial CSV not found: {merged_csv}")
    if not os.path.isfile(ged_mat):
        raise FileNotFoundError(f"GED MAT file not found: {ged_mat}")

    variables = ["GammaFrequency", "GammaPower", "MSRate_bl", "BCEA_bl", "Vel2D_bl"]

    merged = filter_gcp_analysis_cohort(
        _rename_gaze_columns(load_merged_trial_metrics(merged_csv)),
        controls_dir,
    )
    ged = label_condition(load_ged_condition_metrics(ged_mat))
    # Keep GED subjects to the same analysis cohort selected from merged metrics.
    included_ids = set(merged["ID"].astype(str).unique())
    ged = ged.loc[ged["ID"].astype(str).isin(included_ids)].copy()

    required_gaze_cols = ["MSRate_bl", "BCEA_bl", "Vel2D_bl"]
    missing_gaze_cols = [col for col in required_gaze_cols if col not in merged.columns]
    if missing_gaze_cols:
        raise KeyError(f"Merged trial table missing required gaze columns: {missing_gaze_cols}")

    gaze_vars = ["MSRate_bl", "BCEA_bl", "Vel2D_bl"]
    ged_vars = ["GammaFrequency", "GammaPower"]

    gaze_trials = _trial_table_with_outliers(merged, gaze_vars)
    gaze_means = _subject_condition_means(gaze_trials, gaze_vars)

    # Condition-averaged spectral peaks are already one row per subject x condition.
    ged_means = ged.loc[:, ["ID", "Condition", *ged_vars]].copy()
    ged_means["ID"] = ged_means["ID"].astype(str)
    for var in ged_vars:
        dvar = ged_means.loc[ged_means[var].notna(), ["ID", "Condition", var]].copy()
        if dvar.empty:
            continue
        dvar = iqr_outlier_filter(dvar, [var], by="Condition")
        drop_idx = dvar.index[dvar[var].isna()]
        ged_means.loc[drop_idx, var] = np.nan

    gaze_delta = _contrast_delta(gaze_means, gaze_vars)
    ged_delta = _contrast_delta(ged_means, ged_vars)

    delta = pd.merge(gaze_delta, ged_delta, on="ID", how="inner")
    if delta.empty:
        raise RuntimeError("No subjects with complete 25% and 100% means for all variables.")

    print(f"Subjects with complete deltas: N = {len(delta)}")
    for panel in PANELS:
        x = delta[panel["x"]].to_numpy(dtype=float)
        y = delta[panel["y"]].to_numpy(dtype=float)
        n_hyp = _hypothesis_count(x, y, panel["hypothesis"])
        print(
            f"  {panel['title']}: hypothesis quadrant {n_hyp}/{len(delta)}; "
            f"{_correlation_text(x, y).replace(chr(10), '; ')}"
        )

    out_path = plot_gamma_gaze_scatter(delta, output_dir)
    print(f"\nFigure saved to: {out_path}")


if __name__ == "__main__":
    main()
