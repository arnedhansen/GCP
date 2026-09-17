# %% GCP Stats Rainclouds — confirmatory full-window subject x contrast
# One point per subject x condition from GCP_merged_data.csv (Include == 1).
# Asterisks: FDR-BH pairwise MixedLM (follow-up). Slope printed to console (primary).
#
# Prefer running GCP_stats_lmm.py first so brackets can load saved pairwise CSV;
# otherwise pairwise contrasts are fit on the fly via the same helpers.
#
# Run:
#   python GCP_stats_rainclouds.py

# %% Imports
import os
import sys
import warnings

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde
from statsmodels.tools.sm_exceptions import ConvergenceWarning

sys.path.insert(0, os.path.expanduser("/Users/Arne/Documents/GitHub"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from functions.rainclouds_plotting_helpers import add_stat_brackets

from GCP_stats_lmm import (
    CONDITION_ORDER,
    PAIRWISE_COMPARISONS,
    VARIABLES,
    bracket_labels_from_pairwise,
    fit_contrast_slope,
    fit_pairwise_fdr,
    load_merged_subject_data,
    print_slope_row,
)

warnings.filterwarnings("ignore", category=ConvergenceWarning, module="statsmodels")

# %% GCP colours (color_def.m)
PAL = ["#FFE680", "#E69966", "#E66666", "#000000"]

FIGURE_SAVE_DPI = 600
YLABEL_GRID_X = -0.15
YLIM_PAD_FRAC = 0.06
# Primary test is the continuous slope (console); figure shows FDR pairwise asterisks.
SHOW_SLOPE_ANNOTATION = False

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


def data_ylim(yvals: np.ndarray, pad_frac: float = YLIM_PAD_FRAC) -> tuple[float, float]:
    ymin = float(np.nanmin(yvals))
    ymax = float(np.nanmax(yvals))
    yr = ymax - ymin
    if not np.isfinite(yr) or yr == 0:
        yr = max(abs(ymax), abs(ymin), 1.0) * 0.1
    pad = yr * pad_frac
    return ymin - pad, ymax + pad


def format_slope_annotation(stats: dict) -> str:
    beta = stats["beta"]
    ci_lo = stats["CI_low"]
    ci_hi = stats["CI_high"]
    p = stats["p"]
    if not np.isfinite(p):
        p_txt = "p = NA"
    elif p < 0.001:
        p_txt = "p < .001"
    else:
        p_txt = f"p = {p:.3f}".replace("0.", ".")
    return f"β = {beta:.2f}, 95% CI [{ci_lo:.2f}, {ci_hi:.2f}], {p_txt}"


def load_pairwise_csv(stats_dir: str) -> pd.DataFrame | None:
    path = os.path.join(stats_dir, "GCP_pairwise_mixedlm.csv")
    if not os.path.isfile(path):
        return None
    return pd.read_csv(path)


def pairwise_for_var(
    dat: pd.DataFrame,
    var: str,
    pw_all: pd.DataFrame | None,
) -> pd.DataFrame | None:
    """Prefer LMM-script CSV; otherwise fit the same FDR pairwise MixedLM."""
    if pw_all is not None and not pw_all.empty:
        sub = pw_all.loc[pw_all["Variable"] == var].copy()
        if not sub.empty:
            return sub
    return fit_pairwise_fdr(dat, var)


def plot_raincloud(
    dvar: pd.DataFrame,
    var: str,
    ylab: str,
    sname: str,
    pal_dict: dict,
    output_dir: str,
    pw: pd.DataFrame | None,
    slope: dict | None,
) -> None:
    fig, ax = plt.subplots(figsize=(8, 6), facecolor="white")
    ax.set_facecolor("white")

    viol_alpha = 0.60
    n_per_cond = [
        int(dvar.loc[dvar["Condition"] == c, var].notna().sum())
        for c in CONDITION_ORDER
    ]
    few_points = bool(n_per_cond) and max(n_per_cond) <= 20
    dot_alpha = 0.85 if few_points else 0.18
    dot_size = 70 if few_points else 24
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

    range_y = max(ymax_plot - ymin_plot, np.finfo(float).eps)
    n_br = len(PAIRWISE_COMPARISONS)
    step = 0.09 * range_y
    y_positions = [ymax_plot + 0.06 * range_y + i * step for i in range(n_br)]
    labels = bracket_labels_from_pairwise(pw, var, PAIRWISE_COMPARISONS)
    add_stat_brackets(
        ax=ax,
        xcats=CONDITION_ORDER,
        comparisons=PAIRWISE_COMPARISONS,
        y_positions=y_positions,
        labels=labels,
        xmap=xpos,
        fontsize=11,
    )
    ymax_plot = y_positions[-1] + 0.08 * range_y if y_positions else ymax_plot

    if SHOW_SLOPE_ANNOTATION and slope is not None:
        ann = format_slope_annotation(slope)
        ymax_plot = ymax_plot + 0.06 * range_y
        ax.text(
            0.5,
            0.98,
            ann,
            transform=ax.transAxes,
            ha="center",
            va="top",
            fontsize=13,
            color="black",
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
    stats_dir = os.path.join(base_dir, "data", "stats")
    output_dir = os.path.join(base_dir, "figures", "stats", "rainclouds")
    os.makedirs(output_dir, exist_ok=True)

    merged_csv = os.path.join(features_dir, "GCP_merged_data.csv")
    if not os.path.isfile(merged_csv):
        raise FileNotFoundError(f"Merged subject CSV not found: {merged_csv}")

    dat = load_merged_subject_data(merged_csv)
    print(
        f"Loaded subject-level merged data: "
        f"subjects={dat['ID'].nunique()}, rows={len(dat)}"
    )

    pw_all = load_pairwise_csv(stats_dir)
    if pw_all is not None:
        print(f"Loaded pairwise CSV -> {os.path.basename(os.path.join(stats_dir, 'GCP_pairwise_mixedlm.csv'))}")
    else:
        print("No saved pairwise CSV; fitting FDR pairwise MixedLM per DV.")

    pal_dict = dict(zip(CONDITION_ORDER, PAL))

    for var, ylab, sname in VARIABLES:
        if var not in dat.columns:
            print(f"WARNING: {var} missing; skipping.")
            continue

        dvar = dat.loc[dat[var].notna(), ["ID", "Condition", var]].copy()
        if dvar.empty:
            print(f"WARNING: no data for {var}; skipping.")
            continue

        dvar["Condition"] = pd.Categorical(
            dvar["Condition"], categories=CONDITION_ORDER, ordered=True
        )
        dvar = dvar.loc[np.isfinite(dvar[var])].copy()
        if dvar.empty:
            print(f"WARNING: no finite data for {var}; skipping.")
            continue

        print(f"\n{var}: n={len(dvar)}, subjects={dvar['ID'].nunique()}")

        slope = fit_contrast_slope(dat, var)
        if slope is not None:
            print_slope_row(slope)
        else:
            print("  slope: unavailable")

        pw = pairwise_for_var(dat, var, pw_all)
        if pw is not None:
            for _, r in pw.iterrows():
                g1 = r.get("Group1", r.get("group1"))
                g2 = r.get("Group2", r.get("group2"))
                padj = float(r["p_adj"])
                print(f"  {g1} vs {g2}: p_adj={padj:.4g}")
        else:
            print("  pairwise: unavailable")

        plot_raincloud(dvar, var, ylab, sname, pal_dict, output_dir, pw, slope)

    print("\nRaincloud figures saved to:", output_dir)


if __name__ == "__main__":
    main()
