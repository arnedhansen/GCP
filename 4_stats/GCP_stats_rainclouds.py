# %% GCP Stats Rainclouds — GED condition peaks, gaze/behavior trials
# Gamma frequency/power: one point per subject x condition from peaks of
# condition-averaged GED spectra (all_condition_peak_*_full).
# Gaze/behavior: trial-level rainclouds from the merged trial table.
# Raincloud figures with half-kernel densities, boxplots, jittered points.
# Annotation: continuous contrast slope (beta, 95% CI, p), not pairwise brackets.
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
import scipy.io
import statsmodels.formula.api as smf
from scipy.stats import gaussian_kde, norm
from statsmodels.tools.sm_exceptions import ConvergenceWarning

sys.path.insert(0, os.path.expanduser("/Users/Arne/Documents/GitHub"))

from functions.stats_helpers import iqr_outlier_filter

warnings.filterwarnings("ignore", category=ConvergenceWarning, module="statsmodels")

# %% GCP colours (color_def.m)
PAL = ["#FFE680", "#E69966", "#E66666", "#000000"]
CONDITION_ORDER = ["25%", "50%", "75%", "100%"]
COND_NUM_TO_LABEL = {1: "25%", 2: "50%", 3: "75%", 4: "100%"}
COND_LABEL_TO_PCT = {"25%": 25.0, "50%": 50.0, "75%": 75.0, "100%": 100.0}
CONTRAST_LEVELS = np.array([25.0, 50.0, 75.0, 100.0], dtype=float)
CONTRAST_POP_SD = float(np.sqrt(np.mean((CONTRAST_LEVELS - CONTRAST_LEVELS.mean()) ** 2)))

VARIABLES = [
    ("GammaFrequency", "Peak Gamma Frequency [Hz]", "gamma_freq"),
    ("GammaPower", "Peak Gamma Power [dB]", "gamma_power"),
    ("MSRate_bl", "Microsaccade Rate [%]", "ms"),
    ("PupilSize_bl", "Pupil Size [%]", "pupil"),
    ("VelV_bl", "Gaze Velocity Y [%]", "vely"),
    ("ReactionTime", "Reaction Time [s]", "rt"),
]

FIGURE_SAVE_DPI = 600
YLABEL_GRID_X = -0.15
YLIM_PAD_FRAC = 0.06
SHOW_SLOPE_ANNOTATION = True

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


def load_ged_condition_metrics(mat_path: str) -> pd.DataFrame:
    """Subject x condition peaks from condition-averaged GED spectra."""
    ged = scipy.io.loadmat(mat_path, squeeze_me=True, struct_as_record=False)
    freq = np.asarray(ged["all_condition_peak_freq_full"], dtype=float)
    power = np.asarray(ged["all_condition_peak_power_full"], dtype=float)
    subjects = _as_str_list(ged["subjects"])

    if freq.ndim != 2 or power.ndim != 2:
        raise ValueError(
            "Expected all_condition_peak_freq/power_full as [condition x subject] matrices."
        )
    if freq.shape != power.shape:
        raise ValueError(
            f"Frequency/power matrix shape mismatch: {freq.shape} vs {power.shape}"
        )

    n_cond, n_subj = freq.shape
    rows: list[dict] = []
    for c in range(n_cond):
        cond_label = CONDITION_ORDER[c] if c < len(CONDITION_ORDER) else str(c + 1)
        for s in range(n_subj):
            sid = subjects[s] if s < len(subjects) else str(s + 1)
            rows.append(
                {
                    "ID": sid,
                    "Condition": cond_label,
                    "GammaFrequency": float(freq[c, s]) if np.isfinite(freq[c, s]) else np.nan,
                    "GammaPower": float(power[c, s]) if np.isfinite(power[c, s]) else np.nan,
                }
            )
    return pd.DataFrame(rows)


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
            "Gaze_MSRate_bl": "MSRate_bl",
            "Gaze_PupilSize_bl": "PupilSize_bl",
            "Gaze_VelV_bl": "VelV_bl",
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


def contrast_num_c_from_labels(condition: pd.Series) -> np.ndarray:
    """Centered/scaled contrast coding matching section 2.4 / power analysis."""
    pct = condition.map(COND_LABEL_TO_PCT).to_numpy(dtype=float)
    return (pct - float(CONTRAST_LEVELS.mean())) / CONTRAST_POP_SD


def subject_condition_means(dvar: pd.DataFrame, var: str) -> pd.DataFrame:
    """One row per subject x condition (confirmatory analysis unit)."""
    out = (
        dvar.groupby(["ID", "Condition"], observed=True, sort=False)[var]
        .mean()
        .reset_index()
        .rename(columns={var: "value"})
    )
    out = out.loc[np.isfinite(out["value"])].copy()
    out["contrast_num_c"] = contrast_num_c_from_labels(out["Condition"])
    return out


def fit_contrast_slope(dvar: pd.DataFrame, var: str) -> dict | None:
    """
    Fit value ~ contrast_num_c + (1 + contrast_num_c | Subject), with fallback
    to (1 | Subject) on non-convergence / singularity.
    """
    dat = subject_condition_means(dvar, var)
    if dat["ID"].nunique() < 2 or dat["contrast_num_c"].nunique() < 2:
        return None

    formula = "value ~ contrast_num_c"
    fit_kwargs = dict(reml=False, method="lbfgs")

    def _fit(re_formula: str):
        model = smf.mixedlm(
            formula,
            data=dat,
            groups=dat["ID"],
            re_formula=re_formula,
        )
        return model.fit(**fit_kwargs)

    res = None
    used_re = "(1 + contrast_num_c | Subject)"
    try:
        res = _fit("~contrast_num_c")
        if bool(getattr(res, "converged", True)) is False:
            raise RuntimeError("random-slope fit did not converge")
    except Exception:
        used_re = "(1 | Subject)"
        try:
            res = _fit("1")
        except Exception as exc:
            print(f"WARNING: slope model failed for {var}: {exc}")
            return None

    if "contrast_num_c" not in res.params.index:
        print(f"WARNING: contrast_num_c missing from slope fit for {var}")
        return None

    beta = float(res.params["contrast_num_c"])
    # Prefer model SE; fall back to normal approx from bse if needed.
    se = float(res.bse["contrast_num_c"]) if "contrast_num_c" in res.bse.index else np.nan
    if np.isfinite(se) and se > 0:
        z = beta / se
        p = float(2.0 * (1.0 - norm.cdf(abs(z))))
        ci_lo = beta - 1.96 * se
        ci_hi = beta + 1.96 * se
    else:
        z = np.nan
        p = np.nan
        ci_lo = np.nan
        ci_hi = np.nan

    return {
        "beta": beta,
        "se": se,
        "z": z,
        "p": p,
        "ci_lo": ci_lo,
        "ci_hi": ci_hi,
        "re": used_re,
        "n_subjects": int(dat["ID"].nunique()),
        "n_obs": int(len(dat)),
    }


def format_slope_annotation(stats: dict) -> str:
    beta = stats["beta"]
    ci_lo = stats["ci_lo"]
    ci_hi = stats["ci_hi"]
    p = stats["p"]
    if not np.isfinite(p):
        p_txt = "p = NA"
    elif p < 0.001:
        p_txt = "p < .001"
    else:
        p_txt = f"p = {p:.3f}".replace("0.", ".")

    return (
        f"β = {beta:.2f}, 95% CI [{ci_lo:.2f}, {ci_hi:.2f}], {p_txt}"
    )


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
    # Subject-level gamma rainclouds have few points; keep them readable.
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

    range_y = ymax_plot - ymin_plot
    if SHOW_SLOPE_ANNOTATION:
        slope = fit_contrast_slope(dvar, var)
        if slope is not None:
            ann = format_slope_annotation(slope)
            # Leave a little headroom so the annotation does not collide with clouds.
            ymax_plot = ymax_plot + 0.08 * range_y
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
            print(
                f"  slope ({slope['re']}): β={slope['beta']:.3f}, "
                f"p={slope['p']:.4g}, n_subj={slope['n_subjects']}"
            )
        else:
            print(f"  slope annotation unavailable for {var}")

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
        label_condition(load_ged_condition_metrics(ged_mat)),
        controls_dir,
    )

    data_by_var = {
        "GammaFrequency": ged,
        "GammaPower": ged,
        "MSRate_bl": merged,
        "PupilSize_bl": merged,
        "VelV_bl": merged,
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
