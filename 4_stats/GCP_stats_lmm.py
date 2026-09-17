# %% GCP confirmatory LMMs — full-window slope + FDR pairwise
# Primary: value ~ contrast_num_c + (1 + contrast_num_c | ID), fallback (1 | ID).
# Follow-up: categorical MixedLM + all six pairwise contrasts, FDR-BH within DV.
# Analysis unit: one observation per subject x contrast (GCP_merged_data.csv).
#
# Run:
#   python GCP_stats_lmm.py
#
# Outputs (AOC-compatible CSV columns for later table compilation):
#   .../data/stats/GCP_mixedlm_slope.csv
#   .../data/stats/GCP_pairwise_mixedlm.csv
#   .../data/stats/GCP_mixedlm_fixed_<sname>.csv

# %% Imports
import os
import sys
import warnings

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import norm
from statsmodels.tools.sm_exceptions import ConvergenceWarning

sys.path.insert(0, os.path.expanduser("/Users/Arne/Documents/GitHub"))
from functions.stats_helpers import mixedlm_pairwise_contrasts, p_to_signif

warnings.filterwarnings("ignore", category=ConvergenceWarning, module="statsmodels")

# %% Constants
CONDITION_ORDER = ["25%", "50%", "75%", "100%"]
COND_NUM_TO_LABEL = {1: "25%", 2: "50%", 3: "75%", 4: "100%"}
COND_LABEL_TO_PCT = {"25%": 25.0, "50%": 50.0, "75%": 75.0, "100%": 100.0}
CONTRAST_LEVELS = np.array([25.0, 50.0, 75.0, 100.0], dtype=float)
CONTRAST_POP_SD = float(np.sqrt(np.mean((CONTRAST_LEVELS - CONTRAST_LEVELS.mean()) ** 2)))

# Confirmatory full-window DVs (RR): pupil, MS rate, BCEA, Vel2D, GPP, GPF
VARIABLES = [
    ("PupilSize_bl", "Pupil Size [%]", "pupil"),
    ("MSRate_bl", "Microsaccade Rate [%]", "ms"),
    ("BCEA_bl", "BCEA [%]", "bcea"),
    ("Vel2D_bl", "Eye Velocity [%]", "vel2d"),
    ("Power", "Peak Gamma Power [dB]", "gpp"),
    ("Frequency", "Peak Gamma Frequency [Hz]", "gpf"),
]

PAIRWISE_COMPARISONS = [
    (CONDITION_ORDER[i], CONDITION_ORDER[j])
    for i in range(len(CONDITION_ORDER))
    for j in range(i + 1, len(CONDITION_ORDER))
]


def contrast_num_c_from_labels(condition: pd.Series) -> np.ndarray:
    """Centered/scaled contrast coding matching section 2.4 / power analysis."""
    pct = condition.map(COND_LABEL_TO_PCT).to_numpy(dtype=float)
    return (pct - float(CONTRAST_LEVELS.mean())) / CONTRAST_POP_SD


def load_merged_subject_data(csv_path: str) -> pd.DataFrame:
    """Subject-level master matrix; filter Include == 1 when present."""
    dat = pd.read_csv(csv_path)
    if "Include" in dat.columns:
        dat = dat.loc[dat["Include"].astype(bool)].copy()

    dat["ID"] = dat["ID"].astype(str)
    if np.issubdtype(dat["Condition"].dtype, np.number):
        dat["Condition"] = dat["Condition"].map(COND_NUM_TO_LABEL)
    else:
        dat["Condition"] = dat["Condition"].astype(str).str.strip()

    dat["Condition"] = pd.Categorical(
        dat["Condition"], categories=CONDITION_ORDER, ordered=True
    )
    dat["contrast_num_c"] = contrast_num_c_from_labels(dat["Condition"])
    return dat


def subject_condition_frame(dat: pd.DataFrame, var: str) -> pd.DataFrame:
    """One finite observation per subject x condition for `var`."""
    out = dat.loc[dat[var].notna(), ["ID", "Condition", "contrast_num_c", var]].copy()
    out = out.rename(columns={var: "value"})
    out = out.loc[np.isfinite(out["value"])].copy()
    out["Condition"] = pd.Categorical(
        out["Condition"], categories=CONDITION_ORDER, ordered=True
    )
    return out


def fit_contrast_slope(dat: pd.DataFrame, var: str) -> dict | None:
    """
    Primary MixedLM: value ~ contrast_num_c.
    Preferred RE: random intercept + random slope for contrast; fall back to
    random intercept only on non-convergence / singularity.
    """
    d = subject_condition_frame(dat, var)
    if d["ID"].nunique() < 2 or d["contrast_num_c"].nunique() < 2:
        return None

    formula = "value ~ contrast_num_c"
    fit_kwargs = dict(reml=False, method="lbfgs")

    def _fit(re_formula: str):
        model = smf.mixedlm(
            formula,
            data=d,
            groups=d["ID"],
            re_formula=re_formula,
        )
        return model.fit(**fit_kwargs)

    res = None
    used_re = "(1 + contrast_num_c | ID)"
    model_label = f"{var} ~ contrast_num_c + (1 + contrast_num_c | ID)"
    try:
        res = _fit("~contrast_num_c")
        if bool(getattr(res, "converged", True)) is False:
            raise RuntimeError("random-slope fit did not converge")
    except Exception:
        used_re = "(1 | ID)"
        model_label = f"{var} ~ contrast_num_c + (1 | ID)"
        try:
            res = _fit("1")
        except Exception as exc:
            print(f"WARNING: slope model failed for {var}: {exc}")
            return None

    if "contrast_num_c" not in res.params.index:
        print(f"WARNING: contrast_num_c missing from slope fit for {var}")
        return None

    beta = float(res.params["contrast_num_c"])
    se = float(res.bse["contrast_num_c"]) if "contrast_num_c" in res.bse.index else np.nan
    if np.isfinite(se) and se > 0:
        z = beta / se
        p = float(2.0 * (1.0 - norm.cdf(abs(z))))
        ci_lo = beta - 1.96 * se
        ci_hi = beta + 1.96 * se
    else:
        z = p = ci_lo = ci_hi = np.nan

    return {
        "Variable": var,
        "ModelLabel": model_label,
        "Term": "contrast_num_c",
        "beta": beta,
        "SE": se,
        "stat": z,
        "p": p,
        "CI_low": ci_lo,
        "CI_high": ci_hi,
        "N_obs": int(len(d)),
        "N_subjects": int(d["ID"].nunique()),
        "RE": used_re,
        "result": res,
        "data": d,
    }


def slope_to_fixed_df(slope: dict) -> pd.DataFrame:
    """AOC-like fixed-effects row for the primary contrast slope."""
    return pd.DataFrame(
        [
            {
                "DV": slope["Variable"],
                "ModelLabel": slope["ModelLabel"],
                "Term": slope["Term"],
                "beta": slope["beta"],
                "SE": slope["SE"],
                "stat": slope["stat"],
                "p": slope["p"],
                "CI_low": slope["CI_low"],
                "CI_high": slope["CI_high"],
                "N_obs": slope["N_obs"],
                "N_subjects": slope["N_subjects"],
                "RE": slope["RE"],
            }
        ]
    )


def fit_pairwise_fdr(dat: pd.DataFrame, var: str) -> pd.DataFrame | None:
    """Categorical MixedLM + all pairwise contrasts, FDR-BH within DV."""
    d = subject_condition_frame(dat, var)
    if d["ID"].nunique() < 2 or d["Condition"].nunique() < 2:
        return None

    try:
        pw = mixedlm_pairwise_contrasts(
            d,
            value_col="value",
            group_col="Condition",
            id_col="ID",
            p_adjust="fdr_bh",
        )
    except Exception as exc:
        print(f"WARNING: pairwise MixedLM failed for {var}: {exc}")
        return None

    model_label = f"{var} ~ C(Condition) + (1 | ID)"
    out = pd.DataFrame(
        {
            "Variable": var,
            "ModelLabel": model_label,
            "N_obs": int(len(d)),
            "Group1": pw["group1"].astype(str),
            "Group2": pw["group2"].astype(str),
            "Estimate": pw["estimate"].astype(float),
            "SE": pw["se"].astype(float),
            "z": pw["z"].astype(float),
            "p": pw["p"].astype(float),
            "CI95_low": (pw["estimate"] - 1.96 * pw["se"]).astype(float),
            "CI95_high": (pw["estimate"] + 1.96 * pw["se"]).astype(float),
            "p_adj": pw["p_adj"].astype(float),
            "signif": [p_to_signif(float(p)) for p in pw["p_adj"]],
        }
    )
    return out


def bracket_labels_from_pairwise(
    pw: pd.DataFrame | None,
    var: str,
    comparisons: list[tuple[str, str]] | None = None,
) -> list[str]:
    """Map FDR-adjusted pairwise p-values to raincloud bracket labels."""
    if comparisons is None:
        comparisons = PAIRWISE_COMPARISONS
    labels: list[str] = []
    for g1, g2 in comparisons:
        if pw is None or pw.empty:
            labels.append("n.s.")
            continue
        sub = pw.loc[pw["Variable"] == var] if "Variable" in pw.columns else pw
        row = sub.loc[
            ((sub["Group1"] == g1) & (sub["Group2"] == g2))
            | ((sub["Group1"] == g2) & (sub["Group2"] == g1))
        ]
        if row.empty:
            labels.append("n.s.")
        else:
            labels.append(p_to_signif(float(row["p_adj"].iloc[0])))
    return labels


def print_slope_row(slope: dict) -> None:
    p = slope["p"]
    p_txt = "NA" if not np.isfinite(p) else f"{p:.4g}"
    print(
        f"  slope {slope['RE']}: β={slope['beta']:.4g}, SE={slope['SE']:.4g}, "
        f"95% CI [{slope['CI_low']:.4g}, {slope['CI_high']:.4g}], "
        f"p={p_txt}, n_subj={slope['N_subjects']}, n_obs={slope['N_obs']}"
    )


def print_pairwise_rows(pw: pd.DataFrame) -> None:
    for _, r in pw.iterrows():
        print(
            f"  {r['Group1']} vs {r['Group2']}: "
            f"est={r['Estimate']:.4g}, p={r['p']:.4g}, "
            f"p_adj={r['p_adj']:.4g} ({r['signif']})"
        )


def main() -> None:
    base_dir = "/Volumes/g_psyplafor_methlab$/Students/Arne/GCP"
    features_dir = os.path.join(base_dir, "data", "features")
    stats_dir = os.path.join(base_dir, "data", "stats")
    os.makedirs(stats_dir, exist_ok=True)

    merged_csv = os.path.join(features_dir, "GCP_merged_data.csv")
    if not os.path.isfile(merged_csv):
        raise FileNotFoundError(f"Merged subject CSV not found: {merged_csv}")

    dat = load_merged_subject_data(merged_csv)
    print(
        f"Loaded {merged_csv}\n"
        f"  subjects={dat['ID'].nunique()}, rows={len(dat)}"
    )

    slope_rows: list[pd.DataFrame] = []
    pairwise_rows: list[pd.DataFrame] = []

    for var, ylab, sname in VARIABLES:
        if var not in dat.columns:
            print(f"WARNING: {var} missing; skipping.")
            continue

        print(f"\n{var} ({ylab})")
        slope = fit_contrast_slope(dat, var)
        if slope is None:
            print("  slope: unavailable")
        else:
            print_slope_row(slope)
            fixed_df = slope_to_fixed_df(slope)
            slope_rows.append(fixed_df)
            fixed_path = os.path.join(stats_dir, f"GCP_mixedlm_fixed_{sname}.csv")
            fixed_df.to_csv(fixed_path, index=False)
            print(f"  saved {os.path.basename(fixed_path)}")

        pw = fit_pairwise_fdr(dat, var)
        if pw is None:
            print("  pairwise: unavailable")
        else:
            print_pairwise_rows(pw)
            pairwise_rows.append(pw)

    if slope_rows:
        slope_path = os.path.join(stats_dir, "GCP_mixedlm_slope.csv")
        pd.concat(slope_rows, ignore_index=True).to_csv(slope_path, index=False)
        print(f"\nSaved pooled slope table -> {slope_path}")

    if pairwise_rows:
        pw_path = os.path.join(stats_dir, "GCP_pairwise_mixedlm.csv")
        pd.concat(pairwise_rows, ignore_index=True).to_csv(pw_path, index=False)
        print(f"Saved pooled pairwise table -> {pw_path}")

    print("\nDone.")


if __name__ == "__main__":
    main()
