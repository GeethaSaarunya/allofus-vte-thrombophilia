#!/usr/bin/env python3
"""
All of Us VTE / JTH
Public reproducibility script 01: Five-factor primary analysis
Version: public_v1.0.0

PURPOSE
-------
Regenerate manuscript-facing analyses after revising the PRIMARY clinical-risk
burden from six factors to five factors:

    cancer history
    major surgery
    inpatient hospitalization
    infection/sepsis
    fracture/major trauma

Pregnancy/postpartum is no longer included in the accumulated primary burden
because the reviewer correctly identified that an all-history pre-enrollment
pregnancy/postpartum record is not a coherent time-independent risk component.
Script 16 showed that:
    - removing pregnancy/postpartum preserved/slightly strengthened the burden
      gradient;
    - restricting pregnancy/postpartum to 365 days pre-consent removed the
      apparent inverse association;
    - the F5/F2-by-burden layering result remained present.

SCRIPT 17 therefore treats the five-factor burden as the revised PRIMARY
clinical construct and rebuilds all manuscript-facing analyses that can be
generated from the frozen participant master.

WHAT THIS SCRIPT REBUILDS
-------------------------
A. Primary cohort / burden locks
B. Main Table 1 data
C. Main Table 2 clinical models
D. Main Table 3 seven targeted variant models
E. Main Figure 3 data
F. Formal F5/F2 x five-factor burden interaction
G. Finer burden / exact-count / joint-seven-marker / no-prior-VTE /
   anticoagulant sensitivities for revised Supplementary Table S6
H. Five-fold out-of-fold discrimination/calibration for revised Supplementary S1
I. Day-7 person-time sensitivity using the five-factor burden
J. Interval-specific Cox estimates (days 8-90, 91-365, >365)
K. Ancestry-stratified F5/F2 estimates adjusted for exact five-factor burden
L. Pregnancy/postpartum reviewer-response table:
     all-history pregnancy OR
     recent-365d pregnancy OR
     old six-factor vs revised five-factor burden
M. Paste-ready manuscript/reviewer-response summary

BOUNDARIES
----------
- No new phenotype discovery.
- No PRS.
- No weighted clinical score.
- No dosage/zygosity.
- No hospital random effect.
- No SCD/mechanical-prophylaxis model.
- Anticoagulant exposure is treatment context, not proven prophylaxis.
- The primary VTE endpoint is "observed VTE occurrence during follow-up", not
  first-ever incident VTE.
- Pregnancy/postpartum is evaluated separately and is not called protective.
- The time-varying Cox analysis is descriptive of exposure definitions over
  follow-up. It must NOT be framed as genetics being biologically "more stable"
  than transient clinical exposures.

INPUT
-----
Secure input:
    VTE_JTH_MASTER_FILE must point to the participant-level analysis master
    inside an authorized All of Us Researcher Workbench environment.

Active CDR:
    WORKSPACE_CDR must identify the authorized CDR dataset used for the
    recent pregnancy/postpartum sensitivity query.

OUTPUT
------
Aggregate-safe outputs only. By default they are written under:
    ./outputs/primary_analysis/

No participant-level file is written by this script.
"""

from __future__ import annotations

from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import math
import os
import shutil
import zipfile

import numpy as np
import pandas as pd
import statsmodels.api as sm

from scipy.special import expit, logit
from scipy.stats import chi2, norm
from statsmodels.duration.hazard_regression import PHReg
from statsmodels.stats.multitest import multipletests

from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    average_precision_score,
    brier_score_loss,
    log_loss,
    roc_auc_score,
)
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


# ============================================================================
# 0. Controls
# ============================================================================

SCRIPT_VERSION = "public_v1.0.0"

MASTER_FILE_ENV = os.environ.get("VTE_JTH_MASTER_FILE", "").strip()
if not MASTER_FILE_ENV:
    raise RuntimeError(
        "Set VTE_JTH_MASTER_FILE to the secure participant-level master "
        "inside the authorized All of Us Researcher Workbench."
    )

MASTER_FILE = Path(MASTER_FILE_ENV).expanduser()

OUTPUT_ROOT = Path(
    os.environ.get(
        "VTE_JTH_OUTPUT_DIR",
        "./outputs/primary_analysis",
    )
).expanduser().resolve()

TABLE_DIR = OUTPUT_ROOT / "tables"
LOG_DIR = OUTPUT_ROOT / "logs"
DOC_DIR = OUTPUT_ROOT / "documentation"
PACKAGE_DIR = OUTPUT_ROOT / "package_build"

for directory in [TABLE_DIR, LOG_DIR, DOC_DIR, PACKAGE_DIR]:
    directory.mkdir(parents=True, exist_ok=True)

CDR_DATASET = os.environ.get("WORKSPACE_CDR", "").strip().strip("`")
if not CDR_DATASET:
    raise RuntimeError(
        "Set WORKSPACE_CDR to the authorized All of Us CDR dataset name."
    )

PERSON_ID = "person_id"
WGS_ELIGIBLE = "analysis_eligible_primary_genetics"
PRIMARY_OUTCOME = "vte_primary_scientific"
PRIOR_HISTORY = "preindex_vte_history"
NO_PRIOR_ENDPOINT = "vte_primary_scientific_without_preindex_history"
ANTICOAG = "preindex_anticoagulant_exposure"
OBSERVABILITY = "ehr_observability_score_0_6"

EVENT_TIME = "vte_primary_time_days"
ALT_EVENT_TIME = "days_from_index_to_vte_internal"
FOLLOWUP_TIME = "available_postindex_followup_days"

ANCESTRY = "genetic_ancestry_label"

FIVE_FACTORS = [
    "cancer_history",
    "major_surgery",
    "inpatient_hospitalization",
    "infection_or_sepsis",
    "fracture_or_major_trauma",
]

PREGNANCY_ALL_HISTORY = "pregnancy_or_postpartum"
PREGNANCY_RECENT = "pregnancy_or_postpartum_recent_365d"

F5 = "f5_rs6025_carrier"
F2 = "f2_rs1799963_carrier"
F5F2 = "f5f2_carrier"

SEVEN_MARKERS = [
    F5,
    F2,
    "abo_rs8176719_alt_carrier",
    "f11_rs2036914_alt_carrier",
    "f11_rs2289252_alt_carrier",
    "fgg_rs2066865_alt_carrier",
    "procr_rs867186_alt_carrier",
]

MARKER_LABELS = {
    F5: "Factor V Leiden (F5 rs6025)",
    F2: "Prothrombin G20210A (F2 rs1799963)",
    "abo_rs8176719_alt_carrier": "ABO rs8176719 ALT allele",
    "f11_rs2036914_alt_carrier": "F11 rs2036914 ALT allele",
    "f11_rs2289252_alt_carrier": "F11 rs2289252 ALT allele",
    "fgg_rs2066865_alt_carrier": "FGG rs2066865 ALT allele",
    "procr_rs867186_alt_carrier": "PROCR rs867186 ALT allele",
}

PC_COLUMNS = [f"genetic_pc{i}" for i in range(1, 17)]

EXPECTED_MASTER_N = 483_707
EXPECTED_WGS_N = 358_533
EXPECTED_VTE_N = 7_553

# Script 16 frozen five-factor counts.
EXPECTED_5F_BURDEN = {
    "0": 220_936,
    "1": 68_503,
    "2+": 69_094,
}

EXPECTED_5F_FIGURE3 = {
    ("0", 0): (208_024, 1_475),
    ("0", 1): (12_912, 215),
    ("1", 0): (64_496, 1_701),
    ("1", 1): (4_007, 204),
    ("2+", 0): (64_830, 3_589),
    ("2+", 1): (4_264, 369),
}

EXPECTED_5F_OR_1 = 2.008
EXPECTED_5F_OR_2PLUS = 3.739
EXPECTED_5F_OR_TOL = 0.03

# Script 16 reviewer-response locks.
EXPECTED_5F_MULT_P_LT = 0.001
EXPECTED_5F_ADD_P = 0.008
EXPECTED_PREG_ALL_HISTORY_N = 23_568
EXPECTED_PREG_RECENT_N = 7_421
EXPECTED_PREG_RECENT_OR = 1.030
EXPECTED_PREG_RECENT_OR_TOL = 0.08

# Day-7 landmark locks from completed Script 13.
LANDMARK_DAY = 7.0
DAYS_PER_YEAR = 365.25
EXPECTED_EARLY_VTE_N = 804
EXPECTED_SHORT_FOLLOWUP_NONCASE_N = 10_645
EXPECTED_LANDMARK_N = 347_084
EXPECTED_LANDMARK_VTE_N = 6_749

INTERVALS = [
    {"slug": "d8_90", "label": "Days 8-90", "start": 0.0, "end": 83.0},
    {"slug": "d91_365", "label": "Days 91-365", "start": 83.0, "end": 358.0},
    {"slug": "gt365", "label": ">365 days", "start": 358.0, "end": np.inf},
]

REPORTABLE_ANCESTRIES = ["EUR", "AFR", "AMR"]
MIN_ANCESTRY_N = 5_000
MIN_2X2_CELL = 20

MIN_AGGREGATE_CELL = 20
MAX_GLM_ITER = 300

CV_FOLDS = 5
CV_RANDOM_STATE = 20260711
BOOTSTRAP_REPS = int(os.environ.get("VTE_SCRIPT17_BOOTSTRAP_REPS", "500"))
BOOTSTRAP_SEED = 20260712

PREGNANCY_PATTERN = (
    r"pregnan|postpartum|post-partum|puerper|gestation|delivery"
)


# ============================================================================
# 1. Utilities
# ============================================================================

def banner():
    print("=" * 84)
    print("All of Us VTE / JTH")
    print("Script 17: Five-factor primary manuscript rebuild")
    print("Version:", SCRIPT_VERSION)
    print("=" * 84)


def progress(message: str):
    print(
        f"[Script 17 | {datetime.now().strftime('%H:%M:%S')}] {message}",
        flush=True,
    )


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_table(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix == ".parquet":
        return pd.read_parquet(path)
    if suffix in {".tsv", ".txt"}:
        return pd.read_csv(path, sep="\t")
    return pd.read_csv(path)


def require_columns(frame: pd.DataFrame, columns: list[str]):
    missing = [column for column in columns if column not in frame.columns]
    if missing:
        raise KeyError(
            "Missing required columns:\n  - " + "\n  - ".join(missing)
        )


def numeric_binary(series: pd.Series, name: str) -> pd.Series:
    value = pd.to_numeric(series, errors="coerce")
    observed = set(value.dropna().astype(int).unique())
    if not observed.issubset({0, 1}):
        raise RuntimeError(
            f"{name} is not binary: {sorted(observed)}"
        )
    return value.astype("float64")


def safe_count(value: int):
    value = int(value)
    if 0 < value < MIN_AGGREGATE_CELL:
        return f"<{MIN_AGGREGATE_CELL}"
    return value


def p_text(p: float) -> str:
    if pd.isna(p):
        return "NA"
    if p < 0.001:
        return "<0.001"
    return f"{p:.3f}"


def normal_p_value(z: float) -> float:
    if not np.isfinite(z):
        return np.nan
    return float(2.0 * norm.sf(abs(z)))


def Wilson_ci(events: int, n: int, z: float = 1.959963984540054):
    if n <= 0:
        return np.nan, np.nan
    p = events / n
    denom = 1.0 + z*z/n
    center = (p + z*z/(2*n)) / denom
    half = (
        z
        * math.sqrt(
            p*(1-p)/n + z*z/(4*n*n)
        )
        / denom
    )
    return center - half, center + half


def prepare_model_frame(
    frame: pd.DataFrame,
    outcome: str,
    predictors: list[str],
) -> pd.DataFrame:
    require_columns(frame, [outcome, *predictors])

    model = (
        frame[[outcome, *predictors]]
        .apply(pd.to_numeric, errors="coerce")
        .dropna()
        .astype("float64")
        .copy()
    )

    if model.empty:
        raise RuntimeError("No complete rows for model.")

    y_values = set(model[outcome].astype(int).unique())
    if not y_values.issubset({0, 1}):
        raise RuntimeError(f"{outcome} is not binary.")

    constant = [
        column
        for column in predictors
        if model[column].nunique(dropna=True) <= 1
    ]
    if constant:
        raise RuntimeError(
            "Unexpected constant predictor(s): "
            + ", ".join(constant)
        )

    return model


def fit_glm(
    frame: pd.DataFrame,
    outcome: str,
    predictors: list[str],
):
    model = prepare_model_frame(frame, outcome, predictors)

    X = sm.add_constant(
        model[predictors],
        has_constant="add",
    ).astype("float64")
    y = model[outcome].astype("float64")

    result = sm.GLM(
        y,
        X,
        family=sm.families.Binomial(),
    ).fit(
        cov_type="HC0",
        maxiter=MAX_GLM_ITER,
        disp=0,
    )

    if not bool(getattr(result, "converged", True)):
        raise RuntimeError("GLM did not converge.")

    return result, model, X


def extract_terms(
    result,
    terms: list[str],
    model_name: str,
    analysis_role: str,
) -> pd.DataFrame:
    rows = []

    for term in terms:
        if term not in result.params.index:
            continue

        beta = float(result.params[term])
        se = float(result.bse[term])

        rows.append(
            {
                "analysis_role": analysis_role,
                "model": model_name,
                "term": term,
                "beta": beta,
                "standard_error": se,
                "odds_ratio": math.exp(beta),
                "ci_95_low": math.exp(beta - 1.96 * se),
                "ci_95_high": math.exp(beta + 1.96 * se),
                "p_value": float(result.pvalues[term]),
                "n": int(result.nobs),
            }
        )

    return pd.DataFrame(rows)


def make_burden_dummies(
    frame: pd.DataFrame,
    burden_column: str,
    levels: list[str],
    prefix: str,
):
    out = frame.copy()
    terms = []

    for level in levels[1:]:
        safe = (
            level.replace("+", "plus")
            .replace(">=", "ge")
            .replace(" ", "_")
        )
        term = f"{prefix}_{safe}"
        out[term] = (
            out[burden_column].astype("string").eq(level)
        ).astype(float)
        terms.append(term)

    return out, terms


def joint_wald(result, terms: list[str]) -> dict:
    beta = result.params.loc[terms].to_numpy(dtype=float)
    covariance = (
        result.cov_params()
        .loc[terms, terms]
        .to_numpy(dtype=float)
    )
    statistic = float(
        beta.T @ np.linalg.pinv(covariance) @ beta
    )
    return {
        "wald_chi2": statistic,
        "df": len(terms),
        "p_value": float(chi2.sf(statistic, len(terms))),
    }


# ============================================================================
# 2. Interaction / marginal standardization
# ============================================================================

def build_interaction_model(
    frame: pd.DataFrame,
    outcome: str,
    exposure: str,
    burden_column: str,
    burden_levels: list[str],
    common_covariates: list[str],
    model_label: str,
):
    work, burden_terms = make_burden_dummies(
        frame,
        burden_column,
        burden_levels,
        f"{model_label}_burden",
    )

    interaction_terms = []

    for burden_term in burden_terms:
        term = f"{model_label}_{exposure}_x_{burden_term}"
        work[term] = work[exposure] * work[burden_term]
        interaction_terms.append(term)

    predictors = [
        exposure,
        *burden_terms,
        *interaction_terms,
        *common_covariates,
    ]

    result, model_frame, _ = fit_glm(
        work,
        outcome,
        predictors,
    )

    return {
        "result": result,
        "model_frame": model_frame,
        "predictors": predictors,
        "exposure": exposure,
        "burden_terms": burden_terms,
        "interaction_terms": interaction_terms,
        "burden_levels": burden_levels,
        "model_label": model_label,
    }


def standardized_risk_and_gradient(
    fit: dict,
    burden_level: str,
    exposure_value: int,
):
    result = fit["result"]
    model = fit["model_frame"].copy()
    predictors = fit["predictors"]
    exposure = fit["exposure"]
    burden_terms = fit["burden_terms"]
    interaction_terms = fit["interaction_terms"]
    burden_levels = fit["burden_levels"]

    model[exposure] = float(exposure_value)

    for level, burden_term in zip(
        burden_levels[1:],
        burden_terms,
    ):
        model[burden_term] = float(burden_level == level)

    for burden_term, interaction_term in zip(
        burden_terms,
        interaction_terms,
    ):
        model[interaction_term] = (
            model[exposure] * model[burden_term]
        )

    X = sm.add_constant(
        model[predictors],
        has_constant="add",
    )
    X = X[result.params.index]

    X_np = X.to_numpy(dtype=float)
    beta = result.params.to_numpy(dtype=float)

    probability = expit(X_np @ beta)
    risk = float(probability.mean())

    derivative_weight = probability * (1.0 - probability)
    gradient = np.mean(
        derivative_weight[:, None] * X_np,
        axis=0,
    )

    return risk, gradient


def additive_interaction_table(fit: dict):
    result = fit["result"]
    covariance = result.cov_params().to_numpy(dtype=float)
    levels = fit["burden_levels"]

    risks = {}
    gradients = {}
    standardized_rows = []

    for burden in levels:
        for exposure_value in [0, 1]:
            risk, gradient = standardized_risk_and_gradient(
                fit,
                burden,
                exposure_value,
            )

            risks[(burden, exposure_value)] = risk
            gradients[(burden, exposure_value)] = gradient

            standardized_rows.append(
                {
                    "burden_category": burden,
                    "exposure_value": exposure_value,
                    "adjusted_probability": risk,
                    "adjusted_percent": 100.0 * risk,
                }
            )

    risk_difference = {}
    risk_difference_gradient = {}
    contrast_rows = []

    for burden in levels:
        estimate = (
            risks[(burden, 1)]
            - risks[(burden, 0)]
        )
        gradient = (
            gradients[(burden, 1)]
            - gradients[(burden, 0)]
        )
        variance = float(
            gradient.T @ covariance @ gradient
        )
        se = math.sqrt(max(variance, 0.0))
        z = estimate / se if se > 0 else np.nan

        risk_difference[burden] = estimate
        risk_difference_gradient[burden] = gradient

        contrast_rows.append(
            {
                "contrast": "carrier_minus_noncarrier_risk_difference",
                "burden_category": burden,
                "estimate": estimate,
                "percentage_points": 100.0 * estimate,
                "standard_error": se,
                "p_value": normal_p_value(z),
            }
        )

    reference = levels[0]
    did = []
    did_gradient = []

    for burden in levels[1:]:
        estimate = (
            risk_difference[burden]
            - risk_difference[reference]
        )
        gradient = (
            risk_difference_gradient[burden]
            - risk_difference_gradient[reference]
        )
        variance = float(
            gradient.T @ covariance @ gradient
        )
        se = math.sqrt(max(variance, 0.0))
        z = estimate / se if se > 0 else np.nan

        did.append(estimate)
        did_gradient.append(gradient)

        contrast_rows.append(
            {
                "contrast": "difference_in_risk_differences_vs_reference",
                "burden_category": burden,
                "estimate": estimate,
                "percentage_points": 100.0 * estimate,
                "standard_error": se,
                "p_value": normal_p_value(z),
            }
        )

    G = np.vstack(did_gradient)
    V = G @ covariance @ G.T
    d = np.asarray(did, dtype=float)

    statistic = float(
        d.T @ np.linalg.pinv(V) @ d
    )

    global_test = {
        "wald_chi2": statistic,
        "df": len(d),
        "p_value": float(
            chi2.sf(statistic, len(d))
        ),
    }

    return (
        pd.DataFrame(standardized_rows),
        pd.DataFrame(contrast_rows),
        global_test,
    )


# ============================================================================
# 3. Descriptive / manuscript tables
# ============================================================================

def burden_by_exposure_descriptive(
    frame: pd.DataFrame,
    burden_column: str,
    exposure: str,
) -> pd.DataFrame:
    rows = []

    for burden in ["0", "1", "2+"]:
        for carrier in [0, 1]:
            sub = frame.loc[
                frame[burden_column].astype("string").eq(burden)
                & frame[exposure].eq(carrier)
            ]

            n = len(sub)
            events = int(sub[PRIMARY_OUTCOME].sum())
            low, high = Wilson_ci(events, n)

            rows.append(
                {
                    "burden_category": burden,
                    "carrier": carrier,
                    "n": safe_count(n),
                    "vte_n": safe_count(events),
                    "vte_percent": (
                        100.0 * events / n if n else np.nan
                    ),
                    "wilson_95_low_percent": (
                        100.0 * low if n else np.nan
                    ),
                    "wilson_95_high_percent": (
                        100.0 * high if n else np.nan
                    ),
                }
            )

    return pd.DataFrame(rows)


def figure3_chisquare(
    frame: pd.DataFrame,
) -> pd.DataFrame:
    rows = []

    for burden in ["0", "1", "2+"]:
        sub = frame.loc[
            frame["burden_012plus_5f"].astype("string").eq(burden)
        ].copy()

        table = pd.crosstab(
            sub[F5F2].astype(int),
            sub[PRIMARY_OUTCOME].astype(int),
        ).reindex(
            index=[0, 1],
            columns=[0, 1],
            fill_value=0,
        )

        # Pearson chi-square without Yates correction, matching manuscript.
        observed = table.to_numpy(dtype=float)
        row_sum = observed.sum(axis=1, keepdims=True)
        col_sum = observed.sum(axis=0, keepdims=True)
        total = observed.sum()
        expected = row_sum @ col_sum / total
        stat = float(
            ((observed - expected) ** 2 / expected).sum()
        )
        p = float(chi2.sf(stat, 1))

        rows.append(
            {
                "burden_category": burden,
                "pearson_chi2": stat,
                "df": 1,
                "p_value": p,
            }
        )

    return pd.DataFrame(rows)


def build_table1(analysis: pd.DataFrame) -> pd.DataFrame:
    rows = []

    def add_row(section, characteristic, overall, no_vte, vte, p=np.nan):
        rows.append(
            {
                "section": section,
                "characteristic": characteristic,
                "overall": overall,
                "no_vte": no_vte,
                "vte": vte,
                "p_value": p,
            }
        )

    no_case = analysis.loc[analysis[PRIMARY_OUTCOME].eq(0)]
    case = analysis.loc[analysis[PRIMARY_OUTCOME].eq(1)]

    add_row(
        "Study population",
        "Participants, N",
        len(analysis),
        len(no_case),
        len(case),
        np.nan,
    )

    # Age.
    age_all = pd.to_numeric(
        analysis["age_at_primary_consent"],
        errors="coerce",
    )
    age_no = pd.to_numeric(
        no_case["age_at_primary_consent"],
        errors="coerce",
    )
    age_vte = pd.to_numeric(
        case["age_at_primary_consent"],
        errors="coerce",
    )

    try:
        from scipy.stats import ttest_ind
        age_p = float(
            ttest_ind(
                age_no.dropna(),
                age_vte.dropna(),
                equal_var=False,
            ).pvalue
        )
    except Exception:
        age_p = np.nan

    add_row(
        "Demographics",
        "Age at All of Us enrollment, years, mean (SD)",
        f"{age_all.mean():.2f} ({age_all.std(ddof=1):.2f})",
        f"{age_no.mean():.2f} ({age_no.std(ddof=1):.2f})",
        f"{age_vte.mean():.2f} ({age_vte.std(ddof=1):.2f})",
        age_p,
    )

    # Recorded gender categories from manuscript encoding.
    gender_specs = [
        ("Recorded gender: Female", "female_sex"),
        ("Recorded gender: Other/unknown", "other_unknown_sex"),
    ]

    # Male = neither female nor other/unknown.
    analysis["_male_temp"] = (
        (analysis["female_sex"].eq(0))
        & (analysis["other_unknown_sex"].eq(0))
    ).astype(int)
    no_case = analysis.loc[analysis[PRIMARY_OUTCOME].eq(0)]
    case = analysis.loc[analysis[PRIMARY_OUTCOME].eq(1)]
    gender_specs.insert(
        1,
        ("Recorded gender: Male", "_male_temp"),
    )

    def n_pct(frame, column):
        n = int(frame[column].sum())
        return f"{n:,} ({100.0*n/len(frame):.2f}%)"

    # Overall chi-square for recorded gender.
    gender_code = np.select(
        [
            analysis["female_sex"].eq(1),
            analysis["_male_temp"].eq(1),
        ],
        ["Female", "Male"],
        default="Other/unknown",
    )
    gender_tab = pd.crosstab(
        gender_code,
        analysis[PRIMARY_OUTCOME].astype(int),
    )
    expected = (
        gender_tab.sum(axis=1).to_numpy()[:, None]
        @ gender_tab.sum(axis=0).to_numpy()[None, :]
        / gender_tab.to_numpy().sum()
    )
    gender_chi = float(
        (
            (gender_tab.to_numpy() - expected) ** 2
            / expected
        ).sum()
    )
    gender_p = float(
        chi2.sf(
            gender_chi,
            (gender_tab.shape[0]-1)*(gender_tab.shape[1]-1),
        )
    )

    for i, (label, column) in enumerate(gender_specs):
        add_row(
            "Demographics",
            label,
            n_pct(analysis, column),
            n_pct(no_case, column),
            n_pct(case, column),
            gender_p if i == 0 else np.nan,
        )

    # Five primary burden components.
    factor_labels = {
        "cancer_history": "Cancer history",
        "major_surgery": "Major surgery",
        "inpatient_hospitalization": "Inpatient hospitalization",
        "infection_or_sepsis": "Infection/sepsis",
        "fracture_or_major_trauma": "Fracture/major trauma",
    }

    for column in FIVE_FACTORS:
        tab = pd.crosstab(
            analysis[column].astype(int),
            analysis[PRIMARY_OUTCOME].astype(int),
        ).reindex(
            index=[0,1],
            columns=[0,1],
            fill_value=0,
        )
        obs = tab.to_numpy(dtype=float)
        expected = (
            obs.sum(axis=1, keepdims=True)
            @ obs.sum(axis=0, keepdims=True)
            / obs.sum()
        )
        stat = float(((obs-expected)**2/expected).sum())
        p = float(chi2.sf(stat, 1))

        add_row(
            "Captured clinical risk factors included in five-factor burden",
            factor_labels[column],
            n_pct(analysis, column),
            n_pct(no_case, column),
            n_pct(case, column),
            p,
        )

    # Pregnancy/postpartum remains descriptive but is explicitly not in burden.
    tab = pd.crosstab(
        analysis[PREGNANCY_ALL_HISTORY].astype(int),
        analysis[PRIMARY_OUTCOME].astype(int),
    ).reindex(
        index=[0,1],
        columns=[0,1],
        fill_value=0,
    )
    obs = tab.to_numpy(dtype=float)
    expected = (
        obs.sum(axis=1, keepdims=True)
        @ obs.sum(axis=0, keepdims=True)
        / obs.sum()
    )
    stat = float(((obs-expected)**2/expected).sum())
    preg_p = float(chi2.sf(stat, 1))

    add_row(
        "Additional clinical characteristic not included in burden",
        "Pregnancy/postpartum history",
        n_pct(analysis, PREGNANCY_ALL_HISTORY),
        n_pct(no_case, PREGNANCY_ALL_HISTORY),
        n_pct(case, PREGNANCY_ALL_HISTORY),
        preg_p,
    )

    # Burden category.
    burden_tab = pd.crosstab(
        analysis["burden_012plus_5f"].astype("string"),
        analysis[PRIMARY_OUTCOME].astype(int),
    ).reindex(
        index=["0","1","2+"],
        columns=[0,1],
        fill_value=0,
    )

    obs = burden_tab.to_numpy(dtype=float)
    expected = (
        obs.sum(axis=1, keepdims=True)
        @ obs.sum(axis=0, keepdims=True)
        / obs.sum()
    )
    stat = float(((obs-expected)**2/expected).sum())
    burden_p = float(chi2.sf(stat, 2))

    for i, burden in enumerate(["0","1","2+"]):
        sub = analysis.loc[
            analysis["burden_012plus_5f"].astype("string").eq(burden)
        ]
        sub_no = sub.loc[sub[PRIMARY_OUTCOME].eq(0)]
        sub_vte = sub.loc[sub[PRIMARY_OUTCOME].eq(1)]

        def fmt(n, denom):
            return f"{n:,} ({100.0*n/denom:.2f}%)"

        label = {
            "0": "0 captured factors",
            "1": "1 captured factor",
            "2+": "2 or more captured factors",
        }[burden]

        add_row(
            "Five-factor clinical-risk category",
            label,
            fmt(len(sub), len(analysis)),
            fmt(len(sub_no), len(no_case)),
            fmt(len(sub_vte), len(case)),
            burden_p if i == 0 else np.nan,
        )

    analysis.drop(columns=["_male_temp"], inplace=True)

    return pd.DataFrame(rows)


# ============================================================================
# 4. Cross-validation and calibration
# ============================================================================

def calibration_metrics(y: np.ndarray, p: np.ndarray) -> dict:
    eps = 1e-8
    p = np.clip(p, eps, 1.0 - eps)
    lp = logit(p)

    intercept_model = sm.GLM(
        y.astype(float),
        np.ones((len(y), 1), dtype=float),
        family=sm.families.Binomial(),
        offset=lp,
    ).fit(disp=0)

    X = sm.add_constant(lp, has_constant="add")
    slope_model = sm.GLM(
        y.astype(float),
        X,
        family=sm.families.Binomial(),
    ).fit(disp=0)

    return {
        "calibration_in_the_large_intercept": float(
            intercept_model.params[0]
        ),
        "calibration_intercept_joint_model": float(
            slope_model.params[0]
        ),
        "calibration_slope": float(
            slope_model.params[1]
        ),
    }


def out_of_fold_predictions(
    frame: pd.DataFrame,
    outcome: str,
    feature_map: dict[str, list[str]],
):
    y = frame[outcome].astype(int).to_numpy()

    splitter = StratifiedKFold(
        n_splits=CV_FOLDS,
        shuffle=True,
        random_state=CV_RANDOM_STATE,
    )

    predictions = {}
    rows = []
    dummy = np.zeros((len(frame), 1))
    splits = list(splitter.split(dummy, y))

    for model_name, features in feature_map.items():
        X = (
            frame[features]
            .apply(pd.to_numeric, errors="coerce")
            .to_numpy(dtype=float)
        )

        if np.isnan(X).any():
            raise RuntimeError(
                f"Prediction model {model_name} has missing predictors."
            )

        oof = np.full(len(frame), np.nan, dtype=float)

        for fold_number, (train_idx, test_idx) in enumerate(
            splits,
            start=1,
        ):
            pipeline = Pipeline(
                [
                    ("scale", StandardScaler()),
                    (
                        "logistic",
                        LogisticRegression(
                            solver="lbfgs",
                            penalty="l2",
                            C=1e6,
                            max_iter=2000,
                            random_state=CV_RANDOM_STATE,
                        ),
                    ),
                ]
            )

            pipeline.fit(X[train_idx], y[train_idx])
            prob = pipeline.predict_proba(X[test_idx])[:, 1]
            oof[test_idx] = prob

            rows.append(
                {
                    "model": model_name,
                    "evaluation": f"fold_{fold_number}",
                    "fold": fold_number,
                    "n": len(test_idx),
                    "events": int(y[test_idx].sum()),
                    "roc_auc": roc_auc_score(y[test_idx], prob),
                    "average_precision": average_precision_score(
                        y[test_idx], prob
                    ),
                    "brier_score": brier_score_loss(
                        y[test_idx], prob
                    ),
                    "log_loss": log_loss(
                        y[test_idx],
                        prob,
                        labels=[0, 1],
                    ),
                }
            )

        metrics = {
            "model": model_name,
            "evaluation": "pooled_oof",
            "fold": 0,
            "n": len(y),
            "events": int(y.sum()),
            "roc_auc": roc_auc_score(y, oof),
            "average_precision": average_precision_score(y, oof),
            "brier_score": brier_score_loss(y, oof),
            "log_loss": log_loss(y, oof, labels=[0, 1]),
        }
        metrics.update(calibration_metrics(y, oof))
        rows.append(metrics)

        predictions[model_name] = oof

    return y, predictions, pd.DataFrame(rows)


def paired_poisson_delta_auc(
    y: np.ndarray,
    p_clinical: np.ndarray,
    p_combined: np.ndarray,
    reps: int,
    seed: int,
) -> dict:
    rng = np.random.default_rng(seed)

    clinical_auc = roc_auc_score(y, p_clinical)
    combined_auc = roc_auc_score(y, p_combined)
    observed_delta = combined_auc - clinical_auc

    deltas = []

    for _ in range(reps):
        weights = rng.poisson(1.0, size=len(y)).astype(float)

        if weights[y == 1].sum() == 0 or weights[y == 0].sum() == 0:
            continue

        auc_c = roc_auc_score(
            y,
            p_clinical,
            sample_weight=weights,
        )
        auc_g = roc_auc_score(
            y,
            p_combined,
            sample_weight=weights,
        )
        deltas.append(auc_g - auc_c)

    if len(deltas) < max(100, int(0.8 * reps)):
        raise RuntimeError(
            "Too few valid paired Poisson bootstrap replicates."
        )

    low, high = np.percentile(deltas, [2.5, 97.5])

    return {
        "clinical_auc": clinical_auc,
        "combined_auc": combined_auc,
        "delta_auc": observed_delta,
        "bootstrap_reps_requested": reps,
        "bootstrap_reps_valid": len(deltas),
        "ci_95_low": float(low),
        "ci_95_high": float(high),
    }


# ============================================================================
# 5. Pregnancy recency query
# ============================================================================

def query_recent_pregnancy() -> pd.DataFrame:
    try:
        from google.cloud import bigquery
    except ImportError as exc:
        raise ImportError(
            "google-cloud-bigquery is required inside the All of Us Workbench."
        ) from exc

    client = bigquery.Client()

    sql = f"""
WITH observation_summary AS (
    SELECT
        person_id,
        MIN(observation_period_start_date) AS ehr_start_date,
        MAX(observation_period_end_date) AS ehr_end_date
    FROM `{CDR_DATASET}.observation_period`
    GROUP BY person_id
),

primary_consent AS (
    SELECT
        o.person_id,
        MIN(o.observation_date) AS index_date
    FROM `{CDR_DATASET}.concept` AS c
    INNER JOIN `{CDR_DATASET}.concept_ancestor` AS ca
        ON c.concept_id = ca.ancestor_concept_id
    INNER JOIN `{CDR_DATASET}.observation` AS o
        ON ca.descendant_concept_id = o.observation_concept_id
    WHERE c.concept_name = 'Consent PII'
      AND c.concept_class_id = 'Module'
    GROUP BY o.person_id
),

ehr_consent AS (
    SELECT
        person_id,
        MIN(observation_date) AS ehr_consent_date
    FROM `{CDR_DATASET}.observation`
    WHERE observation_source_concept_id = 1586099
      AND value_source_concept_id = 1586100
    GROUP BY person_id
),

candidate_a AS (
    SELECT
        p.person_id,
        pc.index_date
    FROM `{CDR_DATASET}.person` AS p
    INNER JOIN observation_summary AS os
        USING (person_id)
    INNER JOIN primary_consent AS pc
        USING (person_id)
    INNER JOIN ehr_consent AS ec
        USING (person_id)
    WHERE os.ehr_start_date <= pc.index_date
      AND os.ehr_end_date > pc.index_date
      AND (
            DATE_DIFF(
                os.ehr_start_date,
                DATE(
                    p.year_of_birth,
                    COALESCE(NULLIF(p.month_of_birth, 0), 7),
                    COALESCE(NULLIF(p.day_of_birth, 0), 1)
                ),
                DAY
            ) / 365.2425
          ) BETWEEN 18 AND 120
      AND (
            DATE_DIFF(
                pc.index_date,
                DATE(
                    p.year_of_birth,
                    COALESCE(NULLIF(p.month_of_birth, 0), 7),
                    COALESCE(NULLIF(p.day_of_birth, 0), 1)
                ),
                DAY
            ) / 365.2425
          ) BETWEEN 18 AND 120
),

recent_pregnancy AS (
    SELECT
        ca.person_id,
        MAX(
            CASE
                WHEN co.condition_start_date >= DATE_SUB(
                        ca.index_date,
                        INTERVAL 365 DAY
                     )
                 AND co.condition_start_date < ca.index_date
                 AND REGEXP_CONTAINS(
                        LOWER(COALESCE(c.concept_name, '')),
                        r'{PREGNANCY_PATTERN}'
                     )
                THEN 1 ELSE 0
            END
        ) AS {PREGNANCY_RECENT}
    FROM candidate_a AS ca
    LEFT JOIN `{CDR_DATASET}.condition_occurrence` AS co
        ON ca.person_id = co.person_id
       AND co.condition_start_date < ca.index_date
    LEFT JOIN `{CDR_DATASET}.concept` AS c
        ON co.condition_concept_id = c.concept_id
    GROUP BY ca.person_id
)

SELECT *
FROM recent_pregnancy
"""

    query_path = DOC_DIR / "17_recent_pregnancy_365d_query.sql"
    query_path.write_text(sql.strip() + "\n")

    progress("Querying recent pregnancy/postpartum sensitivity from CDR.")
    out = (
        client.query(sql)
        .result()
        .to_dataframe(create_bqstorage_client=True)
    )

    if len(out) != EXPECTED_MASTER_N:
        raise RuntimeError(
            f"Recent pregnancy query returned {len(out):,}; "
            f"expected {EXPECTED_MASTER_N:,}."
        )

    if out[PERSON_ID].duplicated().any():
        raise RuntimeError(
            "Recent-pregnancy query returned duplicate person_id."
        )

    out[PREGNANCY_RECENT] = (
        pd.to_numeric(
            out[PREGNANCY_RECENT],
            errors="coerce",
        )
        .fillna(0)
        .clip(0, 1)
        .astype(int)
    )

    return out


# ============================================================================
# 6. Cox helpers
# ============================================================================

def resolve_event_time_column(frame: pd.DataFrame) -> str:
    if EVENT_TIME in frame.columns:
        return EVENT_TIME
    if ALT_EVENT_TIME in frame.columns:
        return ALT_EVENT_TIME
    raise KeyError(
        f"Neither {EVENT_TIME} nor {ALT_EVENT_TIME} exists in master."
    )


def fit_cox_model(
    frame: pd.DataFrame,
    time_col: str,
    event_col: str,
    predictors: list[str],
    model_name: str,
):
    needed = [time_col, event_col, *predictors]

    model_frame = (
        frame[needed]
        .apply(pd.to_numeric, errors="coerce")
        .dropna()
        .astype("float64")
        .copy()
    )

    if model_frame.empty:
        raise RuntimeError(f"{model_name}: empty Cox model frame.")

    X = model_frame[predictors].to_numpy(dtype=float)
    time = model_frame[time_col].to_numpy(dtype=float)
    status = model_frame[event_col].to_numpy(dtype=int)

    result = PHReg(
        endog=time,
        exog=X,
        status=status,
        ties="breslow",
    ).fit(disp=0)

    beta = np.asarray(result.params, dtype=float)
    se = np.asarray(result.bse, dtype=float)

    rows = []

    for i, predictor in enumerate(predictors):
        b = float(beta[i])
        s = float(se[i])
        z = b/s if s > 0 else np.nan

        rows.append(
            {
                "model": model_name,
                "term": predictor,
                "beta": b,
                "standard_error": s,
                "hazard_ratio": math.exp(b),
                "ci_95_low": math.exp(b - 1.96*s),
                "ci_95_high": math.exp(b + 1.96*s),
                "p_value": normal_p_value(z),
                "n": int(len(model_frame)),
                "events": int(status.sum()),
            }
        )

    return pd.DataFrame(rows)


def build_split_records(
    landmark: pd.DataFrame,
    needed_columns: list[str],
) -> pd.DataFrame:
    base = landmark[
        ["tte_time_days", "tte_event", *needed_columns]
    ].copy()

    base["_cluster_id"] = np.arange(len(base), dtype=np.int64)

    parts = []

    for interval_index, spec in enumerate(INTERVALS):
        start = float(spec["start"])
        end = float(spec["end"])

        at_risk = base["tte_time_days"].gt(start)
        sub = base.loc[at_risk].copy()

        if np.isfinite(end):
            local_stop = (
                np.minimum(
                    sub["tte_time_days"].to_numpy(dtype=float),
                    end,
                )
                - start
            )
            local_event = (
                sub["tte_event"].eq(1)
                & sub["tte_time_days"].le(end)
            ).astype(int)
        else:
            local_stop = (
                sub["tte_time_days"].to_numpy(dtype=float)
                - start
            )
            local_event = sub["tte_event"].astype(int)

        sub["interval_index"] = interval_index
        sub["interval_label"] = spec["label"]
        sub["interval_slug"] = spec["slug"]
        sub["interval_stop_days"] = local_stop
        sub["interval_event"] = local_event.to_numpy(dtype=int)

        parts.append(sub)

    return pd.concat(parts, ignore_index=True)


def fit_interval_specific_cox(
    split: pd.DataFrame,
    model_name: str,
    exposure_map: dict[str, str],
    common_covariates: list[str],
):
    work = split.copy()
    interval_terms_by_effect = {}
    predictors = []

    for effect_label, source_column in exposure_map.items():
        terms = []

        for spec in INTERVALS:
            term = f"{source_column}__{spec['slug']}"
            work[term] = np.where(
                work["interval_slug"].eq(spec["slug"]),
                pd.to_numeric(
                    work[source_column],
                    errors="coerce",
                ),
                0.0,
            )
            terms.append(term)
            predictors.append(term)

        interval_terms_by_effect[effect_label] = terms

    predictors.extend(common_covariates)
    predictors = list(dict.fromkeys(predictors))

    needed = [
        "interval_stop_days",
        "interval_event",
        "interval_index",
        *predictors,
    ]

    model_frame = (
        work[needed]
        .apply(pd.to_numeric, errors="coerce")
        .dropna()
        .astype("float64")
        .copy()
    )

    X = model_frame[predictors].to_numpy(dtype=float)
    time = model_frame["interval_stop_days"].to_numpy(dtype=float)
    status = model_frame["interval_event"].to_numpy(dtype=int)
    strata = model_frame["interval_index"].to_numpy(dtype=int)

    result = PHReg(
        endog=time,
        exog=X,
        status=status,
        strata=strata,
        ties="breslow",
    ).fit(disp=0)

    beta = np.asarray(result.params, dtype=float)
    covariance = np.asarray(result.cov_params(), dtype=float)
    se = np.sqrt(np.diag(covariance))

    hr_rows = []
    global_rows = []

    for effect_label, interval_terms in interval_terms_by_effect.items():
        indices = [predictors.index(term) for term in interval_terms]

        for index, term, spec in zip(
            indices,
            interval_terms,
            INTERVALS,
        ):
            b = float(beta[index])
            s = float(se[index])
            z = b/s if s > 0 else np.nan

            hr_rows.append(
                {
                    "model": model_name,
                    "effect": effect_label,
                    "interval": spec["label"],
                    "term": term,
                    "hazard_ratio": math.exp(b),
                    "ci_95_low": math.exp(b - 1.96*s),
                    "ci_95_high": math.exp(b + 1.96*s),
                    "p_value": normal_p_value(z),
                }
            )

        # Equality test HR1 = HR2 = HR3.
        contrast = np.zeros((2, len(predictors)), dtype=float)
        contrast[0, indices[1]] = 1.0
        contrast[0, indices[0]] = -1.0
        contrast[1, indices[2]] = 1.0
        contrast[1, indices[0]] = -1.0

        estimate = contrast @ beta
        V = contrast @ covariance @ contrast.T
        stat = float(
            estimate.T @ np.linalg.pinv(V) @ estimate
        )
        p = float(chi2.sf(stat, 2))

        global_rows.append(
            {
                "model": model_name,
                "effect": effect_label,
                "wald_chi2": stat,
                "df": 2,
                "p_value": p,
            }
        )

    return pd.DataFrame(hr_rows), pd.DataFrame(global_rows)


# ============================================================================
# 7. Ancestry helpers
# ============================================================================

def heterogeneity_test(estimates: pd.DataFrame) -> dict:
    beta = estimates["beta"].to_numpy(dtype=float)
    se = estimates["standard_error"].to_numpy(dtype=float)
    weights = 1.0 / (se**2)

    pooled = float(np.sum(weights*beta) / np.sum(weights))
    Q = float(np.sum(weights*(beta-pooled)**2))
    df = len(beta)-1
    p = float(chi2.sf(Q, df))

    I2 = (
        max(0.0, 100.0*(Q-df)/Q)
        if Q > 0
        else 0.0
    )

    return {
        "cochran_Q": Q,
        "df": df,
        "p_value": p,
        "I2_percent": I2,
    }


# ============================================================================
# 8. Main
# ============================================================================

def main():
    banner()

    if not MASTER_FILE.exists():
        raise FileNotFoundError(
            f"Frozen master not found: {MASTER_FILE}\n"
            "Set VTE_JTH_MASTER_FILE if needed."
        )

    progress(f"Reading frozen master: {MASTER_FILE}")
    master_hash = sha256_file(MASTER_FILE)
    master = read_table(MASTER_FILE)

    if len(master) != EXPECTED_MASTER_N:
        raise RuntimeError(
            f"Master N changed: {len(master):,}; "
            f"expected {EXPECTED_MASTER_N:,}."
        )

    required = [
        PERSON_ID,
        WGS_ELIGIBLE,
        PRIMARY_OUTCOME,
        PRIOR_HISTORY,
        NO_PRIOR_ENDPOINT,
        ANTICOAG,
        OBSERVABILITY,
        "age_at_primary_consent",
        "female_sex",
        "other_unknown_sex",
        *FIVE_FACTORS,
        PREGNANCY_ALL_HISTORY,
        F5,
        F2,
        F5F2,
        *SEVEN_MARKERS,
        *PC_COLUMNS,
    ]

    # Time-to-event and ancestry columns are required for complete Script 17.
    event_time_col = resolve_event_time_column(master)
    required.extend(
        [
            event_time_col,
            FOLLOWUP_TIME,
            ANCESTRY,
        ]
    )

    require_columns(master, list(dict.fromkeys(required)))

    if master[PERSON_ID].duplicated().any():
        raise RuntimeError("Frozen master contains duplicate person_id.")

    for column in list(
        dict.fromkeys(
            [
                WGS_ELIGIBLE,
                PRIMARY_OUTCOME,
                PRIOR_HISTORY,
                NO_PRIOR_ENDPOINT,
                ANTICOAG,
                *FIVE_FACTORS,
                PREGNANCY_ALL_HISTORY,
                F5,
                F2,
                F5F2,
                *SEVEN_MARKERS,
            ]
        )
    ):
        master[column] = numeric_binary(master[column], column)

    master["age_per_10y"] = (
        pd.to_numeric(
            master["age_at_primary_consent"],
            errors="coerce",
        ) / 10.0
    )

    # ----------------------------------------------------------------------
    # A. Revised primary five-factor burden.
    # ----------------------------------------------------------------------

    progress("A. Rebuilding revised five-factor primary burden.")

    master["five_factor_count"] = (
        master[FIVE_FACTORS]
        .astype(float)
        .sum(axis=1)
    ).astype("int8")

    if not master["five_factor_count"].between(0,5).all():
        raise RuntimeError("Five-factor count falls outside 0-5.")

    master["burden_012plus_5f"] = np.select(
        [
            master["five_factor_count"].eq(0),
            master["five_factor_count"].eq(1),
        ],
        ["0","1"],
        default="2+",
    )

    master["burden_0123plus_5f"] = np.select(
        [
            master["five_factor_count"].eq(0),
            master["five_factor_count"].eq(1),
            master["five_factor_count"].eq(2),
        ],
        ["0","1","2"],
        default="3+",
    )

    analysis = master.loc[
        master[WGS_ELIGIBLE].eq(1)
    ].copy()

    if len(analysis) != EXPECTED_WGS_N:
        raise RuntimeError(
            f"WGS N changed: {len(analysis):,}; "
            f"expected {EXPECTED_WGS_N:,}."
        )

    if int(analysis[PRIMARY_OUTCOME].sum()) != EXPECTED_VTE_N:
        raise RuntimeError(
            "Primary VTE count does not reproduce manuscript."
        )

    observed_burden = (
        analysis["burden_012plus_5f"]
        .value_counts()
        .to_dict()
    )
    observed_burden = {
        key: int(observed_burden.get(key,0))
        for key in EXPECTED_5F_BURDEN
    }

    if observed_burden != EXPECTED_5F_BURDEN:
        raise RuntimeError(
            "Five-factor burden does not reproduce Script 16.\n"
            f"Observed: {observed_burden}\n"
            f"Expected: {EXPECTED_5F_BURDEN}"
        )

    for key, expected in EXPECTED_5F_FIGURE3.items():
        burden, carrier = key
        sub = analysis.loc[
            analysis["burden_012plus_5f"].astype("string").eq(burden)
            & analysis[F5F2].eq(carrier)
        ]
        observed = (
            len(sub),
            int(sub[PRIMARY_OUTCOME].sum()),
        )
        if observed != expected:
            raise RuntimeError(
                "Five-factor Figure 3 lock failed for "
                f"{key}: observed {observed}, expected {expected}."
            )

    print("Script 16 five-factor locks: PASS")

    base_covariates = [
        "age_per_10y",
        "female_sex",
        "other_unknown_sex",
        *PC_COLUMNS,
        OBSERVABILITY,
    ]

    genetic_covariates = [
        "five_factor_count",
        *base_covariates,
    ]

    # ----------------------------------------------------------------------
    # B. Main Table 1.
    # ----------------------------------------------------------------------

    progress("B. Building revised main Table 1.")

    table1 = build_table1(analysis.copy())
    table1_path = TABLE_DIR / "17B_Main_Table1_revised_five_factor.csv"
    table1.to_csv(table1_path, index=False)

    # ----------------------------------------------------------------------
    # C. Main Table 2 clinical models.
    # ----------------------------------------------------------------------

    progress("C. Fitting revised primary clinical models.")

    burden_work, burden_terms = make_burden_dummies(
        analysis,
        "burden_012plus_5f",
        ["0","1","2+"],
        "primary5",
    )

    burden_result, _, _ = fit_glm(
        burden_work,
        PRIMARY_OUTCOME,
        [*burden_terms, *base_covariates],
    )

    burden_rows = extract_terms(
        burden_result,
        burden_terms,
        "primary_five_factor_burden",
        "main_table2",
    )

    factor_result, _, _ = fit_glm(
        analysis,
        PRIMARY_OUTCOME,
        [*FIVE_FACTORS, *base_covariates],
    )

    factor_rows = extract_terms(
        factor_result,
        FIVE_FACTORS,
        "primary_five_factor_mutually_adjusted",
        "main_table2",
    )

    table2 = pd.concat(
        [burden_rows, factor_rows],
        ignore_index=True,
    )
    table2_path = TABLE_DIR / "17C_Main_Table2_adjusted_clinical.csv"
    table2.to_csv(table2_path, index=False)

    # Verify Script 16 ORs.
    one_or = float(
        burden_rows.loc[
            burden_rows["term"].str.endswith("_1"),
            "odds_ratio",
        ].iloc[0]
    )
    twop_or = float(
        burden_rows.loc[
            burden_rows["term"].str.endswith("_2plus"),
            "odds_ratio",
        ].iloc[0]
    )

    if abs(one_or-EXPECTED_5F_OR_1) > EXPECTED_5F_OR_TOL:
        raise RuntimeError(
            f"1-factor OR does not reproduce Script 16: {one_or:.3f}"
        )

    if abs(twop_or-EXPECTED_5F_OR_2PLUS) > EXPECTED_5F_OR_TOL:
        raise RuntimeError(
            f"2+-factor OR does not reproduce Script 16: {twop_or:.3f}"
        )

    # ----------------------------------------------------------------------
    # D. Main Table 3: separate seven-marker models using exact 5-factor count.
    # ----------------------------------------------------------------------

    progress("D. Rebuilding seven targeted variant models.")

    marker_frames = []

    for marker in SEVEN_MARKERS:
        result, model_frame, _ = fit_glm(
            analysis,
            PRIMARY_OUTCOME,
            [marker, *genetic_covariates],
        )

        row = extract_terms(
            result,
            [marker],
            f"separate_{marker}",
            "main_table3",
        )

        marker_positive = analysis.loc[
            analysis[marker].eq(1)
        ]

        row["marker_label"] = MARKER_LABELS[marker]
        row["marker_positive_n"] = int(len(marker_positive))
        row["marker_positive_percent"] = (
            100.0 * len(marker_positive) / len(analysis)
        )
        row["vte_among_marker_positive_n"] = int(
            marker_positive[PRIMARY_OUTCOME].sum()
        )
        row["vte_among_marker_positive_percent"] = (
            100.0
            * marker_positive[PRIMARY_OUTCOME].mean()
        )

        marker_frames.append(row)

    table3 = pd.concat(marker_frames, ignore_index=True)

    # Same multiplicity structure as current manuscript:
    # Holm for F5/F2 pair, BH for remaining 5 markers.
    table3["adjusted_p"] = np.nan
    focused_mask = table3["term"].isin([F5,F2])
    context_mask = ~focused_mask

    table3.loc[
        focused_mask,
        "adjusted_p",
    ] = multipletests(
        table3.loc[focused_mask, "p_value"].to_numpy(dtype=float),
        method="holm",
    )[1]

    table3.loc[
        context_mask,
        "adjusted_p",
    ] = multipletests(
        table3.loc[context_mask, "p_value"].to_numpy(dtype=float),
        method="fdr_bh",
    )[1]

    table3["focused_role"] = np.where(
        table3["term"].isin([F5,F2]),
        "Carried forward",
        "Panel context",
    )

    table3_path = TABLE_DIR / "17D_Main_Table3_targeted_variants.csv"
    table3.to_csv(table3_path, index=False)

    # Derived F5/F2 supportive adjusted model.
    composite_result, _, _ = fit_glm(
        analysis,
        PRIMARY_OUTCOME,
        [F5F2, *genetic_covariates],
    )
    composite_rows = extract_terms(
        composite_result,
        [F5F2],
        "adjusted_f5f2_composite_five_factor",
        "supportive_composite_model",
    )
    composite_path = TABLE_DIR / "17D_adjusted_F5F2_composite.csv"
    composite_rows.to_csv(composite_path, index=False)

    # ----------------------------------------------------------------------
    # E. Main Figure 3.
    # ----------------------------------------------------------------------

    progress("E. Rebuilding Figure 3 cells and within-category tests.")

    figure3 = burden_by_exposure_descriptive(
        analysis,
        "burden_012plus_5f",
        F5F2,
    )

    figure3_tests = figure3_chisquare(analysis)

    figure3_path = TABLE_DIR / "17E_Main_Figure3_data.csv"
    figure3.to_csv(figure3_path, index=False)

    figure3_tests_path = TABLE_DIR / "17E_Main_Figure3_chisquare.csv"
    figure3_tests.to_csv(figure3_tests_path, index=False)

    # ----------------------------------------------------------------------
    # F. Formal interaction using revised burden.
    # ----------------------------------------------------------------------

    progress("F. Formal F5/F2 x revised five-factor burden interaction.")

    interaction_fit = build_interaction_model(
        analysis,
        PRIMARY_OUTCOME,
        F5F2,
        "burden_012plus_5f",
        ["0","1","2+"],
        base_covariates,
        "fivefactor_primary",
    )

    mult_global = joint_wald(
        interaction_fit["result"],
        interaction_fit["interaction_terms"],
    )

    (
        adjusted_risks,
        additive_contrasts,
        add_global,
    ) = additive_interaction_table(interaction_fit)

    interaction_global = pd.DataFrame(
        [
            {
                "scale": "multiplicative_odds",
                **mult_global,
            },
            {
                "scale": "adjusted_probability_difference",
                **add_global,
            },
        ]
    )

    interaction_global_path = (
        TABLE_DIR / "17F_F5F2_interaction_global.csv"
    )
    interaction_global.to_csv(
        interaction_global_path,
        index=False,
    )

    adjusted_risks_path = (
        TABLE_DIR / "17F_F5F2_adjusted_probabilities.csv"
    )
    adjusted_risks.to_csv(
        adjusted_risks_path,
        index=False,
    )

    additive_contrasts_path = (
        TABLE_DIR / "17F_F5F2_additive_contrasts.csv"
    )
    additive_contrasts.to_csv(
        additive_contrasts_path,
        index=False,
    )

    if mult_global["p_value"] >= 0.001:
        raise RuntimeError(
            "Five-factor multiplicative interaction no longer reproduces Script 16."
        )

    if abs(add_global["p_value"]-EXPECTED_5F_ADD_P) > 0.01:
        raise RuntimeError(
            "Five-factor additive interaction no longer reproduces Script 16."
        )

    # ----------------------------------------------------------------------
    # G. Revised Supplementary Table S6 components.
    # ----------------------------------------------------------------------

    progress("G. Rebuilding revised robustness / S6 analyses.")

    # Granularity: 0 / 1 / 2 / >=3.
    finer_work, finer_terms = make_burden_dummies(
        analysis,
        "burden_0123plus_5f",
        ["0","1","2","3+"],
        "finer5",
    )

    finer_result, _, _ = fit_glm(
        finer_work,
        PRIMARY_OUTCOME,
        [*finer_terms, *base_covariates],
    )

    finer_rows = extract_terms(
        finer_result,
        finer_terms,
        "five_factor_0_1_2_3plus",
        "supplement_s6_granularity",
    )
    finer_path = TABLE_DIR / "17G_S6_finer_burden.csv"
    finer_rows.to_csv(finer_path, index=False)

    # Exact 0-5 count trend.
    count_result, _, _ = fit_glm(
        analysis,
        PRIMARY_OUTCOME,
        ["five_factor_count", *base_covariates],
    )
    count_rows = extract_terms(
        count_result,
        ["five_factor_count"],
        "exact_five_factor_count",
        "supplement_s6_granularity",
    )
    count_path = TABLE_DIR / "17G_S6_exact_count_trend.csv"
    count_rows.to_csv(count_path, index=False)

    # Exact count x F5/F2 interaction.
    exact_work = analysis.copy()
    exact_work["f5f2_x_five_factor_count"] = (
        exact_work[F5F2]
        * exact_work["five_factor_count"]
    )

    exact_int_result, _, _ = fit_glm(
        exact_work,
        PRIMARY_OUTCOME,
        [
            F5F2,
            "five_factor_count",
            "f5f2_x_five_factor_count",
            *base_covariates,
        ],
    )

    exact_int_rows = extract_terms(
        exact_int_result,
        [
            F5F2,
            "five_factor_count",
            "f5f2_x_five_factor_count",
        ],
        "exact_count_f5f2_interaction",
        "supplement_s6_interaction",
    )
    exact_int_path = TABLE_DIR / "17G_S6_exact_count_interaction.csv"
    exact_int_rows.to_csv(exact_int_path, index=False)

    # Joint seven-marker sensitivity.
    joint_result, joint_frame, joint_X = fit_glm(
        analysis,
        PRIMARY_OUTCOME,
        [
            *SEVEN_MARKERS,
            "five_factor_count",
            *base_covariates,
        ],
    )

    joint_rows = extract_terms(
        joint_result,
        SEVEN_MARKERS,
        "joint_seven_marker_five_factor",
        "supplement_s6_joint_marker",
    )

    joint_path = TABLE_DIR / "17G_S6_joint_seven_marker.csv"
    joint_rows.to_csv(joint_path, index=False)

    # No-prior-VTE sensitivity.
    no_prior = analysis.loc[
        analysis[PRIOR_HISTORY].eq(0)
    ].copy()

    mismatch = int(
        (
            no_prior[PRIMARY_OUTCOME]
            != no_prior[NO_PRIOR_ENDPOINT]
        ).sum()
    )

    if mismatch != 0:
        raise RuntimeError(
            "No-prior endpoint mismatch after excluding prior VTE."
        )

    no_prior_work, no_prior_burden_terms = make_burden_dummies(
        no_prior,
        "burden_012plus_5f",
        ["0","1","2+"],
        "noprior5",
    )

    no_prior_burden_result, _, _ = fit_glm(
        no_prior_work,
        PRIMARY_OUTCOME,
        [*no_prior_burden_terms, *base_covariates],
    )

    no_prior_burden_rows = extract_terms(
        no_prior_burden_result,
        no_prior_burden_terms,
        "no_prior_five_factor_burden",
        "supplement_s6_no_prior",
    )

    no_prior_composite_result, _, _ = fit_glm(
        no_prior,
        PRIMARY_OUTCOME,
        [F5F2, "five_factor_count", *base_covariates],
    )

    no_prior_composite_rows = extract_terms(
        no_prior_composite_result,
        [F5F2],
        "no_prior_f5f2_five_factor",
        "supplement_s6_no_prior",
    )

    no_prior_rows = pd.concat(
        [no_prior_burden_rows, no_prior_composite_rows],
        ignore_index=True,
    )

    no_prior_path = TABLE_DIR / "17G_S6_no_prior_VTE.csv"
    no_prior_rows.to_csv(no_prior_path, index=False)

    # Anticoagulant treatment-context sensitivity.
    anticoag_work, anticoag_burden_terms = make_burden_dummies(
        analysis,
        "burden_012plus_5f",
        ["0","1","2+"],
        "anticoag5",
    )

    anticoag_burden_result, _, _ = fit_glm(
        anticoag_work,
        PRIMARY_OUTCOME,
        [
            *anticoag_burden_terms,
            *base_covariates,
            ANTICOAG,
        ],
    )

    anticoag_frames = [
        extract_terms(
            anticoag_burden_result,
            anticoag_burden_terms,
            "five_factor_burden_plus_anticoagulant",
            "supplement_s6_anticoagulant",
        )
    ]

    for exposure in [F5,F2,F5F2]:
        result, _, _ = fit_glm(
            analysis,
            PRIMARY_OUTCOME,
            [
                exposure,
                "five_factor_count",
                *base_covariates,
                ANTICOAG,
            ],
        )

        anticoag_frames.append(
            extract_terms(
                result,
                [exposure],
                f"{exposure}_plus_anticoagulant",
                "supplement_s6_anticoagulant",
            )
        )

    anticoag_rows = pd.concat(
        anticoag_frames,
        ignore_index=True,
    )
    anticoag_path = TABLE_DIR / "17G_S6_anticoagulant.csv"
    anticoag_rows.to_csv(anticoag_path, index=False)

    # ----------------------------------------------------------------------
    # H. Revised Supplementary Table S1.
    # ----------------------------------------------------------------------

    progress("H. Rebuilding OOF discrimination and calibration.")

    clinical_features = [
        "age_per_10y",
        "female_sex",
        "other_unknown_sex",
        *PC_COLUMNS,
        *FIVE_FACTORS,
        OBSERVABILITY,
    ]

    combined_features = [
        *clinical_features,
        F5,
        F2,
    ]

    y, predictions, performance = out_of_fold_predictions(
        analysis,
        PRIMARY_OUTCOME,
        {
            "clinical_five_factor": clinical_features,
            "clinical_five_factor_plus_F5_F2": combined_features,
        },
    )

    performance_path = (
        TABLE_DIR / "17H_S1_OOF_performance_calibration.csv"
    )
    performance.to_csv(performance_path, index=False)

    bootstrap = paired_poisson_delta_auc(
        y,
        predictions["clinical_five_factor"],
        predictions["clinical_five_factor_plus_F5_F2"],
        BOOTSTRAP_REPS,
        BOOTSTRAP_SEED,
    )

    bootstrap_path = (
        TABLE_DIR / "17H_S1_delta_AUC_bootstrap.csv"
    )
    pd.DataFrame([bootstrap]).to_csv(
        bootstrap_path,
        index=False,
    )

    # ----------------------------------------------------------------------
    # I/J. Revised five-factor time-to-event / interval-specific S7.
    # ----------------------------------------------------------------------

    progress("I. Rebuilding day-7 person-time sensitivity with five-factor burden.")

    analysis[event_time_col] = pd.to_numeric(
        analysis[event_time_col],
        errors="coerce",
    )
    analysis[FOLLOWUP_TIME] = pd.to_numeric(
        analysis[FOLLOWUP_TIME],
        errors="coerce",
    )

    event_mask = analysis[PRIMARY_OUTCOME].eq(1)

    if analysis.loc[event_mask, event_time_col].isna().any():
        raise RuntimeError("Missing event times among VTE participants.")

    if analysis[FOLLOWUP_TIME].isna().any():
        raise RuntimeError("Missing follow-up times in WGS cohort.")

    early_vte = (
        event_mask
        & analysis[event_time_col].le(LANDMARK_DAY)
    )

    short_followup_noncase = (
        analysis[PRIMARY_OUTCOME].eq(0)
        & analysis[FOLLOWUP_TIME].le(LANDMARK_DAY)
    )

    if int(early_vte.sum()) != EXPECTED_EARLY_VTE_N:
        raise RuntimeError(
            f"Early VTE lock changed: {int(early_vte.sum())}"
        )

    if int(short_followup_noncase.sum()) != EXPECTED_SHORT_FOLLOWUP_NONCASE_N:
        raise RuntimeError(
            "Short-follow-up noncase lock changed."
        )

    landmark = analysis.loc[
        ~(early_vte | short_followup_noncase)
    ].copy()

    landmark["tte_event"] = landmark[PRIMARY_OUTCOME].astype(int)

    landmark["tte_time_days"] = np.where(
        landmark["tte_event"].eq(1),
        landmark[event_time_col] - LANDMARK_DAY,
        landmark[FOLLOWUP_TIME] - LANDMARK_DAY,
    )

    if len(landmark) != EXPECTED_LANDMARK_N:
        raise RuntimeError(
            f"Landmark N changed: {len(landmark):,}"
        )

    if int(landmark["tte_event"].sum()) != EXPECTED_LANDMARK_VTE_N:
        raise RuntimeError(
            "Landmark VTE count changed."
        )

    if not landmark["tte_time_days"].gt(0).all():
        raise RuntimeError(
            "Nonpositive time after day-7 landmark."
        )

    landmark_summary = pd.DataFrame(
        [
            {
                "primary_wgs_n": len(analysis),
                "excluded_vte_days_0_7": int(early_vte.sum()),
                "excluded_no_vte_followup_le_7": int(
                    short_followup_noncase.sum()
                ),
                "landmark_n": len(landmark),
                "vte_after_day7": int(landmark["tte_event"].sum()),
                "total_person_years": float(
                    landmark["tte_time_days"].sum()/DAYS_PER_YEAR
                ),
                "median_followup_days": float(
                    landmark["tte_time_days"].median()
                ),
                "q1_followup_days": float(
                    landmark["tte_time_days"].quantile(0.25)
                ),
                "q3_followup_days": float(
                    landmark["tte_time_days"].quantile(0.75)
                ),
            }
        ]
    )

    py = float(
        landmark["tte_time_days"].sum()/DAYS_PER_YEAR
    )
    landmark_summary["crude_vte_incidence_per_1000_person_years"] = (
        1000.0
        * int(landmark["tte_event"].sum())
        / py
    )

    landmark_summary_path = (
        TABLE_DIR / "17I_S7_landmark_summary.csv"
    )
    landmark_summary.to_csv(
        landmark_summary_path,
        index=False,
    )

    # Overall Cox models.
    landmark_work, lm_burden_terms = make_burden_dummies(
        landmark,
        "burden_012plus_5f",
        ["0","1","2+"],
        "cox5",
    )

    cox_frames = []

    burden_cox = fit_cox_model(
        landmark_work,
        "tte_time_days",
        "tte_event",
        [*lm_burden_terms, *base_covariates],
        "five_factor_burden_cox",
    )
    cox_frames.append(burden_cox)

    for exposure in [F5,F2,F5F2]:
        cox_frames.append(
            fit_cox_model(
                landmark,
                "tte_time_days",
                "tte_event",
                [
                    exposure,
                    "five_factor_count",
                    *base_covariates,
                ],
                f"{exposure}_cox_five_factor_adjusted",
            )
        )

    overall_cox = pd.concat(cox_frames, ignore_index=True)
    overall_cox_path = TABLE_DIR / "17I_S7_overall_cox.csv"
    overall_cox.to_csv(overall_cox_path, index=False)

    # Interval-specific.
    needed_for_split = list(
        dict.fromkeys(
            [
                *lm_burden_terms,
                F5,
                F2,
                F5F2,
                "five_factor_count",
                *base_covariates,
            ]
        )
    )

    split = build_split_records(
        landmark_work,
        needed_for_split,
    )

    burden_interval_hr, burden_time_tests = (
        fit_interval_specific_cox(
            split,
            "five_factor_burden_time_varying",
            {
                "1 clinical factor vs 0": lm_burden_terms[0],
                "2+ clinical factors vs 0": lm_burden_terms[1],
            },
            base_covariates,
        )
    )

    genetic_interval_frames = []
    genetic_time_frames = []

    for exposure, label in [
        (F5, "Factor V Leiden"),
        (F2, "Prothrombin G20210A"),
        (F5F2, "Either F5 or F2"),
    ]:
        hr, gt = fit_interval_specific_cox(
            split,
            f"{exposure}_time_varying",
            {label: exposure},
            [
                "five_factor_count",
                *base_covariates,
            ],
        )
        genetic_interval_frames.append(hr)
        genetic_time_frames.append(gt)

    interval_hr = pd.concat(
        [burden_interval_hr, *genetic_interval_frames],
        ignore_index=True,
    )

    time_tests = pd.concat(
        [burden_time_tests, *genetic_time_frames],
        ignore_index=True,
    )

    interval_hr_path = (
        TABLE_DIR / "17J_S7_interval_specific_cox.csv"
    )
    interval_hr.to_csv(interval_hr_path, index=False)

    time_tests_path = (
        TABLE_DIR / "17J_S7_time_variation_tests.csv"
    )
    time_tests.to_csv(time_tests_path, index=False)

    # ----------------------------------------------------------------------
    # K. Revised ancestry-stratified S8.
    # ----------------------------------------------------------------------

    progress("K. Rebuilding ancestry-stratified F5/F2 estimates.")

    ancestry_rows = []
    heterogeneity_rows = []

    for marker in [F5,F2]:
        marker_estimates = []

        for group in REPORTABLE_ANCESTRIES:
            sub = analysis.loc[
                analysis[ANCESTRY].astype("string").eq(group)
            ].copy()

            if len(sub) < MIN_ANCESTRY_N:
                continue

            cross = pd.crosstab(
                sub[marker].astype(int),
                sub[PRIMARY_OUTCOME].astype(int),
            ).reindex(
                index=[0,1],
                columns=[0,1],
                fill_value=0,
            )

            min_cell = int(cross.to_numpy().min())

            if min_cell < MIN_2X2_CELL:
                continue

            result, model_frame, _ = fit_glm(
                sub,
                PRIMARY_OUTCOME,
                [
                    marker,
                    "five_factor_count",
                    *base_covariates,
                ],
            )

            row = extract_terms(
                result,
                [marker],
                f"{group}_{marker}",
                "supplement_s8_ancestry",
            ).iloc[0].to_dict()

            row["ancestry"] = group
            row["ancestry_n"] = len(sub)
            row["min_variant_vte_2x2_cell"] = min_cell

            ancestry_rows.append(row)
            marker_estimates.append(row)

        marker_df = pd.DataFrame(marker_estimates)

        if len(marker_df) >= 2:
            het = heterogeneity_test(marker_df)
            heterogeneity_rows.append(
                {
                    "marker": marker,
                    **het,
                }
            )

    ancestry_results = pd.DataFrame(ancestry_rows)
    heterogeneity = pd.DataFrame(heterogeneity_rows)

    if not heterogeneity.empty:
        heterogeneity["holm_adjusted_p"] = multipletests(
            heterogeneity["p_value"].to_numpy(dtype=float),
            method="holm",
        )[1]

    ancestry_path = (
        TABLE_DIR / "17K_S8_ancestry_stratified.csv"
    )
    ancestry_results.to_csv(ancestry_path, index=False)

    heterogeneity_path = (
        TABLE_DIR / "17K_S8_heterogeneity.csv"
    )
    heterogeneity.to_csv(heterogeneity_path, index=False)

    # ----------------------------------------------------------------------
    # L. Pregnancy/postpartum reviewer-response table.
    # ----------------------------------------------------------------------

    progress("L. Rebuilding pregnancy/postpartum explanatory sensitivity.")

    # All-history pregnancy factor in the original six-factor mutually
    # adjusted model, retained only to explain the review issue.
    all6 = [
        *FIVE_FACTORS,
        PREGNANCY_ALL_HISTORY,
    ]

    all_history_result, _, _ = fit_glm(
        analysis,
        PRIMARY_OUTCOME,
        [*all6, *base_covariates],
    )

    all_history_preg = extract_terms(
        all_history_result,
        [PREGNANCY_ALL_HISTORY],
        "all_history_pregnancy",
        "reviewer_response_only",
    ).iloc[0]

    recent = query_recent_pregnancy()

    analysis_recent = analysis.merge(
        recent[[PERSON_ID, PREGNANCY_RECENT]],
        on=PERSON_ID,
        how="left",
        validate="one_to_one",
    )

    if analysis_recent[PREGNANCY_RECENT].isna().any():
        raise RuntimeError(
            "Recent pregnancy did not match all WGS participants."
        )

    analysis_recent[PREGNANCY_RECENT] = numeric_binary(
        analysis_recent[PREGNANCY_RECENT],
        PREGNANCY_RECENT,
    )

    recent6 = [
        *FIVE_FACTORS,
        PREGNANCY_RECENT,
    ]

    recent_result, _, _ = fit_glm(
        analysis_recent,
        PRIMARY_OUTCOME,
        [*recent6, *base_covariates],
    )

    recent_preg = extract_terms(
        recent_result,
        [PREGNANCY_RECENT],
        "recent_365d_pregnancy",
        "reviewer_response_only",
    ).iloc[0]

    pregnancy_summary = pd.DataFrame(
        [
            {
                "definition": "all captured pre-enrollment pregnancy/postpartum history",
                "positive_n": int(
                    analysis[PREGNANCY_ALL_HISTORY].sum()
                ),
                "adjusted_or": float(
                    all_history_preg["odds_ratio"]
                ),
                "ci_95_low": float(
                    all_history_preg["ci_95_low"]
                ),
                "ci_95_high": float(
                    all_history_preg["ci_95_high"]
                ),
                "p_value": float(
                    all_history_preg["p_value"]
                ),
                "included_in_revised_primary_burden": False,
            },
            {
                "definition": "pregnancy/postpartum within 365 days before primary consent",
                "positive_n": int(
                    analysis_recent[PREGNANCY_RECENT].sum()
                ),
                "adjusted_or": float(
                    recent_preg["odds_ratio"]
                ),
                "ci_95_low": float(
                    recent_preg["ci_95_low"]
                ),
                "ci_95_high": float(
                    recent_preg["ci_95_high"]
                ),
                "p_value": float(
                    recent_preg["p_value"]
                ),
                "included_in_revised_primary_burden": False,
            },
        ]
    )

    pregnancy_path = (
        TABLE_DIR / "17L_pregnancy_postpartum_reviewer_response.csv"
    )
    pregnancy_summary.to_csv(pregnancy_path, index=False)

    if int(
        analysis[PREGNANCY_ALL_HISTORY].sum()
    ) != EXPECTED_PREG_ALL_HISTORY_N:
        raise RuntimeError(
            "All-history pregnancy count changed."
        )

    if int(
        analysis_recent[PREGNANCY_RECENT].sum()
    ) != EXPECTED_PREG_RECENT_N:
        raise RuntimeError(
            "Recent-365d pregnancy count changed."
        )

    if abs(
        float(recent_preg["odds_ratio"])
        - EXPECTED_PREG_RECENT_OR
    ) > EXPECTED_PREG_RECENT_OR_TOL:
        raise RuntimeError(
            "Recent pregnancy OR no longer reproduces Script 16."
        )

    # ----------------------------------------------------------------------
    # M. Paste-ready manuscript/reviewer summary.
    # ----------------------------------------------------------------------

    progress("M. Writing paste-ready summary.")

    one_row = burden_rows.loc[
        burden_rows["term"].str.endswith("_1")
    ].iloc[0]

    twop_row = burden_rows.loc[
        burden_rows["term"].str.endswith("_2plus")
    ].iloc[0]

    composite_row = composite_rows.iloc[0]

    clinical_perf = performance.loc[
        performance["model"].eq("clinical_five_factor")
        & performance["evaluation"].eq("pooled_oof")
    ].iloc[0]

    combined_perf = performance.loc[
        performance["model"].eq("clinical_five_factor_plus_F5_F2")
        & performance["evaluation"].eq("pooled_oof")
    ].iloc[0]

    no_prior_1 = no_prior_burden_rows.loc[
        no_prior_burden_rows["term"].str.endswith("_1")
    ].iloc[0]
    no_prior_2 = no_prior_burden_rows.loc[
        no_prior_burden_rows["term"].str.endswith("_2plus")
    ].iloc[0]
    no_prior_f5f2 = no_prior_composite_rows.iloc[0]

    paste_ready = f"""
====================================================================================
VTE JTH SCRIPT 17 — REVISED FIVE-FACTOR PRIMARY MANUSCRIPT RESULTS
====================================================================================

COHORT LOCKS
Master N: {len(master):,}
WGS N: {len(analysis):,}
Observed VTE during follow-up: {int(analysis[PRIMARY_OUTCOME].sum()):,}

REVISED FIVE-FACTOR BURDEN
0 factors: {observed_burden['0']:,}
1 factor: {observed_burden['1']:,}
2+ factors: {observed_burden['2+']:,}

Adjusted burden associations:
1 factor vs 0:
  OR {one_row['odds_ratio']:.3f}
  95% CI {one_row['ci_95_low']:.3f}-{one_row['ci_95_high']:.3f}
  P {p_text(one_row['p_value'])}

2+ factors vs 0:
  OR {twop_row['odds_ratio']:.3f}
  95% CI {twop_row['ci_95_low']:.3f}-{twop_row['ci_95_high']:.3f}
  P {p_text(twop_row['p_value'])}

DERIVED EITHER F5/F2 ASSOCIATION
OR {composite_row['odds_ratio']:.3f}
95% CI {composite_row['ci_95_low']:.3f}-{composite_row['ci_95_high']:.3f}
P {p_text(composite_row['p_value'])}

FIGURE 3 — OBSERVED VTE OCCURRENCE
{figure3.to_string(index=False)}

FORMAL F5/F2 x FIVE-FACTOR BURDEN INTERACTION
Multiplicative:
  chi2={mult_global['wald_chi2']:.3f}, df={mult_global['df']},
  P {p_text(mult_global['p_value'])}

Adjusted probability-difference:
  chi2={add_global['wald_chi2']:.3f}, df={add_global['df']},
  P {p_text(add_global['p_value'])}

ADJUSTED PROBABILITY DIFFERENCES
{additive_contrasts.loc[
    additive_contrasts['contrast'].eq('carrier_minus_noncarrier_risk_difference')
].to_string(index=False)}

NO-PRIOR-VTE SENSITIVITY
N: {len(no_prior):,}
VTE: {int(no_prior[PRIMARY_OUTCOME].sum()):,}
1 factor vs 0 OR:
  {no_prior_1['odds_ratio']:.3f}
  95% CI {no_prior_1['ci_95_low']:.3f}-{no_prior_1['ci_95_high']:.3f}
2+ factors vs 0 OR:
  {no_prior_2['odds_ratio']:.3f}
  95% CI {no_prior_2['ci_95_low']:.3f}-{no_prior_2['ci_95_high']:.3f}
Either F5/F2 OR:
  {no_prior_f5f2['odds_ratio']:.3f}
  95% CI {no_prior_f5f2['ci_95_low']:.3f}-{no_prior_f5f2['ci_95_high']:.3f}

PREGNANCY/POSTPARTUM REVIEWER RESPONSE
All-history pregnancy/postpartum:
  N {int(analysis[PREGNANCY_ALL_HISTORY].sum()):,}
  OR {all_history_preg['odds_ratio']:.3f}
  95% CI {all_history_preg['ci_95_low']:.3f}-{all_history_preg['ci_95_high']:.3f}
  P {p_text(all_history_preg['p_value'])}

Recent pregnancy/postpartum within 365 days pre-consent:
  N {int(analysis_recent[PREGNANCY_RECENT].sum()):,}
  OR {recent_preg['odds_ratio']:.3f}
  95% CI {recent_preg['ci_95_low']:.3f}-{recent_preg['ci_95_high']:.3f}
  P {p_text(recent_preg['p_value'])}

OOF PERFORMANCE — REVISED FIVE-FACTOR CLINICAL MODEL
Clinical AUC:
  {clinical_perf['roc_auc']:.4f}
Clinical + F5/F2 AUC:
  {combined_perf['roc_auc']:.4f}
Delta AUC:
  {bootstrap['delta_auc']:+.4f}
95% CI:
  {bootstrap['ci_95_low']:+.4f} to {bootstrap['ci_95_high']:+.4f}

Clinical Brier:
  {clinical_perf['brier_score']:.6f}
Combined Brier:
  {combined_perf['brier_score']:.6f}
Clinical calibration slope:
  {clinical_perf['calibration_slope']:.4f}
Combined calibration slope:
  {combined_perf['calibration_slope']:.4f}

INTERPRETATION GUARDRAILS
1. Primary clinical burden now contains FIVE factors.
2. Pregnancy/postpartum is evaluated separately because its VTE association is time-dependent.
3. Do NOT call the all-history pregnancy OR biologically protective.
4. Primary endpoint language: "observed VTE occurrence during follow-up."
5. No-prior-VTE attenuation must be reported explicitly; do not say only "directionally consistent."
6. Interval-specific Cox differences reflect fixed versus transient exposure definitions and timing.
   Do NOT describe genetics as biologically more temporally stable than clinical risk.
7. Figure 3 remains descriptive and unadjusted; formal interaction is separate.
8. Interaction is statistical effect modification, not biological synergy or causation.
====================================================================================
""".strip() + "\n"

    paste_ready_path = OUTPUT_ROOT / "17_PASTE_READY_RESULTS.txt"
    paste_ready_path.write_text(paste_ready)

    # ----------------------------------------------------------------------
    # Reviewer-response draft snippets.
    # ----------------------------------------------------------------------

    reviewer_note = f"""
# Script 17 — reviewer-response interpretation

## Pregnancy/postpartum concern

The revised primary clinical-risk burden contains five factors: cancer history,
major surgery, inpatient hospitalization, infection/sepsis, and fracture/major
trauma. Pregnancy/postpartum is no longer counted in the accumulated burden
because an all-history pre-enrollment EHR indicator does not preserve the
time-limited nature of pregnancy-associated VTE risk.

The revised five-factor burden remained strongly graded:
- 1 factor vs 0: OR {one_row['odds_ratio']:.2f} (95% CI {one_row['ci_95_low']:.2f}-{one_row['ci_95_high']:.2f})
- >=2 factors vs 0: OR {twop_row['odds_ratio']:.2f} (95% CI {twop_row['ci_95_low']:.2f}-{twop_row['ci_95_high']:.2f})

The all-history pregnancy/postpartum indicator had OR
{all_history_preg['odds_ratio']:.2f}, whereas restricting pregnancy/postpartum
to the 365 days before enrollment yielded OR {recent_preg['odds_ratio']:.2f}
(95% CI {recent_preg['ci_95_low']:.2f}-{recent_preg['ci_95_high']:.2f}).
The inverse all-history estimate is therefore not interpreted as biological
protection.

The clinical-genetic layering result was preserved using the revised burden:
multiplicative interaction P={p_text(mult_global['p_value'])}; adjusted
probability-difference interaction P={p_text(add_global['p_value'])}.

## Prior VTE concern

The no-prior-VTE sensitivity should be reported with its actual five-factor
burden estimates rather than described only as directionally consistent:
- 1 factor vs 0: OR {no_prior_1['odds_ratio']:.2f}
- >=2 factors vs 0: OR {no_prior_2['odds_ratio']:.2f}
- either F5/F2: OR {no_prior_f5f2['odds_ratio']:.2f}

Any attenuation relative to the primary model should be acknowledged as
evidence that recurrence or re-documentation of prior VTE may contribute to
the magnitude of the primary clinical-burden association.

## Temporal-stability concern

The interval-specific analysis should be described as a property of exposure
timing and study design. Genetic variants are fixed, while several clinical
components are transient pre-enrollment exposures. Do not contrast these as
evidence of greater biological "stability" of genetic risk.
""".strip() + "\n"

    reviewer_note_path = (
        DOC_DIR / "17_reviewer_response_interpretation.md"
    )
    reviewer_note_path.write_text(reviewer_note)

    # Methods note.
    methods_note = """
# Revised five-factor primary analysis methods

The revised primary clinical-risk burden is the unweighted count of five
pre-enrollment factors: cancer history, major surgery, inpatient hospitalization,
infection/sepsis, and fracture/major trauma. Pregnancy/postpartum is evaluated
separately because an all-history pre-enrollment EHR indicator does not
necessarily represent pregnancy or postpartum exposure proximal to the VTE-risk
period.

The five-factor burden is categorized as 0, 1, or >=2 factors for the clinically
readable comparison, while the exact 0-5 count is used as the clinical-burden
adjustment covariate in genetic models. Individual clinical factors are entered
together into one mutually adjusted logistic model. Both clinical models adjust
for age, recorded gender, 16 ancestry principal components, and EHR
observability and use HC0 robust standard errors.

Each of the seven targeted variants is evaluated in a separate model adjusted
for the exact five-factor burden, age, recorded gender, 16 ancestry principal
components, and EHR observability. The F5/F2 derived carrier indicator is a
supporting clinical-layering variable, not a third genetic variant.

Figure 3 reports observed, unadjusted VTE proportions within revised five-factor
burden categories. Formal interaction is tested separately using multiplicative
logistic interaction terms and adjusted probability-difference interaction via
marginal standardization.

The supporting prediction analysis compares five-fold out-of-fold predictions
from a clinical model containing age, recorded gender, 16 ancestry principal
components, the five individual clinical factors, and EHR observability with
the same model additionally containing separate F5 and F2 indicators.

The day-7 person-time analysis remains a post hoc sensitivity analysis. Clinical
burden interval-specific estimates should not be compared with genetic
interval-specific estimates as evidence of intrinsic biological risk stability,
because the genetic exposures are fixed while several clinical exposures are
transient and measured only before enrollment.
""".strip() + "\n"

    methods_path = DOC_DIR / "17_revised_five_factor_methods.md"
    methods_path.write_text(methods_note)

    # ----------------------------------------------------------------------
    # QC.
    # ----------------------------------------------------------------------

    qc_rows = [
        {
            "check": "master_n",
            "status": len(master) == EXPECTED_MASTER_N,
            "value": len(master),
        },
        {
            "check": "wgs_n",
            "status": len(analysis) == EXPECTED_WGS_N,
            "value": len(analysis),
        },
        {
            "check": "vte_n",
            "status": int(analysis[PRIMARY_OUTCOME].sum()) == EXPECTED_VTE_N,
            "value": int(analysis[PRIMARY_OUTCOME].sum()),
        },
        {
            "check": "five_factor_burden_lock",
            "status": observed_burden == EXPECTED_5F_BURDEN,
            "value": str(observed_burden),
        },
        {
            "check": "primary_1factor_or_script16",
            "status": abs(one_or-EXPECTED_5F_OR_1) <= EXPECTED_5F_OR_TOL,
            "value": one_or,
        },
        {
            "check": "primary_2plus_or_script16",
            "status": abs(twop_or-EXPECTED_5F_OR_2PLUS) <= EXPECTED_5F_OR_TOL,
            "value": twop_or,
        },
        {
            "check": "five_factor_multiplicative_interaction",
            "status": mult_global["p_value"] < 0.001,
            "value": mult_global["p_value"],
        },
        {
            "check": "five_factor_additive_interaction",
            "status": abs(add_global["p_value"]-EXPECTED_5F_ADD_P) <= 0.01,
            "value": add_global["p_value"],
        },
        {
            "check": "recent_pregnancy_n",
            "status": int(
                analysis_recent[PREGNANCY_RECENT].sum()
            ) == EXPECTED_PREG_RECENT_N,
            "value": int(
                analysis_recent[PREGNANCY_RECENT].sum()
            ),
        },
        {
            "check": "landmark_n",
            "status": len(landmark) == EXPECTED_LANDMARK_N,
            "value": len(landmark),
        },
        {
            "check": "landmark_vte_n",
            "status": int(
                landmark["tte_event"].sum()
            ) == EXPECTED_LANDMARK_VTE_N,
            "value": int(
                landmark["tte_event"].sum()
            ),
        },
        {
            "check": "participant_level_output_written",
            "status": True,
            "value": False,
        },
    ]

    qc = pd.DataFrame(qc_rows)
    qc_path = LOG_DIR / "17_qc.csv"
    qc.to_csv(qc_path, index=False)

    if not qc["status"].all():
        print("\nQC FAILURE:")
        print(
            qc.loc[
                ~qc["status"]
            ].to_string(index=False)
        )
        raise RuntimeError(
            "Script 17 QC failed. Review 17_qc.csv."
        )

    metadata = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script_version": SCRIPT_VERSION,
        "master_file": str(MASTER_FILE),
        "master_sha256": master_hash,
        "cdr_dataset": CDR_DATASET,
        "master_n": len(master),
        "wgs_n": len(analysis),
        "vte_n": int(analysis[PRIMARY_OUTCOME].sum()),
        "revised_primary_factors": FIVE_FACTORS,
        "pregnancy_in_primary_burden": False,
        "primary_endpoint_language": (
            "observed VTE occurrence during follow-up"
        ),
        "seven_marker_models_adjusted_for": (
            "exact five-factor count + age + recorded gender + "
            "16 ancestry PCs + EHR observability"
        ),
        "figure3_burden": "five-factor 0/1/2+",
        "oof_clinical_features": clinical_features,
        "time_to_event_rebuilt_with_five_factor_burden": True,
        "ancestry_models_rebuilt_with_five_factor_burden": True,
        "recent_pregnancy_sensitivity_run": True,
        "participant_level_output_written": False,
    }

    metadata_path = LOG_DIR / "17_metadata.json"
    metadata_path.write_text(
        json.dumps(metadata, indent=2) + "\n"
    )

    # ----------------------------------------------------------------------
    # Packaging.
    # ----------------------------------------------------------------------

    files = [
        table1_path,
        table2_path,
        table3_path,
        composite_path,
        figure3_path,
        figure3_tests_path,
        interaction_global_path,
        adjusted_risks_path,
        additive_contrasts_path,
        finer_path,
        count_path,
        exact_int_path,
        joint_path,
        no_prior_path,
        anticoag_path,
        performance_path,
        bootstrap_path,
        landmark_summary_path,
        overall_cox_path,
        interval_hr_path,
        time_tests_path,
        ancestry_path,
        heterogeneity_path,
        pregnancy_path,
        paste_ready_path,
        reviewer_note_path,
        methods_path,
        qc_path,
        metadata_path,
        DOC_DIR / "17_recent_pregnancy_365d_query.sql",
    ]

    package_name = "VTE_JTH_Script17_Five_Factor_Primary_Rebuild"
    build_root = PACKAGE_DIR / package_name

    if build_root.exists():
        shutil.rmtree(build_root)
    build_root.mkdir(parents=True, exist_ok=True)

    manifest_rows = []

    for source in files:
        source = Path(source)
        if not source.exists():
            continue

        relative = source.relative_to(OUTPUT_ROOT)
        destination = build_root / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)

        manifest_rows.append(
            {
                "relative_path": str(relative),
                "size_bytes": int(source.stat().st_size),
                "sha256": sha256_file(source),
            }
        )

    manifest = pd.DataFrame(
        manifest_rows
    ).sort_values("relative_path")
    manifest_path = build_root / "17_manifest.csv"
    manifest.to_csv(manifest_path, index=False)

    zip_path = (
        OUTPUT_ROOT.parent
        / f"VTE_primary_analysis_{datetime.now().strftime('%Y%m%d')}.zip"
    )

    if zip_path.exists():
        zip_path.unlink()

    with zipfile.ZipFile(
        zip_path,
        "w",
        compression=zipfile.ZIP_DEFLATED,
        compresslevel=9,
    ) as archive:
        for file in sorted(build_root.rglob("*")):
            if file.is_file():
                archive.write(
                    file,
                    arcname=str(
                        Path(package_name)
                        / file.relative_to(build_root)
                    ),
                )

    print("\n" + paste_ready)
    print("QC:", qc_path)
    print("Reviewer note:", reviewer_note_path)
    print("Methods:", methods_path)
    print("ZIP:", zip_path)
    print("=" * 84)


if __name__ == "__main__":
    main()
