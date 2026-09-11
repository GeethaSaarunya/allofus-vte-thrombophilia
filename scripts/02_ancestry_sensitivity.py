#!/usr/bin/env python3
"""
All of Us VTE / JTH
Public reproducibility script 02: Five-factor ancestry-stratified F5/F2 sensitivity
Version: public_v1.0.0

PURPOSE
-------
Run the final ancestry-stratified F5/F2 sensitivity analysis using the same
five-factor clinical burden as the primary manuscript analysis. The script
detects and normalizes the supported ancestry-label representations used in
the secure analysis master.

This script:
1. Loads the same frozen participant-level master used by Scripts 16/17.
2. Detects the actual ancestry column from known historical names.
3. Audits the raw ancestry labels in the WGS cohort.
4. Normalizes only explicit genetic-ancestry labels to canonical codes.
5. Re-fits F5 and F2 ancestry-stratified models using the revised exact
   FIVE-factor clinical burden.
6. Applies the same reportability rules:
      - ancestry stratum N >= 5,000
      - every variant x VTE 2x2 cell >= 20
7. Computes Cochran Q heterogeneity and Holm adjustment for the two tests.
8. Writes aggregate-safe outputs and a paste-ready summary.

SCIENTIFIC BOUNDARIES
---------------------
- This is a repair of Supplementary Table S8 only.
- Primary pooled models remain adjusted for all 16 ancestry PCs.
- Non-significant heterogeneity does NOT establish equivalent effects.
- Do not run ancestry x genotype x clinical-burden three-way interactions.
- No participant-level data are exported.
- Primary endpoint wording remains "observed VTE occurrence during follow-up."

INPUT
-----
Secure input:
    VTE_JTH_MASTER_FILE must point to the participant-level analysis master
    inside an authorized All of Us Researcher Workbench environment.

OUTPUT
------
Aggregate-safe outputs only. By default they are written under:
    ./outputs/ancestry_sensitivity/

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

from scipy.stats import chi2
from statsmodels.stats.multitest import multipletests


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
        "VTE_JTH_ANCESTRY_OUTPUT_DIR",
        "./outputs/ancestry_sensitivity",
    )
).expanduser().resolve()

TABLE_DIR = OUTPUT_ROOT / "tables"
LOG_DIR = OUTPUT_ROOT / "logs"
DOC_DIR = OUTPUT_ROOT / "documentation"
PACKAGE_DIR = OUTPUT_ROOT / "package_build"

for d in [TABLE_DIR, LOG_DIR, DOC_DIR, PACKAGE_DIR]:
    d.mkdir(parents=True, exist_ok=True)

PERSON_ID = "person_id"
WGS_ELIGIBLE = "analysis_eligible_primary_genetics"
PRIMARY_OUTCOME = "vte_primary_scientific"
OBSERVABILITY = "ehr_observability_score_0_6"

F5 = "f5_rs6025_carrier"
F2 = "f2_rs1799963_carrier"

FIVE_FACTORS = [
    "cancer_history",
    "major_surgery",
    "inpatient_hospitalization",
    "infection_or_sepsis",
    "fracture_or_major_trauma",
]

PC_COLUMNS = [f"genetic_pc{i}" for i in range(1, 17)]

ANCESTRY_COLUMN_CANDIDATES = [
    # Earlier successful VTE scripts used ancestry_pred.
    "ancestry_pred",
    # Script 17 attempted this field.
    "genetic_ancestry_label",
    # Other plausible historical master names.
    "genetic_ancestry",
    "ancestry_label",
    "predicted_ancestry",
    "ancestry",
]

REPORTABLE_GROUPS = ["EUR", "AFR", "AMR"]

MIN_ANCESTRY_N = 5_000
MIN_2X2_CELL = 20
MAX_GLM_ITER = 300

EXPECTED_MASTER_N = 483_707
EXPECTED_WGS_N = 358_533
EXPECTED_VTE_N = 7_553
EXPECTED_5F_BURDEN = {
    "0": 220_936,
    "1": 68_503,
    "2+": 69_094,
}

# Prior six-factor ancestry estimates retained only for comparison in output.
PRIOR_SIX_FACTOR_REFERENCE = {
    ("F5", "EUR"): (2.25, 2.03, 2.48),
    ("F5", "AFR"): (1.70, 1.13, 2.56),
    ("F5", "AMR"): (2.00, 1.38, 2.91),
    ("F2", "EUR"): (1.59, 1.36, 1.85),
    ("F2", "AFR"): (2.57, 1.58, 4.18),
    ("F2", "AMR"): (1.80, 1.29, 2.52),
}
PRIOR_HETEROGENEITY_REFERENCE = {
    "F5": {"Q_p": 0.381, "I2_percent": 0.0},
    "F2": {"Q_p": 0.160, "I2_percent": 45.5},
}


# ============================================================================
# 1. Utilities
# ============================================================================

def progress(message: str):
    print(
        f"[Script 18 | {datetime.now().strftime('%H:%M:%S')}] {message}",
        flush=True,
    )


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_table(path: Path) -> pd.DataFrame:
    if path.suffix.lower() == ".parquet":
        return pd.read_parquet(path)
    if path.suffix.lower() in {".tsv", ".txt"}:
        return pd.read_csv(path, sep="\t")
    return pd.read_csv(path)


def require_columns(frame: pd.DataFrame, columns: list[str]):
    missing = [c for c in columns if c not in frame.columns]
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


def detect_ancestry_column(frame: pd.DataFrame) -> str:
    for candidate in ANCESTRY_COLUMN_CANDIDATES:
        if candidate in frame.columns:
            return candidate

    ancestry_like = [
        c for c in frame.columns
        if "ancestry" in c.lower()
    ]

    if len(ancestry_like) == 1:
        return ancestry_like[0]

    raise KeyError(
        "Could not uniquely identify ancestry column.\n"
        f"Known candidates checked: {ANCESTRY_COLUMN_CANDIDATES}\n"
        f"Ancestry-like columns found: {ancestry_like}"
    )


def normalize_ancestry(value):
    """
    Normalize only explicit genetic-ancestry labels.

    We intentionally do not infer ancestry from race/ethnicity variables.
    """
    if pd.isna(value):
        return pd.NA

    text = str(value).strip().lower()

    direct = {
        "eur": "EUR",
        "afr": "AFR",
        "amr": "AMR",
        "eas": "EAS",
        "sas": "SAS",
        "mid": "MID",
    }
    if text in direct:
        return direct[text]

    # Explicit full genetic-ancestry descriptions.
    if "european" in text:
        return "EUR"
    if "african" in text:
        return "AFR"
    if (
        "admixed american" in text
        or "american admixed" in text
        or "admixture american" in text
    ):
        return "AMR"
    if "east asian" in text:
        return "EAS"
    if "south asian" in text:
        return "SAS"
    if "middle eastern" in text:
        return "MID"

    # Preserve unrecognized label for audit, but it will not be modeled.
    return f"UNMAPPED::{text}"


def prepare_model_frame(
    frame: pd.DataFrame,
    outcome: str,
    predictors: list[str],
) -> pd.DataFrame:
    model = (
        frame[[outcome, *predictors]]
        .apply(pd.to_numeric, errors="coerce")
        .dropna()
        .astype("float64")
        .copy()
    )

    if model.empty:
        raise RuntimeError("Empty model frame.")

    constants = [
        c for c in predictors
        if model[c].nunique(dropna=True) <= 1
    ]
    if constants:
        raise RuntimeError(
            "Constant predictors encountered: "
            + ", ".join(constants)
        )

    return model


def fit_glm(
    frame: pd.DataFrame,
    outcome: str,
    predictors: list[str],
):
    model = prepare_model_frame(
        frame,
        outcome,
        predictors,
    )

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
        raise RuntimeError("GLM failed to converge.")

    return result, model


def extract_marker_result(
    result,
    marker: str,
    marker_label: str,
    ancestry: str,
    ancestry_n: int,
    min_cell: int,
):
    beta = float(result.params[marker])
    se = float(result.bse[marker])

    return {
        "marker": marker_label,
        "marker_column": marker,
        "ancestry": ancestry,
        "ancestry_n": int(ancestry_n),
        "min_variant_vte_2x2_cell": int(min_cell),
        "beta": beta,
        "standard_error": se,
        "adjusted_or": math.exp(beta),
        "ci_95_low": math.exp(beta - 1.96*se),
        "ci_95_high": math.exp(beta + 1.96*se),
        "p_value": float(result.pvalues[marker]),
        "model_n": int(result.nobs),
    }


def cochran_q(frame: pd.DataFrame) -> dict:
    beta = frame["beta"].to_numpy(dtype=float)
    se = frame["standard_error"].to_numpy(dtype=float)

    w = 1.0 / np.square(se)
    pooled = float(np.sum(w*beta) / np.sum(w))

    Q = float(np.sum(w*np.square(beta-pooled)))
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
# 2. Main
# ============================================================================

def main():
    print("=" * 84)
    print("All of Us VTE / JTH")
    print("Public reproducibility script 02: Five-factor ancestry-stratified F5/F2 sensitivity")
    print("Version:", SCRIPT_VERSION)
    print("=" * 84)

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

    ancestry_column = detect_ancestry_column(master)
    progress(f"Detected ancestry column: {ancestry_column}")

    required = [
        PERSON_ID,
        WGS_ELIGIBLE,
        PRIMARY_OUTCOME,
        OBSERVABILITY,
        ancestry_column,
        "age_at_primary_consent",
        "female_sex",
        "other_unknown_sex",
        F5,
        F2,
        *FIVE_FACTORS,
        *PC_COLUMNS,
    ]
    require_columns(master, required)

    if master[PERSON_ID].duplicated().any():
        raise RuntimeError(
            "Frozen master contains duplicate person_id."
        )

    for column in [
        WGS_ELIGIBLE,
        PRIMARY_OUTCOME,
        "female_sex",
        "other_unknown_sex",
        F5,
        F2,
        *FIVE_FACTORS,
    ]:
        master[column] = numeric_binary(
            master[column],
            column,
        )

    master["age_per_10y"] = (
        pd.to_numeric(
            master["age_at_primary_consent"],
            errors="coerce",
        ) / 10.0
    )

    master["five_factor_count"] = (
        master[FIVE_FACTORS]
        .astype(float)
        .sum(axis=1)
    ).astype("int8")

    master["burden_012plus_5f"] = np.select(
        [
            master["five_factor_count"].eq(0),
            master["five_factor_count"].eq(1),
        ],
        ["0","1"],
        default="2+",
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
            "VTE count no longer matches manuscript lock."
        )

    burden_counts = (
        analysis["burden_012plus_5f"]
        .value_counts()
        .to_dict()
    )
    burden_counts = {
        k: int(burden_counts.get(k,0))
        for k in EXPECTED_5F_BURDEN
    }

    if burden_counts != EXPECTED_5F_BURDEN:
        raise RuntimeError(
            "Five-factor burden does not reproduce Script 17.\n"
            f"Observed {burden_counts}\n"
            f"Expected {EXPECTED_5F_BURDEN}"
        )

    # ----------------------------------------------------------------------
    # Raw ancestry-label audit.
    # ----------------------------------------------------------------------

    progress("Auditing raw and normalized ancestry labels.")

    raw_counts = (
        analysis[ancestry_column]
        .astype("string")
        .fillna("<MISSING>")
        .value_counts(dropna=False)
        .rename_axis("raw_ancestry_label")
        .reset_index(name="n")
    )
    raw_counts["percent_wgs"] = (
        100.0 * raw_counts["n"] / len(analysis)
    )

    raw_counts_path = (
        TABLE_DIR / "18A_raw_ancestry_label_counts.csv"
    )
    raw_counts.to_csv(raw_counts_path, index=False)

    analysis["ancestry_canonical"] = (
        analysis[ancestry_column]
        .map(normalize_ancestry)
        .astype("string")
    )

    normalized_counts = (
        analysis["ancestry_canonical"]
        .fillna("<MISSING>")
        .value_counts(dropna=False)
        .rename_axis("ancestry_canonical")
        .reset_index(name="n")
    )
    normalized_counts["percent_wgs"] = (
        100.0 * normalized_counts["n"] / len(analysis)
    )

    normalized_counts_path = (
        TABLE_DIR / "18A_normalized_ancestry_counts.csv"
    )
    normalized_counts.to_csv(
        normalized_counts_path,
        index=False,
    )

    print("\nNormalized ancestry counts:")
    print(
        normalized_counts.head(20).to_string(index=False)
    )

    # ----------------------------------------------------------------------
    # Eligibility and models.
    # ----------------------------------------------------------------------

    progress("Evaluating ancestry reportability and fitting F5/F2 models.")

    base_covariates = [
        "age_per_10y",
        "female_sex",
        "other_unknown_sex",
        *PC_COLUMNS,
        OBSERVABILITY,
    ]

    eligibility_rows = []
    result_rows = []

    for marker, marker_label in [
        (F5, "Factor V Leiden"),
        (F2, "Prothrombin G20210A"),
    ]:
        for ancestry in REPORTABLE_GROUPS:
            sub = analysis.loc[
                analysis["ancestry_canonical"].eq(ancestry)
            ].copy()

            ancestry_n = len(sub)

            if ancestry_n == 0:
                eligibility_rows.append(
                    {
                        "marker": marker_label,
                        "ancestry": ancestry,
                        "ancestry_n": 0,
                        "min_variant_vte_2x2_cell": 0,
                        "eligible": False,
                        "reason": "no participants after ancestry-label normalization",
                    }
                )
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

            eligible = (
                ancestry_n >= MIN_ANCESTRY_N
                and min_cell >= MIN_2X2_CELL
            )

            reasons = []
            if ancestry_n < MIN_ANCESTRY_N:
                reasons.append(
                    f"N<{MIN_ANCESTRY_N}"
                )
            if min_cell < MIN_2X2_CELL:
                reasons.append(
                    f"minimum 2x2 cell<{MIN_2X2_CELL}"
                )

            eligibility_rows.append(
                {
                    "marker": marker_label,
                    "ancestry": ancestry,
                    "ancestry_n": ancestry_n,
                    "min_variant_vte_2x2_cell": min_cell,
                    "eligible": eligible,
                    "reason": (
                        "reportable"
                        if eligible
                        else "; ".join(reasons)
                    ),
                    "noncarrier_no_vte": int(cross.loc[0,0]),
                    "noncarrier_vte": int(cross.loc[0,1]),
                    "carrier_no_vte": int(cross.loc[1,0]),
                    "carrier_vte": int(cross.loc[1,1]),
                }
            )

            if not eligible:
                continue

            progress(
                f"Fitting {ancestry} | {marker_label} | N={ancestry_n:,}"
            )

            predictors = [
                marker,
                "five_factor_count",
                *base_covariates,
            ]

            result, model_frame = fit_glm(
                sub,
                PRIMARY_OUTCOME,
                predictors,
            )

            result_rows.append(
                extract_marker_result(
                    result,
                    marker,
                    marker_label,
                    ancestry,
                    ancestry_n,
                    min_cell,
                )
            )

    eligibility = pd.DataFrame(eligibility_rows)
    results = pd.DataFrame(result_rows)

    eligibility_path = (
        TABLE_DIR / "18B_ancestry_model_eligibility.csv"
    )
    eligibility.to_csv(
        eligibility_path,
        index=False,
    )

    results_path = (
        TABLE_DIR / "18C_ancestry_stratified_F5_F2_five_factor.csv"
    )
    results.to_csv(
        results_path,
        index=False,
    )

    # Guardrail: the historically reportable groups should be restored.
    for marker_label in [
        "Factor V Leiden",
        "Prothrombin G20210A",
    ]:
        observed_groups = set(
            results.loc[
                results["marker"].eq(marker_label),
                "ancestry",
            ].tolist()
        )
        expected_groups = set(REPORTABLE_GROUPS)

        if observed_groups != expected_groups:
            raise RuntimeError(
                f"Reportable ancestry groups not restored for {marker_label}.\n"
                f"Observed: {sorted(observed_groups)}\n"
                f"Expected: {sorted(expected_groups)}\n"
                "Review 18A/18B ancestry audit outputs."
            )

    # ----------------------------------------------------------------------
    # Heterogeneity.
    # ----------------------------------------------------------------------

    progress("Computing Cochran Q heterogeneity.")

    heterogeneity_rows = []

    for marker_label in [
        "Factor V Leiden",
        "Prothrombin G20210A",
    ]:
        sub = results.loc[
            results["marker"].eq(marker_label)
        ].copy()

        het = cochran_q(sub)

        heterogeneity_rows.append(
            {
                "marker": marker_label,
                **het,
            }
        )

    heterogeneity = pd.DataFrame(
        heterogeneity_rows
    )

    heterogeneity["holm_adjusted_p"] = multipletests(
        heterogeneity["p_value"].to_numpy(dtype=float),
        method="holm",
    )[1]

    heterogeneity_path = (
        TABLE_DIR / "18D_ancestry_heterogeneity_five_factor.csv"
    )
    heterogeneity.to_csv(
        heterogeneity_path,
        index=False,
    )

    # ----------------------------------------------------------------------
    # Comparison with prior six-factor ancestry estimates.
    # ----------------------------------------------------------------------

    comparison_rows = []

    marker_short = {
        "Factor V Leiden": "F5",
        "Prothrombin G20210A": "F2",
    }

    for _, row in results.iterrows():
        short = marker_short[row["marker"]]
        old = PRIOR_SIX_FACTOR_REFERENCE.get(
            (short, row["ancestry"])
        )

        comparison_rows.append(
            {
                "marker": row["marker"],
                "ancestry": row["ancestry"],
                "old_six_factor_or": old[0] if old else np.nan,
                "old_six_factor_ci_low": old[1] if old else np.nan,
                "old_six_factor_ci_high": old[2] if old else np.nan,
                "revised_five_factor_or": row["adjusted_or"],
                "revised_five_factor_ci_low": row["ci_95_low"],
                "revised_five_factor_ci_high": row["ci_95_high"],
                "absolute_or_change": (
                    row["adjusted_or"] - old[0]
                    if old else np.nan
                ),
            }
        )

    comparison = pd.DataFrame(comparison_rows)

    comparison_path = (
        TABLE_DIR / "18E_old_vs_revised_ancestry_estimates.csv"
    )
    comparison.to_csv(
        comparison_path,
        index=False,
    )

    # ----------------------------------------------------------------------
    # Paste-ready result.
    # ----------------------------------------------------------------------

    progress("Writing paste-ready ancestry summary.")

    def row_for(marker_label, ancestry):
        return results.loc[
            results["marker"].eq(marker_label)
            & results["ancestry"].eq(ancestry)
        ].iloc[0]

    f5_het = heterogeneity.loc[
        heterogeneity["marker"].eq("Factor V Leiden")
    ].iloc[0]
    f2_het = heterogeneity.loc[
        heterogeneity["marker"].eq("Prothrombin G20210A")
    ].iloc[0]

    lines = [
        "=" * 84,
        "VTE JTH SCRIPT 18 — FIVE-FACTOR ANCESTRY REPAIR",
        "=" * 84,
        f"Detected ancestry column: {ancestry_column}",
        f"WGS N: {len(analysis):,}",
        f"Observed VTE during follow-up: {int(analysis[PRIMARY_OUTCOME].sum()):,}",
        "",
        "FACTOR V LEIDEN",
    ]

    for ancestry in REPORTABLE_GROUPS:
        row = row_for("Factor V Leiden", ancestry)
        lines.append(
            f"{ancestry}: OR {row['adjusted_or']:.3f} "
            f"(95% CI {row['ci_95_low']:.3f}-{row['ci_95_high']:.3f}); "
            f"P {('<0.001' if row['p_value'] < 0.001 else f'{row['p_value']:.3f}')}"
        )

    lines += [
        f"Cochran Q P={f5_het['p_value']:.3f}; "
        f"Holm P={f5_het['holm_adjusted_p']:.3f}; "
        f"I2={f5_het['I2_percent']:.1f}%",
        "",
        "PROTHROMBIN G20210A",
    ]

    for ancestry in REPORTABLE_GROUPS:
        row = row_for("Prothrombin G20210A", ancestry)
        lines.append(
            f"{ancestry}: OR {row['adjusted_or']:.3f} "
            f"(95% CI {row['ci_95_low']:.3f}-{row['ci_95_high']:.3f}); "
            f"P {('<0.001' if row['p_value'] < 0.001 else f'{row['p_value']:.3f}')}"
        )

    lines += [
        f"Cochran Q P={f2_het['p_value']:.3f}; "
        f"Holm P={f2_het['holm_adjusted_p']:.3f}; "
        f"I2={f2_het['I2_percent']:.1f}%",
        "",
        "INTERPRETATION",
        "- Primary pooled models remain adjusted for all 16 ancestry PCs.",
        "- EUR, AFR, and AMR meet the prespecified reportability thresholds.",
        "- Non-significant heterogeneity means no statistical evidence of "
        "between-group heterogeneity; it does NOT establish equivalent effects.",
        "- These models use the revised exact FIVE-factor burden.",
        "- No three-way ancestry x genotype x clinical-burden interaction is required.",
        "=" * 84,
    ]

    paste_ready = "\n".join(lines) + "\n"

    paste_path = OUTPUT_ROOT / "18_PASTE_READY_RESULTS.txt"
    paste_path.write_text(paste_ready)

    # ----------------------------------------------------------------------
    # Reviewer/manuscript language.
    # ----------------------------------------------------------------------

    manuscript_note = f"""
# Revised ancestry-stratified manuscript language

## Methods

Population structure in the primary pooled genetic models was controlled using
16 ancestry principal components. Supporting ancestry-stratified analyses were
performed for Factor V Leiden and prothrombin G20210A in ancestry groups with at
least {MIN_ANCESTRY_N:,} participants and at least {MIN_2X2_CELL} participants
in every variant-by-VTE 2x2 cell. Each ancestry-specific model was adjusted for
the exact revised five-factor clinical-risk burden, age, recorded sex, 16
ancestry principal components, and EHR observability. Between-group
heterogeneity was assessed using Cochran's Q, with Holm adjustment across the
two variant-level heterogeneity tests.

## Results

[Use the values in 18_PASTE_READY_RESULTS.txt.]

## Interpretation boundary

A non-significant heterogeneity test indicates that the available data do not
provide statistical evidence of between-ancestry heterogeneity in the
reportable strata. It does not establish equivalent or identical genetic
effects across ancestry groups.
""".strip() + "\n"

    manuscript_note_path = (
        DOC_DIR / "18_ancestry_methods_results_language.md"
    )
    manuscript_note_path.write_text(manuscript_note)

    # ----------------------------------------------------------------------
    # QC and metadata.
    # ----------------------------------------------------------------------

    qc = pd.DataFrame(
        [
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
                "status": burden_counts == EXPECTED_5F_BURDEN,
                "value": str(burden_counts),
            },
            {
                "check": "ancestry_column_detected",
                "status": bool(ancestry_column),
                "value": ancestry_column,
            },
            {
                "check": "F5_reportable_groups",
                "status": set(
                    results.loc[
                        results["marker"].eq("Factor V Leiden"),
                        "ancestry",
                    ]
                ) == set(REPORTABLE_GROUPS),
                "value": ",".join(
                    results.loc[
                        results["marker"].eq("Factor V Leiden"),
                        "ancestry",
                    ].tolist()
                ),
            },
            {
                "check": "F2_reportable_groups",
                "status": set(
                    results.loc[
                        results["marker"].eq("Prothrombin G20210A"),
                        "ancestry",
                    ]
                ) == set(REPORTABLE_GROUPS),
                "value": ",".join(
                    results.loc[
                        results["marker"].eq("Prothrombin G20210A"),
                        "ancestry",
                    ].tolist()
                ),
            },
            {
                "check": "participant_level_output_written",
                "status": True,
                "value": False,
            },
        ]
    )

    qc_path = LOG_DIR / "18_qc.csv"
    qc.to_csv(qc_path, index=False)

    if not qc["status"].all():
        raise RuntimeError(
            "Script 18 QC failed:\n"
            + qc.loc[~qc["status"]].to_string(index=False)
        )

    metadata = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script_version": SCRIPT_VERSION,
        "master_file": str(MASTER_FILE),
        "master_sha256": master_hash,
        "ancestry_column_detected": ancestry_column,
        "wgs_n": len(analysis),
        "vte_n": int(analysis[PRIMARY_OUTCOME].sum()),
        "clinical_adjustment": "exact five-factor burden",
        "population_structure_adjustment": "16 ancestry PCs",
        "reportability_min_group_n": MIN_ANCESTRY_N,
        "reportability_min_2x2_cell": MIN_2X2_CELL,
        "reportable_groups_expected": REPORTABLE_GROUPS,
        "three_way_interaction_run": False,
        "participant_level_output_written": False,
    }

    metadata_path = LOG_DIR / "18_metadata.json"
    metadata_path.write_text(
        json.dumps(metadata, indent=2) + "\n"
    )

    # ----------------------------------------------------------------------
    # Package.
    # ----------------------------------------------------------------------

    files = [
        raw_counts_path,
        normalized_counts_path,
        eligibility_path,
        results_path,
        heterogeneity_path,
        comparison_path,
        paste_path,
        manuscript_note_path,
        qc_path,
        metadata_path,
    ]

    package_name = "VTE_JTH_Script18_Five_Factor_Ancestry_Repair"
    build_root = PACKAGE_DIR / package_name

    if build_root.exists():
        shutil.rmtree(build_root)
    build_root.mkdir(parents=True, exist_ok=True)

    manifest_rows = []

    for source in files:
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

    manifest_path = build_root / "18_manifest.csv"
    manifest.to_csv(manifest_path, index=False)

    zip_path = (
        OUTPUT_ROOT.parent
        / f"VTE_ancestry_sensitivity_{datetime.now().strftime('%Y%m%d')}.zip"
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
    print("Raw ancestry audit:", raw_counts_path)
    print("Eligibility:", eligibility_path)
    print("Results:", results_path)
    print("Heterogeneity:", heterogeneity_path)
    print("ZIP:", zip_path)
    print("=" * 84)


if __name__ == "__main__":
    main()
