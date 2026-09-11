#!/usr/bin/env python3
"""
All of Us VTE / JTH
Public reproducibility script 03: Export VTE phenotype concept specification
Version: 15_v1_3_20260819

PURPOSE
-------
Recover the exact vocabulary-level concept disposition used by the frozen
Script 05B VTE endpoint and export manuscript-ready concept IDs/names.

IMPORTANT CORRECTION IN v1.3
----------------------------
The original Script 05B did NOT require every exclusion category to contain at
least one concept. Its QC required only:
    - primary retained concept set nonempty
    - strict acute/current concept set nonempty
    - primary DVT concept set nonempty
    - primary PE concept set nonempty

Therefore, EXCLUDE_SUPERFICIAL and EXCLUDE_UNUSUAL_SITE are allowed to contain
zero concepts in the reviewed candidate set. If a manuscript category has zero
concepts, this script reports:
    Concept IDs: None in reviewed Script 05B candidate set
rather than inventing additional concepts that were never used by the phenotype.

The candidate-concept query and concept classifier below reproduce the original
Script 05B V2 workflow.

This script does NOT rebuild participant outcomes and does NOT export participant
counts, person IDs, or dates.
"""

from __future__ import annotations

from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import math
import os
import re
import shutil
import zipfile

import pandas as pd


SCRIPT_VERSION = "15_v1_3_20260819"

LOCKED_CDR_DATASET = os.environ.get(
    "VTE_LOCKED_CDR_DATASET",
    "wb-silky-artichoke-2408.C2024Q3R9",
).strip("`")

OUTPUT_ROOT = Path(
    os.environ.get(
        "VTE_SCRIPT15_OUTPUT",
        "./outputs/vte_concept_specification",
    )
).expanduser()
OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)

SEARCH_ROOTS = [
    Path.cwd(),
    Path(os.environ.get("VTE_CONCEPT_SOURCE_DIR", ".")).expanduser(),
]

EXCLUDED_PATH_TOKENS = [
    "VTE_CDRv9_omics_extension",
    "cdrv9_omics_extension",
]

PRIMARY_DISPOSITIONS = {
    "RETAIN_PRIMARY_ACUTE_CURRENT",
    "RETAIN_PRIMARY_RECURRENT",
    "RETAIN_PRIMARY_PREGNANCY",
    "RETAIN_PRIMARY_POSTOPERATIVE",
    "RETAIN_PRIMARY_UNSPECIFIED",
}

STRICT_DISPOSITIONS = {
    "RETAIN_PRIMARY_ACUTE_CURRENT",
    "RETAIN_PRIMARY_RECURRENT",
    "RETAIN_PRIMARY_PREGNANCY",
    "RETAIN_PRIMARY_POSTOPERATIVE",
}

CATEGORY_RULES = [
    (1, "Acute/current DVT", "Included", {"RETAIN_PRIMARY_ACUTE_CURRENT"}, {"DVT"},
     "Retained as current primary VTE."),
    (2, "Acute/current PE", "Included", {"RETAIN_PRIMARY_ACUTE_CURRENT"}, {"PE"},
     "Retained as current primary VTE."),
    (3, "Recurrent DVT/PE", "Included", {"RETAIN_PRIMARY_RECURRENT"}, {"DVT", "PE"},
     "Retained in the primary phenotype."),
    (4, "Postoperative DVT/PE", "Included", {"RETAIN_PRIMARY_POSTOPERATIVE"}, {"DVT", "PE"},
     "Retained in the primary phenotype."),
    (5, "Pregnancy-associated DVT/PE", "Included", {"RETAIN_PRIMARY_PREGNANCY"}, {"DVT", "PE"},
     "Retained in the primary phenotype."),
    (6, "Unspecified-current DVT/PE", "Included", {"RETAIN_PRIMARY_UNSPECIFIED"}, {"DVT", "PE"},
     "Retained after concept review."),
    (7, "History/screening/suspected/rule-out", "Excluded", {"EXCLUDE_NONCASE"}, None,
     "Does not represent a qualifying current VTE record."),
    (8, "Chronic/post-thrombotic VTE", "Excluded", {"EXCLUDE_CHRONIC"}, None,
     "Outside the primary current VTE phenotype."),
    (9, "Superficial thrombosis", "Excluded", {"EXCLUDE_SUPERFICIAL"}, None,
     "Outside the DVT/PE endpoint."),
    (10, "Unusual-site thrombosis", "Excluded", {"EXCLUDE_UNUSUAL_SITE"}, None,
     "Characterized separately and excluded from primary VTE."),
]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def path_is_excluded(path: Path) -> bool:
    text = str(path)
    return any(token.lower() in text.lower() for token in EXCLUDED_PATH_TOKENS)


def normalize_frame(df: pd.DataFrame) -> pd.DataFrame:
    required = ["concept_id", "concept_name", "final_disposition", "final_subtype"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise KeyError("missing columns: " + ", ".join(missing))

    out = df.copy()
    out["concept_id"] = pd.to_numeric(out["concept_id"], errors="raise").astype("int64")
    out["concept_name"] = out["concept_name"].astype("string").fillna("").str.strip()
    out["final_disposition"] = (
        out["final_disposition"].astype("string").fillna("").str.strip().str.upper()
    )
    out["final_subtype"] = (
        out["final_subtype"].astype("string").fillna("").str.strip().str.upper()
    )

    for c in ["final_acuity", "final_context", "final_site_class", "final_rationale"]:
        if c not in out.columns:
            out[c] = ""
        else:
            out[c] = out[c].astype("string").fillna("").str.strip()

    return out


def original_05B_qc_contract(df: pd.DataFrame) -> tuple[bool, str, dict]:
    """
    Match the actual QC contract in the original Script 05B.

    The source script did NOT require chronic, superficial, unusual-site, or
    noncase exclusion sets to be nonempty.
    """
    primary = df["final_disposition"].isin(PRIMARY_DISPOSITIONS)
    strict = df["final_disposition"].isin(STRICT_DISPOSITIONS)

    counts = {
        "primary": int(primary.sum()),
        "strict_acute_current": int(strict.sum()),
        "dvt_primary": int((primary & df["final_subtype"].eq("DVT")).sum()),
        "pe_primary": int((primary & df["final_subtype"].eq("PE")).sum()),
        "chronic_excluded": int(df["final_disposition"].eq("EXCLUDE_CHRONIC").sum()),
        "unusual_site_excluded": int(df["final_disposition"].eq("EXCLUDE_UNUSUAL_SITE").sum()),
        "superficial_excluded": int(df["final_disposition"].eq("EXCLUDE_SUPERFICIAL").sum()),
        "noncase_excluded": int(df["final_disposition"].eq("EXCLUDE_NONCASE").sum()),
        "pregnancy_associated": int(
            (
                primary
                & df["final_context"].astype(str).eq("PREGNANCY_ASSOCIATED")
            ).sum()
        ),
        "postoperative": int(
            (
                primary
                & df["final_context"].astype(str).eq("POSTOPERATIVE")
            ).sum()
        ),
        "recurrent": int(
            (
                primary
                & df["final_acuity"].astype(str).eq("RECURRENT_CURRENT")
            ).sum()
        ),
        "retained_unspecified": int(
            df["final_disposition"].eq("RETAIN_PRIMARY_UNSPECIFIED").sum()
        ),
    }

    if counts["primary"] == 0:
        return False, "primary curated VTE concept set is empty", counts
    if counts["strict_acute_current"] == 0:
        return False, "strict acute/current VTE concept set is empty", counts
    if counts["dvt_primary"] == 0 or counts["pe_primary"] == 0:
        return False, "curated DVT or PE primary concept set is empty", counts

    return True, "PASS", counts


def discover_candidate_files() -> pd.DataFrame:
    rows = []
    seen = set()

    explicit = os.environ.get("VTE_CONCEPT_DISPOSITION_FILE", "").strip()
    if explicit:
        p = Path(explicit).expanduser()
        rows.append(
            {
                "path": str(p),
                "exists": p.exists(),
                "source": "explicit_override",
                "priority": 0,
                "excluded_path": path_is_excluded(p),
            }
        )
        seen.add(str(p))

    known = []
    for root in SEARCH_ROOTS:
        known.extend([
            root / "05B_vte_concept_disposition_internal.csv",
            root / "05B_vte_concept_disposition.csv",
            root / "internal" / "qc" / "05B_vte_concept_disposition_internal.csv",
            root / "outputs" / "aggregate_safe" / "05B_vte_concept_disposition.csv",
        ])

    for p in known:
        if str(p) in seen:
            continue
        rows.append(
            {
                "path": str(p),
                "exists": p.exists(),
                "source": "known_path",
                "priority": 1 if "_internal" in p.name else 2,
                "excluded_path": path_is_excluded(p),
            }
        )
        seen.add(str(p))

    skip_dirs = {
        ".git", ".cache", ".local", "__pycache__", "site-packages",
        "node_modules", "lost+found",
    }

    for root in SEARCH_ROOTS:
        if not root.exists():
            continue
        for current, dirs, files in os.walk(root):
            dirs[:] = [d for d in dirs if d not in skip_dirs]
            for filename in files:
                lower = filename.lower()
                if not (
                    lower.endswith(".csv")
                    and "05b" in lower
                    and "vte" in lower
                    and "concept" in lower
                    and "disposition" in lower
                ):
                    continue

                p = Path(current) / filename
                if str(p) in seen:
                    continue

                rows.append(
                    {
                        "path": str(p),
                        "exists": True,
                        "source": "auto_discovery",
                        "priority": 1 if "_internal" in filename else 2,
                        "excluded_path": path_is_excluded(p),
                    }
                )
                seen.add(str(p))

    return pd.DataFrame(rows)


def choose_valid_frozen_file(discovery: pd.DataFrame):
    checked = discovery.copy()
    checked["schema_valid"] = False
    checked["original_05B_qc_valid"] = False
    checked["validation"] = ""

    for idx, row in checked.sort_values(["priority", "path"]).iterrows():
        p = Path(row["path"])

        if not bool(row["exists"]):
            checked.at[idx, "validation"] = "missing"
            continue

        if bool(row["excluded_path"]):
            checked.at[idx, "validation"] = "rejected unrelated project path"
            continue

        try:
            raw = pd.read_csv(p, low_memory=False)
            norm = normalize_frame(raw)
            checked.at[idx, "schema_valid"] = True
            valid, message, counts = original_05B_qc_contract(norm)
            checked.at[idx, "original_05B_qc_valid"] = valid
            checked.at[idx, "validation"] = (
                message + " | " + json.dumps(counts, sort_keys=True)
            )
        except Exception as exc:
            checked.at[idx, "validation"] = f"read/validation error: {exc}"
            continue

        if valid:
            return p, norm, checked, counts

    return None, None, checked, None


def contains(text: str, pattern: str) -> bool:
    return bool(re.search(pattern, str(text or ""), flags=re.IGNORECASE))


def classify_vte_concept(name: str) -> dict:
    """
    Exact classification logic from original 05B_VTE_V2_curate_endpoint_concepts.py.
    """
    text = str(name or "")

    if contains(text, r"pulmonary embol|pulmonary thromboembol"):
        subtype = "PE"
    elif contains(
        text,
        r"deep vein thromb|deep venous thromb|deep vein embol|"
        r"deep venous embol|phlebothrombosis|thrombophlebitis.*deep",
    ):
        subtype = "DVT"
    else:
        return {
            "final_disposition": "EXCLUDE_NONCASE",
            "final_subtype": "OTHER",
            "final_acuity": "NONCASE",
            "final_context": "NONE",
            "final_site_class": "NONQUALIFYING",
            "final_rationale": "Not an explicit DVT or PE concept.",
        }

    noncase_pattern = (
        r"family history|personal history|history of|screening|"
        r"suspected|rule out|risk of|prophylaxis|encounter for|"
        r"observation for|follow-up|follow up|aftercare"
    )
    superficial_pattern = (
        r"superficial|thrombophlebitis(?!.*deep)|superficial phlebitis"
    )
    unusual_site_pattern = (
        r"cerebral|intracranial|retinal|portal|mesenteric|splanchnic|"
        r"hepatic vein|splenic vein|renal vein|ovarian vein|penile vein|"
        r"vena cava|caval thromb|umbilical vein|placental vein|"
        r"dural sinus|cavernous sinus"
    )
    chronic_pattern = (
        r"chronic|post[- ]thrombotic|postphlebitic|"
        r"post-phlebitic|old thromb|residual thromb|sequela"
    )
    pregnancy_pattern = (
        r"pregnan|antepartum|postpartum|post-partum|puerper|obstetric"
    )
    postoperative_pattern = (
        r"postoperative|post-operative|postprocedural|post-procedural|"
        r"following surgery|following procedure|due to procedure"
    )
    recurrent_pattern = r"recurrent"
    acute_current_pattern = r"acute|current|new onset|newly diagnosed"

    if contains(text, noncase_pattern):
        return {
            "final_disposition": "EXCLUDE_NONCASE",
            "final_subtype": subtype,
            "final_acuity": "NONCASE",
            "final_context": "NONE",
            "final_site_class": "NONQUALIFYING",
            "final_rationale": (
                "History, screening, suspected/rule-out, risk, or aftercare terminology."
            ),
        }

    if contains(text, superficial_pattern):
        return {
            "final_disposition": "EXCLUDE_SUPERFICIAL",
            "final_subtype": subtype,
            "final_acuity": "NONCASE",
            "final_context": "NONE",
            "final_site_class": "NONQUALIFYING",
            "final_rationale": (
                "Superficial thrombosis/thrombophlebitis is outside the DVT/PE endpoint."
            ),
        }

    if contains(text, unusual_site_pattern):
        return {
            "final_disposition": "EXCLUDE_UNUSUAL_SITE",
            "final_subtype": subtype,
            "final_acuity": "CURRENT_OR_UNSPECIFIED",
            "final_context": "NONE",
            "final_site_class": "UNUSUAL_SITE",
            "final_rationale": (
                "Unusual-site thrombosis is characterized separately and excluded "
                "from primary VTE."
            ),
        }

    if contains(text, chronic_pattern):
        context = "NONE"
        if contains(text, pregnancy_pattern):
            context = "PREGNANCY_ASSOCIATED"
        elif contains(text, postoperative_pattern):
            context = "POSTOPERATIVE"

        return {
            "final_disposition": "EXCLUDE_CHRONIC",
            "final_subtype": subtype,
            "final_acuity": "CHRONIC",
            "final_context": context,
            "final_site_class": "TYPICAL_OR_UNSPECIFIED",
            "final_rationale": (
                "Chronic, residual, or post-thrombotic terminology is not a "
                "new/current primary VTE event."
            ),
        }

    context = "NONE"
    if contains(text, pregnancy_pattern):
        context = "PREGNANCY_ASSOCIATED"
    elif contains(text, postoperative_pattern):
        context = "POSTOPERATIVE"

    if contains(text, recurrent_pattern):
        acuity = "RECURRENT_CURRENT"
        disposition = "RETAIN_PRIMARY_RECURRENT"
        rationale = (
            "Recurrent current DVT/PE retained in primary and strict-current definitions."
        )
    elif contains(text, acute_current_pattern):
        acuity = "ACUTE_CURRENT"
        disposition = "RETAIN_PRIMARY_ACUTE_CURRENT"
        rationale = (
            "Explicit acute/current DVT/PE retained in primary and strict-current definitions."
        )
    elif context == "PREGNANCY_ASSOCIATED":
        acuity = "CONTEXT_CURRENT"
        disposition = "RETAIN_PRIMARY_PREGNANCY"
        rationale = (
            "Pregnancy-associated DVT/PE retained as a current qualifying VTE context."
        )
    elif context == "POSTOPERATIVE":
        acuity = "CONTEXT_CURRENT"
        disposition = "RETAIN_PRIMARY_POSTOPERATIVE"
        rationale = (
            "Postoperative/postprocedural DVT/PE retained as a current qualifying VTE context."
        )
    else:
        acuity = "UNSPECIFIED_CURRENT"
        disposition = "RETAIN_PRIMARY_UNSPECIFIED"
        rationale = (
            "Typical or unspecified-site DVT/PE without chronic/history/noncase terminology; "
            "retained in the primary definition and excluded from the strict explicit-current sensitivity."
        )

    return {
        "final_disposition": disposition,
        "final_subtype": subtype,
        "final_acuity": acuity,
        "final_context": context,
        "final_site_class": "TYPICAL_OR_UNSPECIFIED",
        "final_rationale": rationale,
    }


def find_manual_override():
    target = "05B_vte_manual_concept_overrides.csv"
    candidates = []

    for root in SEARCH_ROOTS:
        if not root.exists():
            continue
        for current, dirs, files in os.walk(root):
            if target in files:
                p = Path(current) / target
                if not path_is_excluded(p):
                    candidates.append(p)

    candidates = sorted(
        set(candidates),
        key=lambda p: (
            0 if "VTE_AllOfUs_CDRv8_C2024Q3R9_V2" in str(p) else 1,
            str(p),
        ),
    )
    return candidates[0] if candidates else None


def find_metadata():
    target = "05B_curated_vte_endpoint_metadata.json"
    candidates = []

    for root in SEARCH_ROOTS:
        if not root.exists():
            continue
        for current, dirs, files in os.walk(root):
            if target in files:
                p = Path(current) / target
                if not path_is_excluded(p):
                    candidates.append(p)

    candidates = sorted(
        set(candidates),
        key=lambda p: (
            0 if "VTE_AllOfUs_CDRv8_C2024Q3R9_V2" in str(p) else 1,
            str(p),
        ),
    )
    return candidates[0] if candidates else None


def reconstruct_from_locked_cdr() -> tuple[pd.DataFrame, dict, dict]:
    from google.cloud import bigquery

    client = bigquery.Client()

    print("\nNo valid frozen Script 05B concept-disposition file is mounted.")
    print("Reconstructing the exact vocabulary-level Script 05B candidate set from:")
    print(LOCKED_CDR_DATASET)

    # Exact candidate terminology from original Script 05B.
    candidate_vte_pattern = (
        r"pulmonary embol|pulmonary thromboembol|"
        r"deep vein thromb|deep venous thromb|"
        r"deep vein embol|deep venous embol|"
        r"phlebothrombosis|thrombophlebitis.*deep"
    )

    sql = f"""
WITH observation_summary AS (
    SELECT
        person_id,
        MIN(observation_period_start_date) AS earliest_ehr_start_date,
        MAX(observation_period_end_date) AS latest_ehr_end_date
    FROM `{LOCKED_CDR_DATASET}.observation_period`
    GROUP BY person_id
),
primary_consent AS (
    SELECT
        o.person_id,
        MIN(o.observation_date) AS index_date
    FROM `{LOCKED_CDR_DATASET}.concept` AS c
    INNER JOIN `{LOCKED_CDR_DATASET}.concept_ancestor` AS ca
        ON c.concept_id = ca.ancestor_concept_id
    INNER JOIN `{LOCKED_CDR_DATASET}.observation` AS o
        ON ca.descendant_concept_id = o.observation_concept_id
    WHERE c.concept_name = 'Consent PII'
      AND c.concept_class_id = 'Module'
    GROUP BY o.person_id
),
ehr_consent AS (
    SELECT
        person_id,
        MIN(observation_date) AS ehr_consent_date
    FROM `{LOCKED_CDR_DATASET}.observation`
    WHERE observation_source_concept_id = 1586099
      AND value_source_concept_id = 1586100
    GROUP BY person_id
),
consent_spanning_observation AS (
    SELECT
        pc.person_id,
        pc.index_date
    FROM primary_consent AS pc
    INNER JOIN `{LOCKED_CDR_DATASET}.observation_period` AS op
        ON pc.person_id = op.person_id
       AND op.observation_period_start_date <= pc.index_date
       AND op.observation_period_end_date > pc.index_date
    GROUP BY pc.person_id, pc.index_date
),
candidate_a_exact AS (
    SELECT
        p.person_id,
        cso.index_date,
        os.earliest_ehr_start_date,
        os.latest_ehr_end_date
    FROM `{LOCKED_CDR_DATASET}.person` AS p
    INNER JOIN observation_summary AS os USING (person_id)
    INNER JOIN consent_spanning_observation AS cso USING (person_id)
    INNER JOIN ehr_consent AS ec USING (person_id)
    WHERE (
        DATE_DIFF(
            os.earliest_ehr_start_date,
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
            cso.index_date,
            DATE(
                p.year_of_birth,
                COALESCE(NULLIF(p.month_of_birth, 0), 7),
                COALESCE(NULLIF(p.day_of_birth, 0), 1)
            ),
            DAY
        ) / 365.2425
    ) BETWEEN 18 AND 120
),
candidate_concepts AS (
    SELECT DISTINCT
        co.condition_concept_id AS concept_id,
        c.concept_name
    FROM candidate_a_exact AS ca
    INNER JOIN `{LOCKED_CDR_DATASET}.condition_occurrence` AS co
        ON ca.person_id = co.person_id
    LEFT JOIN `{LOCKED_CDR_DATASET}.concept` AS c
        ON co.condition_concept_id = c.concept_id
    WHERE co.condition_start_date IS NOT NULL
      AND REGEXP_CONTAINS(
            LOWER(COALESCE(c.concept_name, '')),
            r'{candidate_vte_pattern}'
          )
)
SELECT
    concept_id,
    concept_name
FROM candidate_concepts
ORDER BY concept_id
"""

    job = client.query(sql)
    concepts = job.result().to_dataframe(create_bqstorage_client=True)

    if concepts.empty:
        raise RuntimeError("Locked CDR reconstruction returned zero candidate concepts.")

    concepts["concept_id"] = pd.to_numeric(
        concepts["concept_id"], errors="raise"
    ).astype("int64")
    concepts["concept_name"] = concepts["concept_name"].astype("string").fillna("")

    classified = pd.DataFrame(
        [classify_vte_concept(name) for name in concepts["concept_name"]]
    )
    disposition = pd.concat(
        [concepts.reset_index(drop=True), classified],
        axis=1,
    )

    override_path = find_manual_override()
    metadata_path = find_metadata()

    metadata_override_rows = None
    if metadata_path is not None:
        try:
            historical_metadata = json.loads(metadata_path.read_text())
            metadata_override_rows = int(
                historical_metadata.get("manual_override_rows", 0)
            )
        except Exception:
            metadata_override_rows = None

    applied_overrides = 0

    if override_path is not None:
        overrides = pd.read_csv(override_path)

        if len(overrides):
            required_override_columns = {
                "concept_id",
                "final_disposition",
                "final_subtype",
                "final_acuity",
                "final_context",
                "final_site_class",
                "final_rationale",
            }
            missing = required_override_columns - set(overrides.columns)
            if missing:
                raise RuntimeError(
                    "Manual override file exists but is missing required columns: "
                    + ", ".join(sorted(missing))
                )

            overrides["concept_id"] = pd.to_numeric(
                overrides["concept_id"], errors="raise"
            ).astype("int64")

            lookup = overrides.set_index("concept_id")

            for idx, row in disposition.iterrows():
                cid = int(row["concept_id"])
                if cid in lookup.index:
                    o = lookup.loc[cid]
                    for column in [
                        "final_disposition",
                        "final_subtype",
                        "final_acuity",
                        "final_context",
                        "final_site_class",
                        "final_rationale",
                    ]:
                        value = o[column]
                        if pd.notna(value) and str(value).strip():
                            disposition.at[idx, column] = value
                    applied_overrides += 1

    if (
        metadata_override_rows is not None
        and metadata_override_rows > 0
        and override_path is None
    ):
        raise RuntimeError(
            "Historical 05B metadata indicates manual overrides were used "
            f"({metadata_override_rows} rows), but the override file is not mounted. "
            "Do not finalize concept IDs until the override file is recovered."
        )

    disposition = normalize_frame(disposition)
    valid, message, concept_counts = original_05B_qc_contract(disposition)

    if not valid:
        raise RuntimeError(
            "Reconstructed Script 05B concept disposition failed the original "
            "05B QC contract: " + message
        )

    provenance = {
        "mode": "reconstructed_exact_original_05B_vocabulary_workflow",
        "cdr_dataset": LOCKED_CDR_DATASET,
        "bigquery_job_id": job.job_id,
        "manual_override_file": str(override_path) if override_path else "",
        "manual_overrides_applied": applied_overrides,
        "historical_metadata_file": str(metadata_path) if metadata_path else "",
        "historical_metadata_manual_override_rows": metadata_override_rows,
    }

    return disposition, provenance, concept_counts


def build_manuscript_outputs(
    df: pd.DataFrame,
    provenance: dict,
    discovery: pd.DataFrame,
    concept_counts: dict,
):
    long_rows = []
    grouped_rows = []

    for order, component, status, dispositions, subtypes, rationale in CATEGORY_RULES:
        mask = df["final_disposition"].isin(dispositions)
        if subtypes is not None:
            mask &= df["final_subtype"].isin(subtypes)

        subset = df.loc[
            mask,
            [
                "concept_id",
                "concept_name",
                "final_disposition",
                "final_subtype",
                "final_acuity",
                "final_context",
                "final_site_class",
                "final_rationale",
            ],
        ].copy()

        if len(subset):
            for _, row in subset.iterrows():
                long_rows.append(
                    {
                        "order": order,
                        "component": component,
                        "status": status,
                        "concept_id": int(row["concept_id"]),
                        "concept_name": str(row["concept_name"]),
                        "final_disposition": str(row["final_disposition"]),
                        "final_subtype": str(row["final_subtype"]),
                        "final_acuity": str(row["final_acuity"]),
                        "final_context": str(row["final_context"]),
                        "final_site_class": str(row["final_site_class"]),
                        "rationale": rationale,
                    }
                )

            concept_ids = "; ".join(
                subset["concept_id"].astype(int).astype(str).tolist()
            )
            concept_names = "; ".join(subset["concept_name"].astype(str).tolist())
            table_note = ""
        else:
            concept_ids = "None in reviewed Script 05B candidate set"
            concept_names = (
                "Exclusion rule specified in the phenotype workflow; "
                "no candidate concept matched this disposition."
            )
            table_note = (
                "Zero concepts in reviewed Script 05B candidate set; "
                "no concept IDs were added retrospectively."
            )

        grouped_rows.append(
            {
                "order": order,
                "component": component,
                "status": status,
                "n_concepts": int(len(subset)),
                "concept_ids": concept_ids,
                "concept_names": concept_names,
                "rationale": rationale,
                "table_note": table_note,
            }
        )

    long = pd.DataFrame(long_rows)
    if len(long):
        long = (
            long.drop_duplicates(subset=["component", "concept_id"])
            .sort_values(["order", "concept_name", "concept_id"])
            .reset_index(drop=True)
        )

    grouped = pd.DataFrame(grouped_rows)

    long_path = OUTPUT_ROOT / "15_TableS4_VTE_concept_ids_long.csv"
    grouped_path = OUTPUT_ROOT / "15_TableS4_VTE_concept_ids_grouped.csv"
    md_path = OUTPUT_ROOT / "15_TableS4_VTE_concept_ids_grouped.md"
    discovery_path = OUTPUT_ROOT / "15_source_discovery.csv"
    qc_path = OUTPUT_ROOT / "15_concept_export_qc.csv"
    metadata_path = OUTPUT_ROOT / "15_metadata.json"

    long.to_csv(long_path, index=False)
    grouped.to_csv(grouped_path, index=False)
    discovery.to_csv(discovery_path, index=False)

    def esc(value):
        return str(value).replace("|", "\\|")

    md = [
        "# Supplementary Table S4. VTE phenotype definition and concept disposition",
        "",
        "| Component | Status | Concept ID(s) | Concept definition/name | Rationale |",
        "|---|---|---|---|---|",
    ]

    for _, row in grouped.iterrows():
        name_text = row["concept_names"]
        if row["table_note"]:
            name_text += " " + row["table_note"]

        md.append(
            "| "
            + " | ".join(
                [
                    esc(row["component"]),
                    esc(row["status"]),
                    esc(row["concept_ids"]),
                    esc(name_text),
                    esc(row["rationale"]),
                ]
            )
            + " |"
        )

    md.extend(
        [
            "",
            "Concept identifiers and names reflect the reviewed Script 05B "
            "candidate concept set used to define the primary VTE phenotype. "
            "A category with zero concept IDs indicates that the exclusion rule "
            "was prespecified in the phenotype workflow but no concept in the "
            "reviewed candidate DVT/PE set matched that disposition. No additional "
            "concepts were added retrospectively for manuscript presentation.",
            "",
        ]
    )

    md_path.write_text("\n".join(md), encoding="utf-8")

    qc = pd.DataFrame(
        [
            {
                "check": "original_05B_primary_concept_set_nonempty",
                "status": concept_counts["primary"] > 0,
                "value": concept_counts["primary"],
            },
            {
                "check": "original_05B_strict_concept_set_nonempty",
                "status": concept_counts["strict_acute_current"] > 0,
                "value": concept_counts["strict_acute_current"],
            },
            {
                "check": "original_05B_DVT_primary_nonempty",
                "status": concept_counts["dvt_primary"] > 0,
                "value": concept_counts["dvt_primary"],
            },
            {
                "check": "original_05B_PE_primary_nonempty",
                "status": concept_counts["pe_primary"] > 0,
                "value": concept_counts["pe_primary"],
            },
            {
                "check": "empty_exclusion_categories_allowed_by_original_05B",
                "status": True,
                "value": (
                    f"superficial={concept_counts['superficial_excluded']}; "
                    f"unusual_site={concept_counts['unusual_site_excluded']}"
                ),
            },
            {
                "check": "participant_level_output_written",
                "status": True,
                "value": False,
            },
        ]
    )

    qc.to_csv(qc_path, index=False)

    if not qc["status"].all():
        raise RuntimeError(
            "Script 15 QC failed:\n"
            + qc.loc[~qc["status"]].to_string(index=False)
        )

    metadata = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script_version": SCRIPT_VERSION,
        "provenance": provenance,
        "source_unique_concepts": int(df["concept_id"].nunique()),
        "original_05B_concept_set_sizes": concept_counts,
        "table_category_counts": grouped[
            ["component", "n_concepts"]
        ].to_dict(orient="records"),
        "important_note": (
            "The original Script 05B did not require EXCLUDE_SUPERFICIAL or "
            "EXCLUDE_UNUSUAL_SITE to be nonempty. Zero-count categories are "
            "reported faithfully rather than populated with concepts that were "
            "not in the reviewed candidate set."
        ),
        "participant_level_output_written": False,
    }

    metadata_path.write_text(
        json.dumps(metadata, indent=2) + "\n",
        encoding="utf-8",
    )

    package_name = "VTE_JTH_Script15_TableS4_Concept_IDs"
    build_root = OUTPUT_ROOT / "package_build" / package_name

    if build_root.exists():
        shutil.rmtree(build_root)
    build_root.mkdir(parents=True, exist_ok=True)

    for source in [
        long_path,
        grouped_path,
        md_path,
        discovery_path,
        qc_path,
        metadata_path,
    ]:
        shutil.copy2(source, build_root / source.name)

    zip_path = Path("/home/jupyter") / (
        f"VTE_Results_15_{datetime.now().strftime('%Y%m%d')}.zip"
    )

    if zip_path.exists():
        zip_path.unlink()

    with zipfile.ZipFile(
        zip_path,
        "w",
        compression=zipfile.ZIP_DEFLATED,
        compresslevel=9,
    ) as archive:
        for file in sorted(build_root.iterdir()):
            if file.is_file():
                archive.write(
                    file,
                    arcname=str(Path(package_name) / file.name),
                )

    return grouped, long, qc, zip_path


def main():
    print("=" * 80)
    print("VTE / JTH Script 15: Recover exact VTE concept IDs for Table S4")
    print("Version:", SCRIPT_VERSION)
    print("=" * 80)

    discovery = discover_candidate_files()
    source_path, disposition, checked, concept_counts = choose_valid_frozen_file(
        discovery
    )

    if source_path is not None:
        print("Using validated frozen Script 05B concept disposition:")
        print(source_path)

        provenance = {
            "mode": "validated_frozen_Script05B_file",
            "source_file": str(source_path),
            "source_sha256": sha256_file(source_path),
        }
    else:
        print("\nNo mounted CSV passed the original Script 05B QC contract.")
        if not checked.empty:
            columns = [
                "path",
                "excluded_path",
                "schema_valid",
                "original_05B_qc_valid",
                "validation",
            ]
            print(checked[columns].to_string(index=False))

        disposition, provenance, concept_counts = reconstruct_from_locked_cdr()

    print("\nOriginal Script 05B concept-set sizes reproduced:")
    for key, value in concept_counts.items():
        print(f"  {key}: {value}")

    grouped, long, qc, zip_path = build_manuscript_outputs(
        disposition,
        provenance,
        checked,
        concept_counts,
    )

    print("\nSCRIPT 15 COMPLETE")
    print("\nSupplementary Table S4 category counts:")
    print(grouped[["component", "n_concepts"]].to_string(index=False))
    print("\nImportant:")
    print(
        "Zero-count superficial/unusual-site categories are allowed if that is "
        "what the original reviewed Script 05B candidate set produced."
    )
    print("\nGrouped S4:")
    print(OUTPUT_ROOT / "15_TableS4_VTE_concept_ids_grouped.csv")
    print("Long concept list:")
    print(OUTPUT_ROOT / "15_TableS4_VTE_concept_ids_long.csv")
    print("ZIP:")
    print(zip_path)
    print("=" * 80)


if __name__ == "__main__":
    main()
