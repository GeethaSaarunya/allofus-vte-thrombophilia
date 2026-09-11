#!/usr/bin/env python3
"""
All of Us VTE — INTERNAL Workbench utility
Script 04: Reconstruct final clinical-risk concept specification

WHY THIS SCRIPT
---------------
The original VTE V2 project files are no longer mounted in the current
Workbench, so this script does NOT depend on the old project directory,
old config file, or old participant-level parquet.

Instead, it reconstructs the clinical concept-disposition specification
directly from the locked All of Us CDR using the SAME candidate-search
patterns and concept-classification rules from the original:

    Script 05A  — clinical-risk curation
    Script 05A1 — refined infection/sepsis curation

This mirrors the recovery strategy successfully used for the VTE concept
export: reconstruct the reviewed concept set from the locked CDR rather
than relying on stale local paths.

IMPORTANT
---------
This is an INTERNAL Workbench script. It is not intended to be public-facing.

It queries participant records only to determine which vocabulary concepts
actually occurred before enrollment in the eligible cohort. It writes only
COUNT-FREE concept definitions and aggregate QC. It does NOT write person_id,
participant counts, record counts, or participant-level dates.

FINAL MANUSCRIPT CLINICAL FRAMEWORK
-----------------------------------
Primary five-factor burden:
    1. cancer_history
    2. major_surgery
    3. inpatient_hospitalization
    4. infection_or_sepsis
    5. fracture_or_major_trauma

Separate factor:
    6. pregnancy_or_postpartum

Historical estrogen/hormone exposure is deliberately excluded from this
repository export because it is not part of the final manuscript framework.

OUTPUT
------
/home/jupyter/VTE_JTH_Clinical_Concept_Export/
    04_clinical_concept_specification_long.csv
    04_clinical_concept_specification_grouped.csv
    04_clinical_concept_specification_grouped.md
    04_clinical_concept_export_qc.csv
    04_source_discovery.csv
    04_metadata.json

and a ZIP:
    /home/jupyter/VTE_JTH_Clinical_Concept_Export_YYYYMMDD.zip
"""

from pathlib import Path
import os
import re
import json
import datetime
import zipfile
import pandas as pd

from google.cloud import bigquery


# ======================================================================
# 0. Controls
# ======================================================================

SCRIPT_VERSION = "04_reconstruct_clinical_concepts_v1_20260911"

DEFAULT_CDR_DATASET = "wb-silky-artichoke-2408.C2024Q3R9"

CDR_DATASET = (
    os.environ.get("WORKSPACE_CDR", "").strip()
    or DEFAULT_CDR_DATASET
).strip("`")

OUTPUT_DIR = Path(
    "/home/jupyter/VTE_JTH_Clinical_Concept_Export"
)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

LONG_OUT = (
    OUTPUT_DIR
    / "04_clinical_concept_specification_long.csv"
)
GROUPED_OUT = (
    OUTPUT_DIR
    / "04_clinical_concept_specification_grouped.csv"
)
GROUPED_MD_OUT = (
    OUTPUT_DIR
    / "04_clinical_concept_specification_grouped.md"
)
QC_OUT = (
    OUTPUT_DIR
    / "04_clinical_concept_export_qc.csv"
)
SOURCE_OUT = (
    OUTPUT_DIR
    / "04_source_discovery.csv"
)
METADATA_OUT = (
    OUTPUT_DIR
    / "04_metadata.json"
)

FINAL_FEATURES = [
    "cancer_history",
    "major_surgery",
    "inpatient_hospitalization",
    "infection_or_sepsis",
    "fracture_or_major_trauma",
    "pregnancy_or_postpartum",
]

FEATURE_LABELS = {
    "cancer_history": "Cancer history",
    "major_surgery": "Major surgery",
    "inpatient_hospitalization": "Inpatient hospitalization",
    "infection_or_sepsis": "Infection/sepsis",
    "fracture_or_major_trauma": "Fracture/major trauma",
    "pregnancy_or_postpartum": "Pregnancy/postpartum",
}


# ======================================================================
# 1. BigQuery connection
# ======================================================================

print("=" * 78)
print("All of Us VTE — Script 04 clinical concept reconstruction")
print("Version:", SCRIPT_VERSION)
print("CDR:", CDR_DATASET)
print("=" * 78)

project = (
    os.environ.get("GOOGLE_PROJECT", "").strip()
    or os.environ.get("GOOGLE_CLOUD_PROJECT", "").strip()
    or None
)

if project:
    client = bigquery.Client(project=project)
else:
    client = bigquery.Client()

print("BigQuery project:", client.project)


# ======================================================================
# 2. Optional historical manual-override discovery
# ======================================================================

BASE_OVERRIDE_NAME = "05A_manual_concept_overrides.csv"
INFECTION_OVERRIDE_NAME = "05A1_infection_manual_overrides.csv"

def discover_optional_named_file(filename):
    """
    Look for a historical override file without requiring the old project path.
    If more than one exists, prefer a path containing the original V2 project.
    """
    roots = [
        Path.cwd(),
        Path("/home/jupyter"),
        Path("/home/jupyter/workspace"),
        Path("/home/dataproc/workspace"),
    ]

    candidates = []

    for root in roots:
        if not root.exists():
            continue

        # Fast direct candidates.
        candidates.extend([
            root / filename,
            root / "documentation" / "curation" / filename,
        ])

        # Bounded recursive search only under workspace-like roots.
        if root.name in {"workspace"} or str(root).endswith("/workspace"):
            try:
                candidates.extend(root.glob(f"**/{filename}"))
            except Exception:
                pass

    existing = []
    seen = set()

    for p in candidates:
        try:
            p = p.resolve()
        except Exception:
            pass

        if p.exists() and p.is_file() and str(p) not in seen:
            seen.add(str(p))
            existing.append(p)

    if not existing:
        return None, []

    existing = sorted(
        existing,
        key=lambda p: (
            "VTE_AllOfUs_CDRv8_C2024Q3R9_V2" not in str(p),
            len(str(p)),
            str(p),
        )
    )

    return existing[0], existing


base_override_path, base_override_candidates = (
    discover_optional_named_file(BASE_OVERRIDE_NAME)
)

infection_override_path, infection_override_candidates = (
    discover_optional_named_file(INFECTION_OVERRIDE_NAME)
)

source_rows = []

for label, chosen, candidates in [
    (
        "05A_manual_override",
        base_override_path,
        base_override_candidates,
    ),
    (
        "05A1_infection_manual_override",
        infection_override_path,
        infection_override_candidates,
    ),
]:
    if candidates:
        for p in candidates:
            source_rows.append({
                "source_type": label,
                "path": str(p),
                "exists": True,
                "selected": bool(chosen and p == chosen),
            })
    else:
        source_rows.append({
            "source_type": label,
            "path": "",
            "exists": False,
            "selected": False,
        })

pd.DataFrame(source_rows).to_csv(
    SOURCE_OUT,
    index=False,
)

print("\nHistorical manual overrides:")
print("  05A :", base_override_path or "none found")
print(
    "  05A1:",
    infection_override_path or "none found",
)


# ======================================================================
# 3. Exact eligible cohort CTE copied from original Script 05A
# ======================================================================

exact_candidate_cte = f"""
WITH observation_summary AS (
    SELECT
        person_id,
        MIN(observation_period_start_date)
            AS earliest_ehr_start_date,
        MAX(observation_period_end_date)
            AS latest_ehr_end_date
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

consent_spanning_observation AS (
    SELECT
        pc.person_id,
        pc.index_date,
        MIN(op.observation_period_start_date)
            AS spanning_period_start_date,
        MAX(op.observation_period_end_date)
            AS spanning_period_end_date
    FROM primary_consent AS pc
    INNER JOIN `{CDR_DATASET}.observation_period` AS op
        ON pc.person_id = op.person_id
       AND op.observation_period_start_date <= pc.index_date
       AND op.observation_period_end_date > pc.index_date
    GROUP BY
        pc.person_id,
        pc.index_date
),

candidate_a_exact AS (
    SELECT
        p.person_id,
        cso.index_date,
        os.earliest_ehr_start_date,
        os.latest_ehr_end_date,
        cso.spanning_period_start_date,
        cso.spanning_period_end_date
    FROM `{CDR_DATASET}.person` AS p
    INNER JOIN observation_summary AS os
        USING (person_id)
    INNER JOIN consent_spanning_observation AS cso
        USING (person_id)
    INNER JOIN ehr_consent AS ec
        USING (person_id)
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
)
"""


# ======================================================================
# 4. Candidate concept discovery copied from original Script 05A
# ======================================================================

candidate_patterns = {
    "cancer_history": (
        r"cancer|neoplasm|carcinoma|leukemia|leukaemia|lymphoma|"
        r"myeloma|melanoma|sarcoma|metastatic"
    ),
    "infection_or_sepsis": (
        r"sepsis|septic|infectious|infection"
    ),
    "fracture_or_major_trauma": (
        r"fracture|trauma|traumatic|crush injury|dislocation"
    ),
    "pregnancy_or_postpartum": (
        r"pregnan|postpartum|post-partum|puerper|gestation|"
        r"labor|labour|delivery|childbirth|miscarriage|abortion|ectopic"
    ),
    "major_surgery": (
        r"surgery|surgical|resection|excision|transplant|arthroplasty|"
        r"bypass|amputation|colectomy|gastrectomy|pancreatectomy|"
        r"hysterectomy|nephrectomy|lobectomy|pneumonectomy|thoracotomy|"
        r"laparotomy|craniotomy|hepatectomy|splenectomy|esophagectomy|"
        r"prostatectomy|mastectomy|endarterectomy|aneurysm repair|"
        r"valve replacement|valve repair|fusion|fixation|cesarean|"
        r"caesarean|whipple|joint replacement|open reduction"
    ),
}


def condition_candidate_select(feature_name, pattern):
    return f"""
    SELECT DISTINCT
        'condition' AS domain,
        '{feature_name}' AS feature_name,
        co.condition_concept_id AS concept_id,
        c.concept_name
    FROM candidate_a_exact AS ca
    INNER JOIN `{CDR_DATASET}.condition_occurrence` AS co
        ON ca.person_id = co.person_id
       AND co.condition_start_date < ca.index_date
    LEFT JOIN `{CDR_DATASET}.concept` AS c
        ON co.condition_concept_id = c.concept_id
    WHERE REGEXP_CONTAINS(
        LOWER(COALESCE(c.concept_name, '')),
        r'{pattern}'
    )
    """


condition_selects = [
    condition_candidate_select(
        feature,
        candidate_patterns[feature],
    )
    for feature in [
        "cancer_history",
        "infection_or_sepsis",
        "fracture_or_major_trauma",
        "pregnancy_or_postpartum",
    ]
]


procedure_select = f"""
    SELECT DISTINCT
        'procedure' AS domain,
        'major_surgery' AS feature_name,
        po.procedure_concept_id AS concept_id,
        pc.concept_name
    FROM candidate_a_exact AS ca
    INNER JOIN `{CDR_DATASET}.procedure_occurrence` AS po
        ON ca.person_id = po.person_id
       AND po.procedure_date < ca.index_date
    LEFT JOIN `{CDR_DATASET}.concept` AS pc
        ON po.procedure_concept_id = pc.concept_id
    WHERE REGEXP_CONTAINS(
        LOWER(COALESCE(pc.concept_name, '')),
        r'{candidate_patterns["major_surgery"]}'
    )
"""


visit_select = f"""
    SELECT DISTINCT
        'visit' AS domain,
        'inpatient_hospitalization' AS feature_name,
        vo.visit_concept_id AS concept_id,
        vc.concept_name
    FROM candidate_a_exact AS ca
    INNER JOIN `{CDR_DATASET}.visit_occurrence` AS vo
        ON ca.person_id = vo.person_id
       AND vo.visit_start_date < ca.index_date
    LEFT JOIN `{CDR_DATASET}.concept` AS vc
        ON vo.visit_concept_id = vc.concept_id
    WHERE REGEXP_CONTAINS(
        LOWER(COALESCE(vc.concept_name, '')),
        r'inpatient'
    )
"""


candidate_sql = (
    exact_candidate_cte
    + "\n"
    + "\nUNION ALL\n".join(
        condition_selects
        + [procedure_select, visit_select]
    )
)


print("\nQuerying candidate clinical concepts from locked CDR...")

candidate_job = client.query(candidate_sql)

candidate_concepts = (
    candidate_job
    .result()
    .to_dataframe()
)

if len(candidate_concepts) == 0:
    raise RuntimeError(
        "No candidate clinical concepts were returned."
    )

candidate_concepts["concept_id"] = pd.to_numeric(
    candidate_concepts["concept_id"],
    errors="coerce",
).astype("Int64")

candidate_concepts["concept_name"] = (
    candidate_concepts["concept_name"]
    .astype("string")
    .fillna("")
)

candidate_concepts = (
    candidate_concepts
    .dropna(subset=["concept_id"])
    .drop_duplicates(
        subset=[
            "domain",
            "feature_name",
            "concept_id",
        ]
    )
    .reset_index(drop=True)
)

print(
    "Candidate concept rows:",
    f"{len(candidate_concepts):,}",
)


# ======================================================================
# 5. Original Script 05A concept-disposition rules
# ======================================================================

def contains(text, pattern):
    return bool(
        re.search(
            pattern,
            str(text),
            flags=re.IGNORECASE,
        )
    )


def classify_cancer(name):
    exclude_pattern = (
        r"benign|uncertain behavior|in situ|malignant essential hypertension|"
        r"malignant hypertension|screening|suspected|rule out|family history"
    )
    retain_pattern = (
        r"malignant neoplasm|malignant neoplastic|malignant tumor|"
        r"primary malignant|secondary malignant|metastatic|carcinoma|"
        r"leukemia|leukaemia|lymphoma|myeloma|melanoma|sarcoma"
    )

    if contains(name, exclude_pattern):
        return (
            "EXCLUDE",
            "Excluded benign, uncertain, in-situ, non-neoplastic malignant terminology, or noncase context.",
            "",
        )

    if contains(name, retain_pattern):
        return (
            "RETAIN",
            "Explicit malignant neoplasm or hematologic malignancy terminology.",
            "cancer_history",
        )

    return (
        "EXCLUDE",
        "Not on the explicit malignant-neoplasm whitelist.",
        "",
    )


def classify_infection_05a(name):
    exclude_pattern = (
        r"carrier|sequela|family history|screening|suspected|rule out|"
        r"vaccination|immunization|laboratory|antibody|serology|"
        r"cardiomyopathy due to|arthritis .* due to|neutropenia associated"
    )
    retain_pattern = (
        r"sepsis|septic|infectious disease|infection"
    )

    if contains(name, exclude_pattern):
        return (
            "EXCLUDE",
            "Excluded carrier, sequela, screening, or other non-active-infection context.",
            "",
        )

    if contains(name, retain_pattern):
        return (
            "RETAIN",
            "Explicit infection or sepsis terminology.",
            "infection_or_sepsis",
        )

    return (
        "EXCLUDE",
        "Not on the explicit infection/sepsis whitelist.",
        "",
    )


def classify_fracture_trauma(name):
    exclude_pattern = (
        r"pathologic|pathological|stress fracture|fatigue fracture|"
        r"osteoporotic|nonunion|malunion|delayed healing|sequela|"
        r"history of|old fracture|healed fracture|congenital|birth injury|"
        r"screening|suspected|rule out"
    )

    retain_pattern = (
        r"fracture|traumatic injury|major trauma|multiple trauma|"
        r"crush injury|dislocation"
    )

    if contains(name, exclude_pattern):
        return (
            "EXCLUDE",
            "Excluded nontraumatic, historical, healing-complication, or noncase fracture terminology.",
            "",
        )

    if contains(name, retain_pattern):
        return (
            "RETAIN",
            "Explicit traumatic fracture, dislocation, crush injury, or major-trauma terminology.",
            "fracture_or_major_trauma",
        )

    return (
        "EXCLUDE",
        "Not on the explicit traumatic fracture/major-trauma whitelist.",
        "",
    )


def classify_pregnancy(name):
    exclude_pattern = (
        r"pregnancy test|test for pregnancy|pregnancy counseling|"
        r"preconception|family history|trying to conceive|"
        r"desires pregnancy|infertility|contraceptive|screening|"
        r"suspected|rule out"
    )

    retain_pattern = (
        r"pregnan|trimester|gestation|postpartum|post-partum|puerper|"
        r"labor|labour|delivery|childbirth|miscarriage|abortion|ectopic"
    )

    if contains(name, exclude_pattern):
        return (
            "EXCLUDE",
            "Excluded testing, counseling, fertility, contraception, or noncase context.",
            "",
        )

    if contains(name, retain_pattern):
        return (
            "RETAIN",
            "Explicit pregnancy, delivery, pregnancy-loss, or postpartum terminology.",
            "pregnancy_or_postpartum",
        )

    return (
        "EXCLUDE",
        "Not on the explicit pregnancy/postpartum whitelist.",
        "",
    )


def classify_major_surgery(name):
    direct_major_pattern = (
        r"transplant|amputation|arthroplasty|coronary artery bypass|"
        r"gastric bypass|aort.*bypass|femor.*bypass|colectomy|gastrectomy|"
        r"pancreatectomy|hysterectomy|nephrectomy|lobectomy|pneumonectomy|"
        r"thoracotomy|laparotomy|craniotomy|hepatectomy|splenectomy|"
        r"esophagectomy|prostatectomy|mastectomy|endarterectomy|"
        r"aneurysm repair|valve replacement|valve repair|spinal fusion|"
        r"joint replacement|open reduction.*internal fixation|"
        r"internal fixation|cesarean|caesarean|whipple|"
        r"bariatric surgery|radical resection"
    )

    generic_operating_pattern = (
        r"surgery|surgical procedure|operative procedure|resection|"
        r"open repair|excision|fixation|fusion"
    )

    exclusion_pattern = (
        r"pathology|histopath|cytology|specimen|laboratory|"
        r"imaging|ultrasound|ultrasonography|radiograph|x-ray|"
        r"computed tomography|magnetic resonance|angiograph|venograph|"
        r"echocardiograph|myocardial perfusion|infusion|hydration|"
        r"injection|preoperative|pre-operative|pre-surgery|"
        r"evaluation|assessment|consultation|anesthesia|anaesthesia|"
        r"biopsy|endoscopy|colonoscopy|bronchoscopy|cystoscopy|"
        r"arthroscopy|destruction|ablation|cauter|cryotherapy|laser|"
        r"sclerotherapy|catheter|pacemaker guidance|device guidance|"
        r"dressing|wound care|debridement|dental|tooth|nail|"
        r"skin lesion|actinic keratos|cataract|ophthalmic|"
        r"dialysis|radiation|chemotherapy|physical therapy|"
        r"postoperative care|aftercare|follow-up|follow up|history of"
    )

    if contains(name, exclusion_pattern):
        return (
            "EXCLUDE",
            "Excluded diagnostic, pathology, imaging, infusion, minor-procedure, guidance, or perioperative-context terminology.",
            "",
        )

    if contains(name, direct_major_pattern):
        return (
            "RETAIN",
            "Explicit named major operation on the high-specificity whitelist.",
            "major_surgery",
        )

    if contains(name, generic_operating_pattern):
        return (
            "RETAIN_IF_INPATIENT",
            "Generic operative terminology retained only when linked to an inpatient visit.",
            "major_surgery",
        )

    return (
        "EXCLUDE",
        "Not on the explicit major-operation whitelist.",
        "",
    )


def classify_inpatient(name):
    if contains(name, r"inpatient"):
        return (
            "RETAIN",
            "Explicit inpatient visit concept.",
            "inpatient_hospitalization",
        )

    return (
        "EXCLUDE",
        "Not an explicit inpatient visit concept.",
        "",
    )


classifier_map = {
    "cancer_history": classify_cancer,
    "infection_or_sepsis": classify_infection_05a,
    "fracture_or_major_trauma": classify_fracture_trauma,
    "pregnancy_or_postpartum": classify_pregnancy,
    "major_surgery": classify_major_surgery,
    "inpatient_hospitalization": classify_inpatient,
}


rule_rows = []

for row in candidate_concepts.itertuples(index=False):
    feature = str(row.feature_name)
    name = str(row.concept_name or "")

    if feature not in classifier_map:
        raise KeyError(
            f"No classifier defined for feature: {feature}"
        )

    disposition, rationale, target_feature = (
        classifier_map[feature](name)
    )

    rule_rows.append({
        "rule_disposition": disposition,
        "rule_rationale": rationale,
        "rule_target_feature": target_feature,
    })


base_disposition = pd.concat(
    [
        candidate_concepts.reset_index(drop=True),
        pd.DataFrame(rule_rows),
    ],
    axis=1,
)


# ======================================================================
# 6. Apply historical Script 05A manual overrides if recovered
# ======================================================================

if base_override_path is not None:
    base_overrides = pd.read_csv(base_override_path)

    if len(base_overrides):
        required = {
            "domain",
            "feature_name",
            "concept_id",
            "final_disposition",
            "target_feature",
            "final_rationale",
        }

        missing = required - set(base_overrides.columns)

        if missing:
            raise KeyError(
                "Recovered Script 05A override file is missing: "
                + ", ".join(sorted(missing))
            )

        base_overrides["concept_id"] = pd.to_numeric(
            base_overrides["concept_id"],
            errors="raise",
        ).astype("Int64")

        base_disposition = base_disposition.merge(
            base_overrides,
            on=[
                "domain",
                "feature_name",
                "concept_id",
            ],
            how="left",
            validate="many_to_one",
        )
    else:
        base_overrides = pd.DataFrame()
else:
    base_overrides = pd.DataFrame()


for col in [
    "final_disposition",
    "target_feature",
    "final_rationale",
]:
    if col not in base_disposition.columns:
        base_disposition[col] = pd.NA


base_disposition["manual_override_applied_05A"] = (
    base_disposition["final_disposition"].notna()
)

base_disposition["final_disposition"] = (
    base_disposition["final_disposition"]
    .fillna(
        base_disposition["rule_disposition"]
    )
)

base_disposition["target_feature"] = (
    base_disposition["target_feature"]
    .fillna(
        base_disposition["rule_target_feature"]
    )
)

base_disposition["final_rationale"] = (
    base_disposition["final_rationale"]
    .fillna(
        base_disposition["rule_rationale"]
    )
)


# ======================================================================
# 7. Refine infection/sepsis exactly as original Script 05A1
# ======================================================================

DIRECT_SEVERE_PATTERN = (
    r"sepsis|septicemia|septicaemia|septic shock|"
    r"bacteremia|bacteraemia|bloodstream infection|"
    r"fungemia|fungaemia|viremia|viraemia|pyemia|pyaemia|"
    r"systemic infection|disseminated infection|invasive infection|"
    r"infective endocarditis|bacterial endocarditis|"
    r"meningitis|encephalitis|brain abscess|cerebral abscess|"
    r"epidural abscess|spinal abscess|"
    r"peritonitis|empyema|osteomyelitis|"
    r"septic arthritis|infectious arthritis|"
    r"necrotizing fasciitis|necrotising fasciitis|gas gangrene|"
    r"liver abscess|hepatic abscess|splenic abscess|"
    r"intra[- ]abdominal abscess|retroperitoneal abscess|"
    r"deep incisional surgical site infection|"
    r"organ[/ -]space surgical site infection"
)

NONACTIVE_CONTEXT_PATTERN = (
    r"history of|family history|screening|suspected|rule out|"
    r"carrier|colonization|colonisation|exposure to|contact with|"
    r"immunization|immunisation|vaccination|serology|antibody|"
    r"laboratory test|infection status|sequela"
)

SCOPE_PATTERNS = [
    (
        "SYSTEMIC_OR_INVASIVE",
        DIRECT_SEVERE_PATTERN,
    ),
    (
        "RESPIRATORY",
        r"pneumonia|respiratory|bronch|influenza|covid|coronavirus|"
        r"pharyng|tonsill|laryng|sinus|rhinitis",
    ),
    (
        "URINARY_OR_RENAL",
        r"urinary|cystitis|pyeloneph|kidney infection|renal infection|urethrit",
    ),
    (
        "SKIN_SOFT_TISSUE",
        r"cellulitis|skin infection|soft tissue|wound infection|"
        r"cutaneous|subcutaneous|folliculitis|impetigo|abscess",
    ),
    (
        "GASTROINTESTINAL_OR_INTRAABDOMINAL",
        r"gastroenter|enterocol|colitis|cholangitis|cholecystitis|"
        r"diverticulitis|appendicitis|pancreatic infection|"
        r"intra[- ]abdominal|peritonitis",
    ),
    (
        "GENITAL_OR_REPRODUCTIVE",
        r"vaginitis|vulvovag|pelvic inflammatory|endometritis|"
        r"salpingitis|orchitis|epididymitis|genital infection",
    ),
    (
        "ENT_DENTAL_OR_OCULAR",
        r"otitis|ear infection|dental|tooth|gingivitis|periodont|"
        r"conjunctivitis|blepharitis|ocular infection|eye infection",
    ),
    (
        "BONE_JOINT_OR_NEUROLOGIC",
        r"osteomyelitis|septic arthritis|infectious arthritis|"
        r"meningitis|encephalitis|brain abscess|epidural abscess",
    ),
    (
        "VIRAL_FUNGAL_OR_PARASITIC",
        r"viral|virus|fungal|fungus|candid|mycos|parasit|protozo|helminth",
    ),
]


def classify_scope(name):
    for scope, pattern in SCOPE_PATTERNS:
        if contains(name, pattern):
            return scope

    return "OTHER_OR_UNSPECIFIED"


def classify_infection_05a1(name):
    if contains(
        name,
        NONACTIVE_CONTEXT_PATTERN,
    ):
        return (
            "EXCLUDE",
            "Non-active, historical, screening, exposure, carrier, or laboratory context.",
        )

    if contains(
        name,
        DIRECT_SEVERE_PATTERN,
    ):
        return (
            "DIRECT_RETAIN",
            "Explicit severe systemic or invasive infection; retained without requiring visit linkage.",
        )

    return (
        "RETAIN_IF_INPATIENT",
        "Active infection retained only when the condition record is linked to an explicit inpatient visit.",
    )


infection_base = base_disposition.loc[
    base_disposition["feature_name"].eq(
        "infection_or_sepsis"
    )
].copy()

# Original Script 05A1 only refined infection concepts that Script 05A
# had retained for infection_or_sepsis.
infection_to_refine = infection_base.loc[
    infection_base["final_disposition"].eq("RETAIN")
    & infection_base["target_feature"].eq(
        "infection_or_sepsis"
    )
].copy()

infection_refined_rows = []

for row in infection_to_refine.itertuples(
    index=False
):
    disposition, rationale = (
        classify_infection_05a1(
            row.concept_name
        )
    )

    infection_refined_rows.append({
        "concept_id": int(row.concept_id),
        "concept_name": str(row.concept_name),
        "infection_scope": classify_scope(
            row.concept_name
        ),
        "final_disposition_05A1": disposition,
        "final_rationale_05A1": rationale,
    })


infection_refined = pd.DataFrame(
    infection_refined_rows
)


# ======================================================================
# 8. Apply historical Script 05A1 overrides if recovered
# ======================================================================

if infection_override_path is not None:
    infection_overrides = pd.read_csv(
        infection_override_path
    )

    if len(infection_overrides):
        required = {
            "concept_id",
            "final_disposition",
            "final_rationale",
        }

        missing = required - set(
            infection_overrides.columns
        )

        if missing:
            raise KeyError(
                "Recovered Script 05A1 override file is missing: "
                + ", ".join(sorted(missing))
            )

        infection_overrides["concept_id"] = pd.to_numeric(
            infection_overrides["concept_id"],
            errors="raise",
        ).astype("int64")

        infection_refined = (
            infection_refined
            .merge(
                infection_overrides[
                    [
                        "concept_id",
                        "final_disposition",
                        "final_rationale",
                    ]
                ],
                on="concept_id",
                how="left",
                validate="one_to_one",
            )
        )

        infection_refined[
            "manual_override_applied_05A1"
        ] = (
            infection_refined[
                "final_disposition"
            ].notna()
        )

        infection_refined[
            "final_disposition_05A1"
        ] = (
            infection_refined[
                "final_disposition"
            ]
            .fillna(
                infection_refined[
                    "final_disposition_05A1"
                ]
            )
        )

        infection_refined[
            "final_rationale_05A1"
        ] = (
            infection_refined[
                "final_rationale"
            ]
            .fillna(
                infection_refined[
                    "final_rationale_05A1"
                ]
            )
        )

        infection_refined = (
            infection_refined
            .drop(
                columns=[
                    "final_disposition",
                    "final_rationale",
                ]
            )
        )
    else:
        infection_overrides = pd.DataFrame()
else:
    infection_overrides = pd.DataFrame()


if "manual_override_applied_05A1" not in (
    infection_refined.columns
):
    infection_refined[
        "manual_override_applied_05A1"
    ] = False


# ======================================================================
# 9. Construct final manuscript concept specification
# ======================================================================

# Non-infection features come directly from final 05A disposition.
final_noninfection = (
    base_disposition.loc[
        ~base_disposition["feature_name"].eq(
            "infection_or_sepsis"
        )
    ]
    .copy()
)

final_noninfection[
    "definition_source"
] = "05A_original_clinical_curation"

final_noninfection[
    "infection_scope"
] = pd.NA


# Infection concepts excluded by 05A stay excluded.
infection_05a_excluded = (
    infection_base.loc[
        ~(
            infection_base[
                "final_disposition"
            ].eq("RETAIN")
            & infection_base[
                "target_feature"
            ].eq("infection_or_sepsis")
        )
    ]
    .copy()
)

infection_05a_excluded[
    "definition_source"
] = "05A_original_exclusion"

infection_05a_excluded[
    "infection_scope"
] = infection_05a_excluded[
    "concept_name"
].map(classify_scope)


# Concepts retained by 05A take their FINAL definition from 05A1.
infection_05a1_final = (
    infection_to_refine
    .merge(
        infection_refined[
            [
                "concept_id",
                "infection_scope",
                "final_disposition_05A1",
                "final_rationale_05A1",
                "manual_override_applied_05A1",
            ]
        ],
        on="concept_id",
        how="left",
        validate="one_to_one",
    )
)

if (
    infection_05a1_final[
        "final_disposition_05A1"
    ].isna().any()
):
    raise RuntimeError(
        "Some Script 05A retained infection concepts "
        "were not assigned a final Script 05A1 disposition."
    )

infection_05a1_final[
    "final_disposition"
] = (
    infection_05a1_final[
        "final_disposition_05A1"
    ]
)

infection_05a1_final[
    "final_rationale"
] = (
    infection_05a1_final[
        "final_rationale_05A1"
    ]
)

infection_05a1_final[
    "target_feature"
] = (
    "infection_or_sepsis"
)

infection_05a1_final[
    "definition_source"
] = "05A1_refined_infection_curation"


# Common final columns.
PUBLIC_COLUMNS = [
    "domain",
    "feature_name",
    "concept_id",
    "concept_name",
    "final_disposition",
    "target_feature",
    "final_rationale",
    "definition_source",
    "infection_scope",
    "manual_override_applied_05A",
    "manual_override_applied_05A1",
]


frames = [
    final_noninfection,
    infection_05a_excluded,
    infection_05a1_final,
]

for frame in frames:
    for col in PUBLIC_COLUMNS:
        if col not in frame.columns:
            frame[col] = pd.NA


public = pd.concat(
    [
        frame[PUBLIC_COLUMNS]
        for frame in frames
    ],
    ignore_index=True,
)

public = public.loc[
    public["feature_name"].isin(
        FINAL_FEATURES
    )
].copy()

public["feature_label"] = public[
    "feature_name"
].map(FEATURE_LABELS)

public["status"] = (
    public["final_disposition"]
    .map({
        "RETAIN": "Included",
        "DIRECT_RETAIN": "Included",
        "RETAIN_IF_INPATIENT": "Conditional inclusion",
        "RETAIN_CONTEXT": "Context only",
        "CONTEXT_ONLY": "Context only",
        "EXCLUDE": "Excluded",
    })
    .fillna("Other")
)

public = (
    public
    .drop_duplicates(
        subset=[
            "domain",
            "feature_name",
            "concept_id",
            "final_disposition",
            "target_feature",
        ]
    )
    .sort_values(
        [
            "feature_name",
            "status",
            "concept_name",
            "concept_id",
        ]
    )
    .reset_index(drop=True)
)


# ======================================================================
# 10. Guardrails
# ======================================================================

observed_features = set(
    public["feature_name"]
    .dropna()
    .astype(str)
)

missing_final_features = (
    set(FINAL_FEATURES)
    - observed_features
)

if missing_final_features:
    raise RuntimeError(
        "Final concept reconstruction is missing manuscript feature(s): "
        + ", ".join(
            sorted(missing_final_features)
        )
    )

if public["concept_id"].isna().any():
    raise RuntimeError(
        "Final concept specification contains missing concept IDs."
    )

if public["concept_name"].isna().any():
    raise RuntimeError(
        "Final concept specification contains missing concept names."
    )

# Historical hormone must never enter this export.
if public["feature_name"].eq(
    "estrogen_or_hormone_exposure"
).any():
    raise RuntimeError(
        "Historical estrogen/hormone feature entered the export."
    )

# No participant-level or count columns.
for forbidden in [
    "person_id",
    "research_id",
    "participant_count",
    "record_count",
    "earliest_record_date",
    "latest_record_date",
]:
    if forbidden in public.columns:
        raise RuntimeError(
            f"Forbidden output column present: {forbidden}"
        )


# ======================================================================
# 11. Write long specification
# ======================================================================

public.to_csv(
    LONG_OUT,
    index=False,
)


# ======================================================================
# 12. Grouped human-readable specification
# ======================================================================

group_rows = []

for (
    feature_name,
    feature_label,
    final_disposition,
    status,
), sub in public.groupby(
    [
        "feature_name",
        "feature_label",
        "final_disposition",
        "status",
    ],
    dropna=False,
    sort=True,
):
    sub = (
        sub
        .sort_values(
            [
                "concept_name",
                "concept_id",
            ]
        )
    )

    concept_ids = "; ".join(
        str(int(x))
        for x in sub["concept_id"]
        .dropna()
        .astype("int64")
        .tolist()
    )

    concept_names = "; ".join(
        sub["concept_name"]
        .astype(str)
        .tolist()
    )

    rationales = (
        sub["final_rationale"]
        .dropna()
        .astype(str)
        .unique()
        .tolist()
    )

    group_rows.append({
        "feature": feature_label,
        "feature_name": feature_name,
        "status": status,
        "final_disposition": final_disposition,
        "n_concepts": int(
            sub["concept_id"].nunique()
        ),
        "concept_ids": concept_ids,
        "concept_names": concept_names,
        "rationale": " | ".join(rationales),
    })


grouped = pd.DataFrame(group_rows)

grouped.to_csv(
    GROUPED_OUT,
    index=False,
)


md_lines = [
    "# Clinical-risk concept specification",
    "",
    "Reconstructed from the locked All of Us CDR using the original "
    "Script 05A clinical curation rules and Script 05A1 refined "
    "infection/sepsis rules.",
    "",
    "| Feature | Status | Disposition | Concept IDs | Concept names | Rationale |",
    "|---|---|---|---|---|---|",
]

for row in grouped.itertuples(index=False):
    md_lines.append(
        "| "
        + " | ".join([
            str(row.feature).replace("|", "\\|"),
            str(row.status).replace("|", "\\|"),
            str(row.final_disposition).replace("|", "\\|"),
            str(row.concept_ids).replace("|", "\\|"),
            str(row.concept_names).replace("|", "\\|"),
            str(row.rationale).replace("|", "\\|"),
        ])
        + " |"
    )

GROUPED_MD_OUT.write_text(
    "\n".join(md_lines) + "\n",
    encoding="utf-8",
)


# ======================================================================
# 13. QC and metadata
# ======================================================================

qc_rows = [
    {
        "check": "candidate_concept_set_nonempty",
        "status": len(candidate_concepts) > 0,
        "value": len(candidate_concepts),
    },
    {
        "check": "all_final_features_present",
        "status": len(
            missing_final_features
        ) == 0,
        "value": ",".join(
            sorted(observed_features)
        ),
    },
    {
        "check": "infection_05A1_refinement_nonempty",
        "status": len(
            infection_05a1_final
        ) > 0,
        "value": len(
            infection_05a1_final
        ),
    },
    {
        "check": "historical_hormone_excluded",
        "status": not public[
            "feature_name"
        ].eq(
            "estrogen_or_hormone_exposure"
        ).any(),
        "value": True,
    },
    {
        "check": "participant_level_output_written",
        "status": True,
        "value": False,
    },
    {
        "check": "participant_counts_written",
        "status": True,
        "value": False,
    },
    {
        "check": "05A_manual_override_rows",
        "status": True,
        "value": int(
            len(base_overrides)
        ),
    },
    {
        "check": "05A1_manual_override_rows",
        "status": True,
        "value": int(
            len(infection_overrides)
        ),
    },
]

qc = pd.DataFrame(qc_rows)

qc.to_csv(
    QC_OUT,
    index=False,
)

if not qc["status"].all():
    print(qc.to_string(index=False))
    raise RuntimeError(
        "Clinical concept export QC failed."
    )


metadata = {
    "created_utc": (
        datetime.datetime.now(
            datetime.timezone.utc
        ).isoformat()
    ),
    "script_version": SCRIPT_VERSION,
    "cdr_dataset": CDR_DATASET,
    "bigquery_project": client.project,
    "bigquery_job_id": candidate_job.job_id,
    "reconstruction_mode": (
        "locked_CDR_plus_original_05A_and_05A1_rules"
    ),
    "candidate_rows": int(
        len(candidate_concepts)
    ),
    "public_rows": int(
        len(public)
    ),
    "final_features": FINAL_FEATURES,
    "manual_overrides": {
        "05A_path": (
            str(base_override_path)
            if base_override_path
            else ""
        ),
        "05A_rows": int(
            len(base_overrides)
        ),
        "05A1_path": (
            str(infection_override_path)
            if infection_override_path
            else ""
        ),
        "05A1_rows": int(
            len(infection_overrides)
        ),
    },
    "participant_level_output_written": False,
    "participant_counts_written": False,
}

METADATA_OUT.write_text(
    json.dumps(
        metadata,
        indent=2,
        sort_keys=True,
    )
    + "\n",
    encoding="utf-8",
)


# ======================================================================
# 14. ZIP package
# ======================================================================

zip_path = Path(
    "/home/jupyter/"
    f"VTE_JTH_Clinical_Concept_Export_"
    f"{datetime.datetime.now().strftime('%Y%m%d')}.zip"
)

with zipfile.ZipFile(
    zip_path,
    "w",
    compression=zipfile.ZIP_DEFLATED,
) as zf:
    for path in [
        LONG_OUT,
        GROUPED_OUT,
        GROUPED_MD_OUT,
        QC_OUT,
        SOURCE_OUT,
        METADATA_OUT,
    ]:
        zf.write(
            path,
            arcname=(
                "VTE_JTH_Clinical_Concept_Export/"
                + path.name
            ),
        )


# ======================================================================
# 15. Completion
# ======================================================================

print("\n" + "=" * 78)
print("SCRIPT 04 COMPLETE")
print("=" * 78)

print(
    "Final count-free concept rows:",
    f"{len(public):,}",
)

print("\nRows by feature:")
for feature in FINAL_FEATURES:
    print(
        " ",
        feature,
        ":",
        int(
            public["feature_name"]
            .eq(feature)
            .sum()
        ),
    )

print("\nOutputs:")
for path in [
    LONG_OUT,
    GROUPED_OUT,
    GROUPED_MD_OUT,
    QC_OUT,
    SOURCE_OUT,
    METADATA_OUT,
]:
    print(" ", path)

print("\nZIP:")
print(" ", zip_path)

print("\nGuardrails:")
print("  participant-level output: NO")
print("  participant counts: NO")
print("  record counts: NO")
print("  historical hormone feature: NO")
print("  infection final source: Script 05A1 rules")

print("=" * 78)
