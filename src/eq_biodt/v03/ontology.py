"""OWL 2 ontology generation, Kobayashi-MaRust (km) classification and justification extraction."""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import time
import hashlib
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from .evidence import CBB_STEPS, sha256

RAK = "https://rubalkhali.science/kb/"
KEGG = "https://www.kegg.jp/entry/"
OBO = "http://purl.obolibrary.org/obo/"
SIO = "http://semanticscience.org/resource/"

ALL_CBB_KOS = tuple(sorted({ko for alts in CBB_STEPS.values() for ko in alts} | {"K01602"}))

ENDPOINT_CLASSES = (
    "GenomicCapacityEndpoint",
    "MappedReadCompositionEndpoint",
    "DNARecoveryEndpoint",
    "ViableDNAEndpoint",
    "BiomassEndpoint",
    "PhysiologicalActivityEndpoint",
)


def find_km() -> Path:
    """Locate the Kobayashi-MaRust (km) OWL 2 reasoner executable."""
    env_bin = os.environ.get("KM_BIN")
    if env_bin and Path(env_bin).exists():
        return Path(env_bin)
    pinned = Path(__file__).resolve().parents[3] / ".tools" / "km"
    if pinned.exists():
        return pinned
    which_bin = shutil.which("km")
    if which_bin:
        return Path(which_bin)
    local_bin = Path.home() / ".local" / "bin" / "km"
    if local_bin.exists():
        return local_bin
    raise FileNotFoundError("Kobayashi-MaRust executable 'km' not found in KM_BIN, PATH or ~/.local/bin/km")


def _tbox_lines() -> list[str]:
    lines = [
        f"Prefix(:=<{RAK}>)",
        f"Prefix(rak:=<{RAK}>)",
        f"Prefix(kegg:=<{KEGG}>)",
        f"Prefix(obo:=<{OBO}>)",
        f"Prefix(sio:=<{SIO}>)",
        f"Ontology(<{RAK}eq_biodt_v03.owl>",
        "Declaration(ObjectProperty(rak:hasKO))",
        "Declaration(ObjectProperty(rak:hasMarkerGuild))",
        "Declaration(ObjectProperty(rak:hasEndpoint))",
        "Declaration(Class(rak:ReferenceMAG))",
        "Declaration(Class(rak:BacterialMAG))",
        "Declaration(Class(rak:ArchaealMAG))",
        "DisjointClasses(rak:BacterialMAG rak:ArchaealMAG)",
        "Declaration(Class(rak:Quality90_5MAG))",
        "Declaration(Class(rak:QuantityObservation))",
    ]
    for ep in ENDPOINT_CLASSES:
        lines.append(f"Declaration(Class(rak:{ep}))")
    lines.append("DisjointClasses(" + " ".join(f"rak:{ep}" for ep in ENDPOINT_CLASSES) + ")")

    # Verified taxonomic, DNA-partition, power and trace-gas guild hierarchy from rubalkhali_twin_constraints.owl
    lines.extend([
        "Declaration(Class(rak:RAK_TAXON_TOTAL_PROKARYOTA))",
        "Declaration(Class(obo:NCBITaxon_201174))",
        "Declaration(Class(obo:NCBITaxon_2283796))",
        "Declaration(Class(obo:NCBITaxon_3100166))",
        "Declaration(Class(rak:RAK_TAXON_OTHER_PROKARYOTA))",
        "SubClassOf(obo:NCBITaxon_201174 rak:RAK_TAXON_TOTAL_PROKARYOTA)",
        "SubClassOf(obo:NCBITaxon_2283796 rak:RAK_TAXON_TOTAL_PROKARYOTA)",
        "SubClassOf(obo:NCBITaxon_3100166 rak:RAK_TAXON_TOTAL_PROKARYOTA)",
        "SubClassOf(rak:RAK_TAXON_OTHER_PROKARYOTA rak:RAK_TAXON_TOTAL_PROKARYOTA)",
        "DisjointUnion(rak:RAK_TAXON_TOTAL_PROKARYOTA obo:NCBITaxon_201174 obo:NCBITaxon_2283796 obo:NCBITaxon_3100166 rak:RAK_TAXON_OTHER_PROKARYOTA)",
        "Declaration(Class(rak:RAK_DNA_TOTAL_SOIL_FRACTION))",
        "Declaration(Class(rak:RAK_DNA_VIABLE_FRACTION))",
        "Declaration(Class(rak:RAK_DNA_RELIC_FRACTION))",
        "SubClassOf(rak:RAK_DNA_VIABLE_FRACTION rak:RAK_DNA_TOTAL_SOIL_FRACTION)",
        "SubClassOf(rak:RAK_DNA_RELIC_FRACTION rak:RAK_DNA_TOTAL_SOIL_FRACTION)",
        "DisjointUnion(rak:RAK_DNA_TOTAL_SOIL_FRACTION rak:RAK_DNA_VIABLE_FRACTION rak:RAK_DNA_RELIC_FRACTION)",
        "Declaration(Class(rak:RAK_FLUX_TOTAL_TRACE_GAS_POWER))",
        "Declaration(Class(rak:RAK_FLUX_H2_OXIDATION_POWER))",
        "Declaration(Class(rak:RAK_FLUX_CO_OXIDATION_POWER))",
        "SubClassOf(rak:RAK_FLUX_H2_OXIDATION_POWER rak:RAK_FLUX_TOTAL_TRACE_GAS_POWER)",
        "SubClassOf(rak:RAK_FLUX_CO_OXIDATION_POWER rak:RAK_FLUX_TOTAL_TRACE_GAS_POWER)",
        "DisjointUnion(rak:RAK_FLUX_TOTAL_TRACE_GAS_POWER rak:RAK_FLUX_H2_OXIDATION_POWER rak:RAK_FLUX_CO_OXIDATION_POWER)",
        "Declaration(Class(rak:RAK_GUILD_TOTAL_TRACE_GAS))",
        "Declaration(Class(rak:RAK_GUILD_TOTAL_H2ASE))",
        "Declaration(Class(rak:RAK_GUILD_HYDB_1L))",
        "Declaration(Class(rak:RAK_GUILD_HYDB_1H))",
        "Declaration(Class(rak:RAK_GUILD_HYDB_2A))",
        "Declaration(Class(rak:RAK_GUILD_COXL_FORMI))",
        "SubClassOf(rak:RAK_GUILD_HYDB_1L rak:RAK_GUILD_TOTAL_H2ASE)",
        "SubClassOf(rak:RAK_GUILD_HYDB_1H rak:RAK_GUILD_TOTAL_H2ASE)",
        "SubClassOf(rak:RAK_GUILD_HYDB_2A rak:RAK_GUILD_TOTAL_H2ASE)",
        "SubClassOf(rak:RAK_GUILD_TOTAL_H2ASE rak:RAK_GUILD_TOTAL_TRACE_GAS)",
        "SubClassOf(rak:RAK_GUILD_COXL_FORMI rak:RAK_GUILD_TOTAL_TRACE_GAS)",
        "DisjointUnion(rak:RAK_GUILD_TOTAL_H2ASE rak:RAK_GUILD_HYDB_1L rak:RAK_GUILD_HYDB_1H rak:RAK_GUILD_HYDB_2A)",
        "DisjointUnion(rak:RAK_GUILD_TOTAL_TRACE_GAS rak:RAK_GUILD_TOTAL_H2ASE rak:RAK_GUILD_COXL_FORMI)",
        "Declaration(Class(rak:H2CarrierMAG))",
        "EquivalentClasses(rak:H2CarrierMAG ObjectSomeValuesFrom(rak:hasMarkerGuild rak:RAK_GUILD_TOTAL_H2ASE))",
        "Declaration(Class(rak:COCarrierMAG))",
        "EquivalentClasses(rak:COCarrierMAG ObjectSomeValuesFrom(rak:hasMarkerGuild rak:RAK_GUILD_COXL_FORMI))",
    ])
    for ko in ALL_CBB_KOS:
        lines.append(f"Declaration(Class(kegg:{ko}))")
    for step, alts in CBB_STEPS.items():
        lines.append(f"Declaration(Class(rak:StepKO_{step}))")
        lines.append(f"Declaration(Class(rak:CBBStep_{step}))")
        for ko in alts:
            lines.append(f"SubClassOf(kegg:{ko} rak:StepKO_{step})")
        lines.append(f"EquivalentClasses(rak:CBBStep_{step} ObjectSomeValuesFrom(rak:hasKO rak:StepKO_{step}))")

    step_classes = " ".join(f"rak:CBBStep_{s}" for s in CBB_STEPS)
    lines.extend([
        "Declaration(Class(rak:RuBisCOMAG))",
        "EquivalentClasses(rak:RuBisCOMAG rak:CBBStep_carboxylation)",
        "Declaration(Class(rak:ArchaealRuBisCOMAG))",
        "EquivalentClasses(rak:ArchaealRuBisCOMAG ObjectIntersectionOf(rak:ArchaealMAG rak:RuBisCOMAG))",
        "Declaration(Class(rak:RuBisCOPRKMAG))",
        "EquivalentClasses(rak:RuBisCOPRKMAG ObjectIntersectionOf(rak:CBBStep_carboxylation rak:CBBStep_phosphoribulokinase))",
        "Declaration(Class(rak:ThreeGeneCBBMAG))",
        "EquivalentClasses(rak:ThreeGeneCBBMAG ObjectIntersectionOf(rak:CBBStep_carboxylation rak:CBBStep_phosphoribulokinase ObjectSomeValuesFrom(rak:hasKO kegg:K01602)))",
        "Declaration(Class(rak:CBBCompleteMAG))",
        f"EquivalentClasses(rak:CBBCompleteMAG ObjectIntersectionOf({step_classes}))",
        "Declaration(Class(rak:H2AutotrophicCandidateMAG))",
        "EquivalentClasses(rak:H2AutotrophicCandidateMAG ObjectIntersectionOf(rak:H2CarrierMAG rak:CBBCompleteMAG))",
        "Declaration(Class(rak:HighQualityH2AutotrophicCandidateMAG))",
        "EquivalentClasses(rak:HighQualityH2AutotrophicCandidateMAG ObjectIntersectionOf(rak:H2AutotrophicCandidateMAG rak:Quality90_5MAG))",
        "Declaration(Class(rak:DualH2COAutotrophicCandidateMAG))",
        "EquivalentClasses(rak:DualH2COAutotrophicCandidateMAG ObjectIntersectionOf(rak:H2AutotrophicCandidateMAG ObjectSomeValuesFrom(rak:hasMarkerGuild rak:RAK_GUILD_COXL_FORMI)))",
    ])
    return lines


def _mag_axiom_lines(row) -> list[str]:
    m_iri = f"rak:MAG_{row.rep}"
    lines = [
        f"Declaration(Class({m_iri}))",
        f"SubClassOf({m_iri} rak:ReferenceMAG)",
    ]
    domain = getattr(row, "domain", "")
    if domain == "Bacteria":
        lines.append(f"SubClassOf({m_iri} rak:BacterialMAG)")
    elif domain == "Archaea":
        lines.append(f"SubClassOf({m_iri} rak:ArchaealMAG)")
    phylum = getattr(row, "phylum", "")
    if phylum == "Actinomycetota":
        lines.append(f"SubClassOf({m_iri} obo:NCBITaxon_201174)")
    elif phylum == "Thermoplasmatota":
        lines.append(f"SubClassOf({m_iri} obo:NCBITaxon_2283796)")
    elif phylum == "Halobacteriota":
        lines.append(f"SubClassOf({m_iri} obo:NCBITaxon_3100166)")
    elif phylum:
        lines.append(f"SubClassOf({m_iri} rak:RAK_TAXON_OTHER_PROKARYOTA)")
    cbb_kos_str = getattr(row, "cbb_kos", "")
    if isinstance(cbb_kos_str, str) and cbb_kos_str:
        for ko in cbb_kos_str.split(";"):
            lines.append(f"SubClassOf({m_iri} ObjectSomeValuesFrom(rak:hasKO kegg:{ko}))")
    if row.h2_1l:
        lines.append(f"SubClassOf({m_iri} ObjectSomeValuesFrom(rak:hasMarkerGuild rak:RAK_GUILD_HYDB_1L))")
    if row.h2_1h:
        lines.append(f"SubClassOf({m_iri} ObjectSomeValuesFrom(rak:hasMarkerGuild rak:RAK_GUILD_HYDB_1H))")
    if row.h2_2a:
        lines.append(f"SubClassOf({m_iri} ObjectSomeValuesFrom(rak:hasMarkerGuild rak:RAK_GUILD_HYDB_2A))")
    if row.co_form1:
        lines.append(f"SubClassOf({m_iri} ObjectSomeValuesFrom(rak:hasMarkerGuild rak:RAK_GUILD_COXL_FORMI))")
    if row.quality_90_5:
        lines.append(f"SubClassOf({m_iri} rak:Quality90_5MAG)")
    return lines


def build_owl_ontology(mag: pd.DataFrame, output_path: Path) -> Path:
    """Write the self-contained OWL 2 Functional Syntax ontology for all reference MAGs."""
    lines = _tbox_lines()
    for row in mag.itertuples():
        lines.extend(_mag_axiom_lines(row))
    lines.append(")")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines) + "\n")
    return output_path


def classify_ontology(ofn_path: Path) -> dict:
    """Classify an OWL 2 ontology with Kobayashi-MaRust (km)."""
    km = find_km()
    t0 = time.perf_counter()
    res = subprocess.run([str(km), "classify", str(ofn_path)], capture_output=True, text=True, check=True)
    wall_s = time.perf_counter() - t0
    payload = json.loads(res.stdout)
    if not payload.get("consistent", False) or payload.get("unsatisfiable") or payload.get("dropped", 0) != 0:
        raise RuntimeError(f"OWL ontology classification failed consistency/completeness gate: {payload}")
    prof_res = subprocess.run([str(km), "profile", str(ofn_path)], capture_output=True, text=True, check=True)
    profile = json.loads(prof_res.stdout)
    return {"classification": payload, "profile": profile, "wall_seconds": wall_s, "reasoner": str(km)}


def explain_mag_entailment(mag: pd.DataFrame, mag_id: str, target_class: str, tmp_dir: Path) -> dict:
    """Extract a verified, subset-minimal OWL source-axiom justification with km explain."""
    km = find_km()
    row = mag.loc[mag.rep.eq(mag_id)]
    if row.empty:
        raise ValueError(f"Unknown reference MAG: {mag_id}")
    lines = _tbox_lines()
    lines.extend(_mag_axiom_lines(row.iloc[0]))
    lines.append(")")
    focal_path = tmp_dir / f"explain_{mag_id}_{target_class}.ofn"
    focal_path.write_text("\n".join(lines) + "\n")
    res = subprocess.run(
        [
            str(km), "explain", "--max-checks", "2048",
            str(focal_path), "subclass",
            f"{RAK}MAG_{mag_id}", f"{RAK}{target_class}",
        ],
        capture_output=True, text=True, check=True,
    )
    return json.loads(res.stdout)


def explain_endpoint_conflict(endpoint_a: str, endpoint_b: str, tmp_dir: Path) -> dict:
    """Verify via km explain that conflating two distinct observation endpoints is unsatisfiable."""
    km = find_km()
    lines = _tbox_lines()
    lines.extend([
        "Declaration(Class(rak:ConflatedObservationProbe))",
        f"SubClassOf(rak:ConflatedObservationProbe rak:{endpoint_a})",
        f"SubClassOf(rak:ConflatedObservationProbe rak:{endpoint_b})",
        ")",
    ])
    probe_path = tmp_dir / f"conflict_{endpoint_a}_{endpoint_b}.ofn"
    probe_path.write_text("\n".join(lines) + "\n")
    res = subprocess.run(
        [
            str(km), "explain", "--max-checks", "1024",
            str(probe_path), "unsatisfiable",
            f"{RAK}ConflatedObservationProbe",
        ],
        capture_output=True, text=True, check=True,
    )
    return json.loads(res.stdout)


def evaluate_ontology_reasoning(mag: pd.DataFrame, panel: pd.DataFrame, ofn_path: Path) -> dict:
    """Run Kobayashi-MaRust OWL classification, justification checks and semantic ablations."""
    build_owl_ontology(mag, ofn_path)
    run = classify_ontology(ofn_path)
    subsumptions = run["classification"]["subsumptions"]
    profile = run["profile"]

    def entailed_mags(cls_name: str) -> set[str]:
        target_curie = f"rak:{cls_name}"
        target_iri = f"{RAK}{cls_name}"
        prefix_curie = "rak:MAG_"
        prefix_iri = f"{RAK}MAG_"
        out = set()
        for sub, sup in subsumptions:
            if sup in (target_curie, target_iri):
                if sub.startswith(prefix_curie):
                    out.add(sub[len(prefix_curie):])
                elif sub.startswith(prefix_iri):
                    out.add(sub[len(prefix_iri):])
        return out

    entailed = {
        "RuBisCOMAG": entailed_mags("RuBisCOMAG"),
        "ArchaealRuBisCOMAG": entailed_mags("ArchaealRuBisCOMAG"),
        "RuBisCOPRKMAG": entailed_mags("RuBisCOPRKMAG"),
        "ThreeGeneCBBMAG": entailed_mags("ThreeGeneCBBMAG"),
        "CBBCompleteMAG": entailed_mags("CBBCompleteMAG"),
        "H2CarrierMAG": entailed_mags("H2CarrierMAG"),
        "COCarrierMAG": entailed_mags("COCarrierMAG"),
        "H2AutotrophicCandidateMAG": entailed_mags("H2AutotrophicCandidateMAG"),
        "HighQualityH2AutotrophicCandidateMAG": entailed_mags("HighQualityH2AutotrophicCandidateMAG"),
        "DualH2COAutotrophicCandidateMAG": entailed_mags("DualH2COAutotrophicCandidateMAG"),
    }

    expected = {
        "RuBisCOMAG": set(mag.loc[mag.rbcl_detected, "rep"]),
        "ArchaealRuBisCOMAG": set(mag.loc[mag.rbcl_detected & mag.domain.eq("Archaea"), "rep"]),
        "RuBisCOPRKMAG": set(mag.loc[mag.rbcl_prk_detected, "rep"]),
        "ThreeGeneCBBMAG": set(mag.loc[mag.three_gene_screen, "rep"]),
        "CBBCompleteMAG": set(mag.loc[mag.cbb_all_steps_detected, "rep"]),
        "H2CarrierMAG": set(mag.loc[mag.h2_carrier_union, "rep"]),
        "COCarrierMAG": set(mag.loc[mag.co_form1, "rep"]),
        "H2AutotrophicCandidateMAG": set(mag.loc[mag.h2_cbb_all_steps, "rep"]),
        "HighQualityH2AutotrophicCandidateMAG": set(mag.loc[mag.h2_cbb_quality_90_5, "rep"]),
        "DualH2COAutotrophicCandidateMAG": set(mag.loc[mag.h2_cbb_all_steps & mag.co_form1, "rep"]),
    }
    for k in expected:
        if entailed[k] != expected[k]:
            raise RuntimeError(f"OWL reasoner mismatch on {k}: {len(entailed[k])} vs {len(expected[k])}")

    # Extract subset-minimal explanation on the first candidate MAG in lexical order
    first_candidate = sorted(entailed["H2AutotrophicCandidateMAG"])[0]
    with tempfile.TemporaryDirectory(prefix="reason-", dir=ofn_path.parent) as work:
        mag_exp = explain_mag_entailment(mag, first_candidate, "H2AutotrophicCandidateMAG", Path(work))
        conflict_exp = explain_endpoint_conflict("GenomicCapacityEndpoint", "MappedReadCompositionEndpoint", Path(work))

    if mag_exp.get("status") != "entailed" or not mag_exp.get("justifications"):
        raise RuntimeError("Failed to extract subset-minimal OWL justification for candidate MAG")
    if conflict_exp.get("status") != "entailed" or not conflict_exp.get("justifications"):
        raise RuntimeError("Failed to extract unsatisfiability justification for endpoint conflict")

    just = mag_exp["justifications"][0]
    conflict_just = conflict_exp["justifications"][0]

    # Semantic ablations & genome-linkage vs marginal interval analysis on the 551 samples
    rbcl_h2_mags = entailed["RuBisCOMAG"] & entailed["H2CarrierMAG"]
    prk_h2_mags = entailed["RuBisCOPRKMAG"] & entailed["H2CarrierMAG"]
    multi_h2_mags = int((mag[["h2_1l", "h2_1h", "h2_2a"]].sum(axis=1) > 1).sum())

    valid = panel.loc[panel.mapped_read_pct.notna()].copy()
    h = valid["h2_carrier_union_share_mapped"].to_numpy(float)
    c = valid["cbb_all_steps_detected_share_mapped"].to_numpy(float)
    linked = valid["h2_cbb_all_steps_share_mapped"].to_numpy(float)
    lower = np.maximum(0.0, h + c - 1.0)
    upper = np.minimum(h, c)
    width = upper - lower

    cbb_df = mag.loc[mag.rep.isin(entailed["CBBCompleteMAG"])]
    cand_df = mag.loc[mag.rep.isin(entailed["H2AutotrophicCandidateMAG"])]
    arch_rbcl_df = mag.loc[mag.rep.isin(entailed["ArchaealRuBisCOMAG"])]

    compartment_summary = {}
    for comp, grp in panel.groupby("compartment"):
        ge_s = grp["total_h2_marker_ge"].dropna()
        rbcl_s = grp["rbcl_detected_share_mapped"].dropna() * 100.0
        cbb_s = grp["cbb_all_steps_detected_share_mapped"].dropna() * 100.0
        h2_s = grp["h2_carrier_union_share_mapped"].dropna() * 100.0
        h2cbb_s = grp["h2_cbb_all_steps_share_mapped"].dropna() * 100.0
        compartment_summary[str(comp)] = {
            "n_samples": int(len(grp)),
            "total_h2_marker_ge_median": round(float(ge_s.median()), 4),
            "total_h2_marker_ge_iqr": [round(float(ge_s.quantile(0.25)), 4), round(float(ge_s.quantile(0.75)), 4)],
            "h2_carrier_mapped_pct_median": round(float(h2_s.median()), 2),
            "rbcl_single_marker_mapped_pct_median": round(float(rbcl_s.median()), 2),
            "cbb_11_step_mapped_pct_median": round(float(cbb_s.median()), 2),
            "h2_cbb_candidate_mapped_pct_median": round(float(h2cbb_s.median()), 3),
            "h2_cbb_candidate_mapped_pct_max": round(float(h2cbb_s.max()), 2),
            "samples_exceeding_1pct_h2_cbb": int((h2cbb_s > 1.0).sum()),
        }

    return {
        "reasoner": "Kobayashi-MaRust (km)",
        "reasoner_version": mag_exp.get("reasonerVersion"),
        "reasoner_binary_sha256": sha256(find_km()),
        "ontology_sha256": sha256(ofn_path),
        "selected_route": profile.get("selected_route"),
        "expressivity": profile.get("profile", {}).get("expressivity", {}).get("code"),
        "ontology_file": str(ofn_path.name),
        "logical_axioms": profile.get("profile", {}).get("source", {}).get("logical_axioms"),
        "declared_classes": profile.get("profile", {}).get("source", {}).get("declared_classes"),
        "consistent": True,
        "unsatisfiable_classes": 0,
        "dropped_clauses": 0,
        "entailed_class_counts": {k: len(v) for k, v in entailed.items()},
        "exact_agreement_all_5438_mags": True,
        "taxonomic_architecture": {
            "archaeal_rubisco_detected_mags": len(entailed["ArchaealRuBisCOMAG"]),
            "archaeal_rubisco_by_phylum": {k: int(v) for k, v in arch_rbcl_df["phylum"].value_counts().items()},
            "cbb_complete_by_phylum": {k: int(v) for k, v in cbb_df["phylum"].value_counts().items()},
            "h2_autotrophic_candidates_by_phylum": {k: int(v) for k, v in cand_df["phylum"].value_counts().items()},
            "h2_autotrophic_candidates_by_order": {k: int(v) for k, v in cand_df["order"].value_counts().items()},
            "h2_autotrophic_candidates_by_genus": {k: int(v) for k, v in cand_df["genus"].replace("", "unclassified").value_counts().items()},
            "h2_autotrophic_hydrogenase_subtypes": {
                "group_1l": int(cand_df["h2_1l"].sum()),
                "group_1h": int(cand_df["h2_1h"].sum()),
                "group_2a": int(cand_df["h2_2a"].sum()),
                "dual_h2_co_form1": len(entailed["DualH2COAutotrophicCandidateMAG"]),
            },
        },
        "compartment_niche_partitioning": compartment_summary,
        "example_candidate_justification": {
            "mag_id": first_candidate,
            "superClass": f"{RAK}H2AutotrophicCandidateMAG",
            "verified": bool(just["verified"]),
            "subset_minimal": bool(just["subsetMinimal"]),
            "axiom_count": len(just["axioms"]),
            "classification_checks": mag_exp["classificationChecks"],
            "axioms": [a["functionalSyntax"] for a in just["axioms"]],
        },
        "endpoint_disjointness_justification": {
            "probe_class": f"{RAK}ConflatedObservationProbe",
            "endpoints": ["GenomicCapacityEndpoint", "MappedReadCompositionEndpoint"],
            "status": conflict_exp["status"],
            "verified": bool(conflict_just["verified"]),
            "subset_minimal": bool(conflict_just["subsetMinimal"]),
            "axiom_count": len(conflict_just["axioms"]),
            "axioms": [a["functionalSyntax"] for a in conflict_just["axioms"]],
        },
        "ablations": {
            "single_marker_rubisco_without_11_step_conjunction": {
                "entailed_cbb_mags": len(entailed["RuBisCOMAG"]),
                "additional_incomplete_screen_mags": len(entailed["RuBisCOMAG"] - entailed["CBBCompleteMAG"]),
                "entailed_h2_cbb_mags": len(rbcl_h2_mags),
                "additional_h2_incomplete_screen_mags": len(rbcl_h2_mags - entailed["H2AutotrophicCandidateMAG"]),
            },
            "two_gene_rubisco_prk_without_11_step_conjunction": {
                "entailed_cbb_mags": len(entailed["RuBisCOPRKMAG"]),
                "additional_incomplete_screen_mags": len(entailed["RuBisCOPRKMAG"] - entailed["CBBCompleteMAG"]),
                "entailed_h2_cbb_mags": len(prk_h2_mags),
                "additional_h2_incomplete_screen_mags": len(prk_h2_mags - entailed["H2AutotrophicCandidateMAG"]),
            },
            "additive_hydrogenase_guilds_without_carrier_union": {
                "multi_hydrogenase_mags_double_counted": multi_h2_mags,
            },
            "sample_marginals_without_genome_linkage": {
                "usable_mapped_samples": int(len(valid)),
                "missing_mapped_samples": int(panel.mapped_read_pct.isna().sum()),
                "nonzero_ambiguity_samples": int((width > 1e-12).sum()),
                "frechet_width_median": float(np.median(width)),
                "frechet_width_max": float(np.max(width)),
                "linked_intersection_median": float(np.median(linked)),
                "linked_intersection_max": float(np.max(linked)),
                "linked_within_frechet_bounds_all_samples": bool(
                    np.all((linked >= lower - 1e-12) & (linked <= upper + 1e-12))
                ),
            },
        },
        "entailed_candidate_mags": sorted(entailed["H2AutotrophicCandidateMAG"]),
        "entailed_high_quality_candidate_mags": sorted(entailed["HighQualityH2AutotrophicCandidateMAG"]),
        "entailed_memberships": {str(rep): sorted(cls for cls, ids in entailed.items() if rep in ids)
                                  for rep in mag.rep},
    }


def rules_digest() -> str:
    return hashlib.sha256("\n".join(_tbox_lines()).encode()).hexdigest()


def export_classification(report: dict, data: Path) -> None:
    """Publish actual reasoner outputs tied to their exact evidence and TBox."""
    artifact = {"schema_version": 1, "rules_sha256": rules_digest(),
                "evidence_sha256": sha256(data / "mag_evidence.tsv"),
                "ontology_sha256": report["ontology_sha256"],
                "reasoner_binary_sha256": report["reasoner_binary_sha256"],
                "classes_by_mag": report["entailed_memberships"]}
    (data / "entailed_memberships.json").write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    (data / "ontology_report.json").write_text(json.dumps(report, indent=2) + "\n")
    manifest_path = data / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    for name in ("entailed_memberships.json", "ontology_report.json"):
        manifest["artifacts"][name] = sha256(data / name)
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")


def load_classification(data: Path) -> dict:
    manifest = json.loads((data / "manifest.json").read_text())
    for name in ("entailed_memberships.json", "ontology_report.json"):
        if sha256(data / name) != manifest["artifacts"].get(name):
            raise ValueError("OWL artifact checksum mismatch")
    artifact = json.loads((data / "entailed_memberships.json").read_text())
    if artifact["rules_sha256"] != rules_digest() or artifact["evidence_sha256"] != sha256(data / "mag_evidence.tsv"):
        raise ValueError("Stale OWL classification: rebuild evidence and run the reasoner")
    return artifact
