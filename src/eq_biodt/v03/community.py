"""Community factors gamma_g from reasoner memberships and mapped-read weights.

q_s(C) = sum_{g in C} w_sg / sum_g w_sg over mapped reference MAGs, where C is a class whose
members are entailed by the OWL reasoner (data/v03/entailed_memberships.json), not read from
table flags. These are deduced quantities: they have no causal ancestors and no intervention
changes them. Samples with missing mapping weights are excluded, not treated as zero.
"""
from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

from .evidence import sha256, source_directory
from .ontology import load_classification

DATA = Path(__file__).resolve().parents[3] / "data" / "v03"
if not DATA.exists():
    DATA = Path(__file__).parent / "data"
MODES = ("reference", "gate", "proportional")
GAS_CLASS = {"h2": "H2CarrierMAG", "co": "COCarrierMAG"}
CLASSES = ("H2CarrierMAG", "COCarrierMAG", "H2AutotrophicCandidateMAG", "CBBCompleteMAG",
           "HighQualityH2AutotrophicCandidateMAG", "DualH2COAutotrophicCandidateMAG",
           "ThreeGeneCBBMAG", "RuBisCOPRKMAG", "RuBisCOMAG", "ArchaealRuBisCOMAG")
# Entailed subsumptions (sub, super) from the TBox definitions in ontology.py.
SUBSUMPTIONS = (("H2AutotrophicCandidateMAG", "H2CarrierMAG"), ("H2AutotrophicCandidateMAG", "CBBCompleteMAG"),
                ("HighQualityH2AutotrophicCandidateMAG", "H2AutotrophicCandidateMAG"),
                ("DualH2COAutotrophicCandidateMAG", "H2AutotrophicCandidateMAG"),
                ("DualH2COAutotrophicCandidateMAG", "COCarrierMAG"),
                ("CBBCompleteMAG", "RuBisCOPRKMAG"), ("ThreeGeneCBBMAG", "RuBisCOPRKMAG"),
                ("RuBisCOPRKMAG", "RuBisCOMAG"), ("ArchaealRuBisCOMAG", "RuBisCOMAG"))
TOL = 1e-12


@lru_cache(maxsize=1)
def _load(data: Path = DATA) -> dict:
    classified = load_classification(data)
    sources = source_directory()
    entry = __import__("json").loads((sources / "manifest.json").read_text())["sources"]["read_mapping"]
    path = sources / entry["file"]
    digest = sha256(path)
    if digest != entry["sha256"]:
        raise ValueError("Read-mapping source digest mismatch")
    cover = pd.read_csv(path, sep="\t", names=["sample", "rep", "pct"])
    bad = cover.loc[cover.pct.isna(), "sample"].unique()
    cover = cover.loc[~cover["sample"].isin(bad) & cover.rep.ne("unmapped")]
    total = cover.groupby("sample").pct.sum()
    classes = classified["classes_by_mag"]
    shares = {}
    for c in CLASSES:
        members = {m for m, v in classes.items() if c in v}
        shares[c] = cover.pct.where(cover.rep.isin(members), 0.).groupby(cover["sample"]).sum()/total.replace(0, np.nan)
    panel = pd.read_csv(data / "sample_evidence.tsv", sep="\t", usecols=["sample", "site", "compartment", "trip_num"])
    table = pd.DataFrame(shares).dropna()
    table = table.loc[table.index.isin(panel["sample"])]       # the 551 evidence-panel samples
    medians = {c: float(table[c].median()) for c in CLASSES}
    return {"table": table, "panel": panel, "medians": medians, "classes": classes,
            "evidence_sha256": classified["evidence_sha256"], "ontology_sha256": classified["ontology_sha256"],
            "read_mapping_sha256": digest, "excluded_samples": sorted(bad.tolist())}


def medians() -> dict[str, float]:
    return dict(_load()["medians"])


def samples_of(community: dict) -> list[str]:
    """Resolve {"sample": id} or {"site","compartment","trip"} to replicate sample ids."""
    c = _load()
    if not isinstance(community, dict):
        raise ValueError("community must be an object")
    if "sample" in community and set(community) == {"sample"}:
        ids = [str(community["sample"])]
    elif set(community) == {"site", "compartment", "trip"}:
        p = c["panel"]
        hit = p.loc[p.site.astype(str).eq(str(community["site"])) & p.compartment.eq(community["compartment"])
                    & p.trip_num.astype(str).eq(str(community["trip"])), "sample"]
        ids = sorted(hit.tolist())
        if not ids:
            raise ValueError(f"No samples for visit {community}")
    else:
        raise ValueError("community must be {'sample': id} or {'site':..., 'compartment':..., 'trip':...}")
    missing = [i for i in ids if i not in c["table"].index]
    if missing:
        raise ValueError(f"No usable read-mapping weights (unknown sample or excluded for missing weights): {missing}")
    return ids


def carrier_shares(sample_id: str | None = None, *, community: dict | None = None) -> dict:
    """Mean class shares over the replicate samples of a sample or visit, with provenance."""
    ids = samples_of(community if community is not None else {"sample": sample_id})
    c = _load()
    q = c["table"].loc[ids].mean()
    return {"samples": ids, "shares": {k: float(q[k]) for k in CLASSES},
            "evidence_sha256": c["evidence_sha256"], "ontology_sha256": c["ontology_sha256"],
            "read_mapping_sha256": c["read_mapping_sha256"]}


def gamma(shares: dict, mode: str) -> dict[str, float]:
    """gamma_h2, gamma_co from class shares under a coupling mode."""
    if mode not in MODES:
        raise ValueError(f"Unknown coupling mode {mode!r}; choose from {list(MODES)}")
    out = {}
    for gas, cls in GAS_CLASS.items():
        if mode == "reference":
            out[gas] = 1.
            continue
        q = shares[cls]
        if not np.isfinite(q) or q < 0:
            raise ValueError(f"Invalid share for {cls}")
        if mode == "gate":
            out[gas] = 1. if q > 0 else 0.
        else:
            med = medians()[cls]
            if med <= 0:
                raise ValueError(f"Median share of {cls} is zero; proportional coupling undefined")
            out[gas] = q/med
    return out


def lifted_consistency(shares: dict | None = None) -> dict:
    """q(sub) <= q(super) for every entailed subsumption; with no argument, over all samples."""
    c = _load()
    checks = {}
    for sub, sup in SUBSUMPTIONS:
        if shares is None:
            checks[f"{sub}<={sup}"] = bool((c["table"][sub] <= c["table"][sup] + TOL).all())
        else:
            checks[f"{sub}<={sup}"] = bool(shares[sub] <= shares[sup] + TOL)
    if shares is None:
        checks.update(membership_consistency(c["classes"]))
    return checks


def membership_consistency(classes: dict | None = None) -> dict:
    """Every MAG entailed in a subclass is entailed in its superclass (reasoner output only)."""
    classes = classes if classes is not None else _load()["classes"]
    return {f"members({sub})subset({sup})": all(sup in v for v in classes.values() if sub in v)
            for sub, sup in SUBSUMPTIONS}


def closed_gates() -> dict:
    """Counts of samples whose H2 / CO gate is closed (zero carrier share)."""
    t = _load()["table"]
    return {"samples": int(len(t)),
            "h2_closed": t.index[t["H2CarrierMAG"] <= 0].tolist(),
            "co_closed": t.index[t["COCarrierMAG"] <= 0].tolist()}


def describe(sample_id: str) -> dict:
    s = carrier_shares(sample_id)
    return {**s, "gamma": {m: gamma(s["shares"], m) for m in MODES}, "medians": medians(),
            "consistency": lifted_consistency(s["shares"]),
            "mapping_check": "shares from entailed memberships and CoverM mapped weights"}
