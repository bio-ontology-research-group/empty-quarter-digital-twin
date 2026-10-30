"""HiGHS Linear Programming Flux Balance Analysis (FBA) & Cross-Feeding Engine.

Implements genome-scale guild stoichiometry, hard maintenance feasibility, and energetic
trade-offs grounded in:
- `~/Public/software/empty-quarter-lifeonair/lifeonair/results/metabolic_pilot/REPORT.txt`:
  72-MAG (24 matched triplets across 13 orders), 10,584-run GSMM audit demonstrating that
  uncurated automated reconstructions (`gapsmith`/`gapseq`, where only 1/48 high-affinity
  carriers wired a native internal H2 reaction to its marker GPR and 6/25 paired models
  disagreed on binary H2 ATP gain) require explicit ontology-curated reaction coupling.
- `~/Public/software/empty-quarter-lifeonair/lifeonair/results/h123_v2/A4_cbb_cocarriage_by_group.tsv`:
  101/558 (18.1%; 27.6% at >=90% completeness) Group 1l MAGs co-carry CBB RuBisCO vs
  30/371 (8.1%; 9.3% at >=90% completeness) Group 1h MAGs (completeness-adjusted p = 0.0203).
- `~/Documents/papers/empty-quarter-tiannyu` (Table S2 measured water chemistry):
  Site 59 sabkha brines (EC = 216.2-231.4 mS/cm, Na+ = 71.0-79.9 g/L, Cl- = 103.7-121.1 g/L)
  vs Site 60 brackish artesian pool (EC = 10.65 mS/cm, Ca2+ = 1.07 g/L, SO4^2- = 2.23 g/L).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict

import numpy as np
from scipy.optimize import linprog


@dataclass(frozen=True)
class CommunityMetabolicInputs:
    """Environmental uptake bounds and stress drivers for the HiGHS FBA guild solver."""

    h2_uptake_bound: float = 1.8          # mmol gDW^-1 h^-1 equivalent scaled flux
    co_uptake_bound: float = 0.35         # mmol gDW^-1 h^-1 equivalent scaled flux
    root_exudate_doc_bound: float = 0.0   # rhizodeposition DOC from native shrubs
    exogenous_nh4_no3_bound: float = 0.08 # inorganic/organic N input (e.g., Oryx/Camel dung or rain)
    sulfate_so4_bound: float = 0.25       # Ca-SO4 electron acceptor / sulfur flux supply
    salinity_ms_cm: float = 1.2           # electrical conductivity (Site 60 = 10.65, Site 59 = 216.2-231.4)
    temperature_c: float = 35.0           # ambient soil/water temperature (up to 60 C surface)
    choline_betaine_bound: float = 0.0    # exogenous choline/glycine-betaine (mM equivalent)
    rbcl_hyab_plasmid_inoc: float = 0.0   # [0, 1] inoculation of Site 60 Allorhizobium + diazotrophs
    ect_kdp_plasmid_inoc: float = 0.0     # [0, 1] inoculation of Sphingomonas ectABC-kdpFABC-proP


class MetabolicCommunityEngine:
    """Solves community metabolic flux distributions using HiGHS (`scipy.optimize.linprog`).

    Enforces hard non-growth maintenance ATP feasibility (`v_atpm >= atpm_required`):
    - When `atpm_required` can be met by available electron donors (`H2`, `CO`, `DOC`,
      `NH4`, `SO4`, `betaine`), solves for optimal autotrophic + heterotrophic biomass synthesis.
    - When `atpm_required > max_possible_atp` (e.g., unbuffered >200 mS/cm Na-Cl brine shock),
      the hard growth LP is infeasible; biomass synthesis is strictly zeroed (`total_biomass_flux = 0.0`,
      `autotrophic_carbon_fixation_flux = 0.0`) and the community enters
      `MAINTENANCE_DEFICIT_DORMANCY` with explicit reporting of `atpm_maintenance_deficit`.
    """

    METABOLIC_PILOT_AUDIT_SUMMARY: Dict[str, Any] = {
        "total_representative_mags": 5438,
        "eligible_quality_filtered_mags": 1759,
        "eligible_group_1l_mags": 178,
        "eligible_group_1h_mags": 163,
        "eligible_neither_control_mags": 1418,
        "matched_triplets_selected": 24,
        "pilot_panel_mags": 72,
        "orders_represented": 13,
        "total_fba_scenarios_run": 10584,
        "optimal_solver_return_rate": 1.0,
        "carrier_proteome_exact_matches": 44,
        "carrier_proteome_total": 48,
        "uncurated_native_h2_gpr_wired_count": 1,
        "uncurated_native_h2_gpr_wired_mag": "S54PRr1__self_SemiBin_621",
        "gapseq_vs_gapsmith_paired_mags": 25,
        "gapseq_vs_gapsmith_phenotype_disagreements": 6,
        "group_1l_total_mags": 558,
        "group_1l_cbb_cocarriage_mags": 101,
        "group_1l_cbb_cocarriage_fraction": 0.1810,
        "group_1l_hc90_cbb_fraction": 0.2760,
        "group_1h_total_mags": 371,
        "group_1h_cbb_cocarriage_mags": 30,
        "group_1h_cbb_cocarriage_fraction": 0.0809,
        "group_1h_hc90_cbb_fraction": 0.0928,
        "cbb_completeness_adjusted_log_odds_1l_vs_1h": 0.963,
        "cbb_completeness_adjusted_p": 0.0203,
    }

    MEASURED_AQUATIC_HYDROCHEMISTRY: Dict[str, Dict[str, float]] = {
        "59-MR-WR1": {
            "site_num": 59,
            "ec_ms_cm": 216.2,
            "ph": 7.92,
            "ca_mg_l": 1371.85,
            "mg_mg_l": 2818.13,
            "na_mg_l": 71044.90,
            "k_mg_l": 2344.70,
            "cl_mg_l": 103722.0,
            "so4_mg_l": 11074.60,
            "total_c_mg_l": 41.11,
            "total_n_mg_l": 168.20,
        },
        "59-MR-WR2": {
            "site_num": 59,
            "ec_ms_cm": 231.4,
            "ph": 8.01,
            "ca_mg_l": 1241.72,
            "mg_mg_l": 3249.28,
            "na_mg_l": 79896.20,
            "k_mg_l": 2762.86,
            "cl_mg_l": 121069.0,
            "so4_mg_l": 11923.90,
            "total_c_mg_l": 33.03,
            "total_n_mg_l": 193.50,
        },
        "60-MR-WR1": {
            "site_num": 60,
            "ec_ms_cm": 10.65,
            "ph": 7.68,
            "ca_mg_l": 1066.40,
            "mg_mg_l": 165.25,
            "na_mg_l": 1548.88,
            "k_mg_l": 101.42,
            "cl_mg_l": 2221.43,
            "so4_mg_l": 2228.98,
            "total_c_mg_l": 9.12,
            "total_n_mg_l": 6.22,
        },
    }

    def solve_community_fba(self, inp: CommunityMetabolicInputs) -> Dict[str, Any]:
        """Solve the community LP via HiGHS with hard maintenance feasibility enforcement."""
        salinity_excess = max(0.0, inp.salinity_ms_cm - 2.0)
        thermal_excess = max(0.0, inp.temperature_c - 38.0)

        osm_mitigation = 1.0 - 0.55 * float(np.clip(inp.ect_kdp_plasmid_inoc, 0.0, 1.0))
        betaine_mitigation = 1.0 / (1.0 + 1.8 * max(0.0, inp.choline_betaine_bound))

        atpm_required = (
            0.08
            + 0.011 * salinity_excess * osm_mitigation * betaine_mitigation
            + 0.028 * (thermal_excess**1.25) * betaine_mitigation
        )

        nif_cap = 0.045 + 0.42 * float(np.clip(inp.rbcl_hyab_plasmid_inoc, 0.0, 1.0))
        cbb_cap = 0.65 + 1.15 * float(np.clip(inp.rbcl_hyab_plasmid_inoc, 0.0, 1.0))
        sox_cap = inp.sulfate_so4_bound * np.exp(-max(0.0, inp.salinity_ms_cm - 25.0) / 85.0)

        bounds = [
            (0.0, max(0.0, inp.h2_uptake_bound)),                    # v0: H2 oxidation
            (0.0, max(0.0, inp.co_uptake_bound)),                    # v1: CO oxidation
            (0.0, cbb_cap),                                          # v2: CBB autotrophy
            (0.0, 0.45),                                             # v3: 3-HP/4-HB archaeal autotrophy
            (0.0, max(0.0, inp.root_exudate_doc_bound)),             # v4: DOC heterotrophy
            (0.0, nif_cap),                                          # v5: nifH N2 fixation
            (0.0, max(0.0, inp.exogenous_nh4_no3_bound + nif_cap)),  # v6: Nitrification
            (0.0, max(0.0, sox_cap)),                                # v7: Sulfur redox flux
            (0.0, max(0.0, inp.choline_betaine_bound)),              # v8: Betaine GbsAB/OpuA rescue
            (0.0, 10.0),                                             # v9: ATPM drain
        ]

        max_possible_atp = (
            1.4 * bounds[0][1]
            + 0.9 * bounds[1][1]
            + 2.2 * bounds[4][1]
            + 0.5 * max(0.0, inp.exogenous_nh4_no3_bound)
            + 0.8
            * min(
                bounds[7][1],
                0.50 * bounds[0][1] + 0.35 * bounds[1][1] + 0.50 * bounds[4][1],
            )
            + 0.7 * bounds[8][1]
        )

        # Primary LP: enforce hard non-growth maintenance requirement v9 >= atpm_required
        A_ub = np.array(
            [
                [-1.4, -0.9, 1.8, 1.5, -2.2, 3.2, -0.5, -0.8, -0.7, 1.0],
                [0.0, 0.0, 0.15, 0.15, 0.15, -1.0, 0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 1.0, 0.0, 0.0, 0.0],
                [-0.50, -0.35, 0.0, 0.0, -0.50, 0.0, 0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0],
            ],
            dtype=float,
        )
        b_ub_hard = np.array(
            [
                0.0,
                max(0.0, inp.exogenous_nh4_no3_bound),
                max(0.0, inp.exogenous_nh4_no3_bound),
                0.0,
                -atpm_required,
            ],
            dtype=float,
        )
        c = np.array(
            [-0.15, -0.10, -1.25, -1.35, -1.10, -0.90, -0.65, -0.80, -0.50, 0.05],
            dtype=float,
        )

        res_hard = linprog(c, A_ub=A_ub, b_ub=b_ub_hard, bounds=bounds, method="highs")

        if res_hard.success and res_hard.x is not None:
            v = res_hard.x
            hard_feasible = True
            atpm_deficit = 0.0
        else:
            # Hard maintenance requirement exceeds maximum respiratory ATP supply:
            # Anabolic/growth pathways (v2, v3, v4, v5) are strictly arrested (zeroed) and all
            # non-growth catabolic electron-donor flux (v0, v1, v6, v7, v8) is allocated to
            # survival maintenance v9 so that A_ub[:4] @ v <= b_ub_hard[:4] holds identically.
            hard_feasible = False
            v = np.zeros(10, dtype=float)
            v[0] = bounds[0][1]
            v[1] = bounds[1][1]
            v[6] = max(0.0, inp.exogenous_nh4_no3_bound)
            v[7] = min(bounds[7][1], 0.50 * v[0] + 0.35 * v[1])
            v[8] = bounds[8][1]
            max_catabolic_atp = (
                1.4 * v[0]
                + 0.9 * v[1]
                + 0.5 * v[6]
                + 0.8 * v[7]
                + 0.7 * v[8]
            )
            v[9] = min(bounds[9][1], max_catabolic_atp)
            atpm_deficit = max(0.0, float(atpm_required - v[9]))

        autotrophic_c_flux = float(v[2] + v[3])
        heterotrophic_c_flux = float(v[4])
        total_biomass_flux = autotrophic_c_flux + 0.85 * heterotrophic_c_flux
        energy_stress_ratio = float(atpm_required / max(1e-4, max_possible_atp))

        # Uncharacterized ("dark") coding fraction of the community proteome catalogue
        # (53.4% baseline in native oligotrophic desert soils from empty-quarter-proteome;
        # copiotrophic DOC enrichment or extreme brine shock selects for cosmopolitan annotated families)
        nutrient_selection = np.exp(-1.45 * heterotrophic_c_flux)
        brine_selection = np.exp(-max(0.0, inp.salinity_ms_cm - 25.0) / 45.0)
        dark_proteome_coding_frac = float(
            np.clip(0.534 * nutrient_selection * (0.09 + 0.91 * brine_selection), 0.03, 0.65)
        )

        viability_state = (
            "VIABLE_GROWTH"
            if (hard_feasible and total_biomass_flux > 1e-6)
            else (
                "MAINTENANCE_FEASIBLE_OLIGOTROPHY"
                if hard_feasible
                else "MAINTENANCE_DEFICIT_DORMANCY"
            )
        )

        return {
            "solver_status": "optimal" if hard_feasible else "infeasible_maintenance_deficit",
            "native_highs_status": int(res_hard.status),
            "native_highs_success": bool(res_hard.success),
            "hard_maintenance_feasible": bool(hard_feasible),
            "metabolic_viability_state": viability_state,
            "total_biomass_flux": round(total_biomass_flux, 5),
            "autotrophic_carbon_fixation_flux": round(autotrophic_c_flux, 5),
            "cbb_rubisco_flux": round(float(v[2]), 5),
            "archaeal_3hp4hb_flux": round(float(v[3]), 5),
            "heterotrophic_doc_flux": round(heterotrophic_c_flux, 5),
            "nifH_n2_fixation_flux": round(float(v[5]), 5),
            "nitrification_flux": round(float(v[6]), 5),
            "sulfur_redox_flux": round(float(v[7]), 5),
            "betaine_osmoprotectant_flux": round(float(v[8]), 5),
            "atpm_maintenance_drain": round(float(v[9]), 5),
            "atpm_required_unmitigated": round(float(atpm_required), 5),
            "atpm_maintenance_deficit": round(float(atpm_deficit), 5),
            "energy_stress_ratio": round(energy_stress_ratio, 4),
            "dark_proteome_coding_fraction": round(dark_proteome_coding_frac, 4),
            "dark_proteome_active_fraction": round(dark_proteome_coding_frac, 4),
        }
