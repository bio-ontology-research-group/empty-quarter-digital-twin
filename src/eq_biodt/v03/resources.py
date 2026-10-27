"""Element-balanced resource ceilings and explicit carbon accounting."""
from __future__ import annotations

import math
import numpy as np
from scipy.optimize import linprog

# Species: H2, CO, O2, CO2, H2O, CH2O. Columns: H2 oxidation,
# CO oxidation, carbohydrate synthesis. Neutral species, zero charge.
SPECIES = ("H2", "CO", "O2", "CO2", "H2O", "CH2O")
STOICHIOMETRY = np.array([[-1, 0, -2], [0, -1, 0], [-.5, -.5, 0],
                          [0, 1, -1], [1, 0, 1], [0, 0, 1]], dtype=float)
# Rows C,H,O; degree-of-reduction 4C+H-2O also conserves electrons.
ELEMENTS = np.array([[0, 1, 0, 1, 0, 1], [2, 0, 0, 0, 2, 2],
                     [0, 1, 2, 2, 1, 1]])


def resource_ceiling(*, h2: float, co: float, oxygen: float, external_co2: float,
                     capture_fraction: float, maintenance_w_m2: float,
                     synthesis_j_mol_c: float, pathway_supported: bool) -> dict:
    """Maximise CH2O-C with explicit mol m-2 s-1 substrate supplies.

    This is a three-reaction, open-boundary resource ceiling. Pathway support
    permits a ceiling calculation; it provides no enzyme/biomass rate constraint.
    Maintenance is mandatory, independently specified, and energy is charged once.
    """
    vals = (h2, co, oxygen, external_co2, capture_fraction, maintenance_w_m2, synthesis_j_mol_c)
    if any(not math.isfinite(x) or x < 0 for x in vals) or not 0 <= capture_fraction <= 1 or synthesis_j_mol_c == 0:
        raise ValueError("Invalid resource/energy input")
    if not pathway_supported:
        return {"status": "insufficient_pathway_evidence", "fluxes": None,
                "reason": "An incomplete annotation screen cannot establish a zero biological synthesis capacity."}
    # H2, CO, O2, external inorganic C, captured energy including maintenance.
    A = np.array([[1, 0, 2], [0, 1, 0], [.5, .5, 0], [0, -1, 1],
                  [-202000*capture_fraction, -243000*capture_fraction, synthesis_j_mol_c]])
    b = np.array([h2, co, oxygen, external_co2, -maintenance_w_m2])
    # Scale tiny physiological fluxes to O(1) before the numerical LP.
    scale = max(h2, co, external_co2, maintenance_w_m2/202000, 1e-20)
    row_scale = np.maximum(np.max(np.abs(A), axis=1), 1)
    fit = linprog([0, 0, -1], A_ub=A/row_scale[:, None], b_ub=b/scale/row_scale,
                  bounds=[(0, None)]*3, method="highs")
    if fit.status == 2:
        return {"status": "infeasible_resource_or_maintenance_demand", "fluxes": None}
    if not fit.success:
        raise RuntimeError(f"Resource solver failed ({fit.status}): {fit.message}")
    x = fit.x*scale
    exchange = STOICHIOMETRY @ x
    residual = ELEMENTS @ exchange
    if np.max(np.abs(residual)) > 1e-10*scale or np.max((A@x-b)/row_scale) > 1e-7*scale:
        raise RuntimeError("Resource solution failed balance/feasibility certification")
    return {"status": "conditional_carbohydrate_production_ceiling",
            "fluxes": dict(zip(("h2_respiration", "co_respiration", "carbohydrate_c"), x.tolist())),
            "net_species_exchange_mol_m2_s": dict(zip(SPECIES, exchange.tolist())),
            "element_residuals": residual.tolist(),
            "energy_slack_w_m2": float(-A[-1]@x-maintenance_w_m2),
            "net_carbon_removal": None}


def carbon_account(*, stock_change_g_c_m2: float, imported_carbon_g_c_m2: float,
                   geological_carbon_g_c_m2: float, lifecycle_g_co2e_m2: float,
                   counterfactual_stock_change_g_c_m2: float,
                   non_co2_g_co2e_m2: float, boundary: str, interval: str,
                   evidence: dict[str, str]) -> dict:
    """Measured stock-change boundary, including all ecosystem carbon pools.

    Respiration already enters measured stock change; never subtract it twice.
    Carbon imports and geological sources enter explicitly; permanence is external.
    Negative values for stock changes are allowed, source/lifecycle loads are >=0.
    """
    required = {"stock_change", "imports", "geological_sources", "lifecycle", "counterfactual", "non_co2"}
    if not required <= evidence.keys() or not all(evidence[k].strip() for k in required):
        raise ValueError("Every carbon-accounting term requires a measurement/source reference")
    vals = (stock_change_g_c_m2, imported_carbon_g_c_m2, geological_carbon_g_c_m2,
            lifecycle_g_co2e_m2, counterfactual_stock_change_g_c_m2, non_co2_g_co2e_m2)
    if not boundary or not interval or any(not math.isfinite(x) for x in vals):
        raise ValueError("Finite carbon terms, spatial boundary and time interval required")
    if min(imported_carbon_g_c_m2, geological_carbon_g_c_m2, lifecycle_g_co2e_m2, non_co2_g_co2e_m2) < 0:
        raise ValueError("Source and emissions terms must be nonnegative")
    atmospheric_c = stock_change_g_c_m2-imported_carbon_g_c_m2-geological_carbon_g_c_m2-counterfactual_stock_change_g_c_m2
    return {"boundary": boundary, "interval": interval, "evidence": evidence,
            "additional_atmospheric_stock_gain_g_c_m2": atmospheric_c,
            "net_balance_g_co2e_m2": atmospheric_c*44/12-lifecycle_g_co2e_m2-non_co2_g_co2e_m2,
            "status": "supplied_stock_and_lifecycle_account", "permanence": "not_assessed"}
