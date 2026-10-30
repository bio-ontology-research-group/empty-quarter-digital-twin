"""Conditional, one-layer water/salt and diffusion-reaction model (SI units).

Kinetic/retention parameters are scenario assumptions until locally measured.
No metagenomic count is interpreted as biomass, activity, or a rate constant.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, replace
from math import isfinite, log, sqrt, tanh

from . import VERSION


@dataclass(frozen=True)
class Parameters:
    k_h2_s: float                       # bulk-soil first-order reaction coefficient
    k_co_s: float
    porosity: float = .4
    depth_m: float = .1
    psi_entry_mpa: float = -.001         # Campbell retention assumption
    retention_b: float = 4.
    dry_threshold_mpa: float = -85.
    wet_activation_mpa: float = -.1
    activation_shape: float = 1.
    q10: float = 2.
    reference_temp_c: float = 25.
    d_h2_m2_s: float = 7.7e-5
    d_co_m2_s: float = 2.1e-5
    h2_mol_m3: float = 2.13e-5
    co_mol_m3: float = 4.02e-6
    oxygen_supply_mol_m2_s: float = .001
    energy_capture_fraction: float = .25
    maintenance_w_m2: float | None = None
    field_capacity_fraction: float = .5
    drainage_fraction_day: float = .5
    osmotic_coefficient: float = 1.
    osmotic_domain_mol_l: float = .1
    mechanism: str = "activation_diffusion"
    integration_steps_per_day: int = 96
    gamma_h2: float = 1.                 # community factors from ontology memberships (theory.py, community.py)
    gamma_co: float = 1.

    def __post_init__(self):
        for k, v in asdict(self).items():
            if isinstance(v, (int, float)) and not isfinite(v):
                raise ValueError(f"Non-finite parameter: {k}")
        if self.gamma_h2 < 0 or self.gamma_co < 0:
            raise ValueError("Community factors gamma must be non-negative")
        if not 0 < self.porosity < 1 or self.depth_m <= 0:
            raise ValueError("Invalid porosity/depth")
        if not self.dry_threshold_mpa < self.wet_activation_mpa < self.psi_entry_mpa < 0:
            raise ValueError("Require dry < activation < entry < 0 MPa")
        if min(self.retention_b, self.activation_shape, self.q10,
               self.d_h2_m2_s, self.d_co_m2_s, self.osmotic_coefficient,
               self.osmotic_domain_mol_l) <= 0:
            raise ValueError("Invalid positive process parameter")
        for k in ("k_h2_s", "k_co_s", "h2_mol_m3", "co_mol_m3",
                  "oxygen_supply_mol_m2_s"):
            if getattr(self, k) < 0:
                raise ValueError(f"Negative parameter: {k}")
        for k in ("energy_capture_fraction", "field_capacity_fraction", "drainage_fraction_day"):
            if not 0 <= getattr(self, k) <= 1:
                raise ValueError(f"Fraction outside [0,1]: {k}")
        if self.maintenance_w_m2 is not None and self.maintenance_w_m2 <= 0:
            raise ValueError("Maintenance demand must be measured/specified and positive")
        if self.mechanism not in {"activation_diffusion", "diffusion_only"}:
            raise ValueError("Unknown competing mechanism")
        if isinstance(self.integration_steps_per_day, bool) or not isinstance(self.integration_steps_per_day, int) or not 1 <= self.integration_steps_per_day <= 2880:
            raise ValueError("integration_steps_per_day must be an integer in [1,2880]")


def water_potential(theta: float, p: Parameters) -> float | None:
    if not isfinite(theta) or not 0 <= theta <= p.porosity:
        raise ValueError("Water content outside physical pore volume")
    return p.psi_entry_mpa * (theta/p.porosity)**(-p.retention_b) if theta else None


def activation(psi_mpa: float | None, p: Parameters) -> float:
    if p.mechanism == "diffusion_only":
        return 1.                        # explicit falsifiable alternative
    if psi_mpa is None or psi_mpa <= p.dry_threshold_mpa:
        return 0.
    if not isfinite(psi_mpa) or psi_mpa > 0:
        raise ValueError("Invalid water potential")
    if psi_mpa >= p.wet_activation_mpa:
        return 1.
    return (log(abs(p.dry_threshold_mpa/psi_mpa)) /
            log(abs(p.dry_threshold_mpa/p.wet_activation_mpa)))**p.activation_shape


def deposition_flux(concentration: float, diffusivity: float, k: float, depth: float) -> float:
    """Solve D C'' - k C = 0, C(0)=C_atm, C'(L)=0.

    J = C_atm sqrt(D k) tanh(L sqrt(k/D)), mol m-2 s-1.
    k multiplies gas-phase concentration per bulk volume, not GE.
    """
    if any(not isfinite(x) or x < 0 for x in (concentration, diffusivity, k, depth)):
        raise ValueError("Invalid diffusion-reaction input")
    if min(concentration, diffusivity, k, depth) == 0:
        return 0.
    return concentration * sqrt(diffusivity*k) * tanh(depth*sqrt(k/diffusivity))


def gas_exchange(theta: float, temp_c: float, p: Parameters, *,
                 osmotic_mpa: float = 0., measured_matric_mpa: float | None = None) -> dict:
    if not isfinite(temp_c) or not 0 <= temp_c <= 50:
        raise ValueError("Q10 scenario domain is 0–50 C; extrapolation needs a thermal model")
    psi = water_potential(theta, p)
    if measured_matric_mpa is not None:
        if not isfinite(measured_matric_mpa) or measured_matric_mpa > 0:
            raise ValueError("Invalid measured matric potential")
        psi = measured_matric_mpa
    if not isfinite(osmotic_mpa) or osmotic_mpa > 0:
        raise ValueError("Osmotic potential must be non-positive")
    total_psi = psi + osmotic_mpa if psi is not None else None
    a = activation(total_psi, p)
    mq = (p.porosity-theta)**(10/3)/p.porosity**2
    thermal = p.q10**((temp_c-p.reference_temp_c)/10)
    h2 = deposition_flux(p.h2_mol_m3, p.d_h2_m2_s*mq, p.k_h2_s*p.gamma_h2*a*thermal, p.depth_m)
    co = deposition_flux(p.co_mol_m3, p.d_co_m2_s*mq, p.k_co_s*p.gamma_co*a*thermal, p.depth_m)
    demand = .5*(h2+co)
    fraction = min(1., p.oxygen_supply_mol_m2_s/demand) if demand else 1.
    h2 *= fraction
    co *= fraction
    gross = h2*202000 + co*243000
    captured = gross*p.energy_capture_fraction
    return {"matric_potential_mpa": psi, "total_water_potential_mpa": total_psi,
            "activation": a, "gas_diffusivity_ratio": mq,
            "h2_mol_m2_s": h2, "co_mol_m2_s": co,
            "oxygen_used_mol_m2_s": .5*(h2+co),
            "chemical_power_w_m2": gross, "captured_power_w_m2": captured,
            "maintenance_coverage": captured/p.maintenance_w_m2 if p.maintenance_w_m2 else None,
            "status": "conditional_process_simulation"}


@dataclass(frozen=True)
class State:
    water_mm: float
    salt_mol_m2: float = 0.             # neutral NaCl formula units; total inventory

    def __post_init__(self):
        if any(not isfinite(v) or v < 0 for v in asdict(self).values()):
            raise ValueError("Invalid water/salt inventory")


def simulate(initial: State, forcing: list[dict], p: Parameters, state_overrides: dict | None = None) -> dict:
    """Daily, well-mixed bucket; water in mm equals L m-2.

    NaCl is a specified idealised solute, not inferred from EC. At >0.1 M the
    salt ledger continues but uptake abstains pending non-ideal speciation.
    state_overrides maps a day to a function State -> State applied to the state entering
    that day (an intervention on the causal state, do(W_t=w) / do(N_t=n)).
    """
    capacity = 1000*p.depth_m*p.porosity
    if initial.water_mm > capacity:
        raise ValueError("Initial water exceeds pore storage")
    state = initial
    rows = []
    overridden = {}
    for day, f in enumerate(forcing):
        if state_overrides and day in state_overrides:
            state = state_overrides[day](state)
            if state.water_mm > capacity:
                raise ValueError("State intervention exceeds pore storage")
            overridden[day] = asdict(state)
        allowed = {"water_mm", "nacl_mol_l", "evaporation_mm", "temp_c"}
        if set(f)-allowed:
            raise ValueError(f"Unsupported causal levers: {sorted(set(f)-allowed)}")
        if "temp_c" not in f or "evaporation_mm" not in f:
            raise ValueError("Each day requires temperature and prescribed evaporation")
        if any(not isfinite(v) for v in f.values()):
            raise ValueError("Non-finite forcing")
        rain, salt_c, evap = (f.get("water_mm", 0.), f.get("nacl_mol_l", 0.), f["evaporation_mm"])
        if min(rain, salt_c, evap) < 0 or not 0 <= f["temp_c"] <= 50:
            raise ValueError("Invalid forcing")
        mixed_water = state.water_mm + rain
        mixed_salt = state.salt_mol_m2 + rain*salt_c
        overflow = max(0., mixed_water-capacity)
        outwater = overflow
        outsalt = mixed_salt*overflow/mixed_water if mixed_water else 0.
        water, salt = mixed_water-overflow, mixed_salt-outsalt
        evaporation = 0.
        steps = p.integration_steps_per_day
        drain_fraction = 1-(1-p.drainage_fraction_day)**(1/steps)
        gas_sum = None
        undefined = set()
        chemistry_ok = True
        for _ in range(steps):
            previous_water, previous_salt = water, salt
            drainage = max(0., water-capacity*p.field_capacity_fraction)*drain_fraction
            salt_loss = salt*drainage/water if water else 0.
            loss = min(evap/steps, water-drainage)
            water = max(0., water-drainage-loss)
            salt = max(0., salt-salt_loss)
            outwater += drainage
            outsalt += salt_loss
            evaporation += loss
            # Midpoint quadrature resolves within-day activation and drying.
            mid_water, mid_salt = (water+previous_water)/2, (salt+previous_salt)/2
            mid_c = mid_salt/mid_water if mid_water else None
            end_c = salt/water if water else None
            supported = ((mid_salt == 0 or (mid_c is not None and mid_c <= p.osmotic_domain_mol_l))
                         and (salt == 0 or (end_c is not None and end_c <= p.osmotic_domain_mol_l)))
            chemistry_ok &= supported
            if supported:
                osmotic = -2*p.osmotic_coefficient*(mid_c or 0)*1000*8.314462618*(f["temp_c"]+273.15)/1e6
                sample = gas_exchange(mid_water/(1000*p.depth_m), f["temp_c"], p, osmotic_mpa=osmotic)
                if gas_sum is None:
                    gas_sum = {k: 0. for k,v in sample.items() if isinstance(v, (int, float))}
                undefined.update(k for k,v in sample.items() if v is None)
                for key in gas_sum:
                    if sample[key] is not None:
                        gas_sum[key] += sample[key]/steps
        next_state = State(water, salt)
        concentration = next_state.salt_mol_m2/next_state.water_mm if next_state.water_mm else None
        gas = ({**(gas_sum or {}), "status": "conditional_process_simulation",
                "aggregation": "subdaily_midpoint_mean", "maintenance_coverage":
                gas_sum.get("maintenance_coverage") if gas_sum else None} if chemistry_ok else None)
        if gas is not None:
            for key in undefined:
                gas[key] = None
        rows.append({"day": day, **asdict(next_state), "water_input_mm": rain,
                     "salt_input_mol_m2": rain*salt_c, "outflow_mm": outwater,
                     "salt_outflow_mol_m2": outsalt, "evaporation_mm": evaporation,
                     "water_balance_residual_mm": state.water_mm+rain-outwater-evaporation-next_state.water_mm,
                     "salt_balance_residual_mol_m2": state.salt_mol_m2+rain*salt_c-outsalt-next_state.salt_mol_m2,
                     "nacl_mol_l": concentration, "gas": gas,
                     "chemistry_status": "ideal_dilute_assumption" if chemistry_ok else "requires_nonideal_speciation"})
        if day in overridden:
            rows[-1]["state_override"] = overridden[day]
        state = next_state
    valid = all(r["gas"] is not None for r in rows)
    return {"version": VERSION, "parameters": asdict(p), "initial": asdict(initial),
            "forcing": forcing, "trajectory": rows,
            "summary": {"water_input_mm": sum(r["water_input_mm"] for r in rows),
                        "h2_uptake_mol_m2": sum(r["gas"]["h2_mol_m2_s"]*86400 for r in rows) if valid else None,
                        "final_water_mm": state.water_mm, "final_salt_mol_m2": state.salt_mol_m2,
                        "net_co2_removal_g_m2": None,
                        "net_carbon_status": "requires_stock_change_source_and_lifecycle_measurements"},
            "unsupported": ["species_specific_growth", "cloud_seeding_efficacy", "community_turnover",
                            "nitrogen_cycling", "sulfate_reduction", "mineral_precipitation", "net_carbon_removal"]}


def compact(p: Parameters, porosity: float) -> Parameters:
    """Prescribe measured post-disturbance porosity; dose-to-porosity is external."""
    return replace(p, porosity=porosity)
