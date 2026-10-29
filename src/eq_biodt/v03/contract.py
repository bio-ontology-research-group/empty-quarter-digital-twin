"""Complete declared scientific contract, also rendered into the manuscript."""
from dataclasses import MISSING, fields

from .evidence import CBB_STEPS
from .process import Parameters


# IDs are stable audit anchors, not asserted external ontology identifiers.
CLAUSES = [
    ("O1", "ontological", "A quantity has one endpoint, unit, denominator and provenance. Genomic capacity, mapped-read composition, DNA recovery, membrane integrity, biomass and activity are disjoint endpoints in the OWL 2 ontology."),
    ("O2", "ontological", "Disjoint gene-copy categories with the same normalizer add to their union. Gene-copy totals may exceed one genome equivalent. Organisms carrying several genes overlap and are classified via OWL class unions."),
    ("O3", "ontological", "For a specified mapped MAG universe, a carrier set and its exact complement partition mapped-read weight. Their union equals the mapped total; their intersection is empty."),
    ("O4", "ontological", "A nonnegative part cannot exceed its whole under the same measure. Composition has no minimum other-taxon fraction; zero and one are valid boundaries."),
    ("O5", "observation", "Missing measurements remain missing. Aggregates requiring a missing part are unavailable. No projection changes causal leaves or intervention settings."),
    ("M1", "observation", "MAG KO presence comes from corrected post-refinement annotation v2; quality joins use exact MAG identifiers. Missing pathway steps can reflect incompleteness or annotation failure."),
    ("M2", "observation", "CBB support is defined in OWL 2 as the conjunction of the 11 KEGG M00165 module requirements, including the alternative SBPase and transaldolase routes, where each step is satisfied by existential restriction over its alternative KOs and classified by Kobayashi-MaRust (km). K01602 is optional in M00165; a separate three-gene class requires it."),
    ("M3", "observation", "Read percentages are relative mapping weights. Carrier shares divide by mapped reads, while unmapped reads remain explicit. Any nonfinite mapping observation makes that sample's composition unavailable. Shannon uses mapped MAG proportions."),
    ("M4", "observation", "Marker genome equivalents use the 60-percent identity-floor read screen and single-copy-gene normalization. They measure relative gene-copy capacity, with no absolute biomass conversion."),
    ("M5", "observation", "The strict quality sensitivity subset requires completeness at least 90 percent and contamination below 5 percent; it is a declared numerical screen, not a complete MIMAG high-quality assignment."),
    ("C1", "causal assumption", "The soil is a horizontally homogeneous layer of prescribed depth and porosity. Campbell retention uses matric potential = entry potential times (theta/porosity)^(-b); the parameters require local retention measurements."),
    ("C2", "causal assumption", "Activity rises from zero at the dry water-potential threshold to one at the wet activation threshold as a clipped logarithmic ramp raised to a shape exponent. The -70 to -100 MPa threshold range comes from Reji et al.; the ramp, wet endpoint and local transfer are assumptions."),
    ("C3", "causal assumption", "Gas diffusion follows Millington-Quirk: D = D0 (porosity-theta)^(10/3)/porosity^2. Saturation closes the gas pathway. Water-phase transport, crust resistance and gas production are omitted."),
    ("C4", "causal assumption", "A steady, first-order reaction layer obeys D C'' - k C = 0 with atmospheric surface concentration and zero basal gradient. Uptake is C_atm sqrt(D k) tanh(L sqrt(k/D)). The approximation is for trace concentrations; kinetic saturation is unmodelled."),
    ("C5", "causal assumption", "Reaction coefficients k are independent supplied parameters per bulk-soil volume, multiplied by activation and Q10^((T-Tref)/10). They are not derived from marker GE. Temperature is restricted to 0-50 C; diffusion constants are held fixed."),
    ("C6", "causal assumption", "H2 and CO oxidation each consume half a mole of O2 per mole donor. A supplied O2-flux limit scales both uptake demands proportionally when limiting; this is a shared-resource approximation."),
    ("C7", "causal assumption", "Chemical power uses fixed apparent energy yields of 202 and 243 kJ/mol H2 and CO. Captured power is a supplied fraction. Maintenance coverage is unavailable unless independent area-based demand is supplied; no biomass is inferred."),
    ("C8", "causal assumption", "Each day adds and mixes water and NaCl and removes proportional overflow. Within-day steps drain water above field capacity with proportional salt and evaporate water only. The substep drainage fraction is 1-(1-d)^(1/n), where d is the daily fraction and n defaults to 96. Gas uptake integrates midpoint rates over substeps. Pore capacity is 1000 L porosity mm."),
    ("C9", "causal assumption", "Salt is conservative NaCl formula units, not conductivity-derived salinity. Dilute osmotic potential is -2 phi c R T; the specified coefficient is a sensitivity parameter. Above 0.1 mol/L or in dry saline soil, gas outputs abstain pending nonideal speciation. Salt inventory remains conserved."),
    ("C10", "causal assumption", "Compaction prescribes post-disturbance porosity and leaves the retention parameters fixed unless replaced. Temperature and evaporation are prescribed levers. Repeated pulses share water and salt memory; biomass, vegetation and community turnover are unmodelled."),
    ("C11", "causal assumption", "The resource ceiling has three balanced reactions: H2+0.5O2 -> H2O; CO+0.5O2 -> CO2; CO2+2H2 -> CH2O+H2O. It maximizes carbohydrate production with donor, O2, external CO2 and independently supplied energy/maintenance bounds. Checksum-verified OWL entailments admit the conditional ceiling for H2-plus-CBB candidates and permit CO oxidation only with CO marker entailment. Incomplete screens return unavailable; lack of entailment is not biological absence. Nitrogen, growth and sulfate reduction are outside this network."),
    ("C12", "causal assumption", "Gas exchange is a conditional process calculation. The resource ceiling is a separate optimization over supplied potential flux bounds, not additional realized uptake. They must never be added. Solver infeasibility and numerical failure are separate outcomes."),
    ("A1", "accounting", "Net carbon removal is unavailable from scenario gas uptake. The optional measured account uses total ecosystem stock change minus imported C, geological C and counterfactual stock change, multiplied by 44/12, then subtracts lifecycle and non-CO2 emissions. Respiration already contributes to stock change. Every term requires a source, interval and spatial boundary; permanence remains separate."),
    ("F1", "predictive assumption", "Forecast observations are site-compartment-campaign means of marker GE. Log1p is the target scale. Missing actual collection dates, coordinates or targets are excluded. Curated collection dates replace nominal trip dates; no pooled 16S, pH or DNA imputation enters prediction."),
    ("F2", "predictive assumption", "Ridge alpha=10 uses standardized latitude/longitude, annual sine/cosine day-of-year on their natural [-1,1] scale and one-hot compartment; preprocessing fits on training data only. Baselines are global mean, compartment mean and latest site-compartment observation, with compartment/global fallback for unvisited states."),
    ("F3", "predictive assumption", "Backtests hold out whole later campaigns, using observations before the test campaign begins. Nominal 90-percent empirical intervals use 5-fold site-held-out residuals fitted on training visits. Persistence uses compartment fallback on unseen calibration sites. Dependence and temporal shift preclude a distribution-free coverage guarantee. The corrected 2026 reference catalogue is fixed throughout, so catalogue discovery is outside the temporal holdout."),
    ("F4", "predictive assumption", "A forecast freezes issue date, later target date, scoring window (default +/-14 days), source digest, all point predictions and methods. Existing forecasts cannot be overwritten. Scoring verifies the digest and uses only new observations after issuance inside the target window."),
    ("D1", "design assumption", "The experiment crosses six target water potentials (-150,-85,-30,-3,-0.3,-0.03 MPa) with three porosities (0.3,0.4,0.5). Default three randomized complete blocks use seed 42. All -3 MPa units are withheld; repeated time points are not biological replication."),
    ("D2", "design assumption", "Model disagreement compares the water-activation model and a diffusion-only alternative at 25 C, divided by supplied noise (default 1e-9 mol/m2/s). Illustrative theta uses the assumed retention curve. Replication requires pilot variance; the design has no prospective power guarantee."),
    ("R1", "decision assumption", "The 24-member sensitivity grid crosses retention b in {3,5}, dry threshold in {-70,-100} MPa, activation exponent in {0.5,2}, and kinetic-rate multiplier in {0.1,1,10}. Members are assumptions, not posterior probabilities."),
    ("R2", "decision assumption", "Scenario comparison maximizes conditional H2 uptake and minimizes supplied water, residual salt and prescribed disturbed area. It reports per-member Pareto sets, their intersection and uptake rank reversals. No scalar utility, Shannon desirability or vegetation benefit is assumed."),
]


def contract() -> dict:
    defaults = {f.name: ("REQUIRED" if f.default is MISSING else f.default) for f in fields(Parameters)}
    return {"version": "0.5.0", "method": "ontological-causal reasoning with separate observation and prediction layers",
            "clauses": [{"id": i, "kind": k, "statement": s} for i,k,s in CLAUSES],
            "parameter_defaults": defaults, "CBB_step_alternatives": CBB_STEPS,
            "causal_identification": "Mechanisms are assumed structural equations; intervention effects are not identified by observational metagenomes.",
            "calibration": "No local gas-kinetic, retention, active-biomass or growth-rate calibration is claimed.",
            "sources": {"moisture_threshold": "https://doi.org/10.1038/s41467-025-67208-3",
                        "coupled_transport": "https://doi.org/10.1029/2021GB006987",
                        "gas_diffusivity": "https://doi.org/10.1039/TF9615701200",
                        "CBB_reactions": "https://www.kegg.jp/entry/M00165",
                        "forecasting": "https://doi.org/10.1073/pnas.1710231115"}}
