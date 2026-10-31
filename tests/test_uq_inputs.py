from eq_biodt.v03.process import Parameters
from eq_biodt.v03.uq_inputs import Input, inputs_for, load_parameter_ranges


def test_registry_has_required_parameters():
    R = load_parameter_ranges()
    assert {"k_h2_s","k_co_s","psi_entry_mpa","retention_b","wet_activation_mpa","dry_threshold_mpa",
            "activation_shape","q10","d_h2_m2_s","d_co_m2_s","drainage_fraction_day",
            "field_capacity_fraction","porosity"} <= set(R)
    for e in R.values():
        assert e["status"] in {"verified_in_source","modeller_assumption","source_not_checked"}
        assert e["low"] <= e["value"] <= e["high"]
        assert e["value"] == getattr(Parameters(k_h2_s=.001, k_co_s=.001), e["id"])


def test_existing_ensemble_bounds_reused():
    R = load_parameter_ranges()
    assert (R["k_h2_s"]["dist"], R["k_h2_s"]["low"], R["k_h2_s"]["high"]) == ("loguniform", 1e-4, 1e-2)
    assert (R["retention_b"]["low"], R["retention_b"]["high"]) == (2.0, 6.0)
    assert (R["dry_threshold_mpa"]["low"], R["dry_threshold_mpa"]["high"]) == (-100.0, -70.0)
    assert (R["activation_shape"]["low"], R["activation_shape"]["high"]) == (0.5, 2.0)
    assert (R["porosity"]["low"], R["porosity"]["high"]) == (0.35, 0.5)
    assert (R["psi_entry_mpa"]["dist"], R["psi_entry_mpa"]["low"], R["psi_entry_mpa"]["high"]) == ("loguniform", -0.01, -0.0001)
    assert (R["q10"]["low"], R["q10"]["high"]) == (1.0, 2.0)


def test_inputs_for_excludes_intervened_parameter():
    names = [i.name for i in inputs_for(["fao29_ece_per_ecw"], exclude_params={"porosity"})]
    assert "porosity" not in names and "fao29_ece_per_ecw" in names and names == sorted(names)
