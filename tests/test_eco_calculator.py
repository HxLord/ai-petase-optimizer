"""
Unit tests for src.eco_calculator.
"""

import pytest
from src.eco_calculator import calculate_scenario, create_scenario_figure


def test_mass_balance_arithmetic():
    """Verify basic mass balance arithmetic."""
    res = calculate_scenario(
        pet_waste_kg=1000.0,
        baseline_rate_pct=15.0,
        hypothetical_rate_pct=75.0
    )
    assert res["baseline_recovered_kg"] == 150.0
    assert res["baseline_unrecovered_kg"] == 850.0
    assert res["hypothetical_recovered_kg"] == 750.0
    assert res["hypothetical_unrecovered_kg"] == 250.0
    assert res["mass_difference_kg"] == 600.0
    assert res["rate_difference_pct"] == 60.0
    assert res["emissions_scenario"] is None


def test_boundary_percentages():
    """Verify boundary values 0% and 100%."""
    res_zero = calculate_scenario(pet_waste_kg=500.0, baseline_rate_pct=0.0, hypothetical_rate_pct=0.0)
    assert res_zero["baseline_recovered_kg"] == 0.0
    assert res_zero["hypothetical_recovered_kg"] == 0.0

    res_hundred = calculate_scenario(pet_waste_kg=500.0, baseline_rate_pct=100.0, hypothetical_rate_pct=100.0)
    assert res_hundred["baseline_recovered_kg"] == 500.0
    assert res_hundred["hypothetical_recovered_kg"] == 500.0


def test_invalid_percentage_bounds():
    """Verify rates < 0% or > 100% raise ValueError."""
    with pytest.raises(ValueError, match="between 0% and 100%"):
        calculate_scenario(pet_waste_kg=100.0, baseline_rate_pct=-5.0, hypothetical_rate_pct=50.0)

    with pytest.raises(ValueError, match="between 0% and 100%"):
        calculate_scenario(pet_waste_kg=100.0, baseline_rate_pct=20.0, hypothetical_rate_pct=105.0)


def test_invalid_waste_mass():
    """Verify negative waste mass raises ValueError."""
    with pytest.raises(ValueError, match="non-negative number"):
        calculate_scenario(pet_waste_kg=-10.0, baseline_rate_pct=10.0, hypothetical_rate_pct=20.0)


def test_processing_time_scenarios():
    """Verify processing time calculations and errors."""
    res = calculate_scenario(
        pet_waste_kg=1000.0,
        baseline_rate_pct=20.0,
        hypothetical_rate_pct=60.0,
        baseline_time_hours=48.0,
        hypothetical_time_hours=12.0
    )
    time_info = res["time_scenario"]
    assert time_info is not None
    assert time_info["time_difference_hours"] == 36.0
    assert time_info["time_ratio"] == 0.25

    # Non-positive time
    with pytest.raises(ValueError, match="must be > 0"):
        calculate_scenario(pet_waste_kg=100.0, baseline_rate_pct=10.0, hypothetical_rate_pct=20.0, baseline_time_hours=0.0)


def test_emissions_factor_calculation():
    """Verify emissions scenario when factor is provided vs when None."""
    # Factor provided: 1.25 kg CO2e / kg rPET
    res = calculate_scenario(
        pet_waste_kg=1000.0,
        baseline_rate_pct=10.0,   # 100 kg
        hypothetical_rate_pct=50.0, # 500 kg -> +400 kg
        emissions_factor_kg_co2_per_kg=1.25
    )
    emissions = res["emissions_scenario"]
    assert emissions is not None
    assert emissions["net_emissions_difference_kg_co2"] == 500.0  # 400 * 1.25
    assert emissions["baseline_avoided_kg_co2"] == 125.0
    assert emissions["hypothetical_avoided_kg_co2"] == 625.0
    assert emissions["vehicle_miles_offset_equivalent"] > 0


def test_negative_emissions_factor_raises_error():
    """Verify negative emissions factor raises ValueError."""
    with pytest.raises(ValueError, match="Emissions factor must be non-negative"):
        calculate_scenario(
            pet_waste_kg=100.0,
            baseline_rate_pct=10.0,
            hypothetical_rate_pct=20.0,
            emissions_factor_kg_co2_per_kg=-1.0
        )


def test_non_finite_inputs_are_rejected():
    with pytest.raises(ValueError):
        calculate_scenario(float("nan"), 10.0, 20.0)
    with pytest.raises(ValueError):
        calculate_scenario(100.0, float("inf"), 20.0)
    with pytest.raises(ValueError):
        calculate_scenario(100.0, 10.0, 20.0, emissions_factor_kg_co2_per_kg=float("nan"))


def test_negative_ghg_difference_has_signed_vehicle_equivalent():
    result = calculate_scenario(
        pet_waste_kg=100.0,
        baseline_rate_pct=80.0,
        hypothetical_rate_pct=20.0,
        emissions_factor_kg_co2_per_kg=1.0,
    )
    assert result["emissions_scenario"]["net_emissions_difference_kg_co2"] == -60.0
    assert result["emissions_scenario"]["vehicle_miles_offset_equivalent"] < 0


def test_create_scenario_figure():
    """Verify Plotly figure generation for eco scenario."""
    res = calculate_scenario(pet_waste_kg=1000.0, baseline_rate_pct=20.0, hypothetical_rate_pct=80.0)
    fig = create_scenario_figure(res)
    assert fig is not None
    assert len(fig.data) == 2  # Recovered and Unrecovered bars
