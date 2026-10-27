"""
Environmental Scenario Calculator Module for AI-PETase Optimizer.
Models material mass balances, hypothetical recovery differences, and estimated
greenhouse gas emission scenarios based strictly on explicit user assumptions.
"""

import math
from numbers import Real
from typing import Dict, Any, Optional
import pandas as pd
import plotly.graph_objects as go


def calculate_scenario(
    pet_waste_kg: float,
    baseline_rate_pct: float,
    hypothetical_rate_pct: float,
    baseline_time_hours: Optional[float] = None,
    hypothetical_time_hours: Optional[float] = None,
    emissions_factor_kg_co2_per_kg: Optional[float] = None
) -> Dict[str, Any]:
    """
    Compute hypothetical PET recycling and environmental scenario metrics.

    Args:
        pet_waste_kg: Total PET plastic waste mass in kilograms (>= 0).
        baseline_rate_pct: Baseline recycling/recovery percentage (0.0 to 100.0).
        hypothetical_rate_pct: Alternative hypothetical recovery percentage (0.0 to 100.0).
        baseline_time_hours: Optional baseline batch processing time in hours (> 0).
        hypothetical_time_hours: Optional alternative processing time in hours (> 0).
        emissions_factor_kg_co2_per_kg: Optional avoided net GHG emissions per kg recycled PET.

    Returns:
        dict: Mass balances, differences, optional time and emission deltas, and disclaimers.

    Raises:
        ValueError: If input values are negative, percentages outside [0, 100], or times non-positive.
    """
    # 1. Input validations
    if (
        isinstance(pet_waste_kg, bool)
        or not isinstance(pet_waste_kg, Real)
        or not math.isfinite(pet_waste_kg)
        or pet_waste_kg < 0
    ):
        raise ValueError(f"PET waste mass must be a non-negative number, got {pet_waste_kg}.")

    if (
        isinstance(baseline_rate_pct, bool)
        or not isinstance(baseline_rate_pct, Real)
        or not math.isfinite(baseline_rate_pct)
        or not (0.0 <= baseline_rate_pct <= 100.0)
    ):
        raise ValueError(
            f"Baseline recovery rate must be between 0% and 100%, got {baseline_rate_pct}%."
        )

    if (
        isinstance(hypothetical_rate_pct, bool)
        or not isinstance(hypothetical_rate_pct, Real)
        or not math.isfinite(hypothetical_rate_pct)
        or not (0.0 <= hypothetical_rate_pct <= 100.0)
    ):
        raise ValueError(
            f"Hypothetical recovery rate must be between 0% and 100%, got {hypothetical_rate_pct}%."
        )

    # 2. Mass balances
    baseline_recovered_kg = pet_waste_kg * (baseline_rate_pct / 100.0)
    baseline_unrecovered_kg = pet_waste_kg - baseline_recovered_kg

    hypothetical_recovered_kg = pet_waste_kg * (hypothetical_rate_pct / 100.0)
    hypothetical_unrecovered_kg = pet_waste_kg - hypothetical_recovered_kg

    mass_difference_kg = hypothetical_recovered_kg - baseline_recovered_kg

    # 3. Processing time scenario (optional)
    time_data = None
    if baseline_time_hours is not None or hypothetical_time_hours is not None:
        if baseline_time_hours is not None and (
            isinstance(baseline_time_hours, bool)
            or not isinstance(baseline_time_hours, Real)
            or not math.isfinite(baseline_time_hours)
            or baseline_time_hours <= 0
        ):
            raise ValueError(f"Baseline processing time must be > 0, got {baseline_time_hours}.")
        if hypothetical_time_hours is not None and (
            isinstance(hypothetical_time_hours, bool)
            or not isinstance(hypothetical_time_hours, Real)
            or not math.isfinite(hypothetical_time_hours)
            or hypothetical_time_hours <= 0
        ):
            raise ValueError(f"Hypothetical processing time must be > 0, got {hypothetical_time_hours}.")

        if baseline_time_hours is not None and hypothetical_time_hours is not None:
            time_diff_hours = baseline_time_hours - hypothetical_time_hours
            time_ratio = (
                hypothetical_time_hours / baseline_time_hours
                if baseline_time_hours > 0 else None
            )
            time_data = {
                "baseline_time_hours": baseline_time_hours,
                "hypothetical_time_hours": hypothetical_time_hours,
                "time_difference_hours": time_diff_hours,
                "time_ratio": round(time_ratio, 3) if time_ratio else None
            }

    # 4. Emissions scenario (optional)
    emissions_data = None
    if emissions_factor_kg_co2_per_kg is not None:
        if (
            isinstance(emissions_factor_kg_co2_per_kg, bool)
            or not isinstance(emissions_factor_kg_co2_per_kg, Real)
            or not math.isfinite(emissions_factor_kg_co2_per_kg)
            or emissions_factor_kg_co2_per_kg < 0
        ):
            raise ValueError(
                f"Emissions factor must be non-negative, got {emissions_factor_kg_co2_per_kg}."
            )
        baseline_avoided_co2 = baseline_recovered_kg * emissions_factor_kg_co2_per_kg
        hypothetical_avoided_co2 = hypothetical_recovered_kg * emissions_factor_kg_co2_per_kg
        net_emissions_diff_kg_co2 = mass_difference_kg * emissions_factor_kg_co2_per_kg

        # US EPA standard conversion: ~0.404 kg CO2e per average passenger vehicle mile
        vehicle_miles_offset = net_emissions_diff_kg_co2 / 0.404

        emissions_data = {
            "factor_used_kg_co2_per_kg": emissions_factor_kg_co2_per_kg,
            "baseline_avoided_kg_co2": round(baseline_avoided_co2, 2),
            "hypothetical_avoided_kg_co2": round(hypothetical_avoided_co2, 2),
            "net_emissions_difference_kg_co2": round(net_emissions_diff_kg_co2, 2),
            "vehicle_miles_offset_equivalent": round(vehicle_miles_offset, 1)
        }

    # 5. Scientific scenario disclaimer
    disclaimer = (
        "SCENARIO DISCLAIMER: All outputs represent hypothetical mathematical calculations "
        "derived strictly from user-configured parameters. They do NOT represent experimentally "
        "verified degradation performance, commercial yields, or real plant LCA outcomes. "
        "Polymer crystallinity, enzyme inactivation, and byproduct inhibition must be evaluated experimentally."
    )

    return {
        "pet_waste_kg": round(pet_waste_kg, 2),
        "baseline_rate_pct": round(baseline_rate_pct, 2),
        "hypothetical_rate_pct": round(hypothetical_rate_pct, 2),
        "rate_difference_pct": round(hypothetical_rate_pct - baseline_rate_pct, 2),
        "baseline_recovered_kg": round(baseline_recovered_kg, 2),
        "baseline_unrecovered_kg": round(baseline_unrecovered_kg, 2),
        "hypothetical_recovered_kg": round(hypothetical_recovered_kg, 2),
        "hypothetical_unrecovered_kg": round(hypothetical_unrecovered_kg, 2),
        "mass_difference_kg": round(mass_difference_kg, 2),
        "time_scenario": time_data,
        "emissions_scenario": emissions_data,
        "disclaimer": disclaimer,
    }


def create_scenario_figure(results: Dict[str, Any]) -> go.Figure:
    """
    Generate a Plotly stacked bar chart illustrating material mass allocation.
    
    Args:
        results: Dictionary produced by calculate_scenario.
        
    Returns:
        go.Figure: Stacked bar chart comparing Baseline vs Hypothetical recovery.
    """
    scenarios = ["Baseline Practice", "Hypothetical Scenario"]
    recovered = [results["baseline_recovered_kg"], results["hypothetical_recovered_kg"]]
    unrecovered = [results["baseline_unrecovered_kg"], results["hypothetical_unrecovered_kg"]]

    fig = go.Figure()

    # Recovered PET
    fig.add_trace(go.Bar(
        name="Recovered / Recycled PET (kg)",
        x=scenarios,
        y=recovered,
        text=[f"{v:,.1f} kg" for v in recovered],
        textposition="inside",
        marker=dict(color="#168575", line=dict(color="#ffffff", width=1.5)),
        hovertemplate="<b>%{x}</b><br>Recovered: %{y:,.1f} kg<extra></extra>"
    ))

    # Unrecovered PET (Residual)
    fig.add_trace(go.Bar(
        name="Residual / Unrecovered PET (kg)",
        x=scenarios,
        y=unrecovered,
        text=[f"{v:,.1f} kg" for v in unrecovered],
        textposition="inside",
        marker=dict(color="#c27a54", line=dict(color="#ffffff", width=1.5)),
        hovertemplate="<b>%{x}</b><br>Residual: %{y:,.1f} kg<extra></extra>"
    ))

    fig.update_layout(
        barmode="stack",
        title="<b>Mass balance: baseline and hypothetical scenario</b>",
        xaxis_title="Scenario",
        yaxis_title="PET Plastic Mass (kg)",
        template="plotly_white",
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="#ffffff",
        font=dict(family="DM Sans, sans-serif", size=12, color="#243b35"),
        title_font=dict(family="Newsreader, Georgia, serif", size=20, color="#172b29"),
        xaxis=dict(showgrid=False, linecolor="#cbd8d0", tickfont=dict(color="#405851")),
        yaxis=dict(gridcolor="#e5ebe7", zerolinecolor="#cbd8d0", tickfont=dict(color="#405851")),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        margin=dict(l=40, r=40, t=80, b=40)
    )

    return fig
