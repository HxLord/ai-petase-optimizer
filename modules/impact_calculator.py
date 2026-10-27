"""
Eco Impact and Circularity Scenario Calculator Module for AI-PETase Optimizer.
Models material mass balances, hypothetical recovery improvements,
and life cycle greenhouse gas emissions based strictly on user assumptions.
"""

import math
from numbers import Real
from typing import Dict, Any, Optional
import pandas as pd
import plotly.graph_objects as go


def calculate_scenario(
    pet_waste_kg: float,
    baseline_rate_pct: float,
    improvement_pct: float,
    processing_days: int,
    enable_emissions: bool = True,
    virgin_factor_kg_co2: Optional[float] = None,
    recycling_factor_kg_co2: Optional[float] = None
) -> Dict[str, Any]:
    """
    Compute mass balance and emission scenario outcomes based on explicit user assumptions.
    """
    numeric_inputs = {
        "PET waste mass": pet_waste_kg,
        "Baseline recovery rate": baseline_rate_pct,
        "Recovery-rate improvement": improvement_pct,
    }
    for label, value in numeric_inputs.items():
        if isinstance(value, bool) or not isinstance(value, Real) or not math.isfinite(value):
            raise ValueError(f"{label} must be a finite number.")
    if pet_waste_kg < 0:
        raise ValueError("PET waste mass must be non-negative.")
    if not 0 <= baseline_rate_pct <= 100:
        raise ValueError("Baseline recovery rate must be between 0% and 100%.")
    if improvement_pct < 0:
        raise ValueError("Recovery-rate improvement must be non-negative.")
    if isinstance(processing_days, bool) or not isinstance(processing_days, int) or processing_days < 1:
        raise ValueError("Processing days must be a positive integer.")
    if (virgin_factor_kg_co2 is None) != (recycling_factor_kg_co2 is None):
        raise ValueError("Provide both emissions factors or neither.")
    for label, value in (
        ("Virgin PET emissions factor", virgin_factor_kg_co2),
        ("Recycling emissions factor", recycling_factor_kg_co2),
    ):
        if value is not None and (
            isinstance(value, bool)
            or not isinstance(value, Real)
            or not math.isfinite(value)
            or value < 0
        ):
            raise ValueError(f"{label} must be a finite non-negative number.")

    clamped_baseline = baseline_rate_pct
    hypothetical_rate_pct = min(100.0, clamped_baseline + improvement_pct)
    
    # Baseline masses
    baseline_recovered_kg = pet_waste_kg * (clamped_baseline / 100.0)
    baseline_unrecovered_kg = pet_waste_kg - baseline_recovered_kg
    
    # Hypothetical scenario masses
    hypothetical_recovered_kg = pet_waste_kg * (hypothetical_rate_pct / 100.0)
    hypothetical_unrecovered_kg = pet_waste_kg - hypothetical_recovered_kg
    
    # Net delta
    additional_recovered_kg = hypothetical_recovered_kg - baseline_recovered_kg
    pct_increase_in_recovered = (
        (additional_recovered_kg / baseline_recovered_kg * 100.0)
        if baseline_recovered_kg > 0 else 0.0
    )
    
    # Daily processing throughput
    days = processing_days
    daily_recovered_baseline_kg = baseline_recovered_kg / days
    daily_recovered_hypothetical_kg = hypothetical_recovered_kg / days
    
    # Emissions balance
    emissions_data = None
    if enable_emissions and virgin_factor_kg_co2 is not None:
        # Avoided emissions factor per kg rPET displacing virgin PET
        net_avoided_factor = virgin_factor_kg_co2 - recycling_factor_kg_co2
        
        baseline_avoided_kg_co2 = baseline_recovered_kg * net_avoided_factor
        hypo_avoided_kg_co2 = hypothetical_recovered_kg * net_avoided_factor
        net_avoided_difference_kg_co2 = additional_recovered_kg * net_avoided_factor
        
        # Equivalent metrics (US EPA standard: 0.404 kg CO2e per average passenger vehicle mile)
        miles_offset = net_avoided_difference_kg_co2 / 0.404 if net_avoided_difference_kg_co2 > 0 else 0
        barrels_oil_equiv = (net_avoided_difference_kg_co2 / 433.0) if net_avoided_difference_kg_co2 > 0 else 0
        
        emissions_data = {
            "virgin_factor_used": virgin_factor_kg_co2,
            "recycling_factor_used": recycling_factor_kg_co2,
            "net_avoided_factor_per_kg": net_avoided_factor,
            "baseline_avoided_kg_co2": baseline_avoided_kg_co2,
            "hypo_avoided_kg_co2": hypo_avoided_kg_co2,
            "net_avoided_difference_kg_co2": net_avoided_difference_kg_co2,
            "miles_offset_equivalent": miles_offset,
            "barrels_oil_equivalent": barrels_oil_equiv
        }
    
    return {
        "pet_waste_kg": pet_waste_kg,
        "baseline_rate_pct": clamped_baseline,
        "improvement_pct": improvement_pct,
        "hypothetical_rate_pct": hypothetical_rate_pct,
        "processing_days": days,
        "baseline_recovered_kg": baseline_recovered_kg,
        "baseline_unrecovered_kg": baseline_unrecovered_kg,
        "hypothetical_recovered_kg": hypothetical_recovered_kg,
        "hypothetical_unrecovered_kg": hypothetical_unrecovered_kg,
        "additional_recovered_kg": additional_recovered_kg,
        "pct_increase_in_recovered": pct_increase_in_recovered,
        "daily_recovered_baseline_kg": daily_recovered_baseline_kg,
        "daily_recovered_hypothetical_kg": daily_recovered_hypothetical_kg,
        "emissions": emissions_data,
        "disclaimer": (
            "SCENARIO DISCLAIMER: All values represent mathematical projections from user-provided inputs "
            "and illustrative LCA factors. They do not imply that any enzyme variant has achieved this recovery "
            "rate or emissions offset in commercial operations."
        )
    }


def create_scenario_comparison_chart(results: Dict[str, Any]) -> go.Figure:
    """
    Stacked bar chart comparing Baseline vs Hypothetical waste breakdown (Recovered vs Unrecovered).
    """
    scenarios = ["Baseline Practice", "Hypothetical Scenario"]
    recovered = [results["baseline_recovered_kg"], results["hypothetical_recovered_kg"]]
    unrecovered = [results["baseline_unrecovered_kg"], results["hypothetical_unrecovered_kg"]]

    fig = go.Figure()

    fig.add_trace(go.Bar(
        name="Recovered / Recycled PET",
        x=scenarios,
        y=recovered,
        text=[
            f"{r:,.0f} kg ({r / results['pet_waste_kg'] * 100:.1f}%)"
            if results["pet_waste_kg"] > 0 else "0 kg (0.0%)"
            for r in recovered
        ],
        textposition="inside",
        marker_color="#00A896",
        hovertemplate="<b>%{x}</b><br>Recovered: %{y:,.1f} kg<extra></extra>"
    ))

    fig.add_trace(go.Bar(
        name="Unrecovered (Landfill / Thermal Waste)",
        x=scenarios,
        y=unrecovered,
        text=[
            f"{u:,.0f} kg ({u / results['pet_waste_kg'] * 100:.1f}%)"
            if results["pet_waste_kg"] > 0 else "0 kg (0.0%)"
            for u in unrecovered
        ],
        textposition="inside",
        marker_color="#E76F51",
        hovertemplate="<b>%{x}</b><br>Unrecovered: %{y:,.1f} kg<extra></extra>"
    ))

    fig.update_layout(
        barmode="stack",
        title="<b>PET Mass Balance: Baseline vs. Hypothetical Scenario</b>",
        xaxis_title="<b>Scenario</b>",
        yaxis_title="<b>PET Mass (kg)</b>",
        template="plotly_dark",
        paper_bgcolor="#0d1b2a",
        plot_bgcolor="#1b263b",
        font=dict(family="sans-serif", size=12, color="#E0E1DD"),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        margin=dict(l=40, r=40, t=60, b=40),
        height=380
    )

    return fig


def create_timeline_projection_chart(results: Dict[str, Any]) -> go.Figure:
    """
    Area line chart showing cumulative recovery over the designated processing window.
    """
    days = results["processing_days"]
    day_steps = list(range(0, days + 1, max(1, days // 20)))
    if day_steps[-1] != days:
        day_steps.append(days)

    baseline_cum = [step * results["daily_recovered_baseline_kg"] for step in day_steps]
    hypo_cum = [step * results["daily_recovered_hypothetical_kg"] for step in day_steps]

    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=day_steps,
        y=baseline_cum,
        mode="lines",
        name="Baseline Practice",
        line=dict(color="#6C757D", width=2, dash="dash"),
        fill="tozeroy",
        fillcolor="rgba(108, 117, 125, 0.2)",
        hovertemplate="Day %{x}: %{y:,.0f} kg cumulative<extra></extra>"
    ))

    fig.add_trace(go.Scatter(
        x=day_steps,
        y=hypo_cum,
        mode="lines",
        name="Hypothetical Scenario",
        line=dict(color="#02C39A", width=3),
        fill="tonexty",
        fillcolor="rgba(2, 195, 154, 0.25)",
        hovertemplate="Day %{x}: %{y:,.0f} kg cumulative<extra></extra>"
    ))

    fig.update_layout(
        title=f"<b>Cumulative PET Recovered Over {days} Days</b>",
        xaxis_title="<b>Processing Timeline (Days)</b>",
        yaxis_title="<b>Cumulative PET Recovered (kg)</b>",
        template="plotly_dark",
        paper_bgcolor="#0d1b2a",
        plot_bgcolor="#1b263b",
        font=dict(family="sans-serif", size=12, color="#E0E1DD"),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        margin=dict(l=40, r=40, t=60, b=40),
        height=380
    )

    return fig


def create_emissions_comparison_chart(emissions_data: Dict[str, Any]) -> Optional[go.Figure]:
    """
    Plotly chart showing estimated net greenhouse gas emissions avoided.
    """
    if not emissions_data:
        return None

    categories = ["Baseline Avoided GHG", "Hypothetical Avoided GHG", "Net GHG Delta"]
    values = [
        emissions_data["baseline_avoided_kg_co2"],
        emissions_data["hypo_avoided_kg_co2"],
        emissions_data["net_avoided_difference_kg_co2"]
    ]
    colors = ["#415A77", "#00A896", "#02C39A"]

    fig = go.Figure(go.Bar(
        x=categories,
        y=values,
        text=[f"{v:,.0f} kg CO2e" for v in values],
        textposition="outside",
        marker_color=colors,
        hovertemplate="<b>%{x}</b><br>Emissions Offset: %{y:,.1f} kg CO2e<extra></extra>"
    ))

    fig.update_layout(
        title="<b>Estimated Net Avoided GHG Emissions (User Assumptions)</b>",
        xaxis_title="<b>Scope</b>",
        yaxis_title="<b>Avoided Greenhouse Gases (kg CO2-eq)</b>",
        template="plotly_dark",
        paper_bgcolor="#0d1b2a",
        plot_bgcolor="#1b263b",
        font=dict(family="sans-serif", size=12, color="#E0E1DD"),
        margin=dict(l=40, r=40, t=60, b=40),
        height=360
    )

    return fig
