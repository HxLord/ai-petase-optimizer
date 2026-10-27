"""
Variant Research and Analysis Module for AI-PETase Optimizer.
Handles curated variant comparisons, metric sanity checks,
and scientific Plotly visualizations with clear experimental condition grouping.
"""

import os
import json
from typing import Dict, Any, List, Optional
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px

VARIANTS_JSON_PATH = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data", "variants.json")


def load_variants_data() -> Dict[str, Any]:
    """Load curated peer-reviewed variant records from JSON."""
    if os.path.exists(VARIANTS_JSON_PATH):
        with open(VARIANTS_JSON_PATH, "r", encoding="utf-8") as f:
            return json.load(f)
    return {}


def get_variants_dataframe() -> pd.DataFrame:
    """Return tabular representation of variants for scientific comparison."""
    data = load_variants_data()
    rows = []
    for var_id, item in data.items():
        delta = item.get("tm_delta")
        rows.append({
            "id": var_id,
            "Name": item["name"],
            "Scaffold": item["scaffold"],
            "PDB": item["pdb_id"],
            "Mutations": ", ".join(item["mutations_list"]) if item["mutations_list"] else "Wild-Type (0)",
            "Mutation Count": item["mutations_count"],
            "Tm (°C)": item["tm_celsius"],
            "ΔTm (°C)": (
                f"{delta:+.1f}" if pd.notna(delta) else "Not reported"
            ),
            "ΔTm Reference": item.get("tm_delta_reference"),
            "Opt Temp (°C)": item["opt_temp_celsius"],
            "Opt pH": item["opt_ph"],
            "Substrate Tested": item["degradation_substrate"],
            "Assay Conditions": item["assay_conditions"],
            "Reported Metric": item["reported_degradation_metric"],
            "Citation": item["citation"],
            "Is WT": item.get("is_wt_reference", False)
        })
    return pd.DataFrame(rows)


def create_tm_comparison_chart(variants_data: Dict[str, Any], selected_id: Optional[str] = None) -> go.Figure:
    """
    Plotly bar chart of Melting Temperatures (Tm) with scaffold-specific delta references.
    Highlights selected variant and colors by thermostability gain.
    """
    comparable = {
        var_id: item
        for var_id, item in variants_data.items()
        if pd.notna(item.get("tm_celsius"))
    }
    if not comparable:
        fig = go.Figure()
        fig.update_layout(title="No verified Tm measurements available")
        return fig

    names = []
    tms = []
    deltas = []
    colors = []
    borders = []
    
    wt_tm = variants_data.get("wt_petase", {}).get("tm_celsius")
    
    for var_id, item in comparable.items():
        names.append(item["name"])
        tm = item["tm_celsius"]
        tms.append(tm)
        delta = item.get("tm_delta")
        deltas.append(delta)
        
        # Color coding: WT is neutral blue/grey, selected is vibrant teal, others teal gradient
        if var_id == "wt_petase":
            colors.append("#6C757D") # Grey for WT reference
        elif var_id == selected_id:
            colors.append("#00A896") # Vibrant teal for active selection
        else:
            colors.append("#028090") # Steel teal
            
        if var_id == selected_id:
            borders.append("3px solid #02C39A")
        else:
            borders.append("1px solid rgba(255,255,255,0.2)")

    fig = go.Figure()

    # Main Tm bars
    fig.add_trace(go.Bar(
        x=names,
        y=tms,
        text=[
            (
                f"{tm:.1f}°C<br>({delta:+.1f}°C vs {item.get('tm_delta_reference')})"
                if pd.notna(delta) and pd.notna(item.get("tm_delta_reference"))
                else f"{tm:.1f}°C"
            )
            for tm, delta, item in zip(tms, deltas, comparable.values())
        ],
        textposition="outside",
        marker=dict(
            color=colors,
            line=dict(color="#02C39A", width=[2.5 if var_id == selected_id else 0 for var_id in comparable.keys()])
        ),
        hovertemplate="<b>%{x}</b><br>Tm: %{y:.1f}°C<extra></extra>"
    ))

    # Reference line for WT Tm
    if pd.notna(wt_tm):
        fig.add_hline(
            y=wt_tm,
            line_dash="dash",
            line_color="#E63946",
            annotation_text=f"WT IsPETase reported Tm ({wt_tm}°C)",
            annotation_position="bottom right",
            annotation_font_color="#E63946"
        )

    # Reference line for PET Glass Transition Temperature (Tg ~ 65°C - 75°C)
    fig.add_hrect(
        y0=65, y1=75,
        fillcolor="#F4A261", opacity=0.15,
        line_width=0,
        annotation_text="PET Glass Transition Zone (Tg 65-75°C)",
        annotation_position="top left",
        annotation_font_color="#F4A261"
    )

    fig.update_layout(
        title="<b>Reported Tm Values (assays may differ; see sources)</b>",
        xaxis_title="<b>Enzyme Variant</b>",
        yaxis_title="<b>Melting Temperature Tm (°C)</b>",
        yaxis=dict(range=[30, 105]),
        template="plotly_dark",
        paper_bgcolor="#0d1b2a",
        plot_bgcolor="#1b263b",
        font=dict(family="sans-serif", size=12, color="#E0E1DD"),
        margin=dict(l=40, r=40, t=60, b=40),
        height=420
    )

    return fig


def create_temperature_landscape_chart(variants_data: Dict[str, Any], selected_id: Optional[str] = None) -> go.Figure:
    """
    Scatter plot of Tm vs Optimal Operating Temperature (Topt).
    Illustrates the operational gap between unfolding and catalytic turnover.
    """
    names = []
    tms = []
    topts = []
    mut_counts = []
    scaffolds = []
    ids = []

    for var_id, item in variants_data.items():
        if pd.isna(item.get("tm_celsius")) or pd.isna(item.get("opt_temp_celsius")):
            continue
        ids.append(var_id)
        names.append(item["name"])
        tms.append(item["tm_celsius"])
        topts.append(item["opt_temp_celsius"])
        mut_counts.append(item.get("mutations_count") or 0)
        scaffolds.append(item["scaffold"])

    sizes = [max(16, count * 3 + 14) for count in mut_counts]

    fig = go.Figure()

    # Scatter trace
    for i, var_id in enumerate(ids):
        is_sel = (var_id == selected_id)
        is_wt = (var_id == "wt_petase")
        
        color = "#00A896" if is_sel else ("#6C757D" if is_wt else "#48CAE4")
        if var_id == "lcc_iccg":
            color = "#E76F51" # Cutinase distinct color

        fig.add_trace(go.Scatter(
            x=[topts[i]],
            y=[tms[i]],
            mode="markers+text",
            name=names[i],
            text=[names[i]],
            textposition="top center",
            marker=dict(
                size=sizes[i] + (6 if is_sel else 0),
                color=color,
                line=dict(color="#FFFFFF" if is_sel else "rgba(255,255,255,0.4)", width=2 if is_sel else 1)
            ),
            hovertemplate=(
                f"<b>{names[i]}</b><br>"
                f"Opt Operating Temp: %{{x}}°C<br>"
                f"Melting Temp Tm: %{{y}}°C<br>"
                f"Mutations: {mut_counts[i]}<br>"
                f"Scaffold: {scaffolds[i]}<extra></extra>"
            )
        ))

    fig.update_layout(
        title="<b>Reported Topt vs. Tm (measurements depend on assay)</b>",
        xaxis_title="<b>Optimal Operating Temperature Topt (°C)</b>",
        yaxis_title="<b>Melting Temperature Tm (°C)</b>",
        xaxis=dict(range=[25, 85]),
        yaxis=dict(range=[40, 105]),
        template="plotly_dark",
        paper_bgcolor="#0d1b2a",
        plot_bgcolor="#1b263b",
        font=dict(family="sans-serif", size=12, color="#E0E1DD"),
        showlegend=False,
        margin=dict(l=40, r=40, t=60, b=40),
        height=400
    )

    return fig
