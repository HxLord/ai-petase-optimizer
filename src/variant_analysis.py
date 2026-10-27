"""
Variant Analysis and Benchmark Comparison Module for AI-PETase Optimizer.
Manages curated published experimental variants from literature (FAST-PETase, DuraPETase, etc.)
and compares them with user-proposed AI mutation scores without conflating predictions with measurements.
"""

import os
from typing import Dict, Any, List, Optional, Tuple
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px

from src.utils import parse_mutation_code, validate_sequence

DEFAULT_VARIANTS_CSV = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data", "variants.csv")

# Known landmarks use full-length UniProt/PDB residue numbering.
ISPETASE_LANDMARKS = {
    160: {"region": "Catalytic Triad", "role": "Catalytic nucleophile (Ser160) forming covalent acyl-enzyme intermediate."},
    206: {"region": "Catalytic Triad", "role": "Charge relay acid (Asp206) polarizing His237."},
    237: {"region": "Catalytic Triad", "role": "General base (His237) deprotonating Ser160."},
    87: {"region": "Oxyanion Hole", "role": "Hydrogen bond donor (Tyr87) stabilizing tetrahedral oxyanion intermediate."},
    161: {"region": "Oxyanion Hole", "role": "Main-chain amide donor (Met161) stabilizing transition state."},
    203: {"region": "Disulfide Bridge DS1", "role": "Covalent cystine bridge (Cys203-Cys239) anchoring the active-site subsite."},
    239: {"region": "Disulfide Bridge DS1", "role": "Covalent cystine bridge (Cys203-Cys239)."},
    273: {"region": "Disulfide Bridge DS2", "role": "C-terminal cystine bridge (Cys273-Cys289) unique to PETases, critical for stability."},
    289: {"region": "Disulfide Bridge DS2", "role": "C-terminal cystine bridge (Cys273-Cys289)."},
    185: {"region": "Substrate Binding Cleft", "role": "Wobble tryptophan (Trp185) providing aromatic π-stacking with PET benzene rings."},
    159: {"region": "Substrate Binding Cleft", "role": "Hydrophobic clamp (Trp159) guiding polymer chain into catalytic cleft."}
}


def load_published_variants(csv_path: Optional[str] = None) -> pd.DataFrame:
    """
    Load curated experimental variant data from CSV.
    
    Args:
        csv_path: Optional path to variants CSV. Defaults to `data/variants.csv`.
        
    Returns:
        DataFrame containing published experimental variants and provenance fields.
    """
    path = csv_path or DEFAULT_VARIANTS_CSV
    if not os.path.exists(path):
        return pd.DataFrame()

    df = pd.read_csv(path)
    # Ensure numeric columns are properly typed without inventing missing data
    numeric_cols = ["mutation_count", "tm_celsius", "tm_delta", "opt_temp_celsius", "opt_ph", "depolymerization_temp_tested"]
    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    return df


def validate_mutation_against_sequence(mutation_str: str, sequence: str) -> Tuple[str, int, str]:
    """
    Validate that a proposed mutation matches the underlying sequence residue.
    
    Args:
        mutation_str: e.g. 'S160A'
        sequence: Protein sequence string
        
    Returns:
        Tuple of (wt_aa, position, mutant_aa)
        
    Raises:
        ValueError: If mutation syntax is invalid, out of bounds, or does not match sequence.
    """
    clean_seq = validate_sequence(sequence)
    wt_aa, pos, mut_aa = parse_mutation_code(mutation_str)
    
    if pos < 1 or pos > len(clean_seq):
        raise ValueError(
            f"Position {pos} is out of bounds for sequence of length {len(clean_seq)}."
        )
        
    actual_aa = clean_seq[pos - 1]
    if actual_aa != wt_aa:
        raise ValueError(
            f"Residue mismatch at position {pos}: mutation specifies '{wt_aa}', "
            f"but sequence actually contains '{actual_aa}'."
        )
        
    return wt_aa, pos, mut_aa


def inspect_site_features(position: int) -> Optional[Dict[str, str]]:
    """Return structural landmark annotation if position is a known critical site."""
    return ISPETASE_LANDMARKS.get(position)


def build_mutation_history_table(results_list: List[Dict[str, Any]]) -> pd.DataFrame:
    """
    Build a structured, reusable results DataFrame from a list of AI scoring outputs.
    
    Args:
        results_list: List of dictionaries returned by `score_mutation`.
        
    Returns:
        Formatted DataFrame with explicit labeling of predictions vs measurements.
    """
    if not results_list:
        return pd.DataFrame()

    rows = []
    for r in results_list:
        sequence_pos = r.get("sequence_position", r.get("position"))
        reference_pos = r.get("reference_position")
        site_info = (
            inspect_site_features(reference_pos)
            if reference_pos is not None and r.get("reference_accession") == "A0A0K8P6T7"
            else None
        )
        rows.append({
            "Mutation": r.get("mutation"),
            "Position": reference_pos if reference_pos is not None else sequence_pos,
            "Selected Sequence Position": sequence_pos,
            "WT Residue": r.get("wt_aa"),
            "Mutant Residue": r.get("mutant_aa"),
            "WT LogProb": r.get("wt_log_prob"),
            "Mutant LogProb": r.get("mutant_log_prob"),
            "Δ LogProb (Score)": r.get("score_difference"),
            "Active Site Landmark": site_info["region"] if site_info else "Surface/Scaffold",
            "Model": r.get("model_name", "ESM-2"),
            "Experimental Status": "Unvalidated (Computational AI Score Only)"
        })
    return pd.DataFrame(rows)


def create_tm_comparison_figure(df: pd.DataFrame, highlight_name: Optional[str] = None) -> go.Figure:
    """
    Generate an interactive Plotly bar chart comparing published melting temperatures (Tm).
    
    Args:
        df: Variants DataFrame containing 'name', 'tm_celsius', 'tm_delta'.
        highlight_name: Optional variant name to highlight.
        
    Returns:
        go.Figure: Styled bar chart with WT reference line and glass transition zone.
    """
    if df.empty or "tm_celsius" not in df.columns:
        fig = go.Figure()
        fig.update_layout(title="No published variant data available")
        return fig

    # Sort by Tm ascending
    plot_df = df.dropna(subset=["tm_celsius"]).sort_values("tm_celsius")

    colors = []
    for _, row in plot_df.iterrows():
        if row.get("is_wt_reference", False):
            colors.append("#64756f")  # Neutral reference for wild-type
        elif row["name"] == highlight_name:
            colors.append("#c27a35")  # Warm highlight for selected variant
        else:
            colors.append("#168575")  # Primary teal

    fig = go.Figure()

    fig.add_trace(go.Bar(
        x=plot_df["name"],
        y=plot_df["tm_celsius"],
        text=[
            (
                f"{row['tm_celsius']:.1f}°C "
                f"({row['tm_delta']:+.1f}°C vs {row['tm_delta_reference']})"
                if pd.notna(row["tm_delta"]) and pd.notna(row.get("tm_delta_reference"))
                else f"{row['tm_celsius']:.1f}°C"
            )
            for _, row in plot_df.iterrows()
        ],
        textposition="outside",
        marker=dict(color=colors, line=dict(color="#f8faf7", width=1.5)),
        hovertemplate="<b>%{x}</b><br>Tm: %{y:.1f}°C<extra></extra>"
    ))

    # WT reference line
    wt_rows = plot_df[plot_df.get("is_wt_reference", False)]
    if not wt_rows.empty:
        wt_tm = wt_rows.iloc[0]["tm_celsius"]
        fig.add_hline(
            y=wt_tm,
            line_dash="dot",
            line_color="#ad5d44",
            annotation_text=f"WT Baseline Tm ({wt_tm:.1f}°C)",
            annotation_position="bottom right"
        )

    # PET glass transition region (Tg ~ 65-75°C)
    fig.add_hrect(
        y0=65, y1=75,
        fillcolor="#dba257",
        opacity=0.12,
        line_width=0,
        annotation_text="PET Glass Transition Range (Tg: 65–75°C)",
        annotation_position="top left"
    )

    fig.update_layout(
        title="<b>Reported melting temperatures (Tm)</b><br>"
              "<span style='font-size:12px; color:#536763;'>"
              "Source: Peer-reviewed literature (Nature, Science, ACS Catalysis). "
              "Note: Tm measures thermodynamic unfolding, not catalytic reaction speed.</span>",
        xaxis_title="Variant Name",
        yaxis_title="Melting Temperature Tm (°C)",
        template="plotly_white",
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="#ffffff",
        font=dict(family="DM Sans, sans-serif", size=12, color="#243b35"),
        title_font=dict(family="Newsreader, Georgia, serif", size=20, color="#172b29"),
        xaxis=dict(showgrid=False, linecolor="#cbd8d0", tickfont=dict(color="#405851")),
        yaxis=dict(range=[30, 105], gridcolor="#e5ebe7", zerolinecolor="#cbd8d0", tickfont=dict(color="#405851")),
        margin=dict(l=40, r=40, t=70, b=40)
    )

    return fig
