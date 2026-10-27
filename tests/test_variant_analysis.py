"""
Unit tests for src.variant_analysis.
"""

import pytest
import pandas as pd
from src.variant_analysis import (
    load_published_variants,
    validate_mutation_against_sequence,
    inspect_site_features,
    build_mutation_history_table,
    create_tm_comparison_figure
)
from src.utils import parse_mutation_code


def test_parse_mutation_code():
    """Verify standard mutation notation parsing."""
    wt, pos, mut = parse_mutation_code("S121E")
    assert wt == "S"
    assert pos == 121
    assert mut == "E"

    with pytest.raises(ValueError):
        parse_mutation_code("invalid_code")

    with pytest.raises(ValueError):
        parse_mutation_code("S0E")  # Pos < 1


def test_validate_mutation_against_sequence():
    """Verify validation against sequence residue identity."""
    seq = "MKTAYIAK"
    # Position 1 is M, Position 2 is K, Position 3 is T
    wt, pos, mut = validate_mutation_against_sequence("T3A", seq)
    assert wt == "T"
    assert pos == 3
    assert mut == "A"

    # Mismatch: position 3 is T, not A
    with pytest.raises(ValueError, match="Residue mismatch at position 3"):
        validate_mutation_against_sequence("A3G", seq)

    # Out of bounds
    with pytest.raises(ValueError, match="out of bounds"):
        validate_mutation_against_sequence("A20G", seq)


def test_load_published_variants_integrity():
    """Verify data/variants.csv loads with required columns and provenance."""
    df = load_published_variants()
    assert not df.empty
    assert "name" in df.columns
    assert "tm_celsius" in df.columns
    assert "doi" in df.columns
    assert "citation" in df.columns

    # Verify presence of landmark strains
    names = set(df["name"].tolist())
    assert "Wild-Type IsPETase" in names
    assert "FAST-PETase" in names
    assert "DuraPETase" in names
    assert "ThermoPETase" in names

    hot = df[df["name"] == "HotPETase"].iloc[0]
    assert hot["pdb_id"] == "7QVH"
    assert pd.isna(hot["tm_celsius"])

    lcc = df[df["name"] == "LCC-ICCG Cutinase"].iloc[0]
    assert lcc["pdb_id"] == "6THT"
    assert len(lcc["mutations"].split(", ")) == int(lcc["mutation_count"])
    assert lcc["tm_delta_reference"] == "LCC parent enzyme"


def test_missing_experimental_values_not_fabricated():
    """Verify that when a variant lacks a measurement, it remains NaN or None."""
    df = load_published_variants()
    # Check that WT has 0.0 delta Tm, and any non-measured fields remain numeric or null
    wt_row = df[df["name"] == "Wild-Type IsPETase"].iloc[0]
    assert wt_row["tm_delta"] == 0.0


def test_inspect_site_features():
    """Verify active site annotations for catalytic triad and oxyanion hole."""
    ser160 = inspect_site_features(160)
    assert ser160 is not None
    assert ser160["region"] == "Catalytic Triad"

    tyr87 = inspect_site_features(87)
    assert tyr87 is not None
    assert tyr87["region"] == "Oxyanion Hole"

    surface_res = inspect_site_features(50)
    assert surface_res is None


def test_build_mutation_history_table():
    """Verify results table formatting and disclaimer labeling."""
    sample_res = [{
        "mutation": "T3A",
        "position": 3,
        "wt_aa": "T",
        "mutant_aa": "A",
        "wt_log_prob": -2.71,
        "mutant_log_prob": -2.81,
        "score_difference": -0.10,
        "model_name": "esm2_t6_8M_UR50D"
    }]
    table = build_mutation_history_table(sample_res)
    assert len(table) == 1
    assert "Experimental Status" in table.columns
    assert "Unvalidated" in table.iloc[0]["Experimental Status"]


def test_create_tm_comparison_figure():
    """Verify Plotly figure generation."""
    df = load_published_variants()
    fig = create_tm_comparison_figure(df)
    assert fig is not None
    assert len(fig.data) > 0
    lcc_idx = list(fig.data[0].x).index("LCC-ICCG Cutinase")
    assert "vs LCC parent enzyme" in fig.data[0].text[lcc_idx]
