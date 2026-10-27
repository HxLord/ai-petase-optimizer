"""
Unit and integration tests for src.protein_ai.
"""

import pytest
from src.protein_ai import score_mutation, MODEL_IDENTIFIER


def test_invalid_mutation_position():
    """Verify position < 1 or position > len(seq) raises ValueError."""
    seq = "MKTAYIAK"
    with pytest.raises(ValueError, match="out of bounds"):
        score_mutation(sequence=seq, position=0, mutant_aa="A")

    with pytest.raises(ValueError, match="out of bounds"):
        score_mutation(sequence=seq, position=9, mutant_aa="A")


def test_synonymous_mutation_without_change():
    """Verify mutations that do not change the amino acid raise ValueError by default."""
    seq = "MKTAYIAK"
    # Residue at 1-based pos 1 is 'M'
    with pytest.raises(ValueError, match="identical to the wild-type residue"):
        score_mutation(sequence=seq, position=1, mutant_aa="M")


def test_invalid_mutant_symbol():
    """Verify non-standard mutant amino acid symbol raises ValueError."""
    seq = "MKTAYIAK"
    with pytest.raises(ValueError, match="Invalid amino acid"):
        score_mutation(sequence=seq, position=2, mutant_aa="X")


def test_esm_inference_scoring():
    """Verify genuine ESM-2 model inference on small test sequence."""
    seq = "MKTAYIAKQRQISFVKSHFSRQ"
    res = score_mutation(sequence=seq, position=3, mutant_aa="A")
    
    assert res["mutation"] == "T3A"
    assert res["position"] == 3
    assert res["wt_aa"] == "T"
    assert res["mutant_aa"] == "A"
    assert isinstance(res["wt_log_prob"], float)
    assert isinstance(res["mutant_log_prob"], float)
    assert isinstance(res["score_difference"], float)
    assert res["model_name"] == MODEL_IDENTIFIER
    assert "SCIENTIFIC NOTICE" in res["disclaimer"]
    
    # Mathematical delta check
    expected_delta = round(res["mutant_log_prob"] - res["wt_log_prob"], 4)
    assert res["score_difference"] == expected_delta
