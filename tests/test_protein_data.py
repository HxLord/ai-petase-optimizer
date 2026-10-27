"""
Unit tests for src.protein_data and sequence normalization.
"""

import pytest
import requests
from src.protein_data import (
    fetch_uniprot_sequence,
    load_user_fasta,
    parse_fasta_content,
    DEFAULT_UNIPROT_ACCESSION,
    VERIFIED_ISPETASE_FULL,
    VERIFIED_ISPETASE_MATURE
)
from src.utils import validate_sequence, validate_amino_acid, STANDARD_AMINO_ACIDS


def test_sequence_normalization_and_whitespace():
    """Verify sequence whitespace removal and uppercase conversion."""
    raw = "  mkt  ayia \n kqrq  "
    clean = validate_sequence(raw)
    assert clean == "MKTAYIAKQRQ"


def test_invalid_amino_acid_symbols():
    """Verify non-canonical characters (e.g. B, Z, X, numbers) raise ValueError."""
    with pytest.raises(ValueError, match="non-canonical amino acid"):
        validate_sequence("MKTAXB")

    with pytest.raises(ValueError, match="non-canonical amino acid"):
        validate_sequence("MKT123")


def test_empty_sequence_raises_error():
    """Verify empty or blank sequence raises ValueError."""
    with pytest.raises(ValueError, match="cannot be empty"):
        validate_sequence("")

    with pytest.raises(ValueError, match="cannot be empty"):
        validate_sequence("   \n\t  ")


def test_sequence_too_long_raises_error():
    """Verify sequence exceeding max_length raises ValueError."""
    long_seq = "A" * 1023
    with pytest.raises(ValueError, match="exceeds maximum allowed limit"):
        validate_sequence(long_seq, max_length=1022)


def test_validate_single_amino_acid():
    """Verify single amino acid validator."""
    assert validate_amino_acid("a") == "A"
    assert validate_amino_acid("W") == "W"
    
    with pytest.raises(ValueError):
        validate_amino_acid("X")
    with pytest.raises(ValueError):
        validate_amino_acid("ALA")  # Must be single letter


def test_parse_fasta_content():
    """Verify Biopython FASTA parsing logic."""
    fasta_str = ">test_seq Sample Description\nMKTA\nYIAK\n"
    res = parse_fasta_content(fasta_str)
    assert res["id"] == "test_seq"
    assert res["sequence"] == "MKTAYIAK"
    assert res["length"] == 8


def test_parse_fasta_corrupt():
    """Verify corrupt FASTA raises ValueError."""
    with pytest.raises(ValueError):
        parse_fasta_content("not a fasta format at all")


def test_fetch_uniprot_fallback(monkeypatch):
    """Verify offline retrieval is explicitly identified and uses the checked-in sequence."""
    def fail_request(*args, **kwargs):
        raise requests.Timeout("offline test")

    monkeypatch.setattr("src.protein_data.requests.get", fail_request)
    data = fetch_uniprot_sequence(
        accession=DEFAULT_UNIPROT_ACCESSION,
        timeout_sec=1,
        use_cache=False,
    )
    assert data["accession"] == DEFAULT_UNIPROT_ACCESSION
    assert data["is_fallback"]
    assert "offline test" in data["retrieval_warning"]
    assert len(data["sequence"]) == 290
    assert data["sequence"] == VERIFIED_ISPETASE_FULL
    assert len(data["mature_sequence"]) == 261


def test_mature_sequence_is_exact_uniprot_fragment():
    assert len(VERIFIED_ISPETASE_MATURE) == 261
    assert VERIFIED_ISPETASE_MATURE == VERIFIED_ISPETASE_FULL[29:]
    assert VERIFIED_ISPETASE_MATURE[130] == "S"


def test_load_user_fasta_string():
    """Verify loading user sequence from FASTA text string."""
    fasta = ">custom_variant My engineered strain\nMNFPRASRLM\n"
    res = load_user_fasta(fasta)
    assert res["name"] == "custom_variant"
    assert res["sequence"] == "MNFPRASRLM"
    assert res["length"] == 10
    assert not res["is_fallback"]
