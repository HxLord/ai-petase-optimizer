import requests
import pytest

from src import protein_structure
from modules import structure_viewer


def test_fetch_structure_rejects_invalid_pdb_identifiers():
    with pytest.raises(ValueError, match="four-character"):
        protein_structure.fetch_structure_pdb("../5XJH")


def test_non_default_structure_does_not_return_ispetase_fallback(monkeypatch):
    def fail_request(*args, **kwargs):
        raise requests.ConnectionError("offline")

    monkeypatch.setattr(protein_structure.requests, "get", fail_request)
    pdb_text, source, is_fallback = protein_structure.fetch_structure_pdb(
        "1ABC", use_cache=False
    )

    assert pdb_text is None
    assert "RCSB request failed" in source
    assert not is_fallback


def test_default_structure_uses_verified_fallback_with_reason(monkeypatch):
    def fail_request(*args, **kwargs):
        raise requests.ConnectionError("offline")

    monkeypatch.setattr(protein_structure.requests, "get", fail_request)
    pdb_text, source, is_fallback = protein_structure.fetch_structure_pdb(
        "5XJH", use_cache=False
    )

    assert pdb_text and ("ATOM" in pdb_text or "HEADER" in pdb_text)
    assert "Verified Local Fallback" in source
    assert "RCSB request failed" in source
    assert is_fallback


def test_legacy_viewer_does_not_use_5xjh_for_other_pdb_ids(monkeypatch):
    def fail_request(*args, **kwargs):
        raise OSError("offline")

    monkeypatch.setattr(structure_viewer.urllib.request, "urlopen", fail_request)
    pdb_text, source, is_fallback = structure_viewer.fetch_pdb_data(
        "1ABC", timeout_sec=1
    )

    assert pdb_text is None
    assert "RCSB request failed" in source
    assert not is_fallback
