"""
Protein Data Retrieval and Processing Module for AI-PETase Optimizer.
Fetches, parses, and validates the wild-type IsPETase amino acid sequence
from UniProtKB, with Biopython FASTA parsing, local caching, and verified offline fallback.
"""

import io
import os
import re
from typing import Dict, Any, Union, Optional
import requests
from Bio import SeqIO

from src.utils import validate_sequence

# Default target accession for wild-type Ideonella sakaiensis PETase
DEFAULT_UNIPROT_ACCESSION = "A0A0K8P6T7"

# Bundled, verified full-length precursor sequence (290 amino acids, including N-terminal signal peptide 1-29)
# Source: UniProtKB A0A0K8P6T7 (PETH_PISS1 / ISF6_4831)
VERIFIED_ISPETASE_FULL = (
    "MNFPRASRLMQAAVLGGLMAVSAAATAQTNPYARGPNPTAASLEASAGPFTVRSFTVSRPSGYGAGTVYYPTN"
    "AGGTVGAIAIVPGYTARQSSIKWWGPRLASHGFVVITIDTNSTLDQPSSRSSQQMAALRQVASLNGTSSSPIY"
    "GKVDTARMGVMGWSMGGGGSLISAANNPSLKAAAPQAPWDSSTNFSSVTVPTLIFACENDSIAPVNSSALPIY"
    "DSMSRNAKQFLEINGGSHSCANSGNSNQALIGKKGVAWMKRFMDNDTRYSTFACENPNSTRVSDFRTANCS"
)

# UniProt residues 30-290 (261 aa); positions in this fragment map to full-length
# UniProt/PDB numbering by adding 29.
VERIFIED_ISPETASE_MATURE = VERIFIED_ISPETASE_FULL[29:]

CACHE_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data", "cache")


def get_cache_path(accession: str) -> str:
    """Return file path for cached FASTA record."""
    clean_acc = re.sub(r"[^A-Za-z0-9_-]", "", accession.strip())
    return os.path.join(CACHE_DIR, f"{clean_acc}.fasta")


def parse_fasta_content(fasta_text: str) -> Dict[str, Any]:
    """
    Parse FASTA string using Biopython SeqIO.
    
    Args:
        fasta_text: FASTA formatted string.
        
    Returns:
        Dict with sequence metadata and validated amino acid sequence.
        
    Raises:
        ValueError: If FASTA cannot be parsed or contains no records.
    """
    try:
        stream = io.StringIO(fasta_text.strip())
        records = list(SeqIO.parse(stream, "fasta"))
    except Exception as exc:
        raise ValueError(f"Failed to parse FASTA content: {str(exc)}") from exc
    if not records:
        raise ValueError("No valid FASTA records found in input.")
    
    first = records[0]
    seq_str = str(first.seq).strip()
    clean_seq = validate_sequence(seq_str)
    
    return {
        "id": first.id,
        "name": first.name,
        "description": first.description,
        "sequence": clean_seq,
        "length": len(clean_seq),
        "all_records_count": len(records)
    }


def fetch_uniprot_sequence(
    accession: str = DEFAULT_UNIPROT_ACCESSION,
    timeout_sec: int = 5,
    use_cache: bool = True
) -> Dict[str, Any]:
    """
    Retrieve protein sequence by UniProt accession.
    
    Tries:
    1. Local cache in `data/cache/<accession>.fasta` if use_cache is True.
    2. Live network request to UniProt REST API.
    3. Verified bundled fallback sequence if network fails and accession is DEFAULT_UNIPROT_ACCESSION.
    
    Args:
        accession: UniProt accession code (e.g. 'A0A0K8P6T7').
        timeout_sec: Network timeout in seconds.
        use_cache: Whether to inspect/write local cache.
        
    Returns:
        Dict containing accession, sequence, length, description, source, and is_fallback.
        
    Raises:
        RuntimeError: If remote fetch fails and no verified fallback exists.
    """
    if not isinstance(accession, str) or not accession.strip():
        raise ValueError("UniProt accession must be a non-empty string.")
    clean_acc = accession.strip().upper()
    cache_file = get_cache_path(clean_acc)
    retrieval_warning = None

    # 1. Check local cache
    if use_cache and os.path.exists(cache_file):
        try:
            with open(cache_file, "r", encoding="utf-8") as f:
                content = f.read()
            parsed = parse_fasta_content(content)
            # Verify record belongs to requested accession
            if clean_acc in parsed["id"] or clean_acc in parsed["description"]:
                return {
                    "accession": clean_acc,
                    "name": parsed["id"],
                    "description": parsed["description"],
                    "sequence": parsed["sequence"],
                    "length": parsed["length"],
                    "is_fallback": False,
                    "source": f"Local Cache ({os.path.basename(cache_file)})",
                    "retrieval_warning": None,
                    "full_sequence": parsed["sequence"],
                    "mature_sequence": (
                        parsed["sequence"][29:] if len(parsed["sequence"]) == 290 and clean_acc == DEFAULT_UNIPROT_ACCESSION
                        else parsed["sequence"]
                    )
                }
        except (OSError, ValueError) as exc:
            retrieval_warning = f"Cached FASTA could not be used: {exc}"

    # 2. Attempt live network fetch
    url = f"https://rest.uniprot.org/uniprotkb/{clean_acc}.fasta"
    try:
        resp = requests.get(
            url,
            timeout=timeout_sec,
            headers={"User-Agent": "AI-PETase-Optimizer/1.0 (Research App)"}
        )
        resp.raise_for_status()
        if not resp.text:
            raise ValueError("UniProt returned an empty FASTA response.")
        parsed = parse_fasta_content(resp.text)
        if clean_acc not in parsed["id"] and clean_acc not in parsed["description"]:
            raise ValueError(
                f"Returned UniProt record '{parsed['id']}' does not match requested accession '{clean_acc}'."
            )

        if use_cache:
            os.makedirs(CACHE_DIR, exist_ok=True)
            with open(cache_file, "w", encoding="utf-8") as f:
                f.write(resp.text)

        return {
            "accession": clean_acc,
            "name": parsed["id"],
            "description": parsed["description"],
            "sequence": parsed["sequence"],
            "length": parsed["length"],
            "is_fallback": False,
            "source": f"UniProtKB Live API ({clean_acc})",
            "retrieval_warning": None,
            "full_sequence": parsed["sequence"],
            "mature_sequence": (
                parsed["sequence"][29:] if len(parsed["sequence"]) == 290 and clean_acc == DEFAULT_UNIPROT_ACCESSION
                else parsed["sequence"]
            )
        }
    except (requests.RequestException, OSError, ValueError) as exc:
        retrieval_warning = str(exc)

    # 3. Verified bundled fallback for default IsPETase
    if clean_acc == DEFAULT_UNIPROT_ACCESSION:
        return {
            "accession": DEFAULT_UNIPROT_ACCESSION,
            "name": "sp|A0A0K8P6T7|PETH_PISS1",
            "description": "Poly(ethylene terephthalate) hydrolase OS=Piscinibacter sakaiensis (Ideonella sakaiensis 201-F6) [Verified Bundled Fallback]",
            "sequence": VERIFIED_ISPETASE_FULL,
            "length": len(VERIFIED_ISPETASE_FULL),
            "is_fallback": True,
            "source": "Verified Bundled Fallback (A0A0K8P6T7)",
            "retrieval_warning": retrieval_warning,
            "full_sequence": VERIFIED_ISPETASE_FULL,
            "mature_sequence": VERIFIED_ISPETASE_MATURE
        }

    raise RuntimeError(
        f"Unable to retrieve UniProt record for accession '{clean_acc}' (network unreachable or timeout) "
        "and no verified local fallback exists for this protein."
    )


def load_user_fasta(fasta_source: Union[str, bytes, io.IOBase]) -> Dict[str, Any]:
    """
    Parse a user-uploaded FASTA file or string.
    
    Args:
        fasta_source: FASTA string, bytes, or file-like object.
        
    Returns:
        Dict with sequence metadata and validated amino acid sequence.
    """
    if isinstance(fasta_source, bytes):
        text = fasta_source.decode("utf-8", errors="replace")
    elif hasattr(fasta_source, "read"):
        content = fasta_source.read()
        text = content.decode("utf-8", errors="replace") if isinstance(content, bytes) else content
    else:
        text = str(fasta_source)

    parsed = parse_fasta_content(text)
    return {
        "accession": parsed["id"],
        "name": parsed["name"],
        "description": parsed["description"],
        "sequence": parsed["sequence"],
        "length": parsed["length"],
        "is_fallback": False,
        "source": "User-Uploaded FASTA",
        "full_sequence": parsed["sequence"],
        "mature_sequence": parsed["sequence"]
    }
