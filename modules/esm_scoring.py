"""
ESM-2 Protein Language Model Scoring Module for AI-PETase Optimizer.
Mutation scores are produced by pretrained model inference; there is no heuristic fallback.

Scientific limitation: ESM-2 sequence scores do not directly measure free energy (ΔΔG),
melting temperature (Tm), catalytic activity, or depolymerization rate on insoluble polymers.
"""

import re
from typing import Dict, Any, Optional, Tuple, List
from data.reference_data import ISPETASE_PDB_SEQUENCE, STRUCTURAL_FEATURES

# BLOSUM62 utility retained for callers that explicitly request substitution-matrix scores.
# It is not used as a fallback by the ESM-2 scoring path.
BLOSUM62 = {
    ('A','A'): 4, ('A','R'): -1, ('A','N'): -2, ('A','D'): -2, ('A','C'): 0, ('A','Q'): -1, ('A','E'): -1, ('A','G'): 0, ('A','H'): -2, ('A','I'): -1, ('A','L'): -1, ('A','K'): -1, ('A','M'): -1, ('A','F'): -2, ('A','P'): -1, ('A','S'): 1, ('A','T'): 0, ('A','W'): -3, ('A','Y'): -2, ('A','V'): 0,
    ('R','R'): 5, ('R','N'): 0, ('R','D'): -2, ('R','C'): -3, ('R','Q'): 1, ('R','E'): 0, ('R','G'): -2, ('R','H'): 0, ('R','I'): -3, ('R','L'): -2, ('R','K'): 2, ('R','M'): -1, ('R','F'): -3, ('R','P'): -2, ('R','S'): -1, ('R','T'): -1, ('R','W'): -3, ('R','Y'): -2, ('R','V'): -3,
    ('N','N'): 6, ('N','D'): 1, ('N','C'): -3, ('N','Q'): 0, ('N','E'): 0, ('N','G'): 0, ('N','H'): 1, ('N','I'): -3, ('N','L'): -3, ('N','K'): 0, ('N','M'): -2, ('N','F'): -3, ('N','P'): -2, ('N','S'): 1, ('N','T'): 0, ('N','W'): -4, ('N','Y'): -2, ('N','V'): -3,
    ('D','D'): 6, ('D','C'): -3, ('D','Q'): 0, ('D','E'): 2, ('D','G'): -1, ('D','H'): -1, ('D','I'): -3, ('D','L'): -4, ('D','K'): -1, ('D','M'): -3, ('D','F'): -3, ('D','P'): -1, ('D','S'): 0, ('D','T'): -1, ('D','W'): -4, ('D','Y'): -3, ('D','V'): -3,
    ('C','C'): 9, ('C','Q'): -3, ('C','E'): -4, ('C','G'): -3, ('C','H'): -3, ('C','I'): -1, ('C','L'): -1, ('C','K'): -3, ('C','M'): -1, ('C','F'): -2, ('C','P'): -3, ('C','S'): -1, ('C','T'): -1, ('C','W'): -2, ('C','Y'): -2, ('C','V'): -1,
    ('Q','Q'): 5, ('Q','E'): 2, ('Q','G'): -2, ('Q','H'): 0, ('Q','I'): -3, ('Q','L'): -2, ('Q','K'): 1, ('Q','M'): 0, ('Q','F'): -3, ('Q','P'): -1, ('Q','S'): 0, ('Q','T'): -1, ('Q','W'): -2, ('Q','Y'): -1, ('Q','V'): -2,
    ('E','E'): 5, ('E','G'): -2, ('E','H'): 0, ('E','I'): -3, ('E','L'): -3, ('E','K'): 1, ('E','M'): -2, ('E','F'): -3, ('E','P'): -1, ('E','S'): 0, ('E','T'): -1, ('E','W'): -3, ('E','Y'): -2, ('E','V'): -2,
    ('G','G'): 6, ('G','H'): -2, ('G','I'): -4, ('G','L'): -4, ('G','K'): -2, ('G','M'): -3, ('G','F'): -3, ('G','P'): -2, ('G','S'): 0, ('G','T'): -2, ('G','W'): -2, ('G','Y'): -3, ('G','V'): -3,
    ('H','H'): 8, ('H','I'): -3, ('H','L'): -3, ('H','K'): -1, ('H','M'): -2, ('H','F'): -1, ('H','P'): -2, ('H','S'): -1, ('H','T'): -2, ('H','W'): -2, ('H','Y'): 2, ('H','V'): -3,
    ('I','I'): 4, ('I','L'): 2, ('I','K'): -3, ('I','M'): 1, ('I','F'): 0, ('I','P'): -3, ('I','S'): -2, ('I','T'): -1, ('I','W'): -3, ('I','Y'): -1, ('I','V'): 3,
    ('L','L'): 4, ('L','K'): -2, ('L','M'): 2, ('L','F'): 0, ('L','P'): -3, ('L','S'): -2, ('L','T'): -1, ('L','W'): -2, ('L','Y'): -1, ('L','V'): 1,
    ('K','K'): 5, ('K','M'): -1, ('K','F'): -3, ('K','P'): -1, ('K','S'): 0, ('K','T'): -1, ('K','W'): -3, ('K','Y'): -2, ('K','V'): -2,
    ('M','M'): 5, ('M','F'): 0, ('M','P'): -2, ('M','S'): -1, ('M','T'): -1, ('M','W'): -1, ('M','Y'): -1, ('M','V'): 1,
    ('F','F'): 6, ('F','P'): -4, ('F','S'): -2, ('F','T'): -2, ('F','W'): 1, ('F','Y'): 3, ('F','V'): -1,
    ('P','P'): 7, ('P','S'): -1, ('P','T'): -1, ('P','W'): -4, ('P','Y'): -3, ('P','V'): -2,
    ('S','S'): 4, ('S','T'): 1, ('S','W'): -3, ('S','Y'): -2, ('S','V'): -2,
    ('T','T'): 5, ('T','W'): -2, ('T','Y'): -2, ('T','V'): 0,
    ('W','W'): 11, ('W','Y'): 2, ('W','V'): -3,
    ('Y','Y'): 7, ('Y','V'): -1,
    ('V','V'): 4
}


def get_blosum_score(a1: str, a2: str) -> int:
    """Retrieve symmetrical BLOSUM62 substitution score."""
    pair = (a1.upper(), a2.upper())
    rev_pair = (a2.upper(), a1.upper())
    if pair in BLOSUM62:
        return BLOSUM62[pair]
    if rev_pair in BLOSUM62:
        return BLOSUM62[rev_pair]
    return -4


def check_esm_availability() -> Tuple[bool, str]:
    """
    Check that the pretrained ESM-2 model and its weights can actually be loaded.
    Returns (is_available, status_message).
    """
    try:
        from src.protein_ai import MODEL_IDENTIFIER, load_model
        load_model()
        return True, f"Pretrained ESM-2 model '{MODEL_IDENTIFIER}' loaded."
    except (ImportError, RuntimeError) as exc:
        return False, f"Pretrained ESM-2 is unavailable: {exc}"


def parse_mutation_code(mutation_str: str) -> Optional[Tuple[str, int, str]]:
    """
    Parse a mutation code string like 'S121E' or 'N233K'.
    Returns (wt_aa, position, mut_aa) or None.
    """
    clean_str = mutation_str.strip().upper()
    match = re.match(r"^([A-Z])(\d+)([A-Z])$", clean_str)
    if not match:
        return None
    wt_aa = match.group(1)
    pos = int(match.group(2))
    mut_aa = match.group(3)
    return wt_aa, pos, mut_aa


def analyze_mutation_site(pos: int) -> Dict[str, Any]:
    """
    Inspect whether the mutation position resides near critical active-site landmarks.
    """
    # Positions are in full-length UniProt/PDB numbering.
    if pos in [160, 206, 237]:
        return {
            "is_critical": True,
            "region": "Catalytic Triad",
            "warning": "Catalytic triad residue. A substitution may perturb catalysis; its effect is not predicted by this annotation."
        }
    # Oxyanion hole: Tyr87, Met161
    if pos in [87, 161]:
        return {
            "is_critical": True,
            "region": "Oxyanion Hole",
            "warning": "Oxyanion-hole residue. A substitution may affect transition-state interactions; experimental validation is required."
        }
    # Disulfide bonds: 203, 239, 273, 289
    if pos in [203, 239, 273, 289]:
        return {
            "is_critical": True,
            "region": "Covalent Disulfide Bridge",
            "warning": "Disulfide-bridge residue. A substitution may disrupt the bridge; stability effects require structural and experimental assessment."
        }
    # Substrate cleft: 185, 159, 208, 88
    if pos in [185, 159, 208, 88]:
        return {
            "is_critical": False,
            "region": "Substrate Binding Cleft",
            "warning": "SUBSTRATE CONTACT: Modulates aromatic stacking or steric clearance for polymer chains."
        }
    # Surface loop / scaffold
    return {
        "is_critical": False,
        "region": "Surface Scaffold / Loop",
        "warning": "SURFACE POSITION: Commonly targeted by rational design to introduce salt-bridges or reduce local flexibility without distorting the active site."
    }


def score_variant_mutation(mutation_str: str) -> Dict[str, Any]:
    """
    Score a single mutation with the pretrained ESM-2 model.
    No heuristic or zero-score fallback is returned if model inference fails.
    """
    parsed = parse_mutation_code(mutation_str)
    if not parsed:
        return {
            "success": False,
            "error": f"Invalid mutation syntax '{mutation_str}'. Please use standard format like 'S121E' or 'N233K'."
        }

    wt_aa, pos, mut_aa = parsed
    site_info = analyze_mutation_site(pos)
    from src.protein_ai import MODEL_IDENTIFIER, score_mutation

    if pos < 1 or pos > len(ISPETASE_PDB_SEQUENCE):
        return {
            "success": False,
            "error": f"Position {pos} is out of bounds for the {len(ISPETASE_PDB_SEQUENCE)}-residue UniProt sequence.",
        }
    actual_wt = ISPETASE_PDB_SEQUENCE[pos - 1]
    if actual_wt != wt_aa:
        return {
            "success": False,
            "error": (
                f"Wild-type residue mismatch at position {pos}: mutation specifies '{wt_aa}', "
                f"but the UniProt reference contains '{actual_wt}'."
            ),
        }

    try:
        result = score_mutation(ISPETASE_PDB_SEQUENCE, pos, mut_aa)
    except (ValueError, RuntimeError) as exc:
        return {"success": False, "error": str(exc)}

    return {
        **result,
        "success": True,
        "mut_aa": mut_aa,
        "score": result["score_difference"],
        "verdict": "No biological classification; ESM-2 sequence score only",
        "verdict_color": "#8D99AE",
        "confidence": "Not estimated",
        "is_real_ai": True,
        "scoring_engine": f"ESM-2 ({MODEL_IDENTIFIER})",
        "site_info": site_info,
        "disclaimer": result["disclaimer"],
    }
