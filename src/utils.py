"""
Shared utility functions for AI-PETase Optimizer.
Provides input validation, numerical formatting, CSV export helpers,
and error handling utilities across modules and UI.
"""

import re
from typing import Optional, Any, Dict, List, Tuple
import pandas as pd

# 20 canonical proteinogenic amino acids (single-letter codes)
STANDARD_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")

# Three-letter to one-letter mapping
THREE_TO_ONE = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
    "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
    "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
    "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y"
}

# One-letter to full amino acid name
ONE_TO_NAME = {
    "A": "Alanine", "C": "Cysteine", "D": "Aspartate", "E": "Glutamate",
    "F": "Phenylalanine", "G": "Glycine", "H": "Histidine", "I": "Isoleucine",
    "K": "Lysine", "L": "Leucine", "M": "Methionine", "N": "Asparagine",
    "P": "Proline", "Q": "Glutamine", "R": "Arginine", "S": "Serine",
    "T": "Threonine", "V": "Valine", "W": "Tryptophan", "Y": "Tyrosine"
}


def validate_amino_acid(aa: str) -> str:
    """
    Validate that an amino acid string represents a single canonical residue.
    
    Args:
        aa: One-letter amino acid code.
        
    Returns:
        Upper-cased valid single-letter code.
        
    Raises:
        ValueError: If the amino acid is invalid or not in the 20 standard residues.
    """
    if not isinstance(aa, str):
        raise ValueError(f"Amino acid must be a string, got {type(aa).__name__}.")
    clean = aa.strip().upper()
    if len(clean) != 1 or clean not in STANDARD_AMINO_ACIDS:
        raise ValueError(
            f"Invalid amino acid '{aa}'. Must be one of the 20 standard residues: "
            f"{', '.join(sorted(STANDARD_AMINO_ACIDS))}."
        )
    return clean


def validate_sequence(sequence: str, max_length: int = 1022) -> str:
    """
    Normalize and validate a protein amino acid sequence.
    
    Args:
        sequence: Raw protein sequence string.
        max_length: Maximum allowed sequence length (default: 1022 for ESM models).
        
    Returns:
        Clean, uppercase sequence string.
        
    Raises:
        ValueError: If sequence is empty, exceeds max length, or contains non-canonical residues.
    """
    if not isinstance(sequence, str):
        raise ValueError(f"Protein sequence must be a string, got {type(sequence).__name__}.")
    if isinstance(max_length, bool) or not isinstance(max_length, int) or max_length < 1:
        raise ValueError(f"Maximum sequence length must be a positive integer, got {max_length}.")
    clean = "".join(sequence.split()).upper()
    if not clean:
        raise ValueError("Protein sequence cannot be empty.")
    
    invalid = set(clean) - STANDARD_AMINO_ACIDS
    if invalid:
        raise ValueError(
            f"Sequence contains non-canonical amino acid symbol(s): {', '.join(sorted(invalid))}. "
            "Only standard 20 amino acids (ACDEFGHIKLMNPQRSTVWY) are supported."
        )
    if len(clean) > max_length:
        raise ValueError(
            f"Sequence length ({len(clean)}) exceeds maximum allowed limit of {max_length} residues."
        )
    return clean


def parse_mutation_code(mutation_str: str) -> Tuple[str, int, str]:
    """
    Parse a standard mutation string (e.g. 'S121E' or 'N233K').
    
    Args:
        mutation_str: Mutation descriptor (e.g. 'S121E').
        
    Returns:
        Tuple of (wild_type_aa, 1-based position, mutant_aa).
        
    Raises:
        ValueError: If the mutation format is malformed or amino acids are non-canonical.
    """
    if not isinstance(mutation_str, str):
        raise ValueError(f"Mutation notation must be a string, got {type(mutation_str).__name__}.")
    clean = mutation_str.strip().upper()
    match = re.match(r"^([A-Z])(\d+)([A-Z])$", clean)
    if not match:
        raise ValueError(
            f"Invalid mutation notation '{mutation_str}'. "
            "Expected format '<WT><Position><MUT>' (e.g., 'S121E', 'N233K')."
        )
    wt_aa = validate_amino_acid(match.group(1))
    pos = int(match.group(2))
    if pos < 1:
        raise ValueError(f"Residue position must be >= 1, got {pos}.")
    mut_aa = validate_amino_acid(match.group(3))
    return wt_aa, pos, mut_aa


def format_number(val: Optional[float], decimals: int = 2, fallback: str = "N/A") -> str:
    """Safely format a numerical value with fixed decimal precision."""
    if val is None or pd.isna(val):
        return fallback
    try:
        return f"{float(val):.{decimals}f}"
    except (ValueError, TypeError):
        return fallback


def format_log_prob(val: Optional[float], decimals: int = 4, fallback: str = "N/A") -> str:
    """Format log-probability value for scientific presentation."""
    return format_number(val, decimals=decimals, fallback=fallback)


def dataframe_to_csv_bytes(df: pd.DataFrame) -> bytes:
    """Convert pandas DataFrame to UTF-8 CSV bytes for download."""
    return df.to_csv(index=False).encode("utf-8")


def records_to_csv_bytes(records: List[Dict[str, Any]]) -> bytes:
    """Convert list of dictionary records to UTF-8 CSV bytes."""
    df = pd.DataFrame(records)
    return dataframe_to_csv_bytes(df)


def format_error_message(error: Exception) -> str:
    """Generate a clean, informative error message for UI and logs."""
    return f"{type(error).__name__}: {str(error)}"
