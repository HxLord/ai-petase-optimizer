"""
AI Protein Sequence Analysis Module for AI-PETase Optimizer.
Utilizes pretrained ESM-2 (Evolutionary Scale Modeling) embeddings to score
single-point amino acid substitutions via masked marginal log-likelihoods.
"""

from functools import lru_cache
from typing import Dict, Any, Tuple, Optional
import torch
import esm

from src.utils import validate_sequence, validate_amino_acid, STANDARD_AMINO_ACIDS

MODEL_IDENTIFIER = "esm2_t6_8M_UR50D"
MAX_SEQUENCE_LENGTH = 1022


def get_default_device() -> torch.device:
    """Return CUDA device if available, otherwise CPU."""
    return torch.device("cuda" if torch.cuda.is_available() else "cpu")


@lru_cache(maxsize=1)
def load_model(device_str: Optional[str] = None) -> Tuple[torch.nn.Module, Any]:
    """
    Load and cache pretrained ESM-2 model and alphabet.
    
    Args:
        device_str: Optional target device ('cpu' or 'cuda'). Defaults to auto-detection.
        
    Returns:
        Tuple of (model, alphabet).
        
    Raises:
        RuntimeError: If model weights cannot be loaded.
    """
    device = torch.device(device_str) if device_str else get_default_device()
    try:
        model, alphabet = esm.pretrained.esm2_t6_8M_UR50D()
        model.to(device)
        model.eval()
        return model, alphabet
    except Exception as exc:
        raise RuntimeError(
            f"Failed to initialize pretrained model '{MODEL_IDENTIFIER}': {str(exc)}. "
            "Ensure torch and fair-esm are installed and model checkpoint is accessible."
        ) from exc


def score_mutation(
    sequence: str,
    position: int,
    mutant_aa: str,
    allow_synonymous: bool = False
) -> Dict[str, Any]:
    """
    Score a single-site amino acid mutation using ESM-2 masked-marginal log-likelihood.

    Args:
        sequence: Protein sequence (one-letter canonical codes).
        position: 1-based residue position (1 <= position <= len(sequence)).
        mutant_aa: Single-letter target amino acid code.
        allow_synonymous: If False, raises ValueError when mutant_aa == wt_aa.

    Returns:
        dict: Mutation details, log-probabilities, score difference, model ID, and scientific disclaimer.

    Raises:
        ValueError: For invalid sequence, position, amino acid, or synonymous mutation when not allowed.
        RuntimeError: If model weights are not loaded.
    """
    # 1. Sequence normalization and bounds check
    clean_seq = validate_sequence(sequence, max_length=MAX_SEQUENCE_LENGTH)
    seq_len = len(clean_seq)

    if isinstance(position, bool) or not isinstance(position, int) or position < 1 or position > seq_len:
        raise ValueError(
            f"Residue position {position} is out of bounds for sequence of length {seq_len}. "
            f"Valid range is 1 to {seq_len}."
        )

    clean_mut_aa = validate_amino_acid(mutant_aa)
    wt_aa = clean_seq[position - 1]

    if not allow_synonymous and clean_mut_aa == wt_aa:
        raise ValueError(
            f"Proposed substitution '{clean_mut_aa}' at position {position} is identical to the "
            f"wild-type residue '{wt_aa}' (no amino acid change)."
        )

    mutation_str = f"{wt_aa}{position}{clean_mut_aa}"

    # 2. Load model and tokenizer
    model, alphabet = load_model()
    device = next(model.parameters()).device
    batch_converter = alphabet.get_batch_converter()

    # 3. Tokenize sequence and mask target position
    data = [("query_sequence", clean_seq)]
    _, _, batch_tokens = batch_converter(data)
    batch_tokens = batch_tokens.to(device)

    # In ESM: index 0 is <cls> token, index `position` corresponds to residue clean_seq[position - 1]
    batch_tokens[0, position] = alphabet.mask_idx

    # 4. Predict logits using torch.inference_mode()
    with torch.inference_mode():
        results = model(batch_tokens, return_contacts=False)
        logits = results["logits"]

    # Compute log probabilities across alphabet vocabulary
    log_probs = torch.log_softmax(logits[0, position], dim=-1)

    wt_idx = alphabet.get_idx(wt_aa)
    mut_idx = alphabet.get_idx(clean_mut_aa)

    original_log_prob = log_probs[wt_idx].item()
    mutant_log_prob = log_probs[mut_idx].item()
    score_difference = mutant_log_prob - original_log_prob

    # 5. Scientific interpretation and disclaimers
    interpretation = (
        "Positive score difference indicates higher compatibility with the learned sequence "
        "distribution relative to wild-type under masked-marginal scoring. It does NOT prove "
        "improved thermostability (Tm), activity, or faster plastic degradation."
        if score_difference >= 0 else
        "Negative score difference indicates lower compatibility with the learned sequence "
        "distribution relative to wild-type under masked-marginal scoring. It may reflect "
        "evolutionary penalty or structural destabilization, but requires experimental validation."
    )

    disclaimer = (
        "SCIENTIFIC NOTICE: Protein language model log-likelihood differences reflect sequence "
        "context compatibility. They are not direct physical measurements of melting temperature (Tm), "
        "free energy change (ddG), kinetic rate constants, or macroscopic PET depolymerization."
    )

    return {
        "mutation": mutation_str,
        "position": position,
        "wt_aa": wt_aa,
        "original_aa": wt_aa,
        "mutant_aa": clean_mut_aa,
        "original_log_prob": round(original_log_prob, 4),
        "wt_log_prob": round(original_log_prob, 4),
        "mutant_log_prob": round(mutant_log_prob, 4),
        "score_difference": round(score_difference, 4),
        "delta_log_prob": round(score_difference, 4),
        "model_name": MODEL_IDENTIFIER,
        "model_identifier": MODEL_IDENTIFIER,
        "device": str(device),
        "interpretation": interpretation,
        "warning": disclaimer,
        "disclaimer": disclaimer,
    }
