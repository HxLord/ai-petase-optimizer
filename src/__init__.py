"""
AI-PETase Optimizer Core Scientific Package.
Integrates sequence validation, ESM-2 protein language model inference,
experimental variant benchmarking, 3D molecular structure visualization,
and hypothetical circularity scenario calculation.
"""

from src.protein_data import (
    fetch_uniprot_sequence,
    load_user_fasta,
    DEFAULT_UNIPROT_ACCESSION,
    VERIFIED_ISPETASE_FULL,
    VERIFIED_ISPETASE_MATURE,
)
from src.protein_ai import (
    score_mutation,
    load_model,
    MODEL_IDENTIFIER,
)
from src.variant_analysis import (
    load_published_variants,
    validate_mutation_against_sequence,
    build_mutation_history_table,
    create_tm_comparison_figure,
)
from src.protein_structure import (
    fetch_structure_pdb,
    generate_3dmol_html,
    DEFAULT_PDB_ID,
)
from src.eco_calculator import (
    calculate_scenario,
    create_scenario_figure,
)
from src.utils import (
    validate_sequence,
    validate_amino_acid,
    parse_mutation_code,
    format_number,
    format_log_prob,
    dataframe_to_csv_bytes,
    records_to_csv_bytes,
    STANDARD_AMINO_ACIDS,
)

__all__ = [
    "fetch_uniprot_sequence",
    "load_user_fasta",
    "DEFAULT_UNIPROT_ACCESSION",
    "VERIFIED_ISPETASE_FULL",
    "VERIFIED_ISPETASE_MATURE",
    "score_mutation",
    "load_model",
    "MODEL_IDENTIFIER",
    "load_published_variants",
    "validate_mutation_against_sequence",
    "build_mutation_history_table",
    "create_tm_comparison_figure",
    "fetch_structure_pdb",
    "generate_3dmol_html",
    "DEFAULT_PDB_ID",
    "calculate_scenario",
    "create_scenario_figure",
    "validate_sequence",
    "validate_amino_acid",
    "parse_mutation_code",
    "format_number",
    "format_log_prob",
    "dataframe_to_csv_bytes",
    "records_to_csv_bytes",
    "STANDARD_AMINO_ACIDS",
]
