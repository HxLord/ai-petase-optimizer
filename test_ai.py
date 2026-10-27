from src.protein_ai import score_mutation

# Small test sequence to verify the model pipeline.
# This is NOT a PETase sequence.
sequence = "MKTAYIAKQRQISFVKSHFSRQ"
result = score_mutation(
    sequence=sequence,
    position=3,
    mutant_aa="A",
)
for key, value in result.items():
    print(f"{key}: {value}")
