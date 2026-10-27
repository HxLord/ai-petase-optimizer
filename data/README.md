# Data Directory

This directory stores structural files, curated reference datasets, and local caches used by the **AI-PETase Optimizer** application.

## Files and Subdirectories

- **`variants.csv`**: Source-linked literature records for selected PET hydrolases.
  - Columns:
    - `id`: Unique identifier for the variant record.
    - `name`: Common scientific designation (e.g. FAST-PETase, DuraPETase).
    - `scaffold`: Wild-type parent scaffold (*Ideonella sakaiensis* PETase or Leaf-branch compost cutinase).
    - `pdb_id`: Protein Data Bank accession for the experimentally resolved structure.
    - `mutations`: Comma-separated list of substitutions relative to wild-type scaffold.
    - `mutation_count`: Total number of engineered substitutions.
    - `tm_celsius`: Reported melting temperature ($T_m$) in °C when populated.
    - `tm_delta`: Reported $T_m$ difference ($^\circ\text{C}$) when available.
    - `tm_delta_reference`: Scaffold-specific reference used for the reported $T_m$ difference; do not interpret LCC values relative to IsPETase.
    - `opt_temp_celsius`: Reported assay temperature or optimum when available.
    - `opt_ph`: Reported assay pH when available.
    - `reported_degradation_metric`: A short reported observation; not normalized across studies.
    - `degradation_substrate`: Physical polymer substrate tested (e.g. amorphous film vs. semi-crystalline flake).
    - `assay_conditions`: Reaction buffer, ionic strength, pH, and temperature.
    - `depolymerization_temp_tested`: Temperature at which degradation assays were conducted ($^\circ\text{C}$).
    - `citation`: Primary academic citation.
    - `doi`: Digital Object Identifier linking to the publication.
    - `rcsb_url`: Direct link to the RCSB Protein Data Bank entry.
    - `paper_url`: Direct link to the journal article.
    - `is_wt_reference`: Boolean flag indicating whether the entry is the baseline wild-type reference.

- **`structures/5xjh.pdb`**: High-resolution (1.58 Å) X-ray crystal structure of *Ideonella sakaiensis* PETase (Han et al., *Nat Commun* 2017). Provides seamless offline 3D visualization.

- **`cache/`**: Local cache directory for fetched UniProt records and structure coordinates.

## Scientific Integrity & Data Provenance Notice

The table is a manually curated reference, not a normalized or independently reanalyzed meta-analysis. Assay methods, substrates, and conditions differ, so values must be checked against their linked primary sources before scientific use. Blank cells mean no value is supplied; they are not imputed from model output.

PDB identifiers are checked against RCSB records. HotPETase has a verified structure and citation, but measurement fields remain blank until their exact values and assay conditions are confirmed. Environmental scenario outputs elsewhere in the application are arithmetic based on assumptions, not life-cycle assessment results.
