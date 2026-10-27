"""
Scientific reference data, sequences, structural annotations,
and environmental baseline parameters for AI-PETase Optimizer.
"""

from typing import Dict, Any, List

# Wild-Type Ideonella sakaiensis PETase (IsPETase) catalytic domain sequence
# Source: UniProtKB A0A0K8P6T7 / RCSB PDB: 5XJH (Chain A, residues 34-290, with residue 30-33 prefix)
# PDB 5XJH coordinates span residues 30 to 290.
ISPETASE_PDB_SEQUENCE = (
    "MNFPRASRLMQAAVLGGLMAVSAAATAQTNPYARGPNPTAASLEASAGPFTVRSFTVSRPSGYGAGTVYYPTN"
    "AGGTVGAIAIVPGYTARQSSIKWWGPRLASHGFVVITIDTNSTLDQPSSRSSQQMAALRQVASLNGTSSSPIY"
    "GKVDTARMGVMGWSMGGGGSLISAANNPSLKAAAPQAPWDSSTNFSSVTVPTLIFACENDSIAPVNSSALPIY"
    "DSMSRNAKQFLEINGGSHSCANSGNSNQALIGKKGVAWMKRFMDNDTRYSTFACENPNSTRVSDFRTANCS"
)

# UniProt residues 30-290 (261 aa); mature-fragment positions map to UniProt/PDB
# numbering by adding 29.
ISPETASE_MATURE_SEQUENCE = ISPETASE_PDB_SEQUENCE[29:]

# Standard Structural Landmarks in 5XJH (PDB numbering)
STRUCTURAL_FEATURES = {
    "catalytic_triad": {
        "title": "Catalytic Triad",
        "residues": [
            {"resname": "SER", "resseq": 160, "role": "Nucleophile that attacks the ester carbonyl of PET"},
            {"resname": "ASP", "resseq": 206, "role": "Acid/base charge relay stabilizing the catalytic histidine"},
            {"resname": "HIS", "resseq": 237, "role": "General base activating Ser160 for nucleophilic attack"}
        ],
        "description": "The classic catalytic triad of α/β hydrolases. In IsPETase, Asp206 and His237 position Ser160 to cleave ester linkages.",
        "color": "#E63946" # Red/Coral
    },
    "oxyanion_hole": {
        "title": "Oxyanion Hole",
        "residues": [
            {"resname": "TYR", "resseq": 87, "role": "Hydrogen bond donor to the tetrahedral intermediate carbonyl oxygen"},
            {"resname": "MET", "resseq": 161, "role": "Main-chain amide hydrogen bond donor to tetrahedral oxygen"}
        ],
        "description": "Stabilizes the negatively charged oxyanion formed during the transition state of ester bond hydrolysis.",
        "color": "#F4A261" # Orange
    },
    "disulfide_bonds": {
        "title": "Disulfide Bonds",
        "pairs": [
            {"res1": 203, "res2": 239, "name": "DS1 (Cys203-Cys239)", "role": "Conserved disulfide bridge stabilizing the catalytic subsite loop"},
            {"res1": 273, "res2": 289, "name": "DS2 (Cys273-Cys289)", "role": "Unique C-terminal disulfide bridge absent in most cutinases, critical for structural stability"}
        ],
        "residues": [
            {"resname": "CYS", "resseq": 203},
            {"resname": "CYS", "resseq": 239},
            {"resname": "CYS", "resseq": 273},
            {"resname": "CYS", "resseq": 289}
        ],
        "description": "Covalent cystine bridges that anchor the tertiary fold. Absence of DS2 severely reduces IsPETase stability.",
        "color": "#E9C46A" # Gold/Yellow
    },
    "substrate_binding_cleft": {
        "title": "Aromatic Substrate Binding Cleft",
        "residues": [
            {"resname": "TRP", "resseq": 185, "role": "Aromatic stacking with PET terephthalate benzene ring (wobble tryptophan)"},
            {"resname": "TRP", "resseq": 159, "role": "Hydrophobic cleft lining adjacent to catalytic Ser160"},
            {"resname": "ILE", "resseq": 208, "role": "Hydrophobic floor of subsite I"},
            {"resname": "THR", "resseq": 88, "role": "Subsite II channel restriction"}
        ],
        "description": "A shallow, open aromatic surface groove that allows the bulky, insoluble PET polymer chain to bind without severe steric clashes.",
        "color": "#2A9D8F" # Teal
    }
}

# Educational distinction between metrics
METRIC_EXPLANATIONS = {
    "tm": {
        "name": "Melting Temperature (Tm)",
        "unit": "°C",
        "definition": "The temperature at which 50% of the enzyme molecules are denatured (unfolded) at equilibrium in solution.",
        "measurement_method": "Circular Dichroism (CD) spectroscopy or Differential Scanning Fluorimetry (DSF / SYPRO Orange).",
        "interpretation": "A thermodynamic stability metric. While a higher Tm allows the enzyme to survive higher temperatures without unravelling, it does NOT automatically guarantee faster reaction velocity at ambient room temperature.",
        "caution": "Do not treat higher Tm as a direct measurement of faster plastic degradation."
    },
    "opt_temp": {
        "name": "Optimal Operating Temperature (Topt)",
        "unit": "°C",
        "definition": "The temperature at which the kinetic rate of hydrolysis is maximal before thermal inactivation causes activity to drop off.",
        "measurement_method": "Incubation assays across a temperature gradient (e.g., 30°C to 75°C) measuring product release over fixed intervals.",
        "interpretation": "Kinetic parameter reflecting the balance between the Arrhenius increase in catalytic rate and irreversible denaturation.",
        "caution": "PET chain mobility around a glass transition depends on material composition and thermal history. It does not establish that an enzyme will be effective at that temperature."
    },
    "enzyme_activity": {
        "name": "Enzyme Activity (Initial Rate)",
        "unit": "U/mg or μmol/min/mg",
        "definition": "The initial speed of ester bond cleavage, often measured on small soluble synthetic model substrates.",
        "measurement_method": "Spectrophotometric assay using p-nitrophenyl butyrate (pNPB) or bis(2-hydroxyethyl) terephthalate (BHET).",
        "interpretation": "Measures catalytic turnover in homogeneous solution.",
        "caution": "Activity on soluble pNPB or BHET often has VERY POOR correlation with the ability to digest insoluble, crystalline macroscopic PET plastics."
    },
    "pet_degradation": {
        "name": "Solid PET Depolymerization",
        "unit": "mg product / day or % mass loss",
        "definition": "Direct chemical breakdown of solid, insoluble PET film, powder, or post-consumer flakes into monomers (TPA, MHET, BHET, EG).",
        "measurement_method": "High-Performance Liquid Chromatography (HPLC) of released monomers or gravimetric mass loss of the solid polymer.",
        "interpretation": "The true operational benchmark for plastic recycling.",
        "caution": "Rates depend on substrate morphology, pretreatment, enzyme loading, assay duration, and analytical method. Values from different studies are not directly comparable without reconciling these conditions."
    }
}

# Physical Properties of PET & Recycling Lifecycles
PET_PHYSICAL_PROPERTIES = {
    "iupac_name": "Poly(ethylene terephthalate)",
    "repeating_unit": "C10H8O4",
    "repeating_unit_mw_g_mol": 192.17,
    "density_g_cm3": 1.38,
    "glass_transition_temp_range_c": "65°C - 75°C",
    "melting_point_c": "250°C - 260°C",
    "hydrolysis_products": [
        {"name": "Terephthalic Acid (TPA)", "formula": "C8H6O4", "mw": 166.13, "state": "Solid monomer, repolymerizable"},
        {"name": "Mono-(2-hydroxyethyl) terephthalate (MHET)", "formula": "C10H10O5", "mw": 210.18, "state": "Intermediate hydrolysate"},
        {"name": "Bis(2-hydroxyethyl) terephthalate (BHET)", "formula": "C12H14O6", "mw": 254.24, "state": "Soluble oligomer"},
        {"name": "Ethylene Glycol (EG)", "formula": "C2H6O2", "mw": 62.07, "state": "Liquid diol monomer"}
    ]
}

# No project-specific life-cycle inventory has been validated for this prototype.
# Keep these compatibility fields empty rather than presenting assumed values as benchmarks.
LCA_BENCHMARKS = {
    "virgin_pet_emissions_factor_default": None,
    "virgin_pet_emissions_range": None,
    "mechanical_recycling_emissions_default": None,
    "mechanical_recycling_emissions_range": None,
    "enzymatic_recycling_emissions_default": None,
    "enzymatic_recycling_emissions_range": None,
    "global_pet_production_annual_mt": None,
    "ocean_plastic_annual_input_mt": None
}

# Primary Research Citations
KEY_CITATIONS = [
    {
        "title": "A bacterium that degrades and assimilates poly(ethylene terephthalate)",
        "authors": "Yoshida, S., Hiraga, K., Takehana, T., Taniguchi, I., Yamaji, H., Maeda, Y., Toyohara, K., Miyamoto, K., Kimura, Y., & Oda, K.",
        "journal": "Science",
        "volume": "351(6278)",
        "pages": "1196–1199",
        "year": 2016,
        "doi": "10.1126/science.aad6359",
        "url": "https://www.science.org/doi/10.1126/science.aad6359",
        "summary": "Discovered Ideonella sakaiensis 201-F6 and identified the original wild-type PETase enzyme capable of hydrolyzing PET into MHET."
    },
    {
        "title": "Structural insight into catalytic mechanism of PET hydrolase",
        "authors": "Han, X., Liu, W., Huang, J. W., Ma, J., Zheng, Y., Ko, T. P., Xu, L., Cheng, Y. S., Chen, C. C., & Guo, R. T.",
        "journal": "Nature Communications",
        "volume": "8",
        "pages": "2106",
        "year": 2017,
        "doi": "10.1038/s41467-017-02255-z",
        "url": "https://www.nature.com/articles/s41467-017-02255-z",
        "summary": "Published high-resolution crystal structure (PDB 5XJH) of IsPETase, elucidating the catalytic triad and oxyanion hole."
    },
    {
        "title": "Machine learning-aided engineering of hydrolases for PET depolymerization",
        "authors": "Lu, H., Diaz, D. J., Czarnecki, N. J., Zhu, C., Kim, W., Shroff, R., Acosta, D. J., Alexander, B. R., Cole, H. O., Zhang, Y., Lynd, N. A., Ellington, A. D., & Alper, H. S.",
        "journal": "Nature",
        "volume": "604",
        "pages": "662–667",
        "year": 2022,
        "doi": "10.1038/s41586-022-04599-z",
        "url": "https://www.nature.com/articles/s41586-022-04599-z",
        "summary": "Engineered FAST-PETase using the MutCompute 3D convolutional neural network, demonstrating closed-loop circular PET recycling."
    },
    {
        "title": "An engineered PET depolymerase to break down and recycle plastic bottles",
        "authors": "Tournier, V., Topham, C. M., Gilles, A., David, B., Folgoas, C., Moya-Leclair, E., Kamionka, E., Desrousseaux, M. L., Texier, H., Gavalda, S., Cot, M., Guémard, E., Dalibey, M., Nomme, J., Cioci, G., Barbe, S., Chateau, M., André, I., Duquesne, S., & Marty, A.",
        "journal": "Nature",
        "volume": "580",
        "pages": "216–219",
        "year": 2020,
        "doi": "10.1038/s41586-020-2149-4",
        "url": "https://www.nature.com/articles/s41586-020-2149-4",
        "summary": "Engineered LCC-ICCG cutinase with disulfide bridges achieving 90% depolymerization of bottle flakes in 10 hours at 72°C."
    },
    {
        "title": "Robust enzyme design by deep learning",
        "authors": "Cui, Y., Chen, Y., Liu, X., Dong, S., Tian, Y. E., Qiao, Y., Mitra, R., Han, J., Li, C., Han, X., Liu, W., Chen, Q., Wei, W., Wang, X., Du, W., Zheng, S., Liu, S., Ju, H., Chen, C. C., Feng, Y., & Wu, B.",
        "journal": "ACS Catalysis",
        "volume": "11(2)",
        "pages": "738–747",
        "year": 2021,
        "doi": "10.1021/acscatal.0c05126",
        "url": "https://pubs.acs.org/doi/10.1021/acscatal.0c05126",
        "summary": "Computationally redesigned PETase using the GRAPE strategy; variant-level measurements and assay conditions should be interpreted from the linked primary paper."
    },
    {
        "title": "Evolutionary Scale Modeling (ESM) Protein Language Models",
        "authors": "Lin, Z., Akin, H., Rao, R., Hie, B., Zhu, Z., Lu, W., Smetanin, N., Verkuil, R., Kabeli, O., Shmueli, Y., dos Santos Costa, A., Fazel-Zarandi, M., Sercu, T., Candido, S., & Rives, A.",
        "journal": "Science",
        "volume": "379(6637)",
        "pages": "1123–1130",
        "year": 2023,
        "doi": "10.1126/science.ade2574",
        "url": "https://github.com/facebookresearch/esm",
        "summary": "State-of-the-art transformer protein language models (ESM-2, ESM-1v) trained on hundreds of millions of evolutionary sequences."
    }
]
