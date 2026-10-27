# AI-PETase Optimizer 🧬♻️

### AI-Assisted Research Platform for PET Plastic Depolymerization Enzymes

**AI-PETase Optimizer** is a research prototype for exploring IsPETase sequences, calculating protein-language-model mutation scores, viewing a reference structure, and comparing explicitly hypothetical material scenarios.

It integrates optional live UniProt retrieval, pretrained **ESM-2** inference, interactive 3D structure rendering with **3Dmol.js**, a source-linked literature table, and arithmetic scenario modeling. It is not an enzyme-design optimizer, laboratory prediction system, or life-cycle assessment tool.

---

## 🌟 Key Application Capabilities & Pages

1. **Page 1 — Overview**
   - Background on PET hydrolysis products, including MHET intermediates.
   - Contextual reference values with assay-dependent caveats.
   - Transparent research prototype disclaimers.

2. **Page 2 — Protein Explorer**
   - Switch between the 261-residue UniProt fragment, the 290-residue precursor, a live UniProt request, or a user-uploaded FASTA.
   - Mature-fragment mutation indices are translated to full-length UniProt/PDB numbering by adding 29.
   - Client-side interactive 3D WebGL protein visualization of PDB **5XJH** (1.58 Å resolution, Han et al., *Nat Commun* 2017) using 3Dmol.js.
   - Biochemical landmark highlighting:
     - **Catalytic Triad**: Ser160 (nucleophile), Asp206 (charge relay), His237 (general base).
     - **Disulfide Bridges**: Conserved DS1 (Cys203-Cys239) and unique C-terminal DS2 (Cys273-Cys289).
     - **Oxyanion Hole**: Tyr87 and Met161 transition-state stabilizers.
     - **Aromatic Substrate Cleft**: Trp185 (wobble tryptophan) and Trp159 (hydrophobic clamp).
   - Multiple molecular representations: Cartoon Ribbon, Stick, Sphere (Spacefill), and Solvent-Accessible Surface (VDW).

3. **Page 3 — AI Mutation Analysis**
   - Point mutation designer: select any 1-based residue position and propose single amino acid substitutions.
   - Pretrained **ESM-2** (`esm2_t6_8M_UR50D`) protein language model inference using masked-marginal log-likelihood evaluation.
   - Calculates wild-type log-prob $\log p(x_{\text{wt}})$, mutant log-prob $\log p(x_{\text{mut}})$, and delta score $\Delta \log p = \log p(x_{\text{mut}}) - \log p(x_{\text{wt}})$.
   - Active site landmark warning system alerting researchers if proposed substitutions affect catalytic or structural cystine positions.
   - One-click CSV export of session mutation scoring history.

4. **Page 4 — Variant Comparison**
   - Source-linked records for:
     - **FAST-PETase** (*Nature* 2022, Lu et al.)
     - **DuraPETase** (*ACS Catalysis* 2021, Cui et al.)
     - **ThermoPETase** (*ACS Catalysis* 2019, Son et al.)
     - **HotPETase** (*Nature Catalysis* 2022, Bell et al.; unverified measurements are blank)
     - **LCC-ICCG Cutinase** (*Nature* 2020, Tournier et al.)
   - Plotly display of available reported melting temperatures ($T_m$), with missing values left blank.
   - Notes explaining that $T_m$, assay optimum, and PET depolymerization are distinct measurements.

5. **Page 5 — Environmental Impact**
   - Dynamic mass-balance calculation evaluating user-configured PET waste batches (kg), baseline recovery rates, and alternative hypothetical scenarios.
   - Hypothetical mass-balance and optional GHG arithmetic based on displayed inputs.
   - Interactive Plotly chart and downloadable CSV with assumptions and a scenario disclaimer. This is not an LCA.

6. **Page 6 — Research & Methodology**
   - Links to UniProt, RCSB PDB, ESM, and cited literature.
   - Limitations and scientific scope disclosures.

---

## 📁 Project Structure

```
gr/
├── app.py                      # Main Streamlit web application & UI pages
├── requirements.txt            # Python package dependencies
├── README.md                   # Comprehensive documentation & setup guide
├── .gitignore                  # Git exclusions for environments, caches & weights
├── .streamlit/
│   └── config.toml             # Streamlit SaaS dark navy/teal theme configuration
├── src/
│   ├── __init__.py             # Package exports
│   ├── protein_data.py         # UniProt sequence retrieval, Biopython FASTA parser & local cache
│   ├── protein_ai.py           # ESM-2 neural inference, masked marginal scoring & torch runner
│   ├── protein_structure.py    # RCSB PDB coordinate loader & 3Dmol.js HTML generator
│   ├── variant_analysis.py     # Published variant dataset loader & Plotly charts
│   ├── eco_calculator.py       # Material mass balances & GHG scenario calculator
│   └── utils.py                # Validation helpers, safe formatting & CSV export utilities
├── data/
│   ├── README.md               # Dataset documentation & provenance guidelines
│   ├── variants.csv            # Curated published experimental variant measurements
│   ├── cache/                  # Local cache for UniProt FASTA records
│   └── structures/
│       └── 5xjh.pdb            # Offline crystal structure of IsPETase (1.58 Å)
├── tests/
│   ├── test_protein_data.py    # Tests for sequence normalization, FASTA parsing & fallbacks
│   ├── test_protein_ai.py      # Tests for ESM-2 mutation scoring & input validation
│   ├── test_variant_analysis.py# Tests for variant tables, landmarks & Plotly figures
│   └── test_eco_calculator.py  # Tests for mass balances, boundary values & emissions
└── assets/
    └── README.md               # Visual assets and diagram documentation
```

---

## Requirements

- Python 3.10 or newer (tested here with Python 3.13.1).
- Windows 10/11, Linux, or macOS.
- CPU inference is supported; CUDA is used when available.
- Internet is needed for live UniProt retrieval and first-time ESM-2 checkpoint download. The PDB structure has a bundled fallback; 3Dmol.js is loaded from its CDN.

---

## Windows setup (PowerShell)

Open the project folder (`gr`) in PowerShell:

### 1. Create a Python Virtual Environment
```powershell
py -3.10 -m venv .venv
```

If Python 3.10 is not installed, use an installed compatible version, for example `py -3.13 -m venv .venv`. Then activate and install:

```powershell
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

If PowerShell blocks activation, call the virtual-environment interpreter directly:

```powershell
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe -m pytest -q
.\.venv\Scripts\python.exe -m streamlit run app.py
```

Run syntax checks, tests, and the app from the activated environment:

```powershell
python -m compileall -q app.py src tests modules data test_app.py test_ai.py
python -m pytest -q
python -m streamlit run app.py
```

Streamlit opens at `http://localhost:8501`.

### VS Code

1. Open the `gr` folder (**File → Open Folder**).
2. Run **Python: Select Interpreter** from the Command Palette and select `.\.venv\Scripts\python.exe`.
3. Open a new integrated PowerShell terminal and confirm the prompt uses `.venv`, or use the interpreter directly as above.
4. Install dependencies, run `python -m pytest -q`, then launch with `python -m streamlit run app.py`.

---

## 🧠 AI Model Architecture & Checkpoint Caching

- **Model**: `esm2_t6_8M_UR50D` via `fair-esm`.
- ESM-2 weights are downloaded by fair-esm on the first scored mutation and cached by PyTorch. A successful package import alone does not mean the weights are available.
- Scores use masked marginal log-probability difference under `torch.inference_mode()` with the model in evaluation mode.
- If weights cannot be loaded, the application reports an error; it does not substitute random, zero, or heuristic scores.

---

## 🔬 Mandatory Scientific Limitations & Disclosures

1. **Evolutionary Likelihood vs Laboratory Stability**:
   - ESM-2 returns model log-probabilities for a masked residue context. A positive $\Delta \log p$ is not a validated prediction of stability, function, or fitness.
2. **Macroscopic Plastic vs Soluble Esters**:
   - Higher activity on small soluble model substrates (such as p-nitrophenyl butyrate or BHET) often correlates poorly with depolymerization rates on macroscopic, insoluble, semi-crystalline PET bottles.
3. **Polymer Crystallinity**:
   - PET crystallinity, thermal history, and assay conditions affect accessibility. A glass-transition value is not a universal threshold that establishes degradation performance.
4. **Hypothetical Scenarios**:
   - The ecological calculator performs arithmetic on assumed rates and an assumed net GHG factor. Outputs are not experimental yields or a project-specific LCA.

---

## 🛠️ Troubleshooting

- **PyTorch / ESM import error**:
  Ensure PyTorch and `fair-esm` are installed in your active virtual environment:
  ```powershell
  python -m pip install torch fair-esm
  ```
- **Structure does not display in 3D viewer**:
  The 3D viewer loads `3Dmol.js` from `https://3Dmol.org/build/3Dmol-min.js`. Ensure your browser allows WebGL and JavaScript execution. If working completely offline, verify `data/structures/5xjh.pdb` exists.
- **UniProt or model weights cannot be loaded**:
  Check network access. The default IsPETase sequence and reference PDB are bundled, but model scoring requires the pretrained ESM-2 checkpoint to be present or downloadable.
