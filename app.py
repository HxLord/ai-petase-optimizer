"""
AI-PETase Optimizer: AI-Assisted Research Application for PET-Degrading Enzymes.
Integrates UniProt sequence retrieval, ESM-2 protein language model inference,
3D macromolecular structure visualization (3Dmol.js), curated literature benchmarks,
and transparent circularity scenario modeling.
"""

import os
import streamlit as st
import pandas as pd
import plotly.graph_objects as go

# Core scientific modules
from src.protein_data import (
    fetch_uniprot_sequence,
    load_user_fasta,
    DEFAULT_UNIPROT_ACCESSION,
    VERIFIED_ISPETASE_FULL,
    VERIFIED_ISPETASE_MATURE,
)
from src.protein_ai import (
    score_mutation,
    MODEL_IDENTIFIER,
)
from src.protein_structure import (
    fetch_structure_pdb,
    generate_3dmol_html,
    DEFAULT_PDB_ID,
)
from src.variant_analysis import (
    load_published_variants,
    build_mutation_history_table,
    create_tm_comparison_figure,
    inspect_site_features,
)
from src.eco_calculator import (
    calculate_scenario,
    create_scenario_figure,
)
from src.utils import (
    format_number,
    format_log_prob,
    dataframe_to_csv_bytes,
    STANDARD_AMINO_ACIDS,
    ONE_TO_NAME,
)

# ---------------------------------------------------------
# Streamlit Application Configuration
# ---------------------------------------------------------
st.set_page_config(
    page_title="AI-PETase Optimizer | Enzyme Research Dashboard",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)


@st.cache_data(ttl=3600, show_spinner=False)
def fetch_default_uniprot_record():
    return fetch_uniprot_sequence(
        accession=DEFAULT_UNIPROT_ACCESSION,
        timeout_sec=15,
        use_cache=False,
    )


# ---------------------------------------------------------
# Scientific SaaS Dashboard CSS
# ---------------------------------------------------------
CUSTOM_CSS = """
<style>
@import url('https://fonts.googleapis.com/css2?family=DM+Sans:wght@400;500;600;700&family=Newsreader:opsz,wght@6..72,500;6..72,600&family=IBM+Plex+Mono:wght@400;500&display=swap');

:root {
    --paper: #f4f5f0;
    --surface: #ffffff;
    --ink: #172b29;
    --muted: #536763;
    --line: #d9e2dd;
    --teal: #087e70;
    --teal-dark: #075e55;
    --teal-wash: #e2f1ec;
    --amber: #8a5a16;
}
html, body, [class*="css"] { font-family: 'DM Sans', -apple-system, BlinkMacSystemFont, sans-serif; }
code, pre, .stCode, [data-testid="stMetricValue"] { font-family: 'IBM Plex Mono', monospace !important; }
[data-testid="stAppViewContainer"], [data-testid="stMain"], [data-testid="stMainBlockContainer"] { background: var(--paper); color:var(--ink); }
[data-testid="stVerticalBlock"] > [data-testid="stElementContainer"] { color:var(--ink); }
[data-testid="stHeader"] { background: rgba(244,245,240,.92); }
[data-testid="stMain"] > div { padding-top: 2rem; }
[data-testid="stMainBlockContainer"] { max-width: 1480px; padding: 2.4rem clamp(1.25rem, 3vw, 3rem) 4rem; }

/* Research notebook masthead */
.saas-header {
    position: relative; overflow: hidden;
    background: #173c37; color: #f8fbf7;
    border: 1px solid #24544c; border-radius: 18px;
    padding: clamp(26px, 4vw, 44px) clamp(24px, 4vw, 48px);
    margin: 0 0 24px; box-shadow: 0 12px 28px rgba(22,53,47,.10);
}
.saas-header:after { content: ''; position:absolute; width:240px; height:240px; right:-62px; top:-126px; border:1px solid rgba(216,238,224,.20); border-radius:50%; box-shadow:0 0 0 24px rgba(216,238,224,.035),0 0 0 52px rgba(216,238,224,.025); }
.saas-header h1 { position:relative; z-index:1; color:#fbfcf8; font-family:'Newsreader', Georgia, serif; font-size:clamp(2rem,4vw,3rem); font-weight:500; line-height:1.05; margin:0 0 9px; letter-spacing:-.025em; }
.saas-header p { position:relative; z-index:1; max-width:72ch; color:#d2e5dc; font-size:1rem; line-height:1.55; margin:0; font-weight:400; }

h1, h2, h3 { color:var(--ink); letter-spacing:-.025em; }
h2, h3 { font-family:'Newsreader', Georgia, serif; font-weight:600; }
[data-testid="stMarkdownContainer"] p, [data-testid="stMarkdownContainer"] li { color:#344944; line-height:1.68; }
[data-testid="stMarkdownContainer"] a { color:var(--teal-dark); text-underline-offset:3px; }
[data-testid="stMarkdownContainer"] a:hover { color:#064b44; }
[data-testid="stMarkdownContainer"] h2 { font-family:'Newsreader', Georgia, serif; font-size:clamp(1.45rem,2.5vw,1.9rem); line-height:1.18; margin:2rem 0 .65rem; }
[data-testid="stMarkdownContainer"] h3 { font-family:'DM Sans', sans-serif; font-size:1.02rem; line-height:1.35; font-weight:700; margin:1.45rem 0 .45rem; color:#1c3933; }
[data-testid="stMarkdownContainer"] h4 { color:#087e70; font-size:.82rem; font-weight:700; text-transform:uppercase; letter-spacing:.065em; }
[data-testid="stMarkdownContainer"] blockquote { border-left:2px solid #168575; color:#405851; background:#edf4ef; padding:.75rem 1rem; border-radius:0 10px 10px 0; }
[data-testid="stMarkdownContainer"] pre { background:#183a35; border:1px solid #28554a; border-radius:12px; padding:18px; }
[data-testid="stMarkdownContainer"] pre code { color:#e3f2e8; font-size:.84rem; line-height:1.65; }
/* Override Streamlit's inherited dark-theme text tokens on every light surface. */
[data-testid="stCaptionContainer"], [data-testid="stCaptionContainer"] *,
[data-testid="stWidgetLabel"], [data-testid="stWidgetLabel"] *,
[data-testid="stMarkdownContainer"] small,
[data-testid="stFileUploaderDropzone"] *,
[data-testid="stText"] { color:#405851 !important; }
[data-testid="stWidgetLabel"] p, [data-testid="stWidgetLabel"] label { color:#263d37 !important; font-weight:600; }
[data-testid="stHelp"] { color:#536763 !important; }
[data-testid="stTextInput"] input, [data-testid="stNumberInput"] input,
[data-testid="stTextArea"] textarea, [data-baseweb="select"] *,
[data-baseweb="popover"] *, [role="listbox"] *, [role="option"] {
    color:#172b29 !important;
}
[data-testid="stTextArea"] textarea:disabled { color:#405851 !important; -webkit-text-fill-color:#405851; opacity:1; background:#edf1ed; }
[data-testid="stCheckbox"] label, [data-testid="stRadio"] label,
[data-testid="stToggle"] label, [data-testid="stSlider"] label { color:#263d37 !important; }
[data-testid="stAlert"] [data-testid="stMarkdownContainer"],
[data-testid="stAlert"] [data-testid="stMarkdownContainer"] * { color:#263d37 !important; }
[data-testid="stAlert"] svg { fill:currentColor; }
[data-testid="stDataFrame"] *, [data-testid="stTable"] * { color:#243b35; }
::selection { color:#102b27; background:#bde3d4; }
:focus-visible { outline:3px solid #168575 !important; outline-offset:2px !important; }

.metric-box { background:var(--surface); border:1px solid var(--line); border-radius:14px; padding:18px 20px; margin:0 0 12px; box-shadow:0 3px 12px rgba(24,50,43,.035); }
.metric-box-title { color:#465b55; font-size:.75rem; font-weight:700; text-transform:uppercase; letter-spacing:.09em; margin-bottom:9px; }
.metric-box-val { color:var(--ink); font-family:'Newsreader', Georgia, serif; font-size:clamp(1.55rem,2.4vw,2rem); font-weight:600; line-height:1.15; }
.metric-box-sub { color:#405851; font-size:.82rem; line-height:1.5; margin-top:9px; }
.disclaimer-banner { background:#fff5e8; border:1px solid #e5c896; border-radius:12px; padding:15px 18px; margin:16px 0 22px; font-size:.9rem; color:#503814; line-height:1.65; }
.info-banner { background:var(--teal-wash); border:1px solid #afd5c8; border-radius:12px; padding:15px 18px; margin:16px 0 22px; font-size:.9rem; color:#183f37; line-height:1.65; }
.badge-experimental, .badge-ai, .badge-scenario { display:inline-block; border-radius:999px; padding:3px 9px; font-size:.73rem; font-weight:700; letter-spacing:.01em; }
.badge-experimental { background:#dcefe5; color:#1d5c40; border:1px solid #a8cbb5; }
.badge-ai { background:#eee8f7; color:#57427b; border:1px solid #c9b9e2; }
.badge-scenario { background:#fff0d8; color:#70470e; border:1px solid #e2bf82; }

/* Streamlit controls */
section[data-testid="stSidebar"] { background:#e9eee9; border-right:1px solid #d4ded8; }
section[data-testid="stSidebar"] > div { padding:1.35rem 1rem; }
section[data-testid="stSidebar"] h2 { font-family:'Newsreader', Georgia, serif; color:#0c695e !important; font-size:1.42rem !important; line-height:1.1; }
section[data-testid="stSidebar"] [data-testid="stMarkdownContainer"] p { color:#405851; font-size:.82rem; line-height:1.55; }
section[data-testid="stSidebar"] hr { border-color:#cbd8d0; margin:1rem 0; }
section[data-testid="stSidebar"] [role="radiogroup"] { gap:4px; }
section[data-testid="stSidebar"] [role="radiogroup"] > label { min-height:42px; padding:8px 10px; border-radius:9px; transition:background .15s ease; }
section[data-testid="stSidebar"] [role="radiogroup"] > label:hover { background:#dce8e1; }
section[data-testid="stSidebar"] [role="radiogroup"] > label:has(input:checked) { background:#d5e9df; }
section[data-testid="stSidebar"] [role="radiogroup"] > label p { color:#243e37; font-weight:600; }
section[data-testid="stSidebar"] [data-testid="stMarkdownContainer"] code { color:#155a50; background:#dcebe4; }
div[role="radiogroup"] label, [data-testid="stRadio"] label { color:#263d37; }
div[data-baseweb="select"] > div, div[data-baseweb="input"] > div, [data-testid="stFileUploaderDropzone"] { background:#fff; border-color:#bdcdc5; border-radius:10px; }
div[data-baseweb="select"] > div:focus-within, div[data-baseweb="input"] > div:focus-within { border-color:#087e70; box-shadow:0 0 0 2px rgba(8,126,112,.14); }
input, textarea { color:#172b29 !important; }
button[kind="primary"], [data-testid="stBaseButton-primary"] { background:#087e70; border:1px solid #087e70; border-radius:9px; color:white; font-weight:700; min-height:42px; transition:background .16s ease, transform .16s ease; }
button[kind="primary"]:hover, [data-testid="stBaseButton-primary"]:hover { background:#075e55; border-color:#075e55; transform:translateY(-1px); }
button[kind="secondary"], [data-testid="stBaseButton-secondary"] { min-height:40px; border-radius:9px; border-color:#9eb9ae; color:#145b51; background:#fbfcf9; }
[data-testid="stDownloadButton"] button { border-radius:9px; font-weight:700; }
[data-testid="stMetric"] { background:#fff; border:1px solid var(--line); border-radius:13px; padding:14px 16px; box-shadow:0 3px 12px rgba(24,50,43,.035); }
[data-testid="stMetricLabel"] { color:#465b55; }
[data-testid="stMetricValue"] { color:#173c37; }
[data-testid="stDataFrame"], [data-testid="stTable"] { border:1px solid var(--line); border-radius:12px; overflow:hidden; background:#fff; }
[data-testid="stPlotlyChart"] { background:var(--surface); border:1px solid var(--line); border-radius:14px; padding:10px; box-shadow:0 3px 12px rgba(24,50,43,.035); }
[data-testid="stAlert"] { border-radius:11px; border-width:1px; }
[data-testid="stExpander"] { border:1px solid var(--line); background:#fff; border-radius:11px; }
[data-testid="stExpander"] summary { color:#1a4038; font-weight:650; }
[data-testid="stProgressBar"] > div > div { background:#0a8474; }
hr { border-color:#d5dfd9; margin:1.6rem 0; }
::-webkit-scrollbar { width:10px; height:10px; }
::-webkit-scrollbar-thumb { background:#afc1b8; border:2px solid var(--paper); border-radius:8px; }
::-webkit-scrollbar-thumb:hover { background:#829d91; }

@media (max-width: 760px) {
    [data-testid="stMainBlockContainer"] { padding:1.25rem 1rem 2.5rem; }
    .saas-header { border-radius:14px; padding:24px 20px; }
    .saas-header:after { right:-150px; }
    .metric-box { padding:15px 16px; }
    [data-testid="stMarkdownContainer"] h2 { margin-top:1.6rem; }
}
</style>
"""
st.markdown(CUSTOM_CSS, unsafe_allow_html=True)

# ---------------------------------------------------------
# Session State Initialization
# ---------------------------------------------------------
if "sequence_source" not in st.session_state:
    st.session_state["sequence_source"] = "IsPETase Mature (5XJH, 261 aa)"
if "active_sequence" not in st.session_state:
    st.session_state["active_sequence"] = VERIFIED_ISPETASE_MATURE
if "active_sequence_name" not in st.session_state:
    st.session_state["active_sequence_name"] = "IsPETase Mature Catalytic Domain"
if "active_accession" not in st.session_state:
    st.session_state["active_accession"] = DEFAULT_UNIPROT_ACCESSION
if "active_sequence_offset" not in st.session_state:
    st.session_state["active_sequence_offset"] = 29
if "active_sequence_is_custom" not in st.session_state:
    st.session_state["active_sequence_is_custom"] = False
if "active_sequence_source" not in st.session_state:
    st.session_state["active_sequence_source"] = "Bundled UniProt reference sequence"
if "mutation_history" not in st.session_state:
    st.session_state["mutation_history"] = []

# ---------------------------------------------------------
# Sidebar Navigation
# ---------------------------------------------------------
with st.sidebar:
    st.markdown("""
    <div style="padding: 10px 0 16px 0;">
        <div style="width:34px;height:4px;background:#087e70;border-radius:4px;margin:3px 0 12px;"></div>
        <h2 style="margin: 0; font-size: 1.3rem; font-weight: 600;">AI-PETase Optimizer</h2>
        <p style="font-size: 0.8rem; margin: 5px 0 0 0;">PET enzyme research platform</p>
    </div>
    """, unsafe_allow_html=True)

    page = st.radio(
        "Navigation",
        [
            "1. Overview",
            "2. Protein Explorer",
            "3. AI Mutation Analysis",
            "4. Variant Comparison",
            "5. Environmental Impact",
            "6. Research & Methodology"
        ],
        index=0
    )

    st.markdown("---")
    st.markdown("### System status")
    st.markdown(f"""
    - **Active Sequence**: `{st.session_state['active_sequence_name']}`
    - **Length**: `{len(st.session_state['active_sequence'])} aa`
    - **Sequence Source**: `{st.session_state['active_sequence_source']}`
    - **AI Model**: `{MODEL_IDENTIFIER}` (loads on first score)
    - **PDB Reference**: `5XJH` (1.58 Å)
    """)

    st.markdown("---")
    st.markdown("""
    <div style="font-size: 0.75rem; color: #405851; line-height:1.55;">
        <strong>Research Prototype Notice:</strong><br>
        Computational scores represent language model likelihood ratios and are not direct physical measurements.
    </div>
    """, unsafe_allow_html=True)


# =========================================================
# PAGE 1: Overview
# =========================================================
if page == "1. Overview":
    st.markdown("""
    <div class="saas-header">
        <h1>AI-PETase Optimizer</h1>
        <p>AI-assisted research for biocatalytic poly(ethylene terephthalate) (PET) plastic depolymerization</p>
    </div>
    """, unsafe_allow_html=True)

    st.markdown("""
    <div class="disclaimer-banner">
        <strong>Scientific prototype disclosure:</strong> AI-PETase Optimizer is an academic and applied research MVP. 
        Language model likelihood scores reflect evolutionary sequence compatibility, not direct measurements of melting 
        temperature (Tm) or catalytic reaction velocity. Experimental validation is strictly required before laboratory deployment.
    </div>
    """, unsafe_allow_html=True)

    col1, col2 = st.columns([3, 2])

    with col1:
        st.markdown("### The PET plastic challenge & biocatalysis")
        st.markdown("""
        Poly(ethylene terephthalate) (PET) is widely used in packaging and textiles. Conventional mechanical recycling can
        affect polymer properties, depending on the feedstock and process.

        In 2016, the discovery of ***Ideonella sakaiensis* PETase (IsPETase)** revealed a natural biological mechanism 
        to hydrolyze ester bonds in PET. The products depend on reaction conditions and can include the intermediate MHET
        as well as terephthalic acid and ethylene glycol.
        
        A reported IsPETase melting temperature is approximately **48.8°C** under the source study's assay conditions.
        This is not a universal denaturation threshold or a direct measure of PET degradation rate. PET glass-transition
        behavior also depends on polymer composition and thermal history.
        """)

        st.markdown("### Research workflow")
        st.markdown("""
        ```
        [ UniProt Sequence ] ──> [ 3D Structure 5XJH ] ──> [ ESM-2 Language Model ]
                 │                          │                           │
                 ▼                          ▼                           ▼
        Sequence Explorer          Active-Site Landmarks        Masked Marginal Scoring (Δ log p)
                 │                          │                           │
                 └──────────────────────────┴───────────────────────────┘
                                            │
                                            ▼
                           [ Variant Benchmark Comparison ]
                                            │
                                            ▼
                           [ Circularity Scenario Modeling ]
        ```
        """)

    with col2:
        st.markdown("### Reference measurements")
        
        st.markdown("""
        <div class="metric-box">
            <div class="metric-box-title">Reported Wild-Type IsPETase Tm</div>
            <div class="metric-box-val">48.8 °C</div>
            <div class="metric-box-sub"><span class="badge-experimental">Reported</span> Han et al., Nat Commun 2017; assay-dependent</div>
        </div>
        """, unsafe_allow_html=True)

        st.markdown("""
        <div class="metric-box">
            <div class="metric-box-title">Reported FAST-PETase Tm</div>
            <div class="metric-box-val">67.8 °C (+19.0°C)</div>
            <div class="metric-box-sub"><span class="badge-experimental">Reported</span> Lu et al., Nature 2022; assay-dependent</div>
        </div>
        """, unsafe_allow_html=True)

        st.markdown("""
        <div class="metric-box">
            <div class="metric-box-title">PET Glass Transition Range (Tg)</div>
            <div class="metric-box-val">65 – 75 °C</div>
            <div class="metric-box-sub">Approximate range; depends on PET composition and thermal history</div>
        </div>
        """, unsafe_allow_html=True)

    st.markdown("---")
    st.markdown("### Explore the platform")
    nav_c1, nav_c2, nav_c3, nav_c4 = st.columns(4)
    with nav_c1:
        st.markdown("**2. Protein Explorer**")
        st.write("Inspect IsPETase sequences and interactive 3D WebGL crystal structure (PDB 5XJH).")
    with nav_c2:
        st.markdown("**3. AI Mutation Analysis**")
        st.write("Propose amino acid substitutions and compute ESM-2 log-likelihood ratios.")
    with nav_c3:
        st.markdown("**4. Variant Comparison**")
        st.write("Compare engineered benchmarks (FAST-PETase, DuraPETase, LCC-ICCG) from literature.")
    with nav_c4:
        st.markdown("**5. Environmental Impact**")
        st.write("Model material mass balance and carbon emission scenarios.")


# =========================================================
# PAGE 2: Protein Explorer
# =========================================================
elif page == "2. Protein Explorer":
    st.markdown("""
    <div class="saas-header">
        <h1>Protein explorer</h1>
        <p>Verified sequence metadata and interactive 3D macromolecular structure (PDB 5XJH)</p>
    </div>
    """, unsafe_allow_html=True)

    # Sequence selection controls
    st.markdown("### Choose a sequence")
    seq_col1, seq_col2 = st.columns([2, 3])

    with seq_col1:
        seq_choice = st.radio(
            "Available Sequence Profiles:",
            [
                "Mature IsPETase (UniProt residues 30-290, 261 aa; bundled reference)",
                "Full IsPETase precursor (290 aa; bundled reference)",
                "Fetch full precursor from UniProt (live)",
                "Upload Custom FASTA File"
            ],
            index=0
        )

        if seq_choice == "Mature IsPETase (UniProt residues 30-290, 261 aa; bundled reference)":
            st.session_state["active_sequence"] = VERIFIED_ISPETASE_MATURE
            st.session_state["active_sequence_name"] = "IsPETase Mature Catalytic Domain (5XJH)"
            st.session_state["active_accession"] = DEFAULT_UNIPROT_ACCESSION
            st.session_state["active_sequence_offset"] = 29
            st.session_state["active_sequence_is_custom"] = False
            st.session_state["active_sequence_source"] = "Bundled UniProt reference sequence"
            is_custom = False
        elif seq_choice == "Full IsPETase precursor (290 aa; bundled reference)":
            st.session_state["active_sequence"] = VERIFIED_ISPETASE_FULL
            st.session_state["active_sequence_name"] = "IsPETase Full Precursor (UniProt A0A0K8P6T7)"
            st.session_state["active_accession"] = DEFAULT_UNIPROT_ACCESSION
            st.session_state["active_sequence_offset"] = 0
            st.session_state["active_sequence_is_custom"] = False
            st.session_state["active_sequence_source"] = "Bundled UniProt reference sequence"
            is_custom = False
        elif seq_choice == "Fetch full precursor from UniProt (live)":
            is_custom = False
            try:
                record = fetch_default_uniprot_record()
                st.session_state["active_sequence"] = record["sequence"]
                st.session_state["active_sequence_name"] = f"IsPETase ({record['accession']})"
                st.session_state["active_accession"] = record["accession"]
                st.session_state["active_sequence_offset"] = 0
                st.session_state["active_sequence_is_custom"] = False
                st.session_state["active_sequence_source"] = record["source"]
                if record["is_fallback"]:
                    warning = record.get("retrieval_warning") or "Live UniProt response was unavailable."
                    st.warning(f"Using verified bundled sequence; live UniProt fetch failed: {warning}")
                else:
                    st.success(f"Loaded {record['length']} residues from {record['source']}.")
            except (RuntimeError, ValueError, OSError) as e:
                st.error(f"Could not load the UniProt sequence: {e}")
        else:
            is_custom = st.session_state["active_sequence_is_custom"]
            uploaded_file = st.file_uploader("Upload FASTA file (.fasta, .fa, .txt):", type=["fasta", "fa", "txt"])
            if uploaded_file is not None:
                try:
                    parsed_user = load_user_fasta(uploaded_file)
                    st.session_state["active_sequence"] = parsed_user["sequence"]
                    st.session_state["active_sequence_name"] = f"Custom: {parsed_user['name']}"
                    st.session_state["active_accession"] = parsed_user["accession"]
                    st.session_state["active_sequence_offset"] = None
                    st.session_state["active_sequence_is_custom"] = True
                    st.session_state["active_sequence_source"] = "User-uploaded FASTA"
                    is_custom = True
                    st.success(f"Loaded custom sequence '{parsed_user['name']}' ({parsed_user['length']} aa).")
                except (ValueError, OSError) as e:
                    st.error(f"Failed to parse uploaded FASTA: {str(e)}")

    with seq_col2:
        st.markdown(f"**Active Record**: `{st.session_state['active_sequence_name']}`")
        st.markdown(f"**Length**: `{len(st.session_state['active_sequence'])} amino acids`")
        st.caption(f"Sequence source: {st.session_state['active_sequence_source']}")
        
        with st.expander("View Amino Acid Sequence (FASTA)", expanded=True):
            st.text_area(
                "Sequence",
                st.session_state["active_sequence"],
                height=130,
                disabled=True
            )

    st.markdown("---")
    st.markdown("### Crystal structure · PDB 5XJH")

    if is_custom:
        st.markdown("""
        <div class="info-banner">
            <strong>ℹ️ Structural Representation Notice:</strong> You have loaded a custom sequence. 
            The 3D viewer below displays the experimental wild-type crystal structure of <em>Ideonella sakaiensis</em> 
            PETase (PDB 5XJH) as a structural reference. It does <strong>not</strong> represent the folded 3D conformation of your custom sequence.
        </div>
        """, unsafe_allow_html=True)

    # 3D Viewer options
    v_col1, v_col2, v_col3, v_col4 = st.columns(4)
    with v_col1:
        style_opt = st.selectbox("Render Style:", ["cartoon", "stick", "sphere", "surface"], index=0)
    with v_col2:
        color_opt = st.selectbox("Color Scheme:", ["spectrum", "chain", "ss"], index=0)
    with v_col3:
        spin_opt = st.checkbox("Auto-Rotate (Spin)", value=False)
    with v_col4:
        landmark_opts = st.multiselect(
            "Highlight Landmarks:",
            ["catalytic_triad", "disulfide_bonds", "oxyanion_hole", "substrate_cleft"],
            default=["catalytic_triad", "disulfide_bonds"]
        )

    # Fetch PDB
    pdb_text, pdb_src, is_fb = fetch_structure_pdb(DEFAULT_PDB_ID)
    if pdb_text:
        html_code = generate_3dmol_html(
            pdb_text=pdb_text,
            style=style_opt,
            color_scheme=color_opt,
            highlight_features=landmark_opts,
            spin=spin_opt,
            height=540
        )
        st.iframe(html_code, height=560)
        st.caption(f"Structure Source: {pdb_src} • Coordinates: RCSB PDB 5XJH (1.58 Å resolution, X-ray diffraction, Han et al. 2017)")
    else:
        st.error("Unable to load 3D structure coordinates. Please ensure `data/structures/5xjh.pdb` is accessible.")


# =========================================================
# PAGE 3: AI Mutation Analysis
# =========================================================
elif page == "3. AI Mutation Analysis":
    st.markdown("""
    <div class="saas-header">
        <h1>Mutation analysis</h1>
        <p>Pretrained ESM-2 masked-marginal log-likelihood scoring of single-point substitutions</p>
    </div>
    """, unsafe_allow_html=True)

    st.markdown("""
    <div class="disclaimer-banner">
        <strong>Scientific caution:</strong> Protein language model log-odds scores quantify sequence 
        context compatibility learned across millions of natural proteins. A positive Δ log p does <strong>not</strong> 
        guarantee increased melting temperature (Tm), improved catalytic turnover on macroscopic PET, or chemical stability.
    </div>
    """, unsafe_allow_html=True)

    seq = st.session_state["active_sequence"]
    seq_len = len(seq)

    st.markdown(f"**Target Protein**: `{st.session_state['active_sequence_name']}` (Length: `{seq_len}` residues)")
    st.caption("Enter a 1-based index in the selected sequence. IsPETase reference positions are shown separately.")

    col_input, col_info = st.columns([2, 3])

    with col_input:
        st.markdown("#### Configure Point Mutation")
        
        pos = st.number_input(
            "Selected-sequence position (1-based):",
            min_value=1,
            max_value=seq_len,
            value=min(131 if st.session_state["active_sequence_offset"] == 29 else 160, seq_len),
            step=1,
            key=f"mutation_position_{st.session_state['active_sequence_name']}"
        )
        
        current_wt = seq[pos - 1]
        reference_offset = st.session_state["active_sequence_offset"]
        reference_position = pos + reference_offset if reference_offset is not None else None
        st.markdown(f"**Original Residue at sequence position {pos}**: `{current_wt}` ({ONE_TO_NAME.get(current_wt, 'Unknown')})")
        if reference_position is not None and not st.session_state["active_sequence_is_custom"]:
            st.caption(f"Corresponding UniProt/PDB residue: {reference_position}")

        # Check if known landmark
        site_feat = (
            inspect_site_features(reference_position)
            if reference_position is not None and not st.session_state["active_sequence_is_custom"]
            else None
        )
        if site_feat:
            st.markdown(f"""
            <div style="background: #fff1e8; border: 1px solid #d9a47d; border-radius: 10px; padding: 10px 13px; margin: 10px 0; font-size: 0.86rem; color:#553b2a;">
                <strong style="color:#8a4c23;">Active-site landmark:</strong> {site_feat['region']}<br>
                {site_feat['role']}
            </div>
            """, unsafe_allow_html=True)

        # Mutant amino acid selection (canonical 20 standard residues)
        available_mutants = sorted(list(STANDARD_AMINO_ACIDS - {current_wt}))
        mutant_options = [f"{aa} - {ONE_TO_NAME.get(aa, '')}" for aa in available_mutants]
        
        selected_mut_str = st.selectbox(
            "Proposed Mutant Amino Acid:",
            mutant_options,
            index=0
        )
        target_mutant_aa = selected_mut_str.split(" - ")[0]

        run_btn = st.button("Score mutation with ESM-2", type="primary", width="stretch")

    with col_info:
        st.markdown("#### Scientific Scoring Engine")
        st.markdown(f"""
        - **Model Architecture**: Meta `fair-esm` / `{MODEL_IDENTIFIER}`
        - **Scoring Method**: Masked Marginal Log-Likelihood Ratio
        - **Equation**: $\\Delta \\log p = \\log p(x_{{\\text{{pos}}}} = \\text{{mutant}}) - \\log p(x_{{\\text{{pos}}}} = \\text{{WT}})$
        - **Interpretation**: the score compares model compatibility for two amino acids in the masked context; it is not a functional verdict.
        """)

    # Run inference
    if run_btn:
        with st.spinner(f"Computing ESM-2 masked inference for substitution {current_wt}{pos}{target_mutant_aa}..."):
            try:
                res = score_mutation(
                    sequence=seq,
                    position=pos,
                    mutant_aa=target_mutant_aa,
                    allow_synonymous=False
                )
                res["sequence_position"] = pos
                res["reference_position"] = reference_position
                res["reference_accession"] = (
                    st.session_state["active_accession"]
                    if not st.session_state["active_sequence_is_custom"]
                    else None
                )
                if reference_position is not None and not st.session_state["active_sequence_is_custom"]:
                    res["mutation"] = f"{res['wt_aa']}{reference_position}{res['mutant_aa']}"
                # Store in history
                st.session_state["mutation_history"].insert(0, res)
                
                st.success(f"Inference completed for mutation **{res['mutation']}**!")

                # Result cards
                rc1, rc2, rc3 = st.columns(3)
                with rc1:
                    st.markdown(f"""
                    <div class="metric-box">
                        <div class="metric-box-title">Wild-Type Log-Prob</div>
                        <div class="metric-box-val">{format_log_prob(res['wt_log_prob'])}</div>
                        <div class="metric-box-sub">Original reference residue: {res['wt_aa']}{reference_position or pos}</div>
                    </div>
                    """, unsafe_allow_html=True)
                with rc2:
                    st.markdown(f"""
                    <div class="metric-box">
                        <div class="metric-box-title">Mutant Log-Prob</div>
                        <div class="metric-box-val">{format_log_prob(res['mutant_log_prob'])}</div>
                        <div class="metric-box-sub">Mutant: {res['mutant_aa']}</div>
                    </div>
                    """, unsafe_allow_html=True)
                with rc3:
                    score_color = "#087e70" if res['score_difference'] >= 0 else "#a94e3b"
                    st.markdown(f"""
                    <div class="metric-box">
                        <div class="metric-box-title">Score Difference (Δ log p)</div>
                        <div class="metric-box-val" style="color: {score_color};">{res['score_difference']:+.4f}</div>
                        <div class="metric-box-sub">Model: {res['model_name']}</div>
                    </div>
                    """, unsafe_allow_html=True)

                st.markdown(f"""
                <div class="info-banner">
                    <strong>Interpretation:</strong> {res['interpretation']}
                </div>
                """, unsafe_allow_html=True)

            except (ValueError, RuntimeError, OSError) as e:
                st.error(f"ESM-2 inference failed: {str(e)}")

    # History Table
    st.markdown("---")
    st.markdown("### Scoring history")
    
    if st.session_state["mutation_history"]:
        history_df = build_mutation_history_table(st.session_state["mutation_history"])
        st.dataframe(history_df, width="stretch")

        csv_data = dataframe_to_csv_bytes(history_df)
        st.download_button(
            label="Export mutation results (CSV)",
            data=csv_data,
            file_name="esm_mutation_scores.csv",
            mime="text/csv"
        )
    else:
        st.info("No mutations scored in this session yet. Propose a mutation above and click 'Score Mutation with ESM-2'.")


# =========================================================
# PAGE 4: Variant Comparison
# =========================================================
elif page == "4. Variant Comparison":
    st.markdown("""
    <div class="saas-header">
        <h1>Curated variant benchmarks</h1>
        <p>Reported literature benchmarks with source links and assay context</p>
    </div>
    """, unsafe_allow_html=True)

    st.markdown("""
    <div class="info-banner">
        <strong>Data Provenance Notice:</strong> Entries link to their cited primary papers, but assay conditions and
        substrates differ and the values have not been independently reanalyzed here. Compare measurements only with
        their original methods in view. Unverified values are left blank; model estimates are never substituted for them.
    </div>
    """, unsafe_allow_html=True)

    variants_df = load_published_variants()

    if not variants_df.empty:
        # Plotly Tm Bar Chart
        fig_tm = create_tm_comparison_figure(variants_df)
        st.plotly_chart(fig_tm, width="stretch")

        # Tabular View
        st.markdown("### Variant data")
        
        display_cols = [
            "name", "scaffold", "pdb_id", "mutations", "mutation_count",
            "tm_celsius", "tm_delta", "tm_delta_reference", "opt_temp_celsius", "opt_ph",
            "reported_degradation_metric", "citation", "doi"
        ]
        available_display_cols = [c for c in display_cols if c in variants_df.columns]
        
        st.dataframe(
            variants_df[available_display_cols].rename(columns={
                "name": "Variant Name",
                "scaffold": "Parent Scaffold",
                "pdb_id": "PDB",
                "mutations": "Substitutions",
                "mutation_count": "Count",
                "tm_celsius": "Tm (°C)",
                "tm_delta": "ΔTm (°C)",
                "tm_delta_reference": "ΔTm Reference",
                "opt_temp_celsius": "Topt (°C)",
                "opt_ph": "Opt pH",
                "reported_degradation_metric": "Reported Activity / Yield",
                "citation": "Source Citation",
                "doi": "DOI"
            }),
            width="stretch"
        )

        st.markdown("---")
        st.markdown("### How to read these measurements")
        bc1, bc2, bc3 = st.columns(3)
        with bc1:
            st.markdown("**1. Melting Temperature ($T_m$)**")
            st.write("A measured unfolding midpoint under defined conditions; not a universal denaturation threshold and not a direct activity measurement.")
        with bc2:
            st.markdown("**2. Optimal Temperature ($T_{\\text{opt}}$)**")
            st.write("An assay-dependent maximum in a specified activity experiment. It depends on substrate, time, pH, and enzyme inactivation.")
        with bc3:
            st.markdown("**3. Polymer Crystallinity ($T_g$)**")
            st.write("PET chain mobility changes around a composition- and history-dependent glass transition. Proximity to $T_g$ alone does not establish degradation performance.")

    else:
        st.error("No published variant dataset found at `data/variants.csv`.")


# =========================================================
# PAGE 5: Environmental Impact
# =========================================================
elif page == "5. Environmental Impact":
    st.markdown("""
    <div class="saas-header">
        <h1>Environmental scenario calculator</h1>
        <p>Hypothetical material mass balances and life cycle greenhouse gas emission modeling</p>
    </div>
    """, unsafe_allow_html=True)

    st.markdown("""
    <div class="disclaimer-banner">
        <strong>Hypothetical scenario:</strong> The calculations below are mathematical scenarios 
        based strictly on user-specified recovery rates and life-cycle emission factors. They do <strong>not</strong> represent 
        proven industrial yields or laboratory claims for any engineered enzyme.
    </div>
    """, unsafe_allow_html=True)

    # Input parameters
    col_in1, col_in2, col_in3 = st.columns(3)
    with col_in1:
        pet_waste_mass = st.number_input(
            "Total PET Waste Batch Mass (kg):",
            min_value=0.0,
            value=10000.0,
            step=500.0
        )
        base_rate = st.slider(
            "Baseline Recovery Rate (%):",
            min_value=0.0,
            max_value=100.0,
            value=15.0,
            step=1.0
        )

    with col_in2:
        hypo_rate = st.slider(
            "Hypothetical Alternative Recovery Rate (%):",
            min_value=0.0,
            max_value=100.0,
            value=75.0,
            step=1.0
        )
        base_time = st.number_input(
            "Baseline Processing Time (hours, optional):",
            min_value=0.0,
            value=48.0,
            step=4.0
        )

    with col_in3:
        hypo_time = st.number_input(
            "Hypothetical Processing Time (hours, optional):",
            min_value=0.0,
            value=16.0,
            step=2.0
        )
        co2_factor = st.number_input(
            "Assumed net GHG difference (kg CO2e / kg recovered PET):",
            min_value=0.0,
            value=1.25,
            step=0.1
        )

    st.caption(
        "The default recovery rates, processing times, and 1.25 kg CO2e/kg factor are illustrative scenario inputs, "
        "not measured enzyme performance or a project-specific life-cycle assessment."
    )

    # Perform calculation
    try:
        eco_results = calculate_scenario(
            pet_waste_kg=pet_waste_mass,
            baseline_rate_pct=base_rate,
            hypothetical_rate_pct=hypo_rate,
            baseline_time_hours=base_time if base_time > 0 else None,
            hypothetical_time_hours=hypo_time if hypo_time > 0 else None,
            emissions_factor_kg_co2_per_kg=co2_factor
        )

        st.markdown("---")
        st.markdown("### Scenario results")

        res_c1, res_c2, res_c3, res_c4 = st.columns(4)
        with res_c1:
            st.markdown(f"""
            <div class="metric-box">
                <div class="metric-box-title">Baseline Recovered PET</div>
                <div class="metric-box-val">{eco_results['baseline_recovered_kg']:,.1f} kg</div>
                <div class="metric-box-sub">At {eco_results['baseline_rate_pct']}% recovery</div>
            </div>
            """, unsafe_allow_html=True)
        with res_c2:
            st.markdown(f"""
            <div class="metric-box">
                <div class="metric-box-title">Hypothetical Recovered PET</div>
                <div class="metric-box-val">{eco_results['hypothetical_recovered_kg']:,.1f} kg</div>
                <div class="metric-box-sub">At {eco_results['hypothetical_rate_pct']}% recovery</div>
            </div>
            """, unsafe_allow_html=True)
        with res_c3:
            st.markdown(f"""
            <div class="metric-box">
                <div class="metric-box-title">Scenario Mass Difference</div>
                <div class="metric-box-val" style="color:#087e70;">{eco_results['mass_difference_kg']:+,.1f} kg</div>
                <div class="metric-box-sub">Hypothetical minus baseline recovery</div>
            </div>
            """, unsafe_allow_html=True)
        with res_c4:
            em = eco_results.get("emissions_scenario")
            co2_diff = em["net_emissions_difference_kg_co2"] if em else 0.0
            st.markdown(f"""
            <div class="metric-box">
                <div class="metric-box-title">Scenario GHG Difference</div>
                <div class="metric-box-val" style="color:#087e70;">{co2_diff:,.1f} kg</div>
                <div class="metric-box-sub">Based only on the assumed factor above (kg CO2e)</div>
            </div>
            """, unsafe_allow_html=True)

        # Plotly mass balance chart
        fig_eco = create_scenario_figure(eco_results)
        st.plotly_chart(fig_eco, width="stretch")

        # CSV Download
        report_df = pd.DataFrame([{
            "Total_PET_Waste_kg": eco_results["pet_waste_kg"],
            "Baseline_Rate_pct": eco_results["baseline_rate_pct"],
            "Hypothetical_Rate_pct": eco_results["hypothetical_rate_pct"],
            "Baseline_Recovered_kg": eco_results["baseline_recovered_kg"],
            "Hypothetical_Recovered_kg": eco_results["hypothetical_recovered_kg"],
            "Scenario_Mass_Difference_kg": eco_results["mass_difference_kg"],
            "Baseline_Processing_Time_hours": base_time if base_time > 0 else "N/A",
            "Hypothetical_Processing_Time_hours": hypo_time if hypo_time > 0 else "N/A",
            "Assumed_Net_GHG_Factor_kg_CO2e_per_kg": co2_factor if em else "N/A",
            "Scenario_GHG_Difference_kg_CO2e": em["net_emissions_difference_kg_co2"] if em else "N/A",
            "Vehicle_Mile_Equivalent_signed": em["vehicle_miles_offset_equivalent"] if em else "N/A",
            "Scenario_Disclaimer": eco_results["disclaimer"]
        }])
        
        csv_report = dataframe_to_csv_bytes(report_df)
        st.download_button(
            label="Download scenario report (CSV)",
            data=csv_report,
            file_name="pet_circularity_scenario.csv",
            mime="text/csv"
        )

    except ValueError as e:
        st.error(f"Scenario calculation error: {str(e)}")


# =========================================================
# PAGE 6: Research & Methodology
# =========================================================
elif page == "6. Research & Methodology":
    st.markdown("""
    <div class="saas-header">
        <h1>Research & methodology</h1>
        <p>Scientific documentation, data sources, model architectures, and experimental citations</p>
    </div>
    """, unsafe_allow_html=True)

    st.markdown("### Primary references & data repositories")
    st.markdown("""
    - **UniProtKB IsPETase Record**: [UniProt A0A0K8P6T7](https://www.uniprot.org/uniprotkb/A0A0K8P6T7/entry)
    - **RCSB Protein Data Bank 5XJH**: [RCSB PDB Structure 5XJH](https://www.rcsb.org/structure/5XJH) (Han et al., *Nat Commun* 2017)
    - **Evolutionary Scale Modeling (ESM)**: [GitHub facebookresearch/esm](https://github.com/facebookresearch/esm) (Lin et al., *Science* 2023)
    - **FAST-PETase Landmark Engineering**: [Lu et al., *Nature* 604, 662–667 (2022)](https://www.nature.com/articles/s41586-022-04599-z)
    - **Natural Discovery of IsPETase**: [Yoshida et al., *Science* 351, 1196–1199 (2016)](https://www.science.org/doi/10.1126/science.aad6359)
    - **DuraPETase Rigidification (GRAPE)**: [Cui et al., *ACS Catalysis* 11, 738–747 (2021)](https://pubs.acs.org/doi/10.1021/acscatal.0c05126)
    - **HotPETase**: [Bell et al., *Nature Catalysis* 5, 673–681 (2022)](https://doi.org/10.1038/s41929-022-00821-3)
    - **LCC-ICCG Cutinase**: [Tournier et al., *Nature* 580, 216–219 (2020)](https://www.nature.com/articles/s41586-020-2149-4)
    """)

    st.markdown("---")
    st.markdown("### Methods & model details")
    st.markdown("""
    1. **Sequence Language Modeling (ESM-2)**:
       - AI-PETase Optimizer uses `esm2_t6_8M_UR50D` (8 million parameters, 6 layers, 320 embedding dimensions).
       - Log-likelihood differences are computed through **Masked Marginal Scoring**: the query position is replaced with the model's `<mask_idx>` token, and output logits are evaluated via `torch.log_softmax`.
       - Score difference $\\Delta \\log p = \\log p(x_{\\text{pos}} = \\text{mut}) - \\log p(x_{\\text{pos}} = \\text{wt})$.
       - *Caution*: ESM models do not calculate free energy perturbation ($\\Delta\\Delta G$) or kinetic constants ($k_{\\text{cat}}$).

    2. **Structural Visualization**:
       - 3D rendering uses **3Dmol.js WebGL** embedding.
       - Structural coordinates originate from PDB ID `5XJH` (chain A), resolved by X-ray diffraction at 1.58 Å resolution.
       - Mutating sequence positions in the UI updates sequence scores but preserves experimental backbone coordinates to prevent misrepresenting unverified conformations.

    3. **Key Industrial Engineering Bottlenecks**:
       - **Substrate Crystallinity**: Crystallinity and thermal history affect PET accessibility. A glass-transition temperature is not a universal performance threshold; results depend on polymer and assay conditions.
       - **Intermediate Products**: Hydrolysis can produce mono-(2-hydroxyethyl) terephthalate (MHET); product effects depend on the enzyme and reaction conditions.
       - **Process Inputs**: Buffer use, pH control, and heat requirements depend on process scale and operating conditions and require process-specific assessment.
    """)
