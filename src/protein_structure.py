"""
Protein Structure Visualization Module for AI-PETase Optimizer.
Fetches, caches, and renders 3D macromolecular structures using 3Dmol.js WebGL.
Integrates experimentally determined crystal structure of IsPETase (RCSB PDB 5XJH)
with clear scientific separation between experimental coordinates and sequence predictions.
"""

import os
import json
import math
import re
from typing import Dict, Any, Optional, Tuple, List
import requests

STRUCTURE_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data", "structures")
DEFAULT_PDB_ID = "5XJH"
DEFAULT_FALLBACK_FILE = os.path.join(STRUCTURE_DIR, "5xjh.pdb")


def fetch_structure_pdb(
    pdb_id: str = DEFAULT_PDB_ID,
    timeout_sec: int = 5,
    use_cache: bool = True
) -> Tuple[Optional[str], str, bool]:
    """
    Retrieve PDB coordinate file text from local cache or RCSB PDB API.
    
    Args:
        pdb_id: RCSB PDB identifier (default: '5XJH').
        timeout_sec: Network timeout for live fetch.
        use_cache: Whether to use local storage in `data/structures/`.
        
    Returns:
        Tuple of (pdb_text, source_description, is_fallback).
    """
    if not isinstance(pdb_id, str) or not re.fullmatch(r"[A-Za-z0-9]{4}", pdb_id.strip()):
        raise ValueError("PDB ID must be a four-character alphanumeric identifier.")
    if (
        isinstance(timeout_sec, bool)
        or not isinstance(timeout_sec, (int, float))
        or not math.isfinite(timeout_sec)
        or timeout_sec <= 0
    ):
        raise ValueError("Structure request timeout must be a positive finite number.")

    clean_id = pdb_id.strip().upper()
    local_path = os.path.join(STRUCTURE_DIR, f"{clean_id.lower()}.pdb")
    retrieval_warning = None
    
    # 1. Check local cache
    if use_cache and os.path.exists(local_path):
        try:
            with open(local_path, "r", encoding="utf-8", errors="replace") as f:
                content = f.read()
            if "ATOM" in content or "HEADER" in content:
                return content, f"Local Cache ({os.path.basename(local_path)})", False
        except OSError as exc:
            retrieval_warning = f"Local structure cache could not be read: {exc}"

    # 2. Live fetch from RCSB PDB
    url = f"https://files.rcsb.org/download/{clean_id}.pdb"
    try:
        resp = requests.get(
            url,
            timeout=timeout_sec,
            headers={"User-Agent": "AI-PETase-Optimizer/1.0 (Research App)"}
        )
        if resp.status_code == 200 and ("ATOM" in resp.text or "HEADER" in resp.text):
            source = f"RCSB Protein Data Bank ({clean_id})"
            if use_cache:
                try:
                    os.makedirs(STRUCTURE_DIR, exist_ok=True)
                    with open(local_path, "w", encoding="utf-8") as f:
                        f.write(resp.text)
                except OSError as exc:
                    source += f" (cache write failed: {exc})"
            return resp.text, source, False
        retrieval_warning = f"RCSB returned HTTP {resp.status_code} or an invalid PDB response."
    except requests.RequestException as exc:
        retrieval_warning = f"RCSB request failed: {exc}"

    # 3. Verified local fallback for default 5XJH
    if clean_id == DEFAULT_PDB_ID and os.path.exists(DEFAULT_FALLBACK_FILE):
        try:
            with open(DEFAULT_FALLBACK_FILE, "r", encoding="utf-8", errors="replace") as f:
                content = f.read()
            source = "Verified Local Fallback (5xjh.pdb)"
            if retrieval_warning:
                source += f"; {retrieval_warning}"
            return content, source, True
        except OSError as exc:
            retrieval_warning = f"Bundled 5XJH fallback could not be read: {exc}"

    detail = f"Structure {clean_id} unavailable."
    if retrieval_warning:
        detail += f" {retrieval_warning}"
    if clean_id == DEFAULT_PDB_ID:
        detail += " No readable local fallback was found."
    return None, detail, False


def generate_3dmol_html(
    pdb_text: str,
    style: str = "cartoon",
    color_scheme: str = "spectrum",
    highlight_features: Optional[List[str]] = None,
    custom_residue: Optional[int] = None,
    custom_label: Optional[str] = None,
    spin: bool = False,
    height: int = 520
) -> str:
    """
    Generate interactive 3Dmol.js WebGL viewer HTML.
    
    Args:
        pdb_text: Raw PDB file text.
        style: 'cartoon', 'stick', 'sphere', or 'surface'.
        color_scheme: 'spectrum', 'chain', or 'ss'.
        highlight_features: List of landmarks to highlight ('catalytic_triad', 'disulfide_bonds', etc.).
        custom_residue: Optional residue number to highlight with a sphere/stick marker.
        custom_label: Label for the custom residue marker.
        spin: Whether to enable slow rotation.
        height: Height in pixels for the container.
        
    Returns:
        HTML string suitable for an iframe-based Streamlit display.
    """
    if not pdb_text:
        return "<div style='color:#e63946; padding:20px;'>Structure coordinates not available.</div>"

    escaped_pdb = json.dumps(pdb_text)
    highlights = highlight_features or []

    # Style configuration
    if style == "stick":
        style_js = "viewer.setStyle({}, {stick: {colorscheme: 'amino'}});"
    elif style == "sphere":
        style_js = "viewer.setStyle({}, {sphere: {colorscheme: 'amino', scale: 0.8}});"
    elif style == "surface":
        style_js = """
            viewer.setStyle({}, {cartoon: {color: 'spectrum'}});
            viewer.addSurface($3Dmol.SurfaceType.VDW, {opacity: 0.60, color: '#168575'});
        """
    else:  # cartoon
        if color_scheme == "chain":
            style_js = "viewer.setStyle({}, {cartoon: {colorscheme: 'chain'}});"
        elif color_scheme == "ss":
            style_js = "viewer.setStyle({}, {cartoon: {colorscheme: 'ssJmol'}});"
        else:
            style_js = "viewer.setStyle({}, {cartoon: {color: 'spectrum'}});"

    # Feature annotations
    feature_script = ""
    
    # Catalytic Triad: Ser160, Asp206, His237
    if "catalytic_triad" in highlights:
        feature_script += """
        viewer.addStyle({resi: [160, 206, 237]}, {stick: {color: '#a94e3b', radius: 0.30}, sphere: {color: '#a94e3b', radius: 0.65, opacity: 0.7}});
        viewer.addLabel("Ser160 (Nu)", {resi: 160}, {backgroundColor: '#a94e3b', backgroundOpacity: 0.92, fontColor: '#FFFFFF', fontSize: 11, showBackground: true});
        viewer.addLabel("Asp206 (Acid)", {resi: 206}, {backgroundColor: '#a94e3b', backgroundOpacity: 0.92, fontColor: '#FFFFFF', fontSize: 11, showBackground: true});
        viewer.addLabel("His237 (Base)", {resi: 237}, {backgroundColor: '#a94e3b', backgroundOpacity: 0.92, fontColor: '#FFFFFF', fontSize: 11, showBackground: true});
        """

    # Disulfide Bridges: DS1 (203-239), DS2 (273-289)
    if "disulfide_bonds" in highlights:
        feature_script += """
        viewer.addStyle({resi: [203, 239, 273, 289]}, {stick: {color: '#F4A261', radius: 0.32}});
        viewer.addLabel("DS1 (203-239)", {resi: 203}, {backgroundColor: '#F4A261', backgroundOpacity: 0.85, fontColor: '#000000', fontSize: 10, showBackground: true});
        viewer.addLabel("DS2 (273-289)", {resi: 273}, {backgroundColor: '#F4A261', backgroundOpacity: 0.85, fontColor: '#000000', fontSize: 10, showBackground: true});
        """

    # Oxyanion Hole: Tyr87, Met161
    if "oxyanion_hole" in highlights:
        feature_script += """
        viewer.addStyle({resi: [87, 161]}, {stick: {color: '#48CAE4', radius: 0.28}});
        viewer.addLabel("Tyr87", {resi: 87}, {backgroundColor: '#0077B6', backgroundOpacity: 0.85, fontColor: '#FFFFFF', fontSize: 10, showBackground: true});
        viewer.addLabel("Met161", {resi: 161}, {backgroundColor: '#0077B6', backgroundOpacity: 0.85, fontColor: '#FFFFFF', fontSize: 10, showBackground: true});
        """

    # Substrate Cleft: Trp185, Trp159
    if "substrate_cleft" in highlights:
        feature_script += """
        viewer.addStyle({resi: [185, 159]}, {stick: {color: '#52B788', radius: 0.28}});
        viewer.addLabel("Trp185 (Wobble)", {resi: 185}, {backgroundColor: '#2D6A4F', backgroundOpacity: 0.85, fontColor: '#FFFFFF', fontSize: 10, showBackground: true});
        viewer.addLabel("Trp159", {resi: 159}, {backgroundColor: '#2D6A4F', backgroundOpacity: 0.85, fontColor: '#FFFFFF', fontSize: 10, showBackground: true});
        """

    # Optional custom residue highlight (e.g. selected mutation position)
    if custom_residue:
        label_text = json.dumps(custom_label or f"Residue {custom_residue}")
        feature_script += f"""
        viewer.addStyle({{resi: {custom_residue}}}, {{stick: {{color: '#9D4EDD', radius: 0.35}}, sphere: {{color: '#C77DFF', radius: 0.85, opacity: 0.8}}}});
        viewer.addLabel({label_text}, {{resi: {custom_residue}}}, {{backgroundColor: '#7B2CBF', backgroundOpacity: 0.9, fontColor: '#FFFFFF', fontSize: 11, showBackground: true}});
        """

    spin_js = "viewer.spin('y', 0.8);" if spin else "viewer.spin(false);"

    html = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
      <meta charset="utf-8">
      <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
      <style>
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body, html {{ width: 100%; height: 100%; overflow: hidden; background-color: #f4f5f0; font-family: 'DM Sans', sans-serif; }}
        #container {{ width: 100%; height: {height}px; position: relative; border-radius: 14px; border: 1px solid #d9e2dd; }}
        .badge {{
          position: absolute; bottom: 8px; left: 12px;
          background: rgba(255, 255, 255, 0.94); color: #405851;
          font-size: 11px; padding: 6px 10px; border-radius: 7px;
          border: 1px solid #d9e2dd; pointer-events: none;
        }}
        .controls-note {{
          position: absolute; top: 8px; right: 12px;
          background: rgba(255, 255, 255, 0.94); color: #075e55;
          font-size: 11px; padding: 6px 10px; border-radius: 7px;
          border: 1px solid #d9e2dd; pointer-events: none;
        }}
      </style>
    </head>
    <body>
      <div id="container"></div>
      <div class="controls-note">Rotate: Drag | Zoom: Scroll | Pan: Right-Click</div>
      <div class="badge">Experimental Crystal Structure (PDB: 5XJH, 1.58 Å) • Han et al., Nat Commun 2017</div>
      <script>
        var element = document.getElementById('container');
      var config = {{ backgroundColor: '#f4f5f0' }};
        var viewer = $3Dmol.createViewer(element, config);
        var pdbData = {escaped_pdb};
        viewer.addModel(pdbData, "pdb");
        {style_js}
        {feature_script}
        viewer.zoomTo();
        viewer.render();
        {spin_js}
      </script>
    </body>
    </html>
    """
    return html
