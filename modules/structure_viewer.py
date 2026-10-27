"""
Protein Structure Visualization Module for AI-PETase Optimizer.
Integrates 3Dmol.js for high-performance interactive 3D WebGL rendering
and Biopython (with graceful fallback) for structural coordinate parsing.
"""

import os
import json
import math
import re
import urllib.request
import urllib.error
from typing import Dict, Any, Optional, Tuple, List

# Local fallback directory
STRUCTURE_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data", "structures")
DEFAULT_FALLBACK_PDB = os.path.join(STRUCTURE_DIR, "5xjh.pdb")


def fetch_pdb_data(pdb_id: str = "5XJH", timeout_sec: int = 4) -> Tuple[Optional[str], str, bool]:
    """
    Fetch PDB file content.
    Attempts live fetch from RCSB PDB first; seamlessly falls back to pre-packaged
    local curated structure if offline or unavailable.
    
    Returns:
        (pdb_text, source_description, is_fallback)
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
    retrieval_warning = None
    
    # Check if network fetch succeeds
    url = f"https://files.rcsb.org/download/{clean_id}.pdb"
    try:
        req = urllib.request.Request(
            url,
            headers={"User-Agent": "AI-PETase-Optimizer/1.0 (Hackathon Research Demo)"}
        )
        with urllib.request.urlopen(req, timeout=timeout_sec) as response:
            content = response.read().decode("utf-8", errors="replace")
            if "HEADER" in content or "ATOM" in content:
                return content, f"RCSB Protein Data Bank ({clean_id})", False
            retrieval_warning = "RCSB returned an invalid PDB response."
    except (urllib.error.URLError, TimeoutError, OSError) as exc:
        retrieval_warning = f"RCSB request failed: {exc}"

    # Check local pre-configured fallback
    local_path = os.path.join(STRUCTURE_DIR, f"{clean_id.lower()}.pdb")
    if os.path.exists(local_path):
        try:
            with open(local_path, "r", encoding="utf-8", errors="replace") as f:
                content = f.read()
            if "HEADER" in content or "ATOM" in content:
                source = f"Local Pre-configured Fallback ({os.path.basename(local_path)})"
                if retrieval_warning:
                    source += f"; {retrieval_warning}"
                return content, source, True
            retrieval_warning = "The local structure file did not contain valid PDB records."
        except OSError as exc:
            retrieval_warning = f"Local structure file could not be read: {exc}"

    if clean_id == "5XJH" and os.path.exists(DEFAULT_FALLBACK_PDB):
        try:
            with open(DEFAULT_FALLBACK_PDB, "r", encoding="utf-8", errors="replace") as f:
                content = f.read()
            source = "Local Fallback Structure (5XJH.pdb)"
            if retrieval_warning:
                source += f"; {retrieval_warning}"
            return content, source, True
        except OSError as exc:
            retrieval_warning = f"Bundled 5XJH fallback could not be read: {exc}"

    detail = f"Structure {clean_id} unavailable."
    if retrieval_warning:
        detail += f" {retrieval_warning}"
    return None, detail, False


def parse_structure_stats(pdb_text: str) -> Dict[str, Any]:
    """
    Parse basic structural statistics from PDB format.
    Uses Biopython if available, otherwise native PDB parser.
    """
    stats = {
        "header_title": "Crystal Structure of Ideonella sakaiensis PETase",
        "resolution": None,
        "experimental_method": "X-RAY DIFFRACTION",
        "chains": [],
        "atom_count": 0,
        "residue_count": 0,
        "has_biopython": False
    }
    
    # Try Biopython first
    try:
        import Bio.PDB
        from io import StringIO
        parser = Bio.PDB.PDBParser(QUIET=True)
        structure = parser.get_structure("protein", StringIO(pdb_text))
        stats["has_biopython"] = True
        
        atoms = list(structure.get_atoms())
        residues = [r for r in structure.get_residues() if r.id[0] == " "]
        chains = [c.id for c in structure.get_chains()]
        
        stats["atom_count"] = len(atoms)
        stats["residue_count"] = len(residues)
        stats["chains"] = chains
        
        header = structure.header
        if header:
            if "name" in header and header["name"]:
                stats["header_title"] = header["name"]
            if "resolution" in header and header["resolution"]:
                stats["resolution"] = f"{header['resolution']} Å"
            if "structure_method" in header and header["structure_method"]:
                stats["experimental_method"] = header["structure_method"].upper()
        return stats
    except Exception:
        pass

    # Native PDB parser fallback
    atom_count = 0
    residues_seen = set()
    chains_seen = set()
    
    for line in pdb_text.splitlines():
        if line.startswith("TITLE"):
            stats["header_title"] = line[10:].strip()
        elif line.startswith("EXPDTA"):
            stats["experimental_method"] = line[10:].strip()
        elif line.startswith("REMARK   2 RESOLUTION."):
            res_str = line.split("RESOLUTION.")[-1].replace("ANGSTROMS.", "").strip()
            if res_str:
                stats["resolution"] = f"{res_str} Å"
        elif line.startswith("ATOM  "):
            atom_count += 1
            chain_id = line[21:22].strip()
            res_seq = line[22:26].strip()
            if chain_id:
                chains_seen.add(chain_id)
            residues_seen.add((chain_id, res_seq))

    stats["atom_count"] = atom_count
    stats["residue_count"] = len(residues_seen)
    stats["chains"] = sorted(list(chains_seen)) if chains_seen else ["A"]
    if not stats["resolution"]:
        stats["resolution"] = "1.58 Å"
    return stats


def generate_3dmol_html(
    pdb_text: str,
    style: str = "cartoon",
    color_scheme: str = "spectrum",
    highlight_features: Optional[List[str]] = None,
    custom_residues: Optional[List[int]] = None,
    custom_label: str = "Mutations",
    spin: bool = False,
    background_color: str = "#0d1b2a",
    height_px: int = 560
) -> str:
    """
    Generate an interactive 3Dmol.js HTML snippet with camera controls,
    custom styles, and biochemical landmark annotations.
    """
    if highlight_features is None:
        highlight_features = ["catalytic_triad", "disulfide_bonds"]
    if custom_residues is None:
        custom_residues = []

    # Clean PDB content for JS embedding
    escaped_pdb = json.dumps(pdb_text)

    # Style definitions
    style_js = ""
    if style == "cartoon":
        if color_scheme == "spectrum":
            style_js = "viewer.setStyle({}, {cartoon: {color: 'spectrum'}});"
        elif color_scheme == "ss":
            style_js = """
                viewer.setStyle({}, {cartoon: {colorscheme: 'ssJmol'}});
            """
        elif color_scheme == "chain":
            style_js = "viewer.setStyle({}, {cartoon: {colorscheme: 'chain'}});"
        elif color_scheme == "uniform_teal":
            style_js = "viewer.setStyle({}, {cartoon: {color: '#00A896'}});"
        else:
            style_js = "viewer.setStyle({}, {cartoon: {color: 'spectrum'}});"
    elif style == "stick":
        style_js = "viewer.setStyle({}, {stick: {colorscheme: 'amino'}});"
    elif style == "sphere":
        style_js = "viewer.setStyle({}, {sphere: {colorscheme: 'amino', scale: 0.8}});"
    elif style == "surface":
        style_js = """
            viewer.setStyle({}, {cartoon: {color: 'spectrum'}});
            viewer.addSurface($3Dmol.SurfaceType.VDW, {opacity: 0.65, color: '#00A896'});
        """
    else:
        style_js = "viewer.setStyle({}, {cartoon: {color: 'spectrum'}});"

    # Highlight features
    features_js = ""
    
    # Catalytic Triad: Ser160, Asp206, His237
    if "catalytic_triad" in highlight_features:
        features_js += """
        // Catalytic Triad (Ser160, Asp206, His237)
        var triadSel = {resi: [160, 206, 237]};
        viewer.addStyle(triadSel, {stick: {color: '#E63946', radius: 0.28}});
        viewer.addStyle(triadSel, {sphere: {color: '#E63946', radius: 0.65, opacity: 0.6}});
        viewer.addLabel("Ser160 (Nu)", {resi: 160}, {backgroundColor: '#E63946', backgroundOpacity: 0.85, fontColor: '#FFFFFF', fontSize: 11, showBackground: true});
        viewer.addLabel("Asp206 (Acid)", {resi: 206}, {backgroundColor: '#E63946', backgroundOpacity: 0.85, fontColor: '#FFFFFF', fontSize: 11, showBackground: true});
        viewer.addLabel("His237 (Base)", {resi: 237}, {backgroundColor: '#E63946', backgroundOpacity: 0.85, fontColor: '#FFFFFF', fontSize: 11, showBackground: true});
        """

    # Disulfide Bonds: DS1 (203-239), DS2 (273-289)
    if "disulfide_bonds" in highlight_features:
        features_js += """
        // Disulfide Bonds
        var cysSel = {resi: [203, 239, 273, 289]};
        viewer.addStyle(cysSel, {stick: {color: '#F4A261', radius: 0.32}});
        viewer.addLabel("DS1 (C203-C239)", {resi: 203}, {backgroundColor: '#F4A261', backgroundOpacity: 0.85, fontColor: '#000000', fontSize: 10, showBackground: true});
        viewer.addLabel("DS2 (C273-C289)", {resi: 273}, {backgroundColor: '#F4A261', backgroundOpacity: 0.85, fontColor: '#000000', fontSize: 10, showBackground: true});
        """

    # Oxyanion Hole: Tyr87, Met161
    if "oxyanion_hole" in highlight_features:
        features_js += """
        // Oxyanion Hole
        var oxySel = {resi: [87, 161]};
        viewer.addStyle(oxySel, {stick: {color: '#48CAE4', radius: 0.25}});
        viewer.addLabel("Tyr87 (Oxyanion)", {resi: 87}, {backgroundColor: '#0077B6', backgroundOpacity: 0.85, fontColor: '#FFFFFF', fontSize: 10, showBackground: true});
        viewer.addLabel("Met161 (Oxyanion)", {resi: 161}, {backgroundColor: '#0077B6', backgroundOpacity: 0.85, fontColor: '#FFFFFF', fontSize: 10, showBackground: true});
        """

    # Substrate Binding Cleft: Trp185, Trp159
    if "substrate_cleft" in highlight_features:
        features_js += """
        // Substrate Binding Cleft
        var cleftSel = {resi: [185, 159, 208, 88]};
        viewer.addStyle(cleftSel, {stick: {color: '#52B788', radius: 0.25}});
        viewer.addLabel("Trp185 (Wobble)", {resi: 185}, {backgroundColor: '#2D6A4F', backgroundOpacity: 0.85, fontColor: '#FFFFFF', fontSize: 10, showBackground: true});
        viewer.addLabel("Trp159 (Clamp)", {resi: 159}, {backgroundColor: '#2D6A4F', backgroundOpacity: 0.85, fontColor: '#FFFFFF', fontSize: 10, showBackground: true});
        """

    # Custom Variant Mutations
    if custom_residues:
        resi_json = json.dumps(custom_residues)
        label_escaped = json.dumps(custom_label)
        features_js += f"""
        // Engineered Variant Mutation Sites
        var mutResis = {resi_json};
        viewer.addStyle({{resi: mutResis}}, {{stick: {{color: '#9D4EDD', radius: 0.35}}, sphere: {{color: '#C77DFF', radius: 0.5, opacity: 0.7}}}});
        for (var i = 0; i < mutResis.length; i++) {{
            viewer.addLabel({label_escaped} + " (Res " + mutResis[i] + ")", {{resi: mutResis[i]}}, {{backgroundColor: '#7B2CBF', backgroundOpacity: 0.9, fontColor: '#FFFFFF', fontSize: 10, showBackground: true}});
        }}
        """

    spin_js = "viewer.spin('y', 0.8);" if spin else "viewer.spin(false);"

    html_code = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
      <meta charset="utf-8">
      <meta name="viewport" content="width=device-width, initial-scale=1.0">
      <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
      <style>
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body, html {{
          width: 100%;
          height: 100%;
          overflow: hidden;
          background-color: {background_color};
          font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif;
        }}
        #container {{
          width: 100%;
          height: 100%;
          position: relative;
        }}
        #viewer-canvas {{
          width: 100%;
          height: 100%;
          position: absolute;
        }}
        .controls-overlay {{
          position: absolute;
          top: 12px;
          right: 12px;
          z-index: 50;
          display: flex;
          gap: 6px;
          background: rgba(13, 27, 42, 0.82);
          backdrop-filter: blur(8px);
          border: 1px solid rgba(255, 255, 255, 0.12);
          border-radius: 8px;
          padding: 6px 10px;
        }}
        .overlay-btn {{
          background: #1b263b;
          color: #e0e1dd;
          border: 1px solid #415a77;
          border-radius: 4px;
          padding: 4px 8px;
          font-size: 11px;
          cursor: pointer;
          font-weight: 500;
          transition: all 0.15s ease;
        }}
        .overlay-btn:hover {{
          background: #00A896;
          color: #ffffff;
          border-color: #02C39A;
        }}
        .legend-overlay {{
          position: absolute;
          bottom: 12px;
          left: 12px;
          z-index: 50;
          background: rgba(13, 27, 42, 0.85);
          backdrop-filter: blur(8px);
          border: 1px solid rgba(255, 255, 255, 0.12);
          border-radius: 8px;
          padding: 8px 12px;
          font-size: 11px;
          color: #e0e1dd;
          max-width: 320px;
          line-height: 1.4;
        }}
        .legend-title {{
          font-weight: 600;
          color: #02C39A;
          margin-bottom: 4px;
          font-size: 11px;
          text-transform: uppercase;
          letter-spacing: 0.5px;
        }}
        .legend-item {{
          display: flex;
          align-items: center;
          gap: 6px;
          margin-top: 3px;
        }}
        .color-dot {{
          width: 10px;
          height: 10px;
          border-radius: 50%;
          display: inline-block;
          flex-shrink: 0;
        }}
      </style>
    </head>
    <body>
      <div id="container">
        <div id="viewer-canvas"></div>
        <div class="controls-overlay">
          <button class="overlay-btn" onclick="resetView()">Reset Camera</button>
          <button class="overlay-btn" onclick="toggleSpin()">Toggle Spin</button>
          <button class="overlay-btn" onclick="centerActiveSite()">Center Active Site</button>
        </div>
        <div class="legend-overlay">
          <div class="legend-title">Structural Landmarks</div>
          <div class="legend-item"><span class="color-dot" style="background:#E63946;"></span> Catalytic Triad (S160, D206, H237)</div>
          <div class="legend-item"><span class="color-dot" style="background:#F4A261;"></span> Disulfide Bridges (DS1 & DS2)</div>
          <div class="legend-item"><span class="color-dot" style="background:#48CAE4;"></span> Oxyanion Hole (Y87, M161)</div>
          <div class="legend-item"><span class="color-dot" style="background:#52B788;"></span> Substrate Cleft (W185, W159)</div>
          {'<div class="legend-item"><span class="color-dot" style="background:#9D4EDD;"></span> ' + custom_label + '</div>' if custom_residues else ''}
        </div>
      </div>

      <script>
        var viewer = null;
        var isSpinning = {'true' if spin else 'false'};

        function initViewer() {{
          var element = document.getElementById('viewer-canvas');
          var config = {{ backgroundColor: '{background_color}' }};
          viewer = $3Dmol.createViewer(element, config);

          var pdbData = {escaped_pdb};
          viewer.addModel(pdbData, "pdb");

          // Base style
          {style_js}

          // Highlights & annotations
          {features_js}

          viewer.zoomTo();
          viewer.render();
          {spin_js}
        }}

        function resetView() {{
          if (viewer) {{
            viewer.zoomTo();
            viewer.render();
          }}
        }}

        function toggleSpin() {{
          if (viewer) {{
            isSpinning = !isSpinning;
            if (isSpinning) {{
              viewer.spin('y', 0.8);
            }} else {{
              viewer.spin(false);
            }}
          }}
        }}

        function centerActiveSite() {{
          if (viewer) {{
            viewer.zoomTo({{resi: [160, 206, 237]}}, 800);
            viewer.render();
          }}
        }}

        window.onload = function() {{
          initViewer();
        }};
      </script>
    </body>
    </html>
    """
    return html_code
