"""
Comprehensive automated test suite for AI-PETase Optimizer.
Verifies structure fallback, variant dataset integrity, 3Dmol.js generation,
ESM scoring module guardrails, and eco scenario mass balances.
"""

import os
import sys
import json
import unittest

# Ensure root directory is on sys.path
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

from data.reference_data import (
    STRUCTURAL_FEATURES,
    METRIC_EXPLANATIONS,
    PET_PHYSICAL_PROPERTIES,
    LCA_BENCHMARKS,
    KEY_CITATIONS,
    ISPETASE_MATURE_SEQUENCE
)
from modules.structure_viewer import (
    fetch_pdb_data,
    parse_structure_stats,
    generate_3dmol_html
)
from modules.variant_analysis import (
    load_variants_data,
    get_variants_dataframe,
    create_tm_comparison_chart
)
from modules.esm_scoring import (
    check_esm_availability,
    parse_mutation_code,
    analyze_mutation_site,
    score_variant_mutation
)
from modules.impact_calculator import (
    calculate_scenario
)


class TestAIPETaseOptimizer(unittest.TestCase):
    
    def test_reference_data_integrity(self):
        """Verify presence of core biological and chemical constants."""
        self.assertGreater(len(ISPETASE_MATURE_SEQUENCE), 250)
        self.assertIn("catalytic_triad", STRUCTURAL_FEATURES)
        triad_res = [r["resseq"] for r in STRUCTURAL_FEATURES["catalytic_triad"]["residues"]]
        self.assertEqual(triad_res, [160, 206, 237])
        self.assertIn("disulfide_bonds", STRUCTURAL_FEATURES)
        self.assertGreaterEqual(len(KEY_CITATIONS), 5)

    def test_curated_variants_data(self):
        """Verify peer-reviewed variants JSON contains all required benchmark strains."""
        variants = load_variants_data()
        self.assertIn("wt_petase", variants)
        self.assertIn("fast_petase", variants)
        self.assertIn("durapetase", variants)
        self.assertIn("thermopetase", variants)
        self.assertIn("hotpetase", variants)
        self.assertIn("lcc_iccg", variants)

        # Check WT reference flags
        self.assertTrue(variants["wt_petase"]["is_wt_reference"])
        self.assertFalse(variants["fast_petase"]["is_wt_reference"])

        # Check FAST-PETase reported parameters
        fast = variants["fast_petase"]
        self.assertEqual(fast["tm_celsius"], 67.8)
        self.assertEqual(fast["tm_delta"], 19.0)
        self.assertEqual(len(fast["mutations_list"]), 5)
        self.assertIn("N233K", fast["mutations_list"])

        # Check DataFrame construction
        df = get_variants_dataframe()
        self.assertEqual(len(df), len(variants))
        self.assertIn("Name", df.columns)
        self.assertIn("Tm (°C)", df.columns)
        self.assertEqual(
            df.loc[df["Name"] == "LCC-ICCG (Cutinase Reference)", "ΔTm Reference"].iloc[0],
            "LCC parent enzyme"
        )

        chart = create_tm_comparison_chart(variants)
        lcc_index = list(chart.data[0].x).index("LCC-ICCG (Cutinase Reference)")
        self.assertIn("vs LCC parent enzyme", chart.data[0].text[lcc_index])

    def test_pdb_fallback_and_parsing(self):
        """Verify that 5XJH fallback PDB file exists and parses valid coordinates."""
        pdb_content, source_desc, is_fallback = fetch_pdb_data("5XJH", timeout_sec=1)
        self.assertIsNotNone(pdb_content)
        self.assertIn("ATOM", pdb_content)

        stats = parse_structure_stats(pdb_content)
        self.assertGreater(stats["atom_count"], 1500)
        self.assertGreater(stats["residue_count"], 200)
        self.assertIn("A", stats["chains"])

    def test_3dmol_html_generation(self):
        """Verify 3Dmol.js HTML snippet generation with landmarks and styles."""
        pdb_content, _, _ = fetch_pdb_data("5XJH", timeout_sec=1)
        html = generate_3dmol_html(
            pdb_text=pdb_content,
            style="cartoon",
            color_scheme="spectrum",
            highlight_features=["catalytic_triad", "disulfide_bonds"],
            custom_residues=[121, 186, 224, 233, 280],
            custom_label="FAST-PETase",
            spin=False
        )
        self.assertIn("3Dmol-min.js", html)
        self.assertIn("Ser160 (Nu)", html)
        self.assertIn("Asp206 (Acid)", html)
        self.assertIn("His237 (Base)", html)
        self.assertIn("DS1 (C203-C239)", html)
        self.assertIn("FAST-PETase", html)

    def test_esm_scoring_module_guardrails(self):
        """Verify mutation parser, site criticality analyzer, and transparent scoring disclaimers."""
        parsed = parse_mutation_code("N233K")
        self.assertEqual(parsed, ("N", 233, "K"))
        self.assertIsNone(parse_mutation_code("INVALID_INPUT"))

        # Test catalytic triad critical alert
        cat_site = analyze_mutation_site(160)
        self.assertTrue(cat_site["is_critical"])
        self.assertIn("Catalytic Triad", cat_site["region"])

        # Test surface site
        surf_site = analyze_mutation_site(233)
        self.assertFalse(surf_site["is_critical"])

        # Test scoring execution & transparent disclaimer
        result = score_variant_mutation("N233K")
        self.assertTrue(result["success"])
        self.assertIsInstance(result["score"], float)
        self.assertEqual(result["score"], result["score_difference"])
        self.assertIn("esm2_t6_8M_UR50D", result["scoring_engine"])
        self.assertIn("disclaimer", result)
        self.assertIn("No biological classification", result["verdict"])

        # ESM-2 scores are not converted into unsupported functional verdicts.
        cat_result = score_variant_mutation("S160A")
        self.assertTrue(cat_result["success"])
        self.assertIn("No biological classification", cat_result["verdict"])
        self.assertNotIn("Deleterious", cat_result["verdict"])

        mismatch = score_variant_mutation("A233K")
        self.assertFalse(mismatch["success"])
        self.assertIn("Wild-type residue mismatch", mismatch["error"])

    def test_impact_calculator_mass_balance(self):
        """Verify material balance conservation and emissions calculations."""
        # 100,000 kg waste, 15% baseline, +25% improvement, 30 days
        res = calculate_scenario(
            pet_waste_kg=100000.0,
            baseline_rate_pct=15.0,
            improvement_pct=25.0,
            processing_days=30,
            enable_emissions=True,
            virgin_factor_kg_co2=2.50,
            recycling_factor_kg_co2=1.25
        )
        # Baseline = 15,000 kg; Hypothetical = 40,000 kg; Additional = 25,000 kg
        self.assertEqual(res["baseline_recovered_kg"], 15000.0)
        self.assertEqual(res["hypothetical_recovered_kg"], 40000.0)
        self.assertEqual(res["additional_recovered_kg"], 25000.0)
        self.assertEqual(res["hypothetical_rate_pct"], 40.0)
        
        # Mass conservation
        self.assertEqual(res["hypothetical_recovered_kg"] + res["hypothetical_unrecovered_kg"], 100000.0)
        
        # Throughput
        self.assertAlmostEqual(res["daily_recovered_hypothetical_kg"], 40000.0 / 30.0, places=2)

        # Emissions: Net factor = 2.50 - 1.25 = 1.25 kg CO2e / kg
        self.assertIsNotNone(res["emissions"])
        expected_avoided_delta = 25000.0 * 1.25 # 31,250 kg CO2e
        self.assertEqual(res["emissions"]["net_avoided_difference_kg_co2"], expected_avoided_delta)
        self.assertIn("miles_offset_equivalent", res["emissions"])


if __name__ == "__main__":
    unittest.main()
