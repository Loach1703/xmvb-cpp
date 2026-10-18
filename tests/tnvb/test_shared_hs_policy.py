#!/usr/bin/env python3
"""Pure shared-bank safety/routing tests; no Block2 import or computation."""
from pathlib import Path
import importlib.util
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("hs_shared", ROOT / "src/tnvb/worker/hs_shared.py")
SHARED = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(SHARED)


class SharedPolicyTest(unittest.TestCase):
    def test_energy_only_parallel_routing(self):
        self.assertTrue(SHARED.use_shared_hs(None, 2))
        self.assertTrue(SHARED.use_shared_hs(None, 16))
        self.assertFalse(SHARED.use_shared_hs(None, 1))
        for root in (0, 1):
            self.assertFalse(SHARED.use_shared_hs(root, 16))

    def test_file_selection_matches_original(self):
        for tags in (["S1", "S12"], ["S1.P2", "S2"]):
            names = ["F.MPS.S1.0", "F.MPS.S12.0", "F.MPS.S11.0",
                     "S1-mps_info.bin", "S12-mps_info.bin", "S1", "S12.",
                     "F.MPS.S1.P2.0", "S1.P2-mps_info.bin", "F.MPS.S2.0"]
            expected = [name for name in names if any(
                name == f"{tag}-mps_info.bin" or f".{tag}." in name for tag in tags)]
            self.assertEqual(SHARED.selected_names(names, tags), expected)

    def test_bank_integrity_detects_changes(self):
        with tempfile.TemporaryDirectory() as directory:
            file = Path(directory) / "tensor"
            file.write_bytes(b"original")
            before = SHARED.fingerprints(directory)
            self.assertEqual(before, SHARED.fingerprints(directory))
            file.write_bytes(b"modified")
            self.assertNotEqual(before, SHARED.fingerprints(directory))

    def test_invalid_manifests_fail_before_block2(self):
        good = dict(norb=14, nelec=14, bond=200, hs_workers=4, hs_threads=2,
                    symmetry="SU2", state_tags=["S1", "S2"])
        SHARED.validate_manifest(good)
        for change in (dict(hs_threads=0), dict(hs_workers=0), dict(symmetry="SZ"),
                       dict(state_tags=[]), dict(state_tags=["S1", "S1"]),
                       dict(state_tags=["../S1"])):
            with self.assertRaises(ValueError):
                SHARED.validate_manifest(dict(good, **change))

    def test_cutoff_default_and_cache_are_consistent(self):
        source = (ROOT / "src/tnvb/worker/pair_core.py").read_text()
        self.assertIn('os.environ.get("XMVB_BLOCK2_SU2_CUTOFF", "0")', source)
        self.assertIn("('XMVB_BLOCK2_SU2_CUTOFF', '0')", source)


if __name__ == "__main__":
    unittest.main()
