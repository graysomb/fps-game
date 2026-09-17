import importlib.util
from pathlib import Path
import unittest

spec = importlib.util.spec_from_file_location("adaptive_gate", Path(__file__).resolve().parents[1]/"tools/benchmark_adaptive.py")
gate_module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gate_module)

class GateTests(unittest.TestCase):
    def pairs(self, gain=-1):
        result = []
        for i in range(10):
            fine = dict(valid=True, coarsened=True, physics_ms=10+i*.01, frame_ms=12+i*.01,
                        frame_p95_ms=15+i*.01, event_p95_ms=30+i*.01, syncs=60)
            adaptive = dict(fine)
            for metric in ("physics_ms", "frame_ms", "frame_p95_ms", "event_p95_ms"):
                adaptive[metric] += gain
            result.append({"off": fine, "on": adaptive})
        return result

    def test_repeatable_improvement(self):
        self.assertTrue(gate_module.gate(self.pairs())["eligible"])

    def test_insufficient_and_invalid_trials_fail_closed(self):
        pairs = self.pairs()
        self.assertFalse(gate_module.gate(pairs[:9])["eligible"])
        pairs[3]["on"]["valid"] = False
        self.assertFalse(gate_module.gate(pairs)["eligible"])

    def test_noise_is_not_a_speedup(self):
        pairs = self.pairs(0)
        for i, pair in enumerate(pairs):
            pair["on"]["physics_ms"] += .5 if i % 2 else -.5
        self.assertFalse(gate_module.gate(pairs)["eligible"])

    def test_fine_fallback_is_not_an_adaptive_speedup(self):
        pairs = self.pairs()
        pairs[4]["on"]["coarsened"] = False
        result = gate_module.gate(pairs)
        self.assertFalse(result["eligible"])
        self.assertEqual(result["status"], "INELIGIBLE_NO_COARSENING")

    def test_sync_or_event_regression_vetoes_physics_gain(self):
        pairs = self.pairs()
        pairs[1]["on"]["syncs"] += 1
        self.assertFalse(gate_module.gate(pairs)["eligible"])
        pairs = self.pairs()
        for pair in pairs:
            pair["on"]["event_p95_ms"] += 5
        self.assertFalse(gate_module.gate(pairs)["eligible"])

    def test_event_union_uses_same_frames(self):
        pair = {}
        for mode in ("off", "on"):
            pair[mode] = {"rows": [dict(step=i, repacked=i==1, generation=int(mode=="on" and i>=3), frameMs=i) for i in range(1, 5)]}
        gate_module.pair_events(pair)
        self.assertEqual(pair["off"]["event_p95_ms"], pair["on"]["event_p95_ms"])
        self.assertGreater(pair["on"]["event_p95_ms"], 1)

if __name__ == "__main__":
    unittest.main()
