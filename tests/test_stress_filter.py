import unittest

import numpy as np
import pandas as pd

from lincs_microglial.stress_filter import condition_consistency, cosine_similarity_to_centroid


class StressFilterTest(unittest.TestCase):
    def test_condition_consistency_distinguishes_replication(self) -> None:
        targets = pd.DataFrame(
            {"gene_name": ["G1", "G2"], "protective_expression_direction": [1, 1]}
        )
        rows = []
        signatures = [
            ("A1", "A", 1, [2.0, 0.0]),
            ("B1", "B", 1, [2.0, 0.0]),
            ("B2", "B", 10, [2.0, 0.0]),
            ("C1", "C", 1, [2.0, 0.0]),
            ("C2", "C", 10, [-2.0, 0.0]),
        ]
        for sig_id, pert_id, dose, zscores in signatures:
            for gene_name, zscore in zip(targets["gene_name"], zscores):
                rows.append(
                    {
                        "sig_id": sig_id,
                        "pert_id": pert_id,
                        "pert_iname": pert_id,
                        "pert_dose": dose,
                        "pert_dose_unit": "uM",
                        "pert_time": 24,
                        "pert_time_unit": "h",
                        "gene_name": gene_name,
                        "zscore": zscore,
                    }
                )
        result = condition_consistency(pd.DataFrame(rows), targets, 2, 0.6).set_index("pert_id")
        self.assertEqual(result.loc["A", "condition_evidence"], "unreplicated_single_profile")
        self.assertTrue(result.loc["A", "condition_consistency_pass"])
        self.assertEqual(result.loc["B", "condition_evidence"], "replicated_consistent")
        self.assertTrue(result.loc["B", "condition_consistency_pass"])
        self.assertEqual(result.loc["C", "condition_evidence"], "replicated_inconsistent")
        self.assertFalse(result.loc["C", "condition_consistency_pass"])

    def test_cosine_similarity_to_centroid(self) -> None:
        matrix = np.array([[1.0, 0.0], [-1.0, 0.0], [0.0, 1.0]], dtype=np.float32)
        similarity = cosine_similarity_to_centroid(matrix, np.array([1.0, 0.0], dtype=np.float32))
        np.testing.assert_allclose(similarity, [1.0, -1.0, 0.0])


if __name__ == "__main__":
    unittest.main()
