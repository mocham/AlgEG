import unittest

from sage.all import matrix, PolynomialRing, QQ

from algebra import monomial_minor_certificate
from borel import verify_f4_table, closure_axioms, f4_complement_ideals
from cone import grouped_linear_relations, verify_f4_table_cones
from cycles import verify_exceptional_types, verify_type
from bounded_strata import (affine_level_consistent, compute_field_degree_bound,
                             certify_pair_inventory, enumerate_bounded_psi_pairs,
                             load_or_generate_f4_pairs)
from rotation import root_chain, verify_characteristic_uniform_configurations
from roots import positive_roots, structure_constants
from twisted_cycles import _obstructing_subsets


class V4Tests(unittest.TestCase):
    def test_f4_table_and_order_ideals(self):
        result = verify_f4_table()
        self.assertEqual(result["positive_root_count"], 24)
        self.assertEqual(result["complement_ideal_count"], 105)
        self.assertTrue(result["ok"], result["failures"])

    def test_cone_terms_are_grouped_by_target(self):
        beta1, beta2 = (1, 0), (0, 1)
        alpha = (1, 1)
        layers = {1: {beta1: (alpha,), beta2: (alpha,)}}
        constants = {((0, 1), beta1): 2, ((1, 0), beta2): 3}
        coefficient_matrix, targets, variables = grouped_linear_relations(
            layers, (), (), constants)
        self.assertEqual(targets, (alpha,))
        self.assertEqual(coefficient_matrix.nrows(), 1)
        self.assertEqual(set(coefficient_matrix[0]), {2, 3})

    def test_source_column_combinations(self):
        ring = PolynomialRing(QQ, ("x", "y"))
        x, y = ring.gens()
        jacobian = matrix(ring, [[x, x, 0], [0, x, y]])
        certificate = monomial_minor_certificate(jacobian)
        self.assertTrue(certificate["ok"])
        self.assertEqual(certificate["columns"], [0, 1])
        self.assertEqual(certificate["source_combinations_tested"], 1)

    def test_nonobstructing_target_constraint(self):
        subsets = list(_obstructing_subsets(
            (0, 1, 2), ((0, 1),)))
        self.assertEqual(subsets, [(0, 2), (1, 2)])
        self.assertEqual(list(_obstructing_subsets(
            (0, 1, 2), ((0, 1), (2, 2)))), [])

    def test_f4_table_cones(self):
        result = verify_f4_table_cones()
        self.assertTrue(result["ok"])
        self.assertEqual([record["id"] for record in result["records"]],
                         ["I", "II", "III", "IV"])
        for record in result["records"]:
            self.assertEqual(record["deficiency"],
                             record["relation_count"] - record["rank"])
            self.assertEqual(record["dimension"],
                             len(record["variables"]) - record["rank"])

    def test_full_exceptional_layers(self):
        result = verify_exceptional_types()
        self.assertTrue(result["ok"])

    def test_configuration_i_rotation_chain(self):
        constants = structure_constants(("F", 4))
        psi2 = tuple(tuple(int(d) for d in root)
                     for root in ("0001", "0011", "0100", "0110", "0120"))
        delta = (0, 0, 1, 0)
        chain = root_chain(psi2, delta, constants, (0, 1, 0, 0))
        self.assertEqual([(item[0], item[1], item[2]) for item in chain], [
            ((0, 1, 0, 0), 1, 0),
            ((0, 1, 1, 0), 2, 1),
            ((0, 1, 2, 0), -1, 2),
        ])

    def test_f4_rotation_ok(self):
        result = verify_characteristic_uniform_configurations()
        self.assertTrue(result["ok"], result)
        for record in result["records"]:
            self.assertTrue(record["uniform_for_p_gt_12"])
            self.assertTrue(all(prime <= 12
                                for prime in record["spreading_bad_primes"]))

    def test_affine_level_consistency(self):
        self.assertTrue(affine_level_consistent(((1, 0), (0, 1))))
        self.assertFalse(affine_level_consistent(((1, 0), (2, 0))))

    def test_bounded_bfs_on_A2(self):
        result = enumerate_bounded_psi_pairs(("A", 2), degree_bound=1)
        self.assertGreater(result["independent_seed_count"], 0)
        self.assertEqual(result["pair_count"], len(result["pairs"]))

    def test_f4_field_degree_bound(self):
        result = compute_field_degree_bound(("F", 4))
        self.assertEqual(result["automatic_from_degree"], 2)
        self.assertEqual(result["degrees_requiring_enumeration"], [1])
        bad_primes = {prime for record in result["records"]
                      for prime in record["rank_witness"]["bad_primes"]}
        self.assertTrue(all(prime <= 12 for prime in bad_primes))

    def test_reviewed_pair_inventory(self):
        result = load_or_generate_f4_pairs()
        certificate = certify_pair_inventory(result)
        self.assertFalse(result["first_principles_used"])
        self.assertEqual(certificate["terminal_pair_count"], 4862)
        self.assertEqual(certificate["relevant_independent_seed_count"], 3002)
        self.assertEqual(certificate["covered_independent_seed_count"], 3002)
        self.assertTrue(certificate["ok"], certificate["failures"])


if __name__ == "__main__":
    unittest.main()
