import hashlib
import os
from pathlib import Path
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from sage.all import LieAlgebra, matrix, PolynomialRing, QQ

from algebra import jacobian_matrix, monomial_minor_certificate
from bch import (
    bch,
    bch_terms,
    chevalley_positive_lie_algebra,
    herr_bch_obstruction,
    oriented_jacobian,
    root_diagonal_action,
    verify_bch_against_truncated_log,
    verify_bch_formulas,
)
from borel import (
    F4_TABLE,
    closure_axioms,
    f4_borel_component_count,
    f4_complement_ideals,
    verify_f4_table,
)
from cli import process_result
from cone import grouped_linear_relations, verify_f4_table_cones
from cycles import verify_exceptional_types, verify_type
from bounded_strata import (affine_level_consistent, certify_pair_inventory,
                             compare_f4_pair_inventory,
                             compute_field_degree_bound,
                             empty_psi2_epsilon,
                             enumerate_bounded_psi_pairs,
                             field_degree_satisfies_bound, k_height_layers,
                             load_f4_pair_inventory, minimum_positive_degree,
                             required_field_degree, verify_empty_f4_pair)
from graphs import root_poset
from rotation import root_chain, verify_characteristic_uniform_configurations
from roots import (downward_closed_subsets, positive_roots, root_less_than,
                   simple_roots, structure_constants)
from twisted_cycles import _obstructing_subsets
from type_d import verify_type_d


class V4Tests(unittest.TestCase):
    def test_f4_table_and_order_ideals(self):
        result = verify_f4_table()
        self.assertEqual(result["positive_root_count"], 24)
        self.assertEqual(result["complement_ideal_count"], 105)
        self.assertEqual(
            result["additional_borel_components"][
                "configuration_labelled_record_count"], 19)
        self.assertEqual(
            result["additional_borel_components"][
                "distinct_label_tuple_count"], 17)
        self.assertEqual(
            result["additional_borel_components"][
                "configuration_labelled_record_count_formula"],
            "(p-1)+3=p+2")
        self.assertTrue(result["ok"], result["failures"])

    def test_f4_borel_component_labels_at_17(self):
        result = f4_borel_component_count(17)
        self.assertEqual(len(result["row_I_labels"]), 16)
        self.assertEqual(len(result["row_II_to_IV_labels"]), 3)
        self.assertEqual(result["configuration_labelled_record_count"], 19)
        self.assertEqual(result["distinct_label_tuple_count"], 17)
        self.assertEqual(result["row_I_labels"][0]["label"], [0, 1, 0, 1])
        self.assertEqual(result["row_I_labels"][-1]["label"], [15, 1, 0, 1])
        self.assertEqual(
            [item["label"] for item in result["row_II_to_IV_labels"]],
            [[1, 1, 0, 1], [0, 1, 0, 1], [0, 1, 0, 0]])
        self.assertEqual(
            [item["label"] for item in result["repeated_label_tuples"]],
            [[0, 1, 0, 1], [1, 1, 0, 1]])
        self.assertTrue(result["ok"])

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
        choice_check = result["root_choice_check"]
        self.assertEqual([
            (record["id"], record["rotation_root"])
            for record in choice_check["records"]
        ], [
            ("IV", "0010"),
            ("IV", "0011"),
            ("IV", "1000"),
            ("III", "1000"),
            ("III", "0010"),
        ])
        self.assertEqual({
            (record["id"], record["rotation_root"]):
                record["jacobian_determinant_at_point"]
            for record in choice_check["records"]
        }, {
            ("IV", "0010"): "895795200",
            ("IV", "0011"): "-30800",
            ("IV", "1000"): "315/4",
            ("III", "1000"): "15360",
            ("III", "0010"): "-29160",
        })
        self.assertTrue(choice_check["all_choices_succeed"])
        self.assertEqual(choice_check["failed_choices"], [])
        self.assertFalse(
            choice_check["each_configuration_has_exactly_one_successful_tested_root"])
        for record in choice_check["records"]:
            self.assertEqual(record["outcome"], "success")
            self.assertTrue(record["reason"]["all_equations_zero"])
            self.assertTrue(record["reason"]["jacobian_determinant_nonzero"])
            self.assertTrue(record["reason"]["required_at_ap_values_nonzero"])
            self.assertTrue(record["reason"]["all_prime_divisors_at_most_12"])

    def test_affine_level_consistency(self):
        self.assertTrue(affine_level_consistent(()))
        self.assertTrue(affine_level_consistent(((1, 0), (0, 1))))
        self.assertFalse(affine_level_consistent(((1, 0), (2, 0))))

    def test_empty_f4_pair_is_checked_separately(self):
        result = verify_empty_f4_pair()
        self.assertTrue(result["affine_level_consistent"])
        self.assertTrue(result["closure_conditions_hold_for_all_Psi0"])
        self.assertFalse(result["included_in_pair_inventory"])
        self.assertEqual(result["inventory_domain"], "Psi2 != empty")
        self.assertEqual(
            result["Psi0_scope"], "arbitrary Psi0 subset of Phi+")
        self.assertEqual(result["nonempty_K_complement_count"], 104)
        self.assertEqual(
            result["epsilon_formula"],
            "epsilon=|K_complement|+|Psi0 minus Phi_K|")
        self.assertEqual(result["cardinality_cases_checked"], 1048)
        self.assertEqual(result["minimum_epsilon"], 1)
        self.assertEqual(result["zero_epsilon_cardinality_case_count"], 0)
        roots = set(positive_roots(("F", 4)))
        complement = {min(roots)}
        phi_k = roots - complement
        self.assertEqual(
            empty_psi2_epsilon(complement, phi_k, roots), 2)
        self.assertTrue(result["ok"])

    def test_degree_values_do_not_prune_affine_bfs(self):
        roots = ((1, 0), (1, 1), (1, 2))
        with (
            patch("bounded_strata.positive_roots", return_value=roots),
            patch("bounded_strata.minimum_positive_degree", return_value=1),
        ):
            result = enumerate_bounded_psi_pairs(("T", 2))
        self.assertEqual(result["visited_compatible_subsets"], 7)
        self.assertEqual(result["degree_relation_node_count"], 7)
        self.assertEqual(result["pair_count"], 7)
        self.assertIn(["10", "11", "12"], [
            record["Psi2"] for record in result["pairs"]
        ])

    def test_bounded_bfs_on_A2(self):
        result = enumerate_bounded_psi_pairs(("A", 2))
        self.assertGreater(result["independent_seed_count"], 0)
        self.assertEqual(result["pair_count"], len(result["pairs"]))

    def test_exact_minimum_positive_degree(self):
        self.assertEqual(minimum_positive_degree(
            ((1, 1),), ((1, 0), (0, 1))), 2)
        self.assertIsNone(minimum_positive_degree(
            ((1, 0),), ((0, 1),)))
        self.assertIsNone(minimum_positive_degree(((1, 0),), ()))

    def test_corrected_field_degree_bound(self):
        self.assertEqual(required_field_degree(13, 4), 3)
        self.assertEqual(required_field_degree(13, 5), 12)
        self.assertEqual(required_field_degree(13, 6), 2)
        self.assertTrue(field_degree_satisfies_bound(13, 2, 6))
        self.assertFalse(field_degree_satisfies_bound(13, 2, 4))
        self.assertTrue(field_degree_satisfies_bound(13, 1, None))

    def test_f4_field_degree_bound(self):
        result = compute_field_degree_bound(("F", 4))
        self.assertEqual(result["automatic_from_degree"], 2)
        self.assertEqual(result["degrees_requiring_enumeration"], [1])
        bad_primes = {prime for record in result["records"]
                      for prime in record["rank_witness"]["bad_primes"]}
        self.assertTrue(all(prime <= 12 for prime in bad_primes))

    def test_reviewed_pair_inventory(self):
        reviewed = load_f4_pair_inventory()
        recomputed = enumerate_bounded_psi_pairs(("F", 4))
        validation = certify_pair_inventory(reviewed)
        comparison = compare_f4_pair_inventory(recomputed, reviewed)
        self.assertEqual(recomputed["minimum_degree_mode"], "exact")
        self.assertFalse(recomputed["degree_cutoff_used"])
        self.assertFalse(recomputed["empty_subset_visited"])
        self.assertEqual(recomputed["degree_relation_node_count"], 0)
        self.assertEqual(recomputed["visited_compatible_subsets"], 50347)
        self.assertEqual(validation["terminal_pair_count"], 4862)
        self.assertEqual(validation["admissible_terminal_pair_count"], 4862)
        self.assertEqual(validation["finite_minimum_degree_pair_count"], 0)
        self.assertEqual(validation["undefined_minimum_degree_pair_count"], 4862)
        self.assertEqual(
            validation["pair_list_sha256"],
            "e1dc8ced6a02269877dbfbc4a908404657d8a415a3b503f54d429888af29d4be")
        self.assertEqual(validation["relevant_independent_seed_count"], 3002)
        self.assertEqual(validation["covered_independent_seed_count"], 3002)
        self.assertTrue(validation["ok"], validation["failures"])
        self.assertEqual(
            comparison["recomputed_pairs_sha256"],
            "e1dc8ced6a02269877dbfbc4a908404657d8a415a3b503f54d429888af29d4be")
        self.assertEqual(
            comparison["reviewed_pairs_sha256"],
            "e1dc8ced6a02269877dbfbc4a908404657d8a415a3b503f54d429888af29d4be")
        self.assertTrue(comparison["ok"], comparison)

    def test_homogeneous_bch_terms(self):
        algebra = LieAlgebra(QQ, 2, step=4)
        X, Y = algebra.basis(1)
        terms = bch_terms(X, Y, max_class=4)
        self.assertEqual(terms[0], X + Y)
        self.assertEqual(terms[1], QQ(1) / 2 * X.bracket(Y))
        self.assertEqual(
            terms[2],
            QQ(1) / 12 * X.bracket(X.bracket(Y))
            + QQ(1) / 12 * X.bracket(Y).bracket(Y))
        self.assertEqual(
            terms[3], QQ(1) / 24 * X.bracket(X.bracket(Y).bracket(Y)))

    def test_symbolic_a2_bch_obstruction(self):
        ring = PolynomialRing(
            QQ, ("fa", "fb", "ga", "gb", "pa", "pb", "qa", "qb"))
        fa, fb, ga, gb, pa, pb, qa, qb = ring.gens()
        algebra = LieAlgebra(
            ring, {("a", "b"): {"c": 1}},
            names=("a", "b", "c"), nilpotent=True)
        a, b, _ = algebra.lie_algebra_generators()
        obstruction = herr_bch_obstruction(
            fa * a + fb * b,
            ga * a + gb * b,
            pa * a + pb * b,
            qa * a + qb * b,
            max_class=2)
        expected = QQ(1) / 2 * (
            ga * qb - gb * qa - fa * pb + fb * pa)
        self.assertEqual(obstruction["c"], expected)

    def test_symbolic_b2_bch_obstruction(self):
        names = (
            "fa", "fb", "fc", "ga", "gb", "gc",
            "pa", "pb", "pc", "qa", "qb", "qc",
        )
        ring = PolynomialRing(QQ, names)
        fa, fb, fc, ga, gb, gc, pa, pb, pc, qa, qb, qc = ring.gens()
        algebra = LieAlgebra(
            ring,
            {("a", "b"): {"c": 1}, ("a", "c"): {"d": 2}},
            names=("a", "b", "c", "d"), nilpotent=True)
        a, b, c, _ = algebra.lie_algebra_generators()
        obstruction = herr_bch_obstruction(
            fa * a + fb * b + fc * c,
            ga * a + gb * b + gc * c,
            pa * a + pb * b + pc * c,
            qa * a + qb * b + qc * c,
            max_class=3)
        right = (
            ga * qc - gc * qa
            + QQ(1) / 6 * (ga - qa) * (ga * qb - gb * qa))
        left = (
            fa * pc - fc * pa
            + QQ(1) / 6 * (fa - pa) * (fa * pb - fb * pa))
        self.assertEqual(obstruction["d"], right - left)

    def test_b2_top_coordinate(self):
        algebra = LieAlgebra(
            QQ,
            {("a", "b"): {"c": 1}, ("a", "c"): {"d": 2}},
            names=("a", "b", "c", "d"), nilpotent=True)
        a, b, c, d = algebra.lie_algebra_generators()
        basis = {"a": a, "b": b, "c": c, "d": d}
        f_lower = a + b
        g_lower = 2 * a + QQ(4) / 3 * b + QQ(1) / 3 * c
        phi_f_g = root_diagonal_action(
            g_lower, basis, {"a": 2, "b": 4, "c": 8, "d": 16})
        gamma_g_f = root_diagonal_action(
            f_lower, basis, {"a": 3, "b": 5, "c": 15, "d": 45})
        obstruction = herr_bch_obstruction(
            f_lower, g_lower, phi_f_g, gamma_g_f, max_class=3)
        self.assertEqual(obstruction["d"], -4)

        corrected_f = f_lower + QQ(1) / 11 * d
        corrected_gamma_g_f = root_diagonal_action(
            corrected_f, basis, {"a": 3, "b": 5, "c": 15, "d": 45})
        self.assertEqual(
            bch(corrected_f, phi_f_g, max_class=3),
            bch(g_lower, corrected_gamma_g_f, max_class=3))

    def test_f4_bch_height_and_jacobians(self):
        f4 = chevalley_positive_lie_algebra(("F", 4))
        self.assertEqual(max(sum(root) for root in f4.roots), 11)
        self.assertEqual(f4.algebra.step(), 11)
        comparison = verify_bch_against_truncated_log(max_class=11)
        self.assertEqual(comparison["maximum_class"], 11)
        self.assertEqual(
            [record["class"] for record in comparison["classes"]],
            list(range(1, 12)))
        self.assertEqual(
            comparison["classes"][-1]["nonzero_word_count"], 2046)
        self.assertTrue(comparison["ok"], comparison)

        roots = tuple(sorted(positive_roots(("F", 4))))
        constants = structure_constants(("F", 4))
        root_set = set(roots)
        compared = 0
        for complement in downward_closed_subsets(roots, root_less_than):
            phi_k = root_set - set(complement)
            if not phi_k:
                continue
            graph = root_poset(phi_k, simple_roots(phi_k))
            for layer in k_height_layers(graph).values():
                independent, sources, targets = oriented_jacobian(
                    layer, f4.constants)
                local, local_sources, local_targets = jacobian_matrix(
                    layer, constants)
                self.assertEqual(sources, local_sources)
                self.assertEqual(targets, local_targets)
                self.assertEqual(independent, local)
                compared += 1
        self.assertGreater(compared, 0)

    def test_f4_rotation_has_no_higher_bch_contribution(self):
        for record in F4_TABLE:
            psi0 = tuple(
                tuple(int(digit) for digit in root) for root in record["Psi0"])
            psi2 = tuple(
                tuple(int(digit) for digit in root) for root in record["Psi2"])
            for target in psi2:
                for perturbation in psi0:
                    residual = tuple(
                        target[index] - perturbation[index]
                        for index in range(len(target)))
                    if any(coefficient < 0 for coefficient in residual):
                        continue
                    frontier = {(0,) * len(target)}
                    for factor_count in range(1, sum(target) + 1):
                        frontier = {
                            tuple(total[index] + root[index]
                                  for index in range(len(target)))
                            for total in frontier for root in psi2
                            if all(total[index] + root[index] <= residual[index]
                                   for index in range(len(target)))
                        }
                        if factor_count >= 2:
                            self.assertNotIn(
                                residual, frontier,
                                (record["id"], target, perturbation, factor_count))

    def test_bch_result(self):
        self.assertTrue(verify_bch_formulas()["ok"])

    def test_type_d_families(self):
        result = verify_type_d()
        self.assertTrue(result["ok"])
        type_d_path = Path(__file__).resolve().parents[1] / "type_d.py"
        self.assertEqual(
            result["implementation_sha256"]["type_d.py"],
            hashlib.sha256(type_d_path.read_bytes()).hexdigest())
        self.assertEqual(result["full_layers"]["layers_tested"], 928)
        self.assertEqual(result["full_layers"]["nonempty_cores"], 225)
        self.assertEqual(result["target_subsets"]["selections_tested"], 15888)
        self.assertEqual(result["target_subsets"]["distinct_cores"], 25)
        self.assertEqual(result["determinants"]["cores_tested"], 25)
        self.assertEqual(result["determinants"]["coefficients"], [-2, 2])
        self.assertEqual(
            result["determinants"]["d8_height_six"]["determinant"]["variables"],
            ["x00000001", "x00000010", "x00000100", "x00001000", "x10000000"])

    def test_cli_check_does_not_write(self):
        with TemporaryDirectory() as temporary:
            directory = Path(temporary)
            process_result("test-v14", {"ok": True}, directory)
            json_path = directory / "test-v14.json"
            sha_path = directory / "test-v14.sha256"
            before = (
                json_path.read_bytes(), sha_path.read_bytes(),
                json_path.stat().st_mtime_ns, sha_path.stat().st_mtime_ns,
            )
            with patch.object(
                    Path, "write_bytes",
                    side_effect=AssertionError("--check attempted to write")):
                process_result("test-v14", {"ok": True}, directory, check=True)
            self.assertEqual(before, (
                json_path.read_bytes(), sha_path.read_bytes(),
                json_path.stat().st_mtime_ns, sha_path.stat().st_mtime_ns,
            ))

            missing = directory / "missing"
            with self.assertRaises(SystemExit):
                process_result("test-v14", {"ok": True}, missing, check=True)
            self.assertFalse(missing.exists())

    def test_cli_check_clean_copy_writes_nothing(self):
        def snapshot(directory):
            entries = {}
            for path in sorted(directory.rglob("*")):
                relative = str(path.relative_to(directory))
                if path.is_dir():
                    entries[relative] = ("directory",)
                else:
                    entries[relative] = (
                        "file",
                        hashlib.sha256(path.read_bytes()).hexdigest(),
                        path.stat().st_mtime_ns,
                    )
            return entries

        implementation = Path(__file__).resolve().parents[1]
        with TemporaryDirectory() as temporary:
            copied = Path(temporary) / "v4"
            shutil.copytree(
                implementation, copied,
                ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
            before = snapshot(copied)
            environment = os.environ.copy()
            environment.pop("PYTHONDONTWRITEBYTECODE", None)
            environment.pop("PYTHONPYCACHEPREFIX", None)
            completed = subprocess.run(
                [
                    sys.executable, "cli.py", "borel", "--check",
                    "--output-directory", "DATA",
                ],
                cwd=copied,
                env=environment,
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            self.assertEqual(snapshot(copied), before)
            self.assertFalse(any(
                path.name == "__pycache__" for path in copied.rglob("*")))


if __name__ == "__main__":
    unittest.main()
