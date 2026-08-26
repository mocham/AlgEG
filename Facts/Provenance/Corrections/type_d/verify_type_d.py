#!/usr/bin/env python3
"""Verify the omitted stable Type D adjacent-height core families."""

from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import sys

from sage.all import LieAlgebra, QQ, RootSystem
from sage.env import SAGE_VERSION


ROOT = next(
    candidate
    for start in (Path.cwd(), Path(__file__).resolve().parent)
    for candidate in (start, *start.parents)
    if (candidate / "References" / "AlgEG-v4" / "algebra.py").is_file()
)
ALGEBRA_ROOT = ROOT / "References" / "AlgEG-v4"
sys.path.insert(0, str(ALGEBRA_ROOT))

from algebra import jacobian_matrix  # noqa: E402
from graphs import (  # noqa: E402
    bipartite_layers,
    contract_source_leaves,
    induced_bipartite,
    root_poset,
    vertices,
)
from roots import (  # noqa: E402
    positive_roots,
    root_string,
    root_tuple,
    subtract,
)


DEFAULT_FULL_MAX = 32
DEFAULT_SUBSET_MAX = 12
DEFAULT_DETERMINANT_MAX = 12


def minus_root(rank, left, right):
    """Return e_left-e_right in Bourbaki simple-root coordinates."""
    coefficients = [0] * rank
    for index in range(left, right):
        coefficients[index - 1] = 1
    return tuple(coefficients)


def plus_root(rank, left, right):
    """Return e_left+e_right in Bourbaki simple-root coordinates."""
    coefficients = [0] * rank
    if right == rank:
        for index in range(left, rank - 1):
            coefficients[index - 1] = 1
        coefficients[rank - 1] = 1
        return tuple(coefficients)

    for index in range(left, right):
        coefficients[index - 1] = 1
    for index in range(right, rank - 1):
        coefficients[index - 1] = 2
    coefficients[rank - 2] = 1
    coefficients[rank - 1] = 1
    return tuple(coefficients)


def expected_core(rank, height):
    """Return the stable K_(rank,m) core, or the empty graph."""
    if height < 2 or height > rank - 2 or height % 2:
        return {}

    m = height // 2
    anchor = rank - height
    chain_sources = [
        plus_root(rank, anchor + offset, rank - offset)
        for offset in range(m)
    ]
    chain_targets = [
        plus_root(rank, anchor - 1 + offset, rank - offset)
        for offset in range(m + 1)
    ]
    minus_right = minus_root(rank, anchor, rank)
    minus_left = minus_root(rank, anchor - 1, rank - 1)
    minus_target = minus_root(rank, anchor - 1, rank)

    adjacency = {
        source: (chain_targets[offset], chain_targets[offset + 1])
        for offset, source in enumerate(chain_sources)
    }
    adjacency[minus_right] = (minus_target, chain_targets[1])
    adjacency[minus_left] = (minus_target, chain_targets[0])
    return {
        source: tuple(sorted(targets))
        for source, targets in sorted(adjacency.items())
    }


def core_rows(core):
    return [
        {
            "source": root_string(source),
            "targets": [root_string(target) for target in targets],
        }
        for source, targets in sorted(core.items())
    ]


def digest(value):
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return sha256(encoded).hexdigest()


def file_sha256(path):
    return sha256(path.read_bytes()).hexdigest()


def type_d_layers(rank):
    roots = positive_roots(("D", rank))
    basic_roots = tuple(
        tuple(int(position == index) for position in range(rank))
        for index in range(rank)
    )
    graph = root_poset(roots, basic_roots)
    return roots, bipartite_layers(graph)


def verify_full_layers(max_rank):
    observations = []
    layers_tested = 0
    nonempty_cores = 0
    for rank in range(4, max_rank + 1):
        roots, layers = type_d_layers(rank)
        maximum_height = max(sum(root) for root in roots)
        nonempty_heights = []
        for height in range(1, maximum_height):
            layers_tested += 1
            actual, _trace = contract_source_leaves(layers.get(height, {}))
            expected = expected_core(rank, height)
            if actual != expected:
                raise AssertionError(f"D{rank} height {height}: unexpected core")
            if actual:
                nonempty_cores += 1
                nonempty_heights.append(height)
        observations.append(
            {
                "rank": rank,
                "positive_roots": len(roots),
                "nonempty_heights": nonempty_heights,
            }
        )
    return {
        "rank_range": [4, max_rank],
        "layers_tested": layers_tested,
        "nonempty_cores": nonempty_cores,
        "observations_sha256": digest(observations),
    }


def verify_target_subsets(max_rank):
    observations = []
    selections_tested = 0
    nonempty_selections = 0
    distinct_cores = 0
    for rank in range(4, max_rank + 1):
        _roots, layers = type_d_layers(rank)
        rank_tested = 0
        rank_nonempty = 0
        signatures = set()
        for height, layer in sorted(layers.items()):
            sources, targets = vertices(layer)
            expected = expected_core(rank, height)
            expected_targets = set(vertices(expected)[1])
            for size in range(len(targets) + 1):
                for selected in combinations(targets, size):
                    rank_tested += 1
                    core, _trace = contract_source_leaves(
                        induced_bipartite(layer, sources, selected)
                    )
                    predicted_nonempty = bool(expected) and expected_targets <= set(selected)
                    if bool(core) != predicted_nonempty:
                        raise AssertionError(
                            f"D{rank} height {height}: target-subset criterion failed"
                        )
                    if core and core != expected:
                        raise AssertionError(
                            f"D{rank} height {height}: unexpected target-induced core"
                        )
                    if core:
                        rank_nonempty += 1
                        signatures.add((height, tuple(core.items())))
        selections_tested += rank_tested
        nonempty_selections += rank_nonempty
        distinct_cores += len(signatures)
        observations.append(
            {
                "rank": rank,
                "selections_tested": rank_tested,
                "nonempty_selections": rank_nonempty,
                "distinct_cores": len(signatures),
            }
        )
    return {
        "rank_range": [4, max_rank],
        "selections_tested": selections_tested,
        "nonempty_selections": nonempty_selections,
        "distinct_cores": distinct_cores,
        "observations_sha256": digest(observations),
    }


def limited_structure_constants(rank, cores):
    lattice = RootSystem(["D", rank]).root_lattice()
    roots_by_tuple = {root_tuple(root): root for root in lattice.roots()}
    basis = LieAlgebra(QQ, cartan_type=("D", rank)).basis()
    constants = {}
    for core in cores:
        for source, targets in core.items():
            for target in targets:
                delta = subtract(target, source)
                key = (source, delta)
                if key in constants:
                    continue
                left = roots_by_tuple[source]
                right = roots_by_tuple[delta]
                bracket = basis[left].bracket(basis[right])
                coefficient = bracket.monomial_coefficients().get(left + right, 0)
                if not coefficient:
                    raise AssertionError(f"D{rank}: missing structure constant {key}")
                constants[key] = int(coefficient)
    return constants


def expected_determinant_variables(rank, height):
    m = height // 2
    indices = [rank - 2 * m - 1, *range(rank - m, rank + 1)]
    return {
        "x" + root_string(tuple(int(position == index) for position in range(1, rank + 1)))
        for index in indices
    }


def verify_determinants(max_rank):
    observations = []
    cores_tested = 0
    coefficients = set()
    d8_height_six = None
    for rank in range(4, max_rank + 1):
        cores = [
            (height, expected_core(rank, height))
            for height in range(2, rank - 1, 2)
        ]
        constants = limited_structure_constants(rank, [core for _height, core in cores])
        rank_records = []
        for height, core in cores:
            jacobian, sources, targets = jacobian_matrix(core, constants)
            if jacobian.nrows() != jacobian.ncols():
                raise AssertionError(f"D{rank} height {height}: core is not square")
            determinant = jacobian.det()
            terms = determinant.dict()
            if len(terms) != 1:
                raise AssertionError(f"D{rank} height {height}: determinant is not monomial")
            exponents, coefficient = next(iter(terms.items()))
            coefficient = int(coefficient)
            coefficients.add(coefficient)
            if abs(coefficient) != 2:
                raise AssertionError(f"D{rank} height {height}: coefficient is not +/-2")
            actual_variables = {
                name
                for name, exponent in zip(
                    jacobian.base_ring().variable_names(), exponents, strict=True
                )
                if exponent
            }
            if any(exponent != 1 for exponent in exponents if exponent):
                raise AssertionError(f"D{rank} height {height}: repeated variable")
            expected_variables = expected_determinant_variables(rank, height)
            if actual_variables != expected_variables:
                raise AssertionError(f"D{rank} height {height}: wrong determinant variables")
            cores_tested += 1
            record = {
                "height": height,
                "coefficient": coefficient,
                "variables": sorted(actual_variables),
            }
            rank_records.append(record)
            if rank == 8 and height == 6:
                d8_height_six = {
                    "adjacency": core_rows(core),
                    "determinant": record,
                }
        observations.append({"rank": rank, "records": rank_records})
    if d8_height_six is None:
        raise AssertionError("determinant range must include D8 height six")
    return {
        "rank_range": [4, max_rank],
        "cores_tested": cores_tested,
        "coefficients": sorted(coefficients),
        "observations_sha256": digest(observations),
        "d8_height_six": d8_height_six,
    }


def build_certificate(full_max, subset_max, determinant_max):
    if full_max < 8 or subset_max < 8 or determinant_max < 8:
        raise ValueError("all verification ranges must include D8")
    return {
        "schema": "type-d-contracted-layers-v1",
        "sage_version": SAGE_VERSION,
        "algorithms": {
            name: file_sha256(ALGEBRA_ROOT / name)
            for name in ("algebra.py", "graphs.py", "roots.py")
        },
        "stable_family": {
            "nonempty_heights": "h=2m with 1<=m<=floor((n-2)/2)",
            "core": "directed C6 joined at one target to directed P_(2m-1)",
            "other_heights": "empty after source-leaf contraction",
            "determinant": "+/-2*x_(n-2m-1)*product(x_i, i=n-m,...,n)",
        },
        "full_layers": verify_full_layers(full_max),
        "target_subsets": verify_target_subsets(subset_max),
        "determinants": verify_determinants(determinant_max),
        "ok": True,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--full-max", type=int, default=DEFAULT_FULL_MAX)
    parser.add_argument("--subset-max", type=int, default=DEFAULT_SUBSET_MAX)
    parser.add_argument("--determinant-max", type=int, default=DEFAULT_DETERMINANT_MAX)
    parser.add_argument("--check", type=Path)
    args = parser.parse_args()

    certificate = build_certificate(
        args.full_max, args.subset_max, args.determinant_max
    )
    rendered = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    if args.check is None:
        sys.stdout.write(rendered)
        return

    recorded = args.check.read_text(encoding="utf-8")
    if recorded != rendered:
        raise SystemExit(f"certificate mismatch: {args.check}")
    print(
        f"verified {args.check}: "
        f"{certificate['full_layers']['layers_tested']} full layers, "
        f"{certificate['target_subsets']['selections_tested']} target selections, "
        f"{certificate['determinants']['cores_tested']} determinants"
    )


if __name__ == "__main__":
    main()
