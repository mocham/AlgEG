"""Self-contained bounded fine-strata generation.

Candidates are not obtained by scanning all subsets of positive roots.  The
algorithm starts from linearly independent seeds for the equations f(alpha)=1
and expands compatible subsets by deterministic breadth-first search.  At each
node it computes the integral difference lattice, the forced Psi0 set, closure,
and the requested degree data.
"""

from collections import deque
import hashlib
import json
from itertools import combinations
from math import gcd
from pathlib import Path

from sage.all import FreeModule, matrix, QQ, vector, ZZ

from cone import linear_cone_certificate
from algebra import jacobian_matrix, polynomial_rank_certificate
from graphs import compactify_graph, contract_source_leaves, root_poset, vertices
from roots import (add, downward_closed_subsets, positive_roots, root_less_than,
                   simple_roots, structure_constants)


PAIR_INVENTORY = Path(__file__).resolve().parent / "DATA" / "psi-pairs-f4-v12.json"
EXPECTED_F4_PAIR_DIGEST = "e1dc8ced6a02269877dbfbc4a908404657d8a415a3b503f54d429888af29d4be"


def affine_level_consistent(subset):
    """Return whether one linear functional takes value one on the subset."""
    if not subset:
        return True
    coefficient = matrix(QQ, subset)
    augmented = coefficient.augment(matrix(QQ, len(subset), 1, [1] * len(subset)))
    return coefficient.rank() == augmented.rank()


def is_independent_seed(subset):
    if not subset:
        return False
    return matrix(QQ, subset).rank() == len(subset)


def forced_psi0(subset, roots):
    """Roots in the integral span of differences of roots in Psi2."""
    if len(subset) <= 1:
        return ()
    base = vector(ZZ, subset[0])
    differences = [vector(ZZ, root) - base for root in subset[1:]]
    lattice = FreeModule(ZZ, len(base)).submodule(differences)
    return tuple(root for root in roots if vector(ZZ, root) in lattice)


def closure_satisfied(psi0, psi2, root_set):
    psi2 = set(psi2)
    return all(add(alpha, beta) not in root_set or add(alpha, beta) in psi2
               for alpha in psi2 for beta in psi0)


def minimum_positive_degree(psi0, psi2):
    """Find the least positive Psi2-sum in Psi0, if one exists."""
    targets = set(psi0)
    if not targets or not psi2:
        return None

    # Positive-root height is additive, so a sum beyond the largest target
    # height cannot return to Psi0.  This makes the search exact and finite.
    maximum = max(sum(target) for target in targets)
    frontier = {(0,) * len(psi2[0])}
    for degree in range(1, maximum + 1):
        frontier = {add(total, root) for total in frontier for root in psi2}
        if frontier & targets:
            return degree
    return None


def required_field_degree(prime, minimum_degree):
    """Return the corrected degree lower bound attached to a finite minimum."""
    if prime <= 1:
        raise ValueError("prime must be greater than one")
    if minimum_degree is None:
        return None
    if minimum_degree <= 0:
        raise ValueError("minimum degree must be positive")
    return (prime - 1) // gcd(prime - 1, minimum_degree)


def field_degree_satisfies_bound(prime, field_degree, minimum_degree):
    """Test d_F >= (p-1)/gcd(p-1, minimum_degree)."""
    if field_degree <= 0:
        raise ValueError("field degree must be positive")
    required = required_field_degree(prime, minimum_degree)
    return required is None or field_degree >= required


def independent_seeds(roots):
    """All nonempty linearly independent subsets, in deterministic order."""
    rank = len(roots[0])
    for size in range(1, rank + 1):
        for subset in combinations(roots, size):
            if is_independent_seed(subset):
                yield tuple(subset)


def k_height_layers(graph):
    """Bipartite layers for the intrinsic K-height, not ambient root height."""
    parents = {root: set() for root in graph}
    for source, targets in graph.items():
        for target in targets:
            parents[target].add(source)
    heights = {}
    for root in sorted(graph, key=lambda item: (sum(item), item)):
        heights[root] = 1 + max((heights[parent] for parent in parents[root]),
                                default=0)
    layers = {}
    for height in sorted(set(heights.values())):
        sources = tuple(root for root in graph if heights[root] == height)
        layer = {
            source: tuple(target for target in graph[source]
                          if heights[target] == height + 1)
            for source in sources
        }
        if any(layer.values()):
            layers[height] = layer
    return layers


def compute_field_degree_bound(cartan_type=("F", 4)):
    """Compute the degree at which the maximal cone estimate is automatic."""
    roots = tuple(sorted(positive_roots(cartan_type)))
    root_set = set(roots)
    constants = structure_constants(cartan_type)
    ideals = downward_closed_subsets(roots, root_less_than)
    records = []
    automatic_degree = 1
    for complement in ideals:
        phi_k = root_set - set(complement)
        if not phi_k:
            continue
        graph = root_poset(phi_k, simple_roots(phi_k))
        certificate = linear_cone_certificate(
            k_height_layers(graph), roots, roots, constants)
        needed = (certificate["dimension"] + len(roots) - 1) // len(roots)
        automatic_degree = max(automatic_degree, needed)
        records.append({
            "K_complement": ["".join(map(str, root))
                             for root in sorted(complement)],
            "residual_cone_dimension": certificate["dimension"],
            "automatic_degree": needed,
            "rank_witness": certificate["rank_witness"],
        })
    return {
        "cartan_type": f"{cartan_type[0]}{cartan_type[1]}",
        "positive_root_count": len(roots),
        "K_count": len(ideals),
        "automatic_from_degree": automatic_degree,
        "degrees_requiring_enumeration": list(range(1, automatic_degree)),
        "records": records,
    }


def _generic_obstruction(layers, constants):
    """Sum generic Jacobian coranks, recording monomial uniformity witnesses."""
    obstruction = 0
    witnesses = []
    for height, layer in sorted(layers.items()):
        core, trace = contract_source_leaves(compactify_graph(layer))
        if not core:
            witnesses.append({"height": height, "contracted": len(trace),
                              "corank": 0})
            continue
        jacobian, sources, targets = jacobian_matrix(core, constants)
        rank_certificate = polynomial_rank_certificate(jacobian)
        rank = rank_certificate["rank"]
        corank = jacobian.nrows() - rank
        obstruction += corank
        witnesses.append({
            "height": height,
            "sources": ["".join(map(str, root)) for root in sources],
            "targets": ["".join(map(str, root)) for root in targets],
            "generic_rank": rank,
            "corank": corank,
            "rank_certificate": rank_certificate,
            "uniform_for_p_gt_12": all(
                prime <= 12 for prime in rank_certificate["bad_primes"]),
        })
    return obstruction, witnesses


def classify_degree_one_f4(psi_result, prime=13):
    """Classify every generated F4 pair against every intrinsic K-poset."""
    cartan_type = ("F", 4)
    roots = tuple(sorted(positive_roots(cartan_type)))
    root_set = set(roots)
    constants = structure_constants(cartan_type)
    ideals = downward_closed_subsets(roots, root_less_than)
    pairs = [
        (tuple(tuple(int(digit) for digit in root) for root in item["Psi0"]),
         tuple(tuple(int(digit) for digit in root) for root in item["Psi2"]),
         item["affine_solution_dimension"])
        for item in psi_result["pairs"]
        if field_degree_satisfies_bound(
            prime, 1, item.get("minimum_positive_degree"))
    ]
    epsilon_zero = []
    epsilon_one = []
    checked = 0
    prefiltered = 0
    uniform = True
    digest = hashlib.sha256()
    for complement in ideals:
        if not complement:
            continue
        phi_k = root_set - set(complement)
        delta_k = set(simple_roots(phi_k))
        graph = root_poset(phi_k, tuple(sorted(delta_k)))
        maximal_layers = k_height_layers(graph)
        maximal_target_count = sum(len(vertices(layer)[1])
                                   for layer in maximal_layers.values())
        for psi0, psi2, solution_dimension in pairs:
            relations = (len(complement) + (len(roots[0]) - solution_dimension)
                         - len(set(psi2) & delta_k)
                         + len(set(psi0) - phi_k))
            # Corank is at most the number of target rows.  This sound bound
            # discards configurations that cannot have epsilon <= 1.
            if relations - maximal_target_count >= 2:
                prefiltered += 1
                digest.update(json.dumps([
                    ["".join(map(str, root)) for root in sorted(complement)],
                    ["".join(map(str, root)) for root in psi0],
                    ["".join(map(str, root)) for root in psi2],
                    "prefiltered", relations, maximal_target_count,
                ], separators=(",", ":")).encode())
                continue
            filtered_layers = {}
            for height, layer in maximal_layers.items():
                filtered = {
                    source: tuple(target for target in targets if target in psi2)
                    for source, targets in layer.items()
                }
                if any(filtered.values()):
                    filtered_layers[height] = filtered
            obstruction, witnesses = _generic_obstruction(filtered_layers, constants)
            uniform &= all(witness.get("uniform_for_p_gt_12", True)
                           for witness in witnesses)
            epsilon = relations - obstruction
            checked += 1
            digest.update(json.dumps([
                ["".join(map(str, root)) for root in sorted(complement)],
                ["".join(map(str, root)) for root in psi0],
                ["".join(map(str, root)) for root in psi2],
                "checked", relations, obstruction, epsilon,
            ], separators=(",", ":")).encode())
            if epsilon <= 1:
                item = {
                    "K_complement": ["".join(map(str, root))
                                     for root in sorted(complement)],
                    "Psi0": ["".join(map(str, root)) for root in psi0],
                    "Psi2": ["".join(map(str, root)) for root in psi2],
                    "relations": relations,
                    "obstruction": obstruction,
                    "epsilon": epsilon,
                    "type": ("D" if set(psi2).issubset(delta_k)
                             else "U" if obstruction else "O"),
                    "witnesses": witnesses,
                }
                (epsilon_zero if epsilon <= 0 else epsilon_one).append(item)
    return {
        "K_count": len(ideals),
        "admissible_nonempty_K_count": len(ideals) - 1,
        "field_degree": 1,
        "prime_used_for_degree_filter": prime,
        "psi_pair_count": len(pairs),
        "cartesian_product_size": len(ideals) * len(pairs),
        "excluded_empty_complement_pairs": len(pairs),
        "classified_pair_count": prefiltered + checked,
        "prefiltered": prefiltered,
        "jacobian_checked": checked,
        "epsilon_zero_count": len(epsilon_zero),
        "epsilon_one_count": len(epsilon_one),
        "classification_sha256": digest.hexdigest(),
        "uniform_for_p_gt_12": uniform,
        "ok": len(epsilon_zero) == 4 and uniform,
        "epsilon_zero": epsilon_zero,
        "epsilon_one": epsilon_one,
    }


def _pair_key(item):
    return (
        tuple(tuple(int(digit) for digit in root) for root in item["Psi0"]),
        tuple(tuple(int(digit) for digit in root) for root in item["Psi2"]),
    )


def certify_pair_inventory(psi_result, prime=13, field_degree=1):
    """Check the stored terminal pairs and all relevant independent seeds."""
    roots = tuple(sorted(positive_roots(("F", 4))))
    root_set = set(roots)
    pair_digest = hashlib.sha256(json.dumps(
        psi_result["pairs"], sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    pair_keys = set()
    finite_minimum_degrees = []
    failures = []
    if psi_result.get("pair_count") != len(psi_result["pairs"]):
        failures.append({"kind": "pair_count"})
    if pair_digest != EXPECTED_F4_PAIR_DIGEST:
        failures.append({"kind": "pair_digest", "actual": pair_digest})
    for item in psi_result["pairs"]:
        psi0, psi2 = _pair_key(item)
        if not affine_level_consistent(psi2):
            failures.append({"kind": "affine", "Psi2": item["Psi2"]})
        if forced_psi0(psi2, roots) != psi0:
            failures.append({"kind": "Psi0", "Psi2": item["Psi2"]})
        if not closure_satisfied(psi0, psi2, root_set):
            failures.append({"kind": "closure", "Psi2": item["Psi2"]})
        degree = minimum_positive_degree(psi0, psi2)
        if item.get("minimum_positive_degree") != degree:
            failures.append({"kind": "minimum_degree", "Psi2": item["Psi2"]})
        if degree is not None:
            finite_minimum_degrees.append(degree)
        if field_degree_satisfies_bound(prime, field_degree, degree):
            pair_keys.add((psi0, psi2))

    relevant_seeds = 0
    covered_seeds = 0
    for seed in independent_seeds(roots):
        psi0 = forced_psi0(seed, roots)
        if not closure_satisfied(psi0, seed, root_set):
            continue
        degree = minimum_positive_degree(psi0, seed)
        if not field_degree_satisfies_bound(prime, field_degree, degree):
            continue
        relevant_seeds += 1
        if (psi0, seed) in pair_keys:
            covered_seeds += 1
        else:
            failures.append({
                "kind": "independent_seed",
                "Psi2": ["".join(map(str, root)) for root in seed],
            })
    return {
        "terminal_pair_count": len(psi_result["pairs"]),
        "pair_list_sha256": pair_digest,
        "matches_first_principles_pair_list": pair_digest == EXPECTED_F4_PAIR_DIGEST,
        "admissible_terminal_pair_count": len(pair_keys),
        "finite_minimum_degree_pair_count": len(finite_minimum_degrees),
        "undefined_minimum_degree_pair_count": (
            len(psi_result["pairs"]) - len(finite_minimum_degrees)),
        "finite_minimum_degrees": sorted(set(finite_minimum_degrees)),
        "degree_filter": {
            "prime": prime,
            "field_degree": field_degree,
            "required_field_degree": "(p-1)/gcd(p-1,minimum_degree)",
        },
        "relevant_independent_seed_count": relevant_seeds,
        "covered_independent_seed_count": covered_seeds,
        "failures": failures,
        "ok": not failures,
    }


def load_or_generate_f4_pairs(progress=None):
    """Use the reviewed inventory, regenerating it only when absent."""
    if PAIR_INVENTORY.exists():
        document = json.loads(PAIR_INVENTORY.read_text())
        result = document["payload"]
        result["inventory_source"] = "reviewed-data"
        result["first_principles_used"] = False
        return result

    result = enumerate_bounded_psi_pairs(("F", 4), progress=progress)
    document = {
        "schema": "grob-pair-inventory-v12",
        "name": "psi-pairs-f4-v12",
        "payload": result,
    }
    encoded = (json.dumps(document, sort_keys=True, indent=2) + "\n").encode()
    PAIR_INVENTORY.parent.mkdir(exist_ok=True)
    PAIR_INVENTORY.write_bytes(encoded)
    digest = hashlib.sha256(encoded).hexdigest()
    PAIR_INVENTORY.with_suffix(".sha256").write_text(
        f"{digest}  {PAIR_INVENTORY.name}\n")
    result["inventory_source"] = "generated-first-principles"
    result["first_principles_used"] = True
    return result


def verify_f4_bounded_enumeration(progress=None):
    """Run the complete v12 proof computation from generated root data."""
    degree = compute_field_degree_bound(("F", 4))
    degrees = degree["degrees_requiring_enumeration"]
    if degrees != [1]:
        return {"ok": False, "error": f"unexpected degree range {degrees}"}
    pairs = load_or_generate_f4_pairs(progress=progress)
    inventory = certify_pair_inventory(pairs, prime=13, field_degree=1)
    pair_digest = hashlib.sha256(json.dumps(
        pairs["pairs"], sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    classification = classify_degree_one_f4(pairs, prime=13)
    return {
        "degree_bound": {
            key: value for key, value in degree.items() if key != "records"
        },
        "degree_rank_bad_primes": sorted({
            prime for record in degree["records"]
            for prime in record["rank_witness"]["bad_primes"]
        }),
        "generation": {
            key: value for key, value in pairs.items() if key != "pairs"
        },
        "psi_pairs_sha256": pair_digest,
        "inventory_certification": inventory,
        "classification": {
            key: value for key, value in classification.items()
            if key not in ("epsilon_zero", "epsilon_one")
        },
        "epsilon_zero": classification["epsilon_zero"],
        "ok": (classification["ok"] and inventory["ok"]
               and pairs["pair_count"] == 4862
               and inventory["finite_minimum_degree_pair_count"] == 0
               and degree["automatic_from_degree"] == 2),
    }


def enumerate_bounded_psi_pairs(cartan_type=("F", 4), progress=None):
    """Enumerate all affine-compatible saturated Psi pairs by BFS.

    Every compatible set contains an independent subset with the same affine
    span.  We seed by those independent subsets and add only roots preserving
    f(alpha)=1.  Canonical frozenset keys remove duplicates reached from
    different seeds.  Exact minimum degree is computed at every accepted node
    and included in the output.  Degree is not used for branch pruning because
    the corrected field-degree condition is not monotone in minimum degree.
    """
    roots = tuple(sorted(positive_roots(cartan_type)))
    root_set = set(roots)
    queue = deque(frozenset(seed) for seed in independent_seeds(roots))
    queued = set(queue)
    visited = set()
    records = {}
    while queue:
        current = queue.popleft()
        queued.discard(current)
        if current in visited:
            continue
        visited.add(current)
        subset = tuple(sorted(current))
        if not affine_level_consistent(subset):
            continue
        psi0 = forced_psi0(subset, roots)
        degree = minimum_positive_degree(psi0, subset)
        if closure_satisfied(psi0, subset, root_set):
            key = (psi0, subset)
            records[key] = {
                "Psi0": ["".join(map(str, root)) for root in psi0],
                "Psi2": ["".join(map(str, root)) for root in subset],
                "minimum_positive_degree": degree,
                "affine_solution_dimension": len(roots[0]) - matrix(QQ, subset).rank(),
            }
        for root in roots:
            if root in current:
                continue
            expanded = frozenset((*current, root))
            if expanded in visited or expanded in queued:
                continue
            if affine_level_consistent(tuple(sorted(expanded))):
                queue.append(expanded)
                queued.add(expanded)
        if progress and len(visited) % progress == 0:
            print(f"visited={len(visited)} queue={len(queue)} records={len(records)}")
    return {
        "cartan_type": f"{cartan_type[0]}{cartan_type[1]}",
        "minimum_degree_mode": "exact",
        "independent_seed_count": sum(1 for _ in independent_seeds(roots)),
        "visited_compatible_subsets": len(visited),
        "pair_count": len(records),
        "pairs": [records[key] for key in sorted(records)],
    }
