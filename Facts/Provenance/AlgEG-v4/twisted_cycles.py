"""Twisted-cycle certificates with the non-obstructing constraint.

Input data are generated in-process by the local ``elliptic_weyl`` module.
"""

from itertools import combinations

from algebra import jacobian_matrix, monomial_minor_certificate
from graphs import contract_source_leaves, signature, vertices
from roots import positive_roots, root_string, structure_constants, subtract


EXCEPTIONAL_TYPES = ("F4", "E6", "E7", "E8")


def _absolute_roots(relative_root):
    return tuple(tuple(int(digit) for digit in root)
                 for root in relative_root.split("-"))


def _obstructing_subsets(target_ids, nonobstructing_pairs):
    """Yield target subsets containing no non-obstructing singleton or pair."""
    target_ids = tuple(sorted(target_ids))
    forbidden = tuple(frozenset(pair) for pair in nonobstructing_pairs)
    for size in range(2, len(target_ids) + 1):
        for subset in combinations(target_ids, size):
            selected = frozenset(subset)
            if all(not pair.issubset(selected) for pair in forbidden):
                yield subset


def _restrict_targets(graph, selected_targets):
    selected_targets = set(selected_targets)
    return {
        source: tuple(target for target in targets if target in selected_targets)
        for source, targets in sorted(graph.items())
    }


def _absolute_graph(relative_graph, roots):
    relative_sources, relative_targets = vertices(relative_graph)
    sources = tuple(sorted({root for orbit in relative_sources
                            for root in _absolute_roots(orbit)}))
    targets = tuple(sorted({root for orbit in relative_targets
                            for root in _absolute_roots(orbit)}))
    return {
        source: tuple(target for target in targets
                      if subtract(target, source) in roots)
        for source in sources
    }


def verify_twisted_records(records):
    """Verify obstructing relative cycles via their absolute source minors.

    records: dict mapping CartanType -> { Levi -> { word -> { ... } } }
    """
    type_records = []
    overall_ok = True
    for type_name in EXCEPTIONAL_TYPES:
        if type_name not in records:
            type_records.append({
                "cartan_type": type_name,
                "configurations": 0,
                "ok": True,
                "failures": [],
                "note": "no data for this type",
            })
            continue
        roots = set(positive_roots((type_name[0], int(type_name[1:]))))
        constants = structure_constants((type_name[0], int(type_name[1:])))
        cache = {}
        failures = []
        configurations = relative_cores = target_subsets = 0
        contracted = absolute_cores = source_combinations = 0
        for levi, words in sorted(records[type_name].items()):
            for word, record in sorted(words.items()):
                configurations += 1
                orbit_ids = {root: int(index)
                             for root, index in record["orbit_dict"].items()}
                orbit_strings = {int(index): root
                                 for index, root in record["orbit_str"].items()}
                nonobstructing_pairs = tuple(
                    tuple(sorted(pair))
                    for pair in record["nonobstructing_orbit_pairs"])
                for height_str, layer in sorted(
                        record["bipartites"].items(),
                        key=lambda item: int(item[0])):
                    graph = {source: tuple(sorted(targets))
                             for source, targets in layer.items()}
                    relative_core, _ = contract_source_leaves(graph)
                    if not relative_core:
                        continue
                    relative_cores += 1
                    _, targets = vertices(relative_core)
                    target_ids = tuple(orbit_ids[target] for target in targets)
                    for selected_ids in _obstructing_subsets(
                            target_ids, nonobstructing_pairs):
                        target_subsets += 1
                        selected = tuple(orbit_strings[index]
                                         for index in selected_ids)
                        restricted = _restrict_targets(relative_core, selected)
                        restricted_core, _ = contract_source_leaves(restricted)
                        if not restricted_core:
                            contracted += 1
                            continue
                        absolute = _absolute_graph(restricted_core, roots)
                        absolute_core, _ = contract_source_leaves(absolute)
                        if not absolute_core:
                            contracted += 1
                            continue
                        absolute_cores += 1
                        core_signature = signature(absolute_core)
                        if core_signature not in cache:
                            jacobian, sources, core_targets = jacobian_matrix(
                                absolute_core, constants)
                            certificate = monomial_minor_certificate(jacobian)
                            certificate["sources"] = [root_string(root)
                                                       for root in sources]
                            certificate["targets"] = [root_string(root)
                                                       for root in core_targets]
                            certificate["selected_sources"] = (
                                [root_string(sources[index])
                                 for index in certificate["columns"]]
                                if certificate["columns"] is not None else None)
                            cache[core_signature] = certificate
                            source_combinations += certificate[
                                "source_combinations_tested"]
                        certificate = cache[core_signature]
                        if not certificate["ok"]:
                            failures.append({
                                "levi": levi,
                                "word": word,
                                "height": int(height_str),
                                "relative_targets": list(selected),
                                "absolute_sources": certificate["sources"],
                                "absolute_targets": certificate["targets"],
                            })
        ok = not failures
        overall_ok &= ok
        type_records.append({
            "cartan_type": type_name,
            "configurations": configurations,
            "relative_cores": relative_cores,
            "obstructing_target_subsets": target_subsets,
            "contracted": contracted,
            "absolute_cores": absolute_cores,
            "unique_absolute_cores": len(cache),
            "source_combinations_tested": source_combinations,
            "certificates": list(cache.values()),
            "failures": failures,
            "ok": ok,
        })
    return {"records": type_records, "ok": overall_ok}
