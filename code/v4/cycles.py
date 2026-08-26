"""Checks for bipartite root-poset graphs of exceptional Chevalley type."""

from itertools import combinations

from algebra import jacobian_matrix, monomial_minor_certificate
from graphs import (bipartite_layers, contract_source_leaves, induced_bipartite,
                    root_poset, signature, vertices)
from roots import positive_roots, root_string, simple_roots, structure_constants


EXCEPTIONAL_TYPES = (("F", 4), ("E", 6), ("E", 7), ("E", 8))


def _target_subsets(targets):
    for size in range(len(targets) + 1):
        yield from combinations(targets, size)


def verify_type(cartan_type):
    roots = positive_roots(cartan_type)
    graph = root_poset(roots, simple_roots(roots))
    constants = structure_constants(cartan_type)
    layer_records = []
    failures = []
    cache = {}
    for layer_height, layer in bipartite_layers(graph).items():
        sources, targets = vertices(layer)
        tested = contracted = nontrivial = source_combinations = 0
        representatives = {}
        for target_subset in _target_subsets(targets):
            tested += 1
            induced = induced_bipartite(layer, sources, target_subset)
            core, _ = contract_source_leaves(induced)
            if not core:
                contracted += 1
                continue
            nontrivial += 1
            core_signature = signature(core)
            if core_signature not in cache:
                jacobian, core_sources, core_targets = jacobian_matrix(core, constants)
                certificate = monomial_minor_certificate(jacobian)
                cache[core_signature] = (certificate, core_sources, core_targets)
                source_combinations += certificate["source_combinations_tested"]
            certificate, core_sources, core_targets = cache[core_signature]
            representatives.setdefault(str(core_signature), {
                "sources": [root_string(root) for root in core_sources],
                "targets": [root_string(root) for root in core_targets],
                "selected_sources": (
                    [root_string(core_sources[index])
                     for index in certificate["columns"]]
                if certificate["columns"] is not None else None),
                **certificate,
            })
            if not certificate["ok"]:
                failures.append({
                    "height": layer_height,
                    "targets": [root_string(root) for root in target_subset],
                })
        layer_records.append({
            "height": layer_height,
            "source_count": len(sources),
            "target_count": len(targets),
            "tested": tested,
            "contracted": contracted,
            "nontrivial": nontrivial,
            "source_combinations_tested": source_combinations,
            "cores": list(representatives.values()),
        })
    return {
        "cartan_type": f"{cartan_type[0]}{cartan_type[1]}",
        "layers": layer_records,
        "failures": failures,
        "ok": not failures,
    }


def verify_exceptional_types():
    records = [verify_type(cartan_type) for cartan_type in EXCEPTIONAL_TYPES]
    return {"records": records, "ok": all(record["ok"] for record in records)}
