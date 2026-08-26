"""Small deterministic bipartite-graph operations."""

from collections import defaultdict

from roots import add, height


def root_poset(roots, basic_roots):
    roots = set(roots)
    return {
        root: tuple(sorted(candidate for delta in basic_roots
                           if (candidate := add(root, delta)) in roots))
        for root in sorted(roots)
    }


def bipartite_layers(graph):
    by_height = defaultdict(list)
    for root in graph:
        by_height[height(root)].append(root)
    layers = {}
    for current_height in sorted(by_height):
        sources = tuple(sorted(by_height[current_height]))
        targets = set(by_height.get(current_height + 1, ()))
        if not targets:
            continue
        layers[current_height] = {
            source: tuple(target for target in graph[source] if target in targets)
            for source in sources
        }
    return layers


def induced_bipartite(graph, source_subset, target_subset):
    sources = set(source_subset)
    targets = set(target_subset)
    return {
        source: tuple(target for target in graph.get(source, ()) if target in targets)
        for source in sorted(sources)
    }


def contract_source_leaves(graph):
    """Remove degree-one graph sources and record each removed edge."""
    adjacency = {source: set(targets) for source, targets in graph.items()}
    trace = []
    while True:
        adjacency = {source: targets for source, targets in adjacency.items() if targets}
        leaves = [(source, next(iter(targets)))
                  for source, targets in sorted(adjacency.items()) if len(targets) == 1]
        if not leaves:
            break
        source, target = leaves[0]
        trace.append((source, target))
        adjacency.pop(source)
        for other in adjacency:
            adjacency[other].discard(target)
    core = {source: tuple(sorted(targets)) for source, targets in sorted(adjacency.items())}
    return core, tuple(trace)


def compactify_graph(graph):
    return {source: targets for source, targets in graph.items() if targets}


def vertices(graph):
    sources = tuple(sorted(graph))
    targets = tuple(sorted(set().union(*(set(v) for v in graph.values())))) if graph else ()
    return sources, targets


def signature(graph):
    return tuple((source, tuple(targets)) for source, targets in sorted(graph.items()))


def height_dict(graph):
    """Assign each graph vertex a height, increasing along each edge."""
    layers = bipartite_layers(graph)
    if not layers:
        return {}
    ht = {}
    for layer_height, layer in sorted(layers.items()):
        for source in layer:
            ht[source] = layer_height
        seen_targets = set(ht.values())
        for target_set in layer.values():
            for target in target_set:
                if target not in seen_targets:
                    ht[target] = layer_height + 1
    return ht
