"""Deterministic root-system helpers used by the current computations."""

from itertools import combinations

from sage.all import LieAlgebra, QQ, RootSystem


def root_tuple(root):
    return tuple(int(c) for c in root.to_vector())


def root_string(root):
    return "".join(str(c) for c in root)


def positive_roots(cartan_type):
    lattice = RootSystem(cartan_type).root_lattice()
    return tuple(sorted(root_tuple(root) for root in lattice.positive_roots()))


def positive_roots_set(cartan_type):
    return set(positive_roots(cartan_type))


def simple_roots(roots):
    roots = set(roots)
    return tuple(sorted(
        root for root in roots
        if not any(tuple(a + b for a, b in zip(left, right)) == root
                   for left in roots for right in roots)
    ))


def add(left, right):
    return tuple(a + b for a, b in zip(left, right))


def subtract(left, right):
    return tuple(a - b for a, b in zip(left, right))


def height(root):
    return sum(root)


def structure_constants(cartan_type):
    """Return N_(a,b) in [e_a,e_b]=N_(a,b)e_(a+b)."""
    lattice = RootSystem(cartan_type).root_lattice()
    lie_algebra = LieAlgebra(QQ, cartan_type=cartan_type)
    basis = lie_algebra.basis()
    constants = {}
    for left in lattice.roots():
        for right in lattice.roots():
            if left == -right:
                continue
            bracket = basis[left].bracket(basis[right])
            if bracket == 0:
                continue
            coefficient = bracket.monomial_coefficients().get(left + right, 0)
            if coefficient:
                constants[(root_tuple(left), root_tuple(right))] = int(coefficient)
    return constants


def root_less_than(left, right):
    return all(a <= b for a, b in zip(left, right))


def downward_closed_subsets(items, less_than):
    """Enumerate order ideals once, in deterministic order."""
    items = tuple(sorted(items))
    closures = []
    for item in items:
        closures.append(frozenset(
            other for other in items
            if other == item or less_than(other, item)
        ))
    ideals = set()
    for mask in range(1 << len(items)):
        ideal = frozenset().union(*(
            closures[index] for index in range(len(items)) if mask & (1 << index)
        ))
        ideals.add(ideal)
    return tuple(sorted(ideals, key=lambda ideal: (len(ideal), tuple(sorted(ideal)))))


def subsets(items):
    items = tuple(sorted(items))
    for size in range(len(items) + 1):
        for subset in combinations(items, size):
            yield subset
