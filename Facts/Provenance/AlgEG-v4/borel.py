"""Deterministic F4 Borel-stratum helpers, table certificates, and full enumeration.

The bounded exhaustive enumeration is implemented independently in
``bounded_strata.py``.
"""

from graphs import bipartite_layers, root_poset
from roots import downward_closed_subsets, positive_roots, root_string, simple_roots


F4_TABLE = (
    {
        "id": "I",
        "K_complement": ("0010",),
        "Psi0": ("0010",),
        "Psi2": ("0001", "0011", "0100", "0110", "0120"),
    },
    {
        "id": "II",
        "K_complement": ("0010",),
        "Psi0": ("0010",),
        "Psi2": ("0001", "0011", "0100", "0110", "0120", "1000"),
    },
    {
        "id": "III",
        "K_complement": ("0010", "1000"),
        "Psi0": ("0010", "1000"),
        "Psi2": ("0001", "0011", "0100", "0110", "0120", "1100", "1110", "1120"),
    },
    {
        "id": "IV",
        "K_complement": ("0001", "0010", "0011", "1000"),
        "Psi0": ("0001", "0010", "0011", "1000"),
        "Psi2": ("0100", "0110", "0111", "0120", "0121", "0122",
                 "1100", "1110", "1111", "1120", "1121", "1122"),
    },
)


def _parse(root):
    return tuple(int(digit) for digit in root)


def root_less_than(left, right):
    return all(a <= b for a, b in zip(left, right))


def f4_complement_ideals():
    """Return the complement of Phi_K in Phi+, as downward-closed sets."""
    return downward_closed_subsets(positive_roots(("F", 4)), root_less_than)


def closure_axioms(psi0, psi2, roots):
    roots = set(roots)
    psi0, psi2 = set(psi0), set(psi2)
    first = all(tuple(a - b for a, b in zip(alpha, beta)) in psi0
                for alpha in psi2 for beta in psi2
                if tuple(a - b for a, b in zip(alpha, beta)) in roots)
    second = all(tuple(a + b for a, b in zip(alpha, beta)) in psi2
                 for alpha in psi2 for beta in psi0
                 if tuple(a + b for a, b in zip(alpha, beta)) in roots)
    return first and second


def verify_f4_table():
    roots = positive_roots(("F", 4))
    root_set = set(roots)
    complement_ideals = set(f4_complement_ideals())
    records = []
    failures = []
    for record in F4_TABLE:
        complement = frozenset(_parse(root) for root in record["K_complement"])
        phi_k = root_set - complement
        delta_k = set(simple_roots(phi_k))
        psi0 = tuple(_parse(root) for root in record["Psi0"])
        psi2 = tuple(_parse(root) for root in record["Psi2"])
        checks = {
            "complement_is_order_ideal": complement in complement_ideals,
            "fine_closure": closure_axioms(psi0, psi2, roots),
            "psi2_is_K_simple": set(psi2).issubset(delta_k),
        }
        if not all(checks.values()):
            failures.append({"id": record["id"], "checks": checks})
        records.append({
            **record,
            "DeltaK": [root_string(root) for root in sorted(delta_k)],
            "checks": checks,
        })
    return {
        "positive_root_count": len(roots),
        "complement_ideal_count": len(complement_ideals),
        "records": records,
        "failures": failures,
        "ok": not failures,
    }
