"""Deterministic F4 Borel-stratum helpers and table checks.

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


def f4_borel_component_count(prime):
    """Return configuration-labelled records and tuple counts for split F4."""
    if prime <= 12:
        raise ValueError("the computation assumes p > 12")
    row_i = [
        {"configuration": "I", "a": a, "label": [a, 1, 0, 1]}
        for a in range(prime - 1)
    ]
    other_rows = [
        {"configuration": "II", "label": [1, 1, 0, 1]},
        {"configuration": "III", "label": [0, 1, 0, 1]},
        {"configuration": "IV", "label": [0, 1, 0, 0]},
    ]
    records = row_i + other_rows
    records_by_tuple = {}
    for record in records:
        label = tuple(record["label"])
        records_by_tuple.setdefault(label, []).append({
            key: value for key, value in record.items() if key != "label"
        })
    repeated = [
        {"label": list(label), "records": records_by_tuple[label]}
        for label in sorted(records_by_tuple)
        if len(records_by_tuple[label]) > 1
    ]
    record_count = len(records)
    tuple_count = len(records_by_tuple)
    return {
        "prime": prime,
        "row_I_formula": "(a,1,0,1), a in Z/(p-1)Z",
        "row_I_count_formula": "p-1",
        "row_II_formula": "(1,1,0,1)",
        "row_III_formula": "(0,1,0,1)",
        "row_IV_formula": "(0,1,0,0)",
        "configuration_labelled_record_count_formula": "(p-1)+3=p+2",
        "distinct_label_tuple_count_formula": "(p-1)+1=p",
        "row_I_labels": row_i,
        "row_II_to_IV_labels": other_rows,
        "configuration_labelled_record_count": record_count,
        "distinct_label_tuple_count": tuple_count,
        "repeated_label_tuples": repeated,
        "ok": record_count == prime + 2 and tuple_count == prime,
    }


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
    component_count = f4_borel_component_count(17)
    return {
        "positive_root_count": len(roots),
        "complement_ideal_count": len(complement_ideals),
        "additional_borel_components": component_count,
        "records": records,
        "failures": failures,
        "ok": (
            not failures
            and component_count["ok"]
            and component_count["configuration_labelled_record_count"] == 19
            and component_count["distinct_label_tuple_count"] == 17
        ),
    }
