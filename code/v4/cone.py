"""Linear cone equations after fixing simple-root coordinates.

Supports both the four F4 table configurations and full
fine strata enumeration.
"""

from collections import defaultdict

from sage.all import matrix, QQ

from algebra import rational_rank_certificate
from borel import F4_TABLE
from graphs import root_poset
from roots import positive_roots, root_string, simple_roots, structure_constants, subtract


def grouped_linear_relations(layers, psi0, psi2, constants):
    """Return one coefficient row per target, never one generator per edge.

    The simple-root t-coordinate is specialized to 1.  Its s-coordinate is
    specialized to 1 when that coordinate occurs.  Remaining variables are
    t_beta and, when present, s_beta.  This is the linear fiber used by the
    dimension comparison in the manuscript.
    """
    psi_with_s = set(psi0) | set(psi2)
    rows = defaultdict(lambda: defaultdict(QQ))
    variables = set()
    for layer in layers.values():
        for beta, targets in layer.items():
            t_beta = f"t{root_string(beta)}"
            variables.add(t_beta)
            if beta in psi_with_s:
                variables.add(f"s{root_string(beta)}")
            for alpha in targets:
                delta = subtract(alpha, beta)
                coefficient = constants.get((delta, beta))
                if coefficient is None:
                    raise ValueError(f"missing structure constant for {(delta, beta)}")
                rows[alpha][t_beta] += coefficient
                if delta in psi_with_s and beta in psi_with_s:
                    rows[alpha][f"s{root_string(beta)}"] += coefficient
    variable_order = tuple(sorted(variables))
    target_order = tuple(sorted(rows))
    coefficient_matrix = matrix(QQ, [
        [rows[target].get(variable, 0) for variable in variable_order]
        for target in target_order
    ])
    return coefficient_matrix, target_order, variable_order


def linear_cone_certificate(layers, psi0, psi2, constants):
    coefficient_matrix, targets, variables = grouped_linear_relations(
        layers, psi0, psi2, constants)
    rank_witness = rational_rank_certificate(coefficient_matrix)
    rank = rank_witness["rank"]
    relation_count = coefficient_matrix.nrows()
    return {
        "targets": [root_string(root) for root in targets],
        "variables": list(variables),
        "matrix": [[int(entry) for entry in row] for row in coefficient_matrix.rows()],
        "relation_count": relation_count,
        "rank": rank,
        "deficiency": relation_count - rank,
        "dimension": len(variables) - rank,
        "rank_witness": rank_witness,
        "uniform_for_p_gt_12": all(prime <= 12 for prime in rank_witness["bad_primes"]),
    }


def verify_f4_table_cones():
    """Certify the grouped linear cone systems for Configurations I--IV."""
    roots = set(positive_roots(("F", 4)))
    constants = structure_constants(("F", 4))
    records = []
    for table_record in F4_TABLE:
        complement = {tuple(int(digit) for digit in root)
                      for root in table_record["K_complement"]}
        phi_k = roots - complement
        graph = root_poset(phi_k, simple_roots(phi_k))
        psi0 = tuple(tuple(int(digit) for digit in root)
                     for root in table_record["Psi0"])
        psi2 = tuple(tuple(int(digit) for digit in root)
                     for root in table_record["Psi2"])
        records.append({
            "id": table_record["id"],
            **linear_cone_certificate({0: graph}, psi0, psi2, constants),
        })
    return {
        "records": records,
        "ok": all(record["uniform_for_p_gt_12"] for record in records),
    }


def verify_fine_strata_cones(strata_result, cartan_type=("F", 4)):
    """Verify cone models for all fine strata enumerated by strata module."""
    if "type_D" not in strata_result:
        return {"ok": False, "error": "invalid strata result format"}
    constants = structure_constants(cartan_type)
    roots = set(positive_roots(cartan_type))

    def _parse_roots(strs):
        return tuple(tuple(int(d) for d in r) for r in strs)

    all_ok = True
    cone_records = []
    for category, items in [("D", strata_result["type_D"]),
                              ("U", strata_result["type_U"]),
                              ("O", strata_result["type_O"])]:
        for item in items:
            K_strs = item["K"]
            complement = {_parse(r) for r in K_strs}
            phi_k = roots - complement
            if not phi_k:
                continue
            psi0 = _parse_roots(item["Psi0"])
            psi2 = _parse_roots(item["Psi2"])
            graph = root_poset(phi_k, simple_roots(phi_k))
            cert = linear_cone_certificate({0: graph}, psi0, psi2, constants)
            cert["category"] = category
            cert["K"] = K_strs
            cone_records.append(cert)
            if cert["dimension"] > 24 * 25:  # safety: huge dimension = likely bug
                all_ok = False

    return {"records": cone_records, "ok": all_ok}


def _parse(root):
    return tuple(int(digit) for digit in root)
