"""Polynomial-matrix checks."""

from itertools import combinations

from sage.all import factor, gcd, matrix, PolynomialRing, QQ, ZZ

from roots import root_string, subtract


def jacobian_matrix(graph, constants, base_ring=QQ):
    sources = tuple(sorted(graph))
    targets = tuple(sorted(set().union(*(set(v) for v in graph.values()))))
    deltas = tuple(sorted({subtract(target, source)
                            for source in sources for target in graph[source]}))
    ring = PolynomialRing(base_ring, [f"x{root_string(delta)}" for delta in deltas])
    variables = {delta: ring(f"x{root_string(delta)}") for delta in deltas}
    rows = []
    for target in targets:
        row = []
        for source in sources:
            if target not in graph[source]:
                row.append(ring.zero())
                continue
            delta = subtract(target, source)
            key = (source, delta)
            if key not in constants:
                raise ValueError(f"missing structure constant for {key}")
            row.append(constants[key] * variables[delta])
        rows.append(row)
    return matrix(ring, rows), sources, targets


def monomial_minor_certificate(jacobian):
    rows, columns = jacobian.nrows(), jacobian.ncols()
    if rows == 0:
        return {"ok": True, "columns": [], "determinant": "1",
                "source_combinations_tested": 0}
    if columns < rows:
        return {"ok": False, "columns": None, "determinant": None,
                "source_combinations_tested": 0}
    tested = 0
    for selected in combinations(range(columns), rows):
        tested += 1
        determinant = jacobian.matrix_from_columns(selected).det()
        if determinant and len(determinant.dict()) == 1:
            return {
                "ok": True,
                "columns": list(selected),
                "determinant": str(determinant),
                "source_combinations_tested": tested,
            }
    return {"ok": False, "columns": None, "determinant": None,
            "source_combinations_tested": tested}


def rational_rank_certificate(coefficient_matrix):
    """Certify a rational rank by an explicit nonzero integral minor.

    The returned primes are exactly the prime divisors of the witness
    determinant.  At every other prime the rank cannot decrease after reduction.
    """
    rational = coefficient_matrix.change_ring(QQ)
    rank = rational.rank()
    if rank == 0:
        return {
            "rank": 0,
            "rows": [],
            "columns": [],
            "determinant": 1,
            "bad_primes": [],
        }
    rows = list(rational.transpose().pivots())
    columns = list(rational.matrix_from_rows(rows).pivots())
    witness = rational.matrix_from_rows_and_columns(rows, columns)
    determinant = witness.det()
    if not determinant:
        raise AssertionError("pivot minor is unexpectedly singular")
    numerator = ZZ(determinant.numerator())
    denominator = ZZ(determinant.denominator())
    bad_primes = sorted({int(prime) for value in (numerator, denominator)
                         for prime, _ in factor(abs(value))})
    return {
        "rank": rank,
        "rows": rows,
        "columns": columns,
        "determinant": str(determinant),
        "bad_primes": bad_primes,
    }


def polynomial_rank_certificate(polynomial_matrix):
    """Compute rank over the fraction field and its integral witness."""
    fraction = polynomial_matrix.base_ring().fraction_field()
    generic = polynomial_matrix.change_ring(fraction)
    rank = generic.rank()
    if rank == 0:
        return {"rank": 0, "rows": [], "columns": [], "determinant": "1",
                "content": 1, "bad_primes": []}
    rows = list(generic.transpose().pivots())
    columns = list(generic.matrix_from_rows(rows).pivots())
    determinant = polynomial_matrix.matrix_from_rows_and_columns(rows, columns).det()
    if not determinant:
        raise AssertionError("generic pivot minor is unexpectedly zero")
    coefficients = [ZZ(coefficient) for coefficient in determinant.coefficients()]
    content = abs(gcd(coefficients))
    bad_primes = sorted(int(prime) for prime, _ in factor(content)) if content else []
    return {
        "rank": rank,
        "rows": rows,
        "columns": columns,
        "determinant": str(determinant),
        "content": int(content),
        "bad_primes": bad_primes,
    }
