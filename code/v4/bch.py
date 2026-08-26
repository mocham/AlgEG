"""Exact BCH and root-graded checks for the Herr obstruction."""

from dataclasses import dataclass
import hashlib
import json

from sage.all import (
    bernoulli,
    factorial,
    LieAlgebra,
    matrix,
    PolynomialRing,
    QQ,
    RootSystem,
)


@dataclass(frozen=True)
class RootLieData:
    algebra: object
    roots: tuple
    basis: dict
    constants: dict


def _positive_compositions(total, length):
    """Yield positive compositions in deterministic lexicographic order."""
    if length == 1:
        if total >= 1:
            yield (total,)
        return
    for first in range(1, total - length + 2):
        for rest in _positive_compositions(total - first, length - 1):
            yield (first, *rest)


def bch_terms(X, Y, max_class=None):
    r"""Return homogeneous terms of ``log(exp(X) exp(Y))``.

    If ``Z_n`` has bracket length ``n``, the recursion is

    .. math::

       Z_1=X+Y,

       m Z_m={1\over2}[X-Y,Z_{m-1}]
       +\sum_{1\le p\le (m-1)/2}{B_{2p}\over(2p)!}
        \sum_{k_1+\cdots+k_{2p}=m-1}
        [Z_{k_1},[\cdots,[Z_{k_{2p}},X+Y]]\cdots]].

    Every ``k_i`` is positive.  In positive characteristic, the relevant
    denominators must be invertible.
    """
    parent = X.parent()
    if Y.parent() != parent:
        Y = parent(Y)
    if max_class is None:
        try:
            max_class = int(parent.step())
        except (AttributeError, TypeError, ValueError) as exc:
            raise ValueError(
                "max_class is required for a non-nilpotent parent") from exc
    max_class = int(max_class)
    if max_class < 0:
        raise ValueError("max_class must be nonnegative")
    if max_class == 0:
        return ()

    terms = [parent.zero(), X + Y]
    for degree in range(2, max_class + 1):
        term = QQ(1) / (2 * degree) * (X - Y).bracket(terms[degree - 1])
        for p in range(1, (degree - 1) // 2 + 1):
            coefficient = bernoulli(2 * p) / (factorial(2 * p) * degree)
            for parts in _positive_compositions(degree - 1, 2 * p):
                nested = terms[1]
                for part in parts:
                    nested = terms[part].bracket(nested)
                term += coefficient * nested
        terms.append(term)
    return tuple(terms[1:])


def bch(X, Y, max_class=None):
    """Return the BCH sum through ``max_class``."""
    terms = bch_terms(X, Y, max_class=max_class)
    return sum(terms, X.parent().zero())


def _multiply_word_polynomials(left, right, max_class):
    result = {}
    for left_word, left_coefficient in left.items():
        for right_word, right_coefficient in right.items():
            word = left_word + right_word
            if len(word) > max_class:
                continue
            coefficient = (
                result.get(word, QQ.zero())
                + left_coefficient * right_coefficient)
            if coefficient:
                result[word] = coefficient
            elif word in result:
                del result[word]
    return result


def truncated_log_exp_product(max_class):
    r"""Return the word coefficients of
    ``log(exp(tX) exp(tY)) mod t^(max_class+1)``.

    A word is a tuple of zeros and ones, denoting ``X`` and ``Y``.  This uses
    multiplication in the free associative algebra and does not call the BCH
    recursion.
    """
    max_class = int(max_class)
    if max_class < 0:
        raise ValueError("max_class must be nonnegative")
    if max_class == 0:
        return ()

    product_minus_one = {}
    for x_power in range(max_class + 1):
        for y_power in range(max_class - x_power + 1):
            degree = x_power + y_power
            if degree == 0:
                continue
            word = (0,) * x_power + (1,) * y_power
            product_minus_one[word] = (
                QQ.one() / (factorial(x_power) * factorial(y_power)))

    logarithm = {}
    power = product_minus_one
    for factor_count in range(1, max_class + 1):
        coefficient = QQ((-1) ** (factor_count + 1)) / factor_count
        for word, value in power.items():
            total = logarithm.get(word, QQ.zero()) + coefficient * value
            if total:
                logarithm[word] = total
            elif word in logarithm:
                del logarithm[word]
        power = _multiply_word_polynomials(
            power, product_minus_one, max_class)

    return tuple({
        word: logarithm[word]
        for word in sorted(logarithm)
        if len(word) == degree
    } for degree in range(1, max_class + 1))


def verify_bch_against_truncated_log(max_class=11):
    """Compare every BCH class with an independent associative expansion."""
    max_class = int(max_class)
    free = LieAlgebra(QQ, names=("X", "Y"))
    X, Y = free.lie_algebra_generators()
    terms = bch_terms(X, Y, max_class=max_class)
    expected_terms = truncated_log_exp_product(max_class)
    lifted_X, lifted_Y = X.lift(), Y.lift()
    parent = lifted_X.parent()
    monomials = {(): parent.one()}

    def monomial(word):
        if word not in monomials:
            monomials[word] = monomial(word[:-1]) * (
                lifted_X if word[-1] == 0 else lifted_Y)
        return monomials[word]

    comparisons = []
    rows = []
    for degree, (actual, expected) in enumerate(
            zip(terms, expected_terms, strict=True), start=1):
        expected_value = sum(
            (coefficient * monomial(word)
             for word, coefficient in expected.items()),
            parent.zero())
        class_ok = actual.lift() == expected_value
        comparisons.append({
            "class": degree,
            "nonzero_word_count": len(expected),
            "ok": class_ok,
        })
        rows.extend([
            degree,
            "".join("X" if letter == 0 else "Y" for letter in word),
            str(coefficient),
        ] for word, coefficient in expected.items())

    encoded = json.dumps(rows, separators=(",", ":")).encode()
    return {
        "maximum_class": max_class,
        "method": "log(exp(tX)exp(tY)) mod t^(maximum_class+1)",
        "word_coefficients_sha256": hashlib.sha256(encoded).hexdigest(),
        "classes": comparisons,
        "ok": all(record["ok"] for record in comparisons),
    }


def herr_bch_obstruction(
    f_lower,
    g_lower,
    phi_f_g_lower,
    gamma_g_f_lower,
    max_class=None,
):
    r"""Return the lower-root obstruction with its required sign.

    Here ``phi_f_g_lower`` is ``Ad(f_ss)(phi(g_lower))`` and
    ``gamma_g_f_lower`` is ``Ad(g_ss)(gamma(f_lower))``.  The equation is

    ``(f_alpha-gamma_g(f_alpha))-(g_alpha-phi_f(g_alpha)) = D_alpha``.
    """
    right = bch(g_lower, gamma_g_f_lower, max_class=max_class)
    left = bch(f_lower, phi_f_g_lower, max_class=max_class)
    return right - left


def root_diagonal_action(element, basis, eigenvalues, coefficient_map=None):
    r"""Apply a diagonal torus action, optionally after a coefficient map."""
    coefficient_map = coefficient_map or (lambda value: value)
    result = element.parent().zero()
    for root, basis_element in basis.items():
        support = tuple(basis_element.monomial_coefficients())
        if len(support) != 1:
            raise ValueError("basis entries must be root vectors")
        coefficient = element[support[0]]
        if coefficient:
            result += eigenvalues[root] * coefficient_map(coefficient) * basis_element
    return result


def _root_tuple(root):
    return tuple(int(coefficient) for coefficient in root.to_vector())


def _root_name(root):
    return "r_" + "_".join(str(coefficient) for coefficient in root)


def chevalley_positive_lie_algebra(cartan_type, base_ring=QQ):
    """Build the positive Chevalley Lie algebra from its root system."""
    lattice = RootSystem(cartan_type).root_lattice()
    lattice_roots = {_root_tuple(root): root for root in lattice.positive_roots()}
    roots = tuple(sorted(lattice_roots))
    chevalley = LieAlgebra(QQ, cartan_type=cartan_type)
    chevalley_basis = chevalley.basis()
    constants = {}
    for left in roots:
        for right in roots:
            result_root = tuple(a + b for a, b in zip(left, right))
            if result_root not in lattice_roots:
                continue
            bracket = chevalley_basis[lattice_roots[left]].bracket(
                chevalley_basis[lattice_roots[right]])
            coefficient = bracket.monomial_coefficients().get(
                lattice_roots[result_root], 0)
            if coefficient:
                constants[(left, right)] = int(coefficient)

    names = {_root_name(root): root for root in roots}
    name_by_root = {root: name for name, root in names.items()}
    structure = {}
    for left_index, left in enumerate(roots):
        for right in roots[left_index + 1:]:
            coefficient = constants.get((left, right))
            if coefficient:
                result_root = tuple(a + b for a, b in zip(left, right))
                structure[(name_by_root[left], name_by_root[right])] = {
                    name_by_root[result_root]: coefficient
                }
    algebra = LieAlgebra(
        base_ring,
        structure,
        names=tuple(name_by_root[root] for root in roots),
        nilpotent=True,
    )
    algebra_basis = algebra.basis()
    basis = {root: algebra_basis[name_by_root[root]] for root in roots}
    return RootLieData(algebra, roots, basis, constants)


def _root_string(root):
    return "".join(str(coefficient) for coefficient in root)


def oriented_jacobian(graph, constants, base_ring=QQ):
    """Compute the leading root-height Jacobian from structure constants."""
    sources = tuple(sorted(graph))
    targets = tuple(sorted({value for values in graph.values() for value in values}))
    deltas = tuple(sorted({
        tuple(target[index] - source[index] for index in range(len(source)))
        for source in sources for target in graph[source]
    }))
    ring = PolynomialRing(
        base_ring, [f"x{_root_string(delta)}" for delta in deltas])
    variables = {delta: ring(f"x{_root_string(delta)}") for delta in deltas}
    rows = []
    for target in targets:
        row = []
        for source in sources:
            if target not in graph[source]:
                row.append(ring.zero())
                continue
            delta = tuple(
                target[index] - source[index] for index in range(len(source)))
            row.append(constants[(source, delta)] * variables[delta])
        rows.append(row)
    return matrix(ring, rows), sources, targets


def _b2(base_ring):
    structure = {("a", "b"): {"c": 1}, ("a", "c"): {"d": 2}}
    algebra = LieAlgebra(
        base_ring, structure, names=("a", "b", "c", "d"), nilpotent=True)
    a, b, c, d = algebra.lie_algebra_generators()
    return algebra, {"a": a, "b": b, "c": c, "d": d}


def verify_bch_formulas():
    """Check the homogeneous recursion and the A2, B2, and F4 cases."""
    free = LieAlgebra(QQ, 2, step=4)
    X, Y = free.basis(1)
    terms = bch_terms(X, Y, max_class=4)
    homogeneous_ok = (
        terms[0] == X + Y
        and terms[1] == QQ(1) / 2 * X.bracket(Y)
        and terms[2] == (
            QQ(1) / 12 * X.bracket(X.bracket(Y))
            + QQ(1) / 12 * X.bracket(Y).bracket(Y)
        )
        and terms[3] == QQ(1) / 24 * X.bracket(X.bracket(Y).bracket(Y))
    )

    a2_ring = PolynomialRing(
        QQ, ("fa", "fb", "ga", "gb", "pa", "pb", "qa", "qb"))
    fa, fb, ga, gb, pa, pb, qa, qb = a2_ring.gens()
    a2 = LieAlgebra(
        a2_ring, {("a", "b"): {"c": 1}},
        names=("a", "b", "c"), nilpotent=True)
    a, b, _ = a2.lie_algebra_generators()
    a2_obstruction = herr_bch_obstruction(
        fa * a + fb * b,
        ga * a + gb * b,
        pa * a + pb * b,
        qa * a + qb * b,
        max_class=2,
    )
    a2_expected = QQ(1) / 2 * (ga * qb - gb * qa - fa * pb + fb * pa)
    a2_ok = a2_obstruction["c"] == a2_expected

    names = (
        "fa", "fb", "fc", "ga", "gb", "gc",
        "pa", "pb", "pc", "qa", "qb", "qc",
    )
    b2_ring = PolynomialRing(QQ, names)
    fa, fb, fc, ga, gb, gc, pa, pb, pc, qa, qb, qc = b2_ring.gens()
    _, basis = _b2(b2_ring)
    a, b, c = basis["a"], basis["b"], basis["c"]
    b2_obstruction = herr_bch_obstruction(
        fa * a + fb * b + fc * c,
        ga * a + gb * b + gc * c,
        pa * a + pb * b + pc * c,
        qa * a + qb * b + qc * c,
        max_class=3,
    )
    right = ga * qc - gc * qa + QQ(1) / 6 * (ga - qa) * (ga * qb - gb * qa)
    left = fa * pc - fc * pa + QQ(1) / 6 * (fa - pa) * (fa * pb - fb * pa)
    b2_ok = b2_obstruction["d"] == right - left

    independent_check = verify_bch_against_truncated_log(max_class=11)
    f4 = chevalley_positive_lie_algebra(("F", 4))
    f4_ok = (
        max(sum(root) for root in f4.roots) == 11
        and int(f4.algebra.step()) == 11
        and independent_check["ok"]
    )

    return {
        "homogeneous_terms": {
            "formulas": [
                "Z_1=X+Y",
                "Z_2=1/2[X,Y]",
                "Z_3=1/12[X,[X,Y]]+1/12[[X,Y],Y]",
                "Z_4=1/24[X,[[X,Y],Y]]",
            ],
            "ok": homogeneous_ok,
        },
        "A2": {
            "root_c_formula": "1/2*(ga*qb-gb*qa-fa*pb+fb*pa)",
            "ok": a2_ok,
        },
        "B2": {
            "root_d_right_formula": (
                "ga*qc-gc*qa+1/6*(ga-qa)*(ga*qb-gb*qa)"),
            "root_d_left_formula": (
                "fa*pc-fc*pa+1/6*(fa-pa)*(fa*pb-fb*pa)"),
            "ok": b2_ok,
        },
        "F4": {
            "maximum_positive_root_height": 11,
            "bch_classes_checked": independent_check["maximum_class"],
            "ok": f4_ok,
        },
        "free_lie_comparison": independent_check,
        "ok": homogeneous_ok and a2_ok and b2_ok and f4_ok,
    }
