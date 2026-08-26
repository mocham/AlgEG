"""Deterministic rotational-deformation certificates for Configurations I--IV."""

from random import Random

from sage.all import factor, GF, matrix, PolynomialRing, QQ, ZZ

from borel import F4_TABLE
from roots import add, height, root_string, structure_constants, subtract


ROTATION_ROOTS = {"I": "0010", "II": "0010", "III": "1000", "IV": "0010"}


def _parse(root):
    return tuple(int(digit) for digit in root)


def _negate(root):
    return tuple(-digit for digit in root)


def root_chain(psi2, rotation_root, constants, alpha):
    """Coefficients of c_(alpha+n*delta) in the rotated c_alpha."""
    psi2 = set(psi2)
    negative_rotation = _negate(rotation_root)
    chain = [(alpha, QQ.one(), 0)]
    beta = alpha
    coefficient = QQ.one()
    step = 0
    while True:
        beta = add(beta, rotation_root)
        if beta not in psi2:
            break
        step += 1
        key = (negative_rotation, beta)
        if key not in constants:
            raise ValueError(f"missing rotation structure constant {key}")
        coefficient = coefficient * constants[key] / step
        chain.append((beta, coefficient, step))
    return tuple(chain)


def substitution_certificate(psi2, rotation_root, constants):
    return {
        root_string(alpha): [
            {"root": root_string(beta), "coefficient": str(coefficient), "mu_power": power}
            for beta, coefficient, power in root_chain(psi2, rotation_root, constants, alpha)
        ]
        for alpha in sorted(psi2)
    }


def _multiplicity(simple_root, root):
    nonzero = [index for index, coefficient in enumerate(simple_root) if coefficient]
    if len(nonzero) != 1 or simple_root[nonzero[0]] != 1:
        raise ValueError(f"expected an absolute simple root, got {simple_root}")
    return root[nonzero[0]]


def _rotation_system(record, rotation_root, specialization, constants, prime=11,
                     base_ring=None):
    ring, equations, determinant = _rotation_polynomials(
        record, rotation_root, specialization, constants, prime, base_ring)
    ideal = ring.ideal(equations)
    saturation, _ = ideal.saturation(ring.ideal(determinant))
    return ring, equations, ideal, determinant, saturation


def _rotation_polynomials(record, rotation_root, specialization, constants,
                          prime=11, base_ring=None):
    field = base_ring if base_ring is not None else GF(prime)
    psi0 = tuple(_parse(root) for root in record["Psi0"])
    psi2 = tuple(_parse(root) for root in record["Psi2"])
    rotation_root = _parse(rotation_root)
    variable_roots = tuple(root for root in psi0 if height(root) == 1)
    names = (["mu"] + [f"lambda{root_string(root)}" for root in variable_roots]
             + [f"bu{root_string(root)}" for root in psi0]
             + [f"br{root_string(root)}" for root in psi0])
    ring = PolynomialRing(field, names, order="degrevlex")
    variables = {name: ring(name) for name in names}
    mu = variables["mu"]

    def rotated(prefix, alpha):
        value = ring.zero()
        for beta, coefficient, power in root_chain(psi2, rotation_root, constants, alpha):
            value += field(coefficient) * mu**power * field(specialization[f"{prefix}{root_string(beta)}"])
        return value

    equations = [ring.one() - sum(variables[f"lambda{root_string(root)}"]
                                  for root in variable_roots)]
    for alpha in psi2:
        if height(alpha) < 2:
            continue
        lhs = sum(variables[f"lambda{root_string(delta)}"]
                  * _multiplicity(delta, alpha) * rotated("at", alpha)
                  for delta in variable_roots)
        rhs = ring.zero()
        for delta in psi0:
            beta = subtract(alpha, delta)
            if beta not in psi2:
                continue
            coefficient = constants.get((beta, delta))
            if coefficient is None:
                raise ValueError(f"missing equation structure constant {(beta, delta)}")
            if height(beta) == 1:
                coefficient *= 2
            rhs += field(coefficient) * (
                rotated("at", beta) * variables[f"bu{root_string(delta)}"]
                + rotated("ap", beta) * variables[f"br{root_string(delta)}"]
            )
        equations.append(rhs - lhs)
    jacobian = matrix(ring, [[equation.derivative(variable) for variable in ring.gens()]
                              for equation in equations])
    if jacobian.nrows() != jacobian.ncols():
        raise ValueError(f"non-square Jacobian {jacobian.nrows()}x{jacobian.ncols()}")
    determinant = jacobian.det()
    return ring, equations, determinant


def find_certificate(record, rotation_root, prime=11, attempts=200, seed=20260801):
    constants = structure_constants(("F", 4))
    psi2 = tuple(_parse(root) for root in record["Psi2"])
    random = Random(seed + int(record["id"].encode().hex(), 16) % 100000)
    for attempt in range(1, attempts + 1):
        specialization = {}
        for root in sorted(psi2):
            specialization[f"at{root_string(root)}"] = random.randrange(1, prime)
            specialization[f"ap{root_string(root)}"] = random.randrange(1, prime)
        ring, equations, ideal, determinant, saturation = _rotation_system(
            record, rotation_root, specialization, constants, prime)
        if ideal.is_one() or saturation.is_one():
            continue
        basis = ideal.groebner_basis()
        remainder = determinant.reduce(basis)
        return {
            "id": record["id"],
            "rotation_root": rotation_root,
            "prime": prime,
            "attempt": attempt,
            "variable_order": [str(variable) for variable in ring.gens()],
            "specialization": specialization,
            "equations": [str(equation) for equation in equations],
            "groebner_basis": [str(element) for element in basis],
            "jacobian_determinant": str(determinant),
            "determinant_remainder": str(remainder),
            "saturation_dimension": int(saturation.dimension()),
            "ok": True,
        }
    return {
        "id": record["id"],
        "rotation_root": rotation_root,
        "prime": prime,
        "attempts": attempts,
        "ok": False,
    }


def _bad_primes(values):
    primes = set()
    for value in values:
        rational = QQ(value)
        for integer in (rational.numerator(), rational.denominator()):
            if integer == 0:
                continue
            for prime, _ in factor(abs(ZZ(integer))):
                primes.add(int(prime))
    return sorted(primes)


def find_characteristic_uniform_certificate(record, rotation_root, attempts=200,
                                              seed=20260801):
    """Produce an exact point with invertible Jacobian over Z[1/12!]."""
    constants = structure_constants(("F", 4))
    psi2 = tuple(_parse(root) for root in record["Psi2"])
    psi0 = tuple(_parse(root) for root in record["Psi0"])
    parameter_names = tuple(
        f"{prefix}{root_string(root)}"
        for root in sorted(psi2) for prefix in ("at", "ap"))
    random = Random(seed + int(record["id"].encode().hex(), 16) % 100000)
    for attempt in range(1, attempts + 1):
        zero_parameters = {name: 0 for name in parameter_names}
        ring, zero_equations, _ = _rotation_polynomials(
            record, rotation_root, zero_parameters, constants, base_ring=QQ)
        point = {"mu": QQ(random.randrange(-2, 3))}
        variable_roots = tuple(root for root in psi0 if height(root) == 1)
        for index, root in enumerate(variable_roots):
            point[f"lambda{root_string(root)}"] = QQ(1 if index == 0 else 0)
        for root in psi0:
            point[f"bu{root_string(root)}"] = QQ(random.randrange(-2, 3))
            point[f"br{root_string(root)}"] = QQ(random.randrange(-2, 3))
        substitution = {ring(name): value for name, value in point.items()}
        if zero_equations[0].subs(substitution) != 0:
            raise AssertionError("normalization point is invalid")

        columns = []
        for parameter_name in parameter_names:
            parameters = dict(zero_parameters)
            parameters[parameter_name] = 1
            _, equations, _ = _rotation_polynomials(
                record, rotation_root, parameters, constants, base_ring=QQ)
            columns.append([equation.subs(substitution) for equation in equations[1:]])
        coefficient_matrix = matrix(QQ, len(zero_equations) - 1,
                                    len(parameter_names),
                                    lambda i, j: columns[j][i])
        kernel = coefficient_matrix.right_kernel().basis()
        if not kernel:
            continue
        for _ in range(200):
            parameter_vector = sum(
                (QQ(random.randrange(-3, 4)) * vector for vector in kernel),
                kernel[0].parent().zero())
            if not parameter_vector or any(value == 0 for value in parameter_vector):
                continue
            specialization = dict(zip(parameter_names, parameter_vector))
            ring, equations, determinant = _rotation_polynomials(
                record, rotation_root, specialization, constants, base_ring=QQ)
            substitution = {ring(name): value for name, value in point.items()}
            values = [equation.subs(substitution) for equation in equations]
            determinant_value = determinant.subs(substitution)
            if any(values) or determinant_value == 0:
                continue
            bad_primes = _bad_primes(
                list(point.values()) + list(parameter_vector) + [determinant_value])
            if any(prime > 12 for prime in bad_primes):
                continue
            return {
                "id": record["id"],
                "rotation_root": rotation_root,
                "attempt": attempt,
                "base_ring": "QQ",
                "variable_order": [str(variable) for variable in ring.gens()],
                "point": {name: str(value) for name, value in point.items()},
                "specialization": {name: str(value)
                                   for name, value in specialization.items()},
                "equation_values": [str(value) for value in values],
                "jacobian_determinant_at_point": str(determinant_value),
                "spreading_bad_primes": bad_primes,
                "uniform_for_p_gt_12": True,
                "ok": True,
            }
    return {
        "id": record["id"],
        "rotation_root": rotation_root,
        "attempts": attempts,
        "base_ring": "QQ",
        "ok": False,
    }


def verify_configurations():
    constants = structure_constants(("F", 4))
    records = []
    for record in F4_TABLE:
        rotation_root = ROTATION_ROOTS[record["id"]]
        certificate = find_certificate(record, rotation_root)
        certificate["substitutions"] = substitution_certificate(
            tuple(_parse(root) for root in record["Psi2"]),
            _parse(rotation_root), constants)
        records.append(certificate)
    return {"records": records, "ok": all(record["ok"] for record in records)}


def verify_characteristic_uniform_configurations():
    constants = structure_constants(("F", 4))
    records = []
    for record in F4_TABLE:
        rotation_root = ROTATION_ROOTS[record["id"]]
        certificate = find_characteristic_uniform_certificate(record, rotation_root)
        certificate["substitutions"] = substitution_certificate(
            tuple(_parse(root) for root in record["Psi2"]),
            _parse(rotation_root), constants)
        records.append(certificate)
    return {
        "records": records,
        "ok": all(record.get("ok") and record.get("uniform_for_p_gt_12")
                  for record in records),
    }
