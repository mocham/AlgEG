"""Command-line entry point for the self-contained V14 computations."""

import sys

# This must precede local imports so --check cannot create __pycache__ files.
if "--check" in sys.argv:
    sys.dont_write_bytecode = True

import argparse
import hashlib
import json
from pathlib import Path

from sage.env import SAGE_VERSION

from bch import verify_bch_formulas
from borel import verify_f4_table
from bounded_strata import verify_f4_bounded_enumeration
from cone import verify_f4_table_cones
from cycles import verify_exceptional_types
from rotation import verify_characteristic_uniform_configurations
from twisted_cycles import verify_twisted_records
from elliptic_weyl import compute_all_elliptic_weyl_orbits, compute_non_borel_strata
from type_d import verify_type_d


DATA_DIRECTORY = Path(__file__).resolve().parent / "DATA"


def render_result(name, payload):
    document = {
        "schema": "grob-computation-v14",
        "sage_version": SAGE_VERSION,
        "name": name,
        "payload": payload,
    }
    encoded = (json.dumps(document, sort_keys=True, indent=2) + "\n").encode()
    digest = hashlib.sha256(encoded).hexdigest()
    return document, encoded, digest


def process_result(name, payload, output_directory=DATA_DIRECTORY, check=False):
    """Write a result, or compare it without changing any file."""
    document, encoded, digest = render_result(name, payload)
    output = output_directory / f"{name}.json"
    sidecar = output_directory / f"{name}.sha256"
    sidecar_bytes = f"{digest}  {output.name}\n".encode("ascii")
    if check:
        if not output.is_file() or output.read_bytes() != encoded:
            raise SystemExit(f"result mismatch: {output}")
        if not sidecar.is_file() or sidecar.read_bytes() != sidecar_bytes:
            raise SystemExit(f"checksum mismatch: {sidecar}")
        print(f"{output}: verified sha256={digest}")
        return document

    output_directory.mkdir(parents=True, exist_ok=True)
    output.write_bytes(encoded)
    sidecar.write_bytes(sidecar_bytes)
    print(f"{output}: sha256={digest}")
    return document


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=(
        "bch", "borel", "cone", "cycles", "twisted-cycles", "rotation",
        "type-d", "f4-exhaustive", "all", "quick"))
    parser.add_argument(
        "--check", action="store_true",
        help="recompute and compare with existing files; do not write")
    parser.add_argument(
        "--output-directory", type=Path, default=DATA_DIRECTORY,
        help="directory to write, or to read when --check is used")
    args = parser.parse_args()

    def output(name, result):
        return process_result(
            name, result, output_directory=args.output_directory, check=args.check)

    ok = True
    if args.command in ("bch", "all", "quick"):
        result = verify_bch_formulas()
        output("bch-root-graded-v14", result)
        ok &= result["ok"]
    if args.command in ("borel", "all", "quick"):
        result = verify_f4_table()
        output("borel-f4-table-v14", result)
        ok &= result["ok"]
    if args.command in ("cone", "all", "quick"):
        result = verify_f4_table_cones()
        output("cone-f4-table-v14", result)
        ok &= result["ok"]
    if args.command in ("cycles", "all", "quick"):
        result = verify_exceptional_types()
        output("cycles-exceptional-v14", result)
        ok &= result["ok"]
    if args.command in ("rotation", "all", "quick"):
        result = verify_characteristic_uniform_configurations()
        output("rotation-f4-uniform-v14", result)
        ok &= result["ok"]
    if args.command in ("type-d", "all"):
        result = verify_type_d()
        output("cycles-type-d-v14", result)
        ok &= result["ok"]
    if args.command in ("f4-exhaustive", "all"):
        result = verify_f4_bounded_enumeration(progress=10000)
        output("fine-strata-f4-exhaustive-v14", result)
        ok &= result["ok"]
    if args.command in ("twisted-cycles", "all"):
        print("Generating elliptic Weyl orbits (this may take a while)...")
        orbits = compute_all_elliptic_weyl_orbits()
        output("elliptic-weyl-orbits-v14", {
            "types": sorted(orbits.keys()),
            "li_count_per_type": {k: len(v) for k, v in orbits.items()},
        })
        print("Computing non-Borel strata...")
        non_borel_data = compute_non_borel_strata(orbits)
        result = verify_twisted_records(non_borel_data)
        output("cycles-twisted-exceptional-v14", result)
        ok &= result["ok"]
    raise SystemExit(0 if ok else 1)


if __name__ == "__main__":
    main()
