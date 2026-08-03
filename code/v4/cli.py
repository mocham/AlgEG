"""Command-line entry point for self-contained v12 certificates."""

import argparse
import hashlib
import json
from pathlib import Path

from sage.env import SAGE_VERSION

from borel import verify_f4_table
from bounded_strata import verify_f4_bounded_enumeration
from cone import verify_f4_table_cones
from cycles import verify_exceptional_types
from rotation import verify_characteristic_uniform_configurations
from twisted_cycles import verify_twisted_records
from elliptic_weyl import compute_all_elliptic_weyl_orbits, compute_non_borel_strata


DATA_DIRECTORY = Path(__file__).resolve().parent / "DATA"


def write_certificate(name, payload):
    document = {
        "schema": "grob-computation-v12",
        "sage_version": SAGE_VERSION,
        "name": name,
        "payload": payload,
    }
    encoded = (json.dumps(document, sort_keys=True, indent=2) + "\n").encode()
    digest = hashlib.sha256(encoded).hexdigest()
    DATA_DIRECTORY.mkdir(exist_ok=True)
    output = DATA_DIRECTORY / f"{name}.json"
    output.write_bytes(encoded)
    (DATA_DIRECTORY / f"{name}.sha256").write_text(f"{digest}  {output.name}\n")
    print(f"{output}: sha256={digest}")
    return document


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=(
        "borel", "cone", "cycles", "twisted-cycles", "rotation",
        "f4-exhaustive", "all", "quick"))
    args = parser.parse_args()

    ok = True
    if args.command in ("borel", "all", "quick"):
        result = verify_f4_table()
        write_certificate("borel-f4-table-v12", result)
        ok &= result["ok"]
    if args.command in ("cone", "all", "quick"):
        result = verify_f4_table_cones()
        write_certificate("cone-f4-table-v12", result)
        ok &= result["ok"]
    if args.command in ("cycles", "all", "quick"):
        result = verify_exceptional_types()
        write_certificate("cycles-exceptional-v12", result)
        ok &= result["ok"]
    if args.command in ("rotation", "all", "quick"):
        result = verify_characteristic_uniform_configurations()
        write_certificate("rotation-f4-uniform-v12", result)
        ok &= result["ok"]
    if args.command in ("f4-exhaustive", "all"):
        result = verify_f4_bounded_enumeration(progress=10000)
        write_certificate("fine-strata-f4-exhaustive-v12", result)
        ok &= result["ok"]
    if args.command in ("twisted-cycles", "all"):
        print("Generating elliptic Weyl orbits (this may take a while)...")
        orbits = compute_all_elliptic_weyl_orbits()
        write_certificate("elliptic-weyl-orbits-v12", {
            "types": sorted(orbits.keys()),
            "li_count_per_type": {k: len(v) for k, v in orbits.items()},
        })
        print("Computing non-Borel strata...")
        non_borel_data = compute_non_borel_strata(orbits)
        result = verify_twisted_records(non_borel_data)
        write_certificate("cycles-twisted-exceptional-v12", result)
        ok &= result["ok"]
    raise SystemExit(0 if ok else 1)


if __name__ == "__main__":
    main()
