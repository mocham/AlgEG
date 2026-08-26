# Algorithmic Emerton-Gee Stacks

This repository computes irreducible components of reduced Emerton-Gee stacks.

# Computational Checks

The reviewed implementation is
[`code/v4/`](https://github.com/mocham/AlgEG/tree/grob-v14/code/v4).
This file records the algorithms, commands, output formats, checksums, and rough
runtimes. The paper contains only the mathematical pseudocode and numerical
conclusions used in its proofs.

## Rules

1. `code/v4/` is self-contained and does not import archived implementations.
2. Root systems, order ideals, Chevalley constants, elliptic Weyl orbits, and
   twisted-cycle data are generated internally.
3. The reviewed F4 pair inventory is comparison data. The no-cutoff
   recomputation is decisive and is run whenever the F4 computation is run.
4. Generated V14 JSON files are outputs, not mathematical inputs.
5. Each machine-assisted theorem records its range, counts, exact rank
   witnesses, SageMath version, and SHA-256 digest.
6. Minimum positive degree is computed exactly, using target-root height as the
   intrinsic finite search bound. It is not used for BFS pruning because the
   corrected field-degree condition is not monotone in minimum degree.
7. The F4 pair inventory has `Psi2 != empty`. The affine-consistent empty subset
   is checked separately.
8. The exhaustive method computes the field-degree bound before considering
   pair data. It never scans all `2^|Phi+|` subsets.

## F4 Pair Computation

The comparison file is
[`psi-pairs-f4-v12.json`](https://github.com/mocham/AlgEG/blob/grob-v14/code/v4/DATA/psi-pairs-f4-v12.json).
It contains 4,862 reviewed pairs `(Psi0, Psi2)`.

The pair inventory is defined for `Psi2 != empty`. The empty subset satisfies
the affine equations, and both closure conditions hold for every
`Psi0 subset Phi+`. For each of the 104 nonempty `K` complements,

```text
Psi2=empty,
obstruction=0,
epsilon=|K_complement|+|Psi0 minus Phi_K|>=1.
```

The second cardinality is nonnegative, so the minimum over arbitrary `Psi0` is
the minimum nonzero `|K_complement|`, namely 1. Thus empty `Psi2` contributes no
epsilon-zero case. The program checks all 1,048 pairs consisting of a nonempty
`K` complement and a possible value of `|Psi0 minus Phi_K|`. The full
nonempty-subset computation does the following.

1. Generate all 9,780 linearly independent subsets of the 24 positive roots.
2. Adjoin roots by breadth-first search whenever the equations
   `f(alpha)=1` remain consistent. No degree test stops a branch.
3. At every one of the 50,347 nonempty affine-consistent subsets, recompute the
   integral difference lattice, `Psi0`, the closure conditions, and the degree
   value.
4. Obtain 4,862 pairs and compare the ordered pair-record list with the reviewed
   file. Both record lists have digest
   `e1dc8ced6a02269877dbfbc4a908404657d8a415a3b503f54d429888af29d4be`.
5. Check all 3,002 relevant independent subsets directly.
6. Apply each pair to the 104 nonempty `K` complements.

There is no positive-degree relation at an affine-consistent node: if
`f(alpha)=1` on `Psi2`, then `f` vanishes on the difference lattice, whereas a
sum of `d>0` roots in `Psi2` has value `d`. The degree computation is retained
as an assertion, not as pruning logic.

For nonempty `Psi2`, the formal product has `105 * 4862 = 510510` entries. The
4,862 empty-complement entries do not define nontrivial `K`-strata. Of the
remaining 505,648 entries, 294,388 are removed by the graph-target row bound and
211,260 receive exact polynomial Jacobian-rank checks. Exactly four have epsilon
zero. The separate empty-`Psi2` calculation has no epsilon-zero case.

## F4 Borel Count

For split simple F4, Configuration I has labels

```text
(a,1,0,1),  a in Z/(p-1)Z,
```

and Configurations II, III, and IV have labels

```text
(1,1,0,1), (0,1,0,1), (0,1,0,0).
```

Thus the number of configuration-labelled component records is

```text
(p-1)+3=p+2.
```

For `p=17`, Configuration I contributes 16 records and Configurations II-IV
contribute 3, for a total of 19. There are 17 distinct label tuples: the tuple
`(1,1,0,1)` occurs in Configuration I with `a=1` and in Configuration II, while
`(0,1,0,1)` occurs in Configuration I with `a=0` and in Configuration III.

## BCH Checks

`bch.py` computes every homogeneous class of `log(exp(X)exp(Y))`. Independently,
it multiplies noncommutative words in

```text
log(exp(tX)exp(tY)) mod t^12
```

and compares the coefficient of `t^m` with the recursive result for every
`1<=m<=11` in the enveloping algebra of the free Lie algebra. The checks also
include the exact A2 and B2 root formulas. Other tests compare the F4 Chevalley
constants and leading Jacobians with the local root-system implementation and
show that terms containing at least two `Psi2` factors do not enter the
first-order F4 rotation equations.

## Rotation-Root Checks

`verify_rotation_root_choices()` does not read or write generated files. It
uses the same equations as the four-row rotation computation and evaluates
exactly the following five choices.

| Configuration | Rotation root | `det(J)` | Prime divisors | Result |
|---|---:|---:|---|---|
| IV | `0010` | `895795200` | 2, 3, 5 | success |
| IV | `0011` | `-30800` | 2, 3, 5, 7, 11 | success |
| IV | `1000` | `315/4` | 2, 3, 5, 7 | success |
| III | `1000` | `15360` | 2, 3, 5 | success |
| III | `0010` | `-29160` | 2, 3, 5 | success |

For each row, the recorded exact values satisfy

```text
f_1=...=f_n=0, det(J)!=0, product(at_alpha*ap_alpha)!=0,
and every prime divisor is at most 12.
```

Thus all five choices work for every `p>12`. Within this list, neither
Configuration III nor Configuration IV has a unique successful root. There is
no failed choice to classify as having no solution, an everywhere-zero
Jacobian determinant, or a failure of the required nonvanishing condition.

## Type D Checks

For `D_n`, the program verifies

```text
h=2m, 1<=m<=floor((n-2)/2): C6 joined at one graph target to P_(2m-1),
all other h: empty after repeated degree-one graph-source removal,
det(J_(n,m))=+/-2*x_(n-2m-1)*product(x_i, i=n-m,...,n).
```

The finite checks use `4<=n<=32` for full adjacent-height graphs and
`4<=n<=12` for all graph-target subsets and determinants. They examine 928
full graphs, 15,888 graph-target selections, and 25 determinants. The JSON also
records the first longer case, `D_8` at height 6, and the SHA-256 values of
`algebra.py`, `graphs.py`, `roots.py`, and `type_d.py`.

## Effect of the Corrected Minimum-Degree Bound

The corrected necessary condition for a finite minimum degree `delta` is
`d_F >= (p-1)/gcd(p-1, delta)`, not `d_F >= delta`. This correction does not
change any previously reported computational result. Exact recomputation shows
that the minimum positive degree is undefined for every one of the 4,862 F4
pairs, so the corrected finite-degree filter removes no pair. The exhaustive
rerun reproduces the same pair count, classification digest, 91 epsilon-one
configurations, and four epsilon-zero configurations. The automatic result for
`d_F >= 2` comes from the independent maximal-cone rank calculation and does
not use the minimum-degree lemma.

## Reproduction

Run with SageMath 10.9 from `code/v4/`:

```sh
python -B -c 'from rotation import verify_rotation_root_choices as check; assert check()["ok"]'
python -B -m unittest discover -s tests -v
python cli.py quick --check
python cli.py f4-exhaustive --check
python cli.py type-d --check
python cli.py twisted-cycles --check
```

`python cli.py all --output-directory DATA` writes all V14 JSON and SHA-256
files. Without `--output-directory`, the default is the `DATA` directory next
to `cli.py`. With `--check`, the command reads and compares files in that
directory and does not create a directory or write any file. The reviewed F4
pair file is never written by the CLI. Check mode disables Python bytecode
before importing local modules, so the commands above do not create
`__pycache__` or `.pyc` files.

Build the consolidated computation document from the repository root with:

```sh
SOURCE_DATE_EPOCH=0 pdflatex -interaction=nonstopmode -halt-on-error computation-details.tex
SOURCE_DATE_EPOCH=0 pdflatex -interaction=nonstopmode -halt-on-error computation-details.tex
SOURCE_DATE_EPOCH=0 pdflatex -interaction=nonstopmode -halt-on-error computation-details.tex
```

The nine-page PDF has SHA-256
`fe20e5de03df2cc1afc6df1c227f91257e123871f942e35a0d84804e9474eaa5`.

V14 JSON files use schema `grob-computation-v14`; the reviewed pair file keeps
its existing machine schema. All JSON files use sorted indentation and one
terminal newline. Check every manifest from `code/v4/DATA/` with:

```sh
sha256sum -c -- *.sha256
```

## Rough Runtimes

These wall times were observed in the August 2026 workspace and are planning
estimates.

| Task | Rough time | Remarks |
|---|---:|---|
| No-cutoff F4 pair recomputation | 2m | Visits 50,347 subsets and compares 4,862 pairs |
| Unit tests | 9m | 27 tests, including F4, BCH, Type D, and CLI checks |
| Full non-writing recomputation | 10m20s | Recomputes and compares all nine V14 JSON files |
| Elliptic Weyl and twisted cycles | 2m | Run separately from `quick` |

## Module Map

| Module | Responsibility |
|---|---|
| `roots.py` | Cartan-generated roots, root order, structure constants |
| `graphs.py` | Adjacent-height bipartite graphs and deterministic edge removals |
| `algebra.py` | Exact polynomial and integral rank witnesses |
| `bounded_strata.py` | Degree bound, no-cutoff F4 pairs, and F4 classification |
| `borel.py`, `cone.py` | F4 table, Borel count, and grouped cone checks |
| `bch.py` | Homogeneous BCH recursion and root-graded formulas |
| `cycles.py` | Absolute checks for exceptional Chevalley types |
| `type_d.py` | Stable `D_n` graph and determinant formulas |
| `elliptic_weyl.py` | Elliptic Weyl orbits and non-Borel data generation |
| `twisted_cycles.py` | Relative and associated absolute cycle checks |
| `rotation.py` | Q-points and integral reduction data |
| `cli.py` | Deterministic JSON output and non-writing comparison mode |

## Artifacts

All files are in
[`code/v4/DATA/`](https://github.com/mocham/AlgEG/tree/grob-v14/code/v4/DATA).

| Artifact | SHA-256 |
|---|---|
| `psi-pairs-f4-v12.json` | `c92c488abfccab63f1f61ab0a34cb15d38fcd8feb4a0c87a46043d228710cc2f` |
| `bch-root-graded-v14.json` | `26d7f40bb9e34e3c3cfe65ee3bb22f07c14ec07f13a305b95d0894018710c1d5` |
| `borel-f4-table-v14.json` | `be156ca676e7f26bb81a1e8720e5f6e0ebfc87c7e95e30e642617ea9da813481` |
| `cone-f4-table-v14.json` | `abc597a0da91a7f2ff78a83a6767413ac30cc03c572ff0136c5f07f7f72ca6e7` |
| `cycles-exceptional-v14.json` | `59c41cd00d68e01d33e7b560ffa8d6ada149404f3cea5975fb1c6d9e0bdb8c70` |
| `cycles-type-d-v14.json` | `f951439a6657dd407305330aaf9e261c845e6b369bf2903eee51acf3aa65bfdd` |
| `fine-strata-f4-exhaustive-v14.json` | `693bdabb5c3233dd9f620866e7a00264d63f48405cf9b36e9f00ee3004351e28` |
| `rotation-f4-uniform-v14.json` | `34d56e56781b6bbbb7faca1ff273e6d82a2c37c7693cacc30f92ac20253230c8` |
| `elliptic-weyl-orbits-v14.json` | `2d9f51d0098e4c9aecca79830e3d763f6857dc35f492dbc7406f587bba078ef3` |
| `cycles-twisted-exceptional-v14.json` | `f0c8f809430b16a88828933d6b7af3f9744515f6e345bf221d92f290db01da57` |
