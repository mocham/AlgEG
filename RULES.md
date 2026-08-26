# Current Revision Rules

1. Files under `LEGACY/` and `EXTERNAL/` are immutable.
2. The current theorem statements retain their established strength unless an
   actual counterexample is produced or concretely suggested. Repair proof or
   computational gaps instead of silently weakening a statement.
3. The current manuscript must not refer to files under `LEGACY/`.
4. Use `EXTERNAL/unobs/unobs.tex` as the principal external reference for the
   integral parabolic-extension and Weil--Deligne results proved there. Cite
   those results and explain only the crystalline or twisted variants needed
   here.
5. Keep computation out of the main exposition. The main paper may state a few
   useful numerical conclusions; algorithms belong in the appendices, and
   implementation rules, artifacts, checksums, and runtimes belong in
   `code/README.md`.
6. The reviewed F4 pair file is comparison data for
   `empty != Psi2 <= Phi+`. The empty subset is affine-consistent and must be
   checked separately for arbitrary `Psi0` using
   `epsilon=|K_complement|+|Psi0 minus Phi_K|>=1` for every nonempty
   `K_complement`. The breadth-first computation must examine every nonempty
   affine-consistent subset and must not stop a branch on the basis of a degree
   value.
7. `cli.py --check` must recompute and compare outputs without writing. It must
   disable Python bytecode before importing local modules. Its output directory
   is selected by `--output-directory`, with `DATA/` as the documented default.
8. `computation-details.tex` and `computation-details.pdf` must build
   deterministically with all fonts embedded.
