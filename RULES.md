# Current Revision Rules

1. Files under `LEGACY/` and `EXTERNAL/` are immutable.
2. The current theorem statements retain their established strength unless an
   actual counterexample is produced or concretely suggested. Repair proof or
   certification gaps instead of silently weakening a statement.
3. The current manuscript must not refer to files under `LEGACY/`.
4. Use `EXTERNAL/unobs/unobs.tex` as the principal external source for the
   integral fixed-type, parabolic-extension, and Weil--Deligne results proved
   there. Cite those results and explain only the crystalline or twisted
   variants needed here.
5. Keep computation out of the main exposition. The main paper may state a few
   useful numerical conclusions; algorithms belong in the appendices, and
   implementation rules, artifacts, checksums, and runtimes belong in
   `code/README.md`.
