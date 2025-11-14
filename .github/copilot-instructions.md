Repository-specific guidance for AI coding agents working on the BIH / ArchAngel project

Summary
- This repository provides an ArchAngel community-detection implementation (igraph-based) in `angel/alg/` and a higher-level BIH (Based Importance Hiding) pipeline implemented in `BihHiding.py` that: (1) runs ArchAngel on a dynamic .ncol edge list, (2) applies BIH rewiring operations per-snapshot to attempt to "hide" overlapping nodes, and (3) re-runs ArchAngel to compare results.

Immediate goals for edits
- Preserve the behavior of `BihHiding.process_all_snapshots` unless the change explicitly targets the BIH algorithm or output layout.
- Keep the `angel` package API stable: main entrypoints are `angel.Angel(...)` and `angel.ArchAngel(...)` (see `angel/__init__.py`).

Key files and what they teach you
- `BihHiding.py` — pipeline orchestration, I/O conventions (expects `merged_snapshots.ncol` in repo root), output under `BIH_Results/`. Look here for snapshot parsing, community file naming conventions (prefixes `original_detection` and `modified_detection`) and CSV output format `BIH_Results/comparison/hiding_comparison.csv`.
- `angel/alg/iArchAngel.py` — reads dynamic .ncol files by slicing on the 3rd column (snapshot id). It writes per-slice files named `<outfile_path>ArchAngel_coms_<snapshot>.txt` and a `*_ct_matches.csv` file. Use `ArchAngel(...).execute()` as shown in `BihHiding.run_archangel` and the tests.
- `angel/alg/iAngel.py` — core community detection (igraph). Important parameters: `threshold`, `match_threshold` (ArchAngel), and `min_comsize`. Outputs are written when `save=True`.
- `angel/test/angel_test.py` — shows minimal test usage: how to construct Angel/ArchAngel and expected outputs (tests remove generated files). Follow the same patterns for quick runnable checks.

Project conventions and gotchas for code edits
- File formats: dynamic networks use 3-column tab-separated lines: nodeA\tnodeB\tsnapshotID. Per-snapshot output files from ArchAngel follow the pattern `<prefix>ArchAngel_coms_<snapshot>.txt` and a matches CSV named `<prefix>ArchAngel_coms_ct_matches.csv`.
- Naming/renaming: `BihHiding.run_archangel` expects ArchAngel to generate filenames prefixed with `angel` and then renames them to either `original_detection...` or `modified_detection...`. When changing ArchAngel filename logic, update `run_archangel` accordingly.
- Verbosity: many classes accept `verbose` / `save` flags; tests exercise non-verbose usage. Keep `verbose=False` in programmatic calls where progress logging would otherwise flood CI logs.
- Dependencies: this project uses `igraph`, `networkx`, `numpy`, and `tqdm`. Ensure they are installed in the environment used for running tests or scripts.

Developer workflows (how to build, run and test)
- Run the main BIH pipeline locally (assumes repo root contains `merged_snapshots.ncol`):
  - Edit `BihHiding.py` `__main__` base_path if needed, then run `python BihHiding.py` from repo root.
- Run unit tests (simple): run the test file `angel/test/angel_test.py` with `python -m unittest angel.test.angel_test`.
- Quick checks: the tests write and delete files named `angels_coms.txt` and `ArchAngel*` — review file removal logic if you change output filenames.

What to change (and what not to change) in PRs
- Safe to change: algorithm internals inside `BihHiding.apply_bih_to_snapshot`, improving clarity, adding unit tests for edge cases, or making ArchAngel output path configurable — but preserve on-disk naming unless you update callers in `BihHiding.py` and tests.
- Risky to change without coordination: the slice-reading logic in `iArchAngel.__read_snapshot` (it affects how snapshots are discovered), the output file naming conventions, and the CSV fields written by `process_all_snapshots` — these are relied upon by the pipeline.

Useful examples to copy
- How to construct ArchAngel programmatically (from tests and `BihHiding.run_archangel`):
  - aa = angel.ArchAngel(network_filename, threshold=0.35, match_threshold=0.35, min_comsize=3, save=True, outfile_path=output_dir)
  - aa.execute()  # produces per-snapshot files and a matches CSV
- How snapshots are parsed from `.ncol` in `iArchAngel.__read_snapshot`: lines are split on tabs and the third token is used as the snapshot id.

If unsure, ask about
- Whether renaming of `angel`-prefixed output files (done in `BihHiding.run_archangel`) should be preserved or simplified.
- Whether the repository should add a small `requirements.txt` (recommended) listing `networkx, igraph, numpy, tqdm` for reproducible runs.

End of file — after updates, run tests and a short BIH run on a small `merged_snapshots.ncol` to validate behavior.
