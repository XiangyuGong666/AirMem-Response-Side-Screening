# Release history

## v0.1.2 — Calculation and figure provenance audit

- Remove Figure 1; retain original numbering for Figures 3–23.
- Require a fresh calculation with verified intermediate files for figures, and record actual per-figure inputs and output hashes.
- Remove unused legacy damping group targets and the obsolete 300-window cap.
- Screen the end-aligned window explicitly and record the decision in the window ledger.

## v0.1.1 — Windows single executable

- Distribute one Windows x64 EXE with its runtime, Tk desktop interface, PDF/PNG renderers and bare FE reference inputs.
- Update threadpoolctl to handle OpenBLAS builds that do not report a version.
- Keep the V19.0 response-side calculation and figure definitions from v0.1.0. Monitoring records remain a separate input.

## v0.1.0 — initial standalone source and Windows folder bundle

- Extracted the V19.0 response-side workflow into its own project and validated the 2022Q4 results.
