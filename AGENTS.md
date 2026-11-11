# start-ss

## Project

Early-stage bioinformatics research on whether RNA secondary structure upstream
of or around translation start codons (AUG) influences translation initiation in
mammalian transcripts, especially human genes. Questions include conservation,
local structure and initiation, and associations with translation efficiency.

## Repository structure

- Human and mammalian analyses are separate pipelines.
- The existing human scripts, data, and results remain under
  `scripts/human/`, `data/human/`, and `results/human/`.
- Planned mammalian scripts and data belong under `scripts/mammals/` and
  `data/mammals/`. Mammalian results may later be stored under
  `results/mammals/`.
- Research documentation is organized by dated folders under `docs/2026/`.
- Standalone, dated research experiments belong under `work/`; keep each
  experiment's scripts, README, raw data, and results together in its own
  subdirectory. Do not move experiment-specific files into shared pipeline
  directories unless the task explicitly calls for it.
- No shared `constants.py` or mammalian scripts exist yet. Do not assume the
  planned mammalian workflow has been implemented.
- `django/` is a generated Django project scaffold; keep research analysis
  separate from Django configuration and interfaces.

## Development practices

- Keep Python readable and modular; use clear names, small functions, and
  comments for non-obvious scientific choices. Follow existing project style.
- Add dependencies to the relevant environment file only when needed and
  document why.
- Do not commit raw or restricted datasets, credentials, or machine-specific
  paths. Record data source, release/version, and access terms.

## Execution policy

- Never execute scripts, Python files, shell commands, or other project
  commands unless the current prompt explicitly asks you to do so.
- This prohibition includes running tests, linters, formatters, installation
  commands, build commands, and pipeline scripts.
- Do not download data, trigger network requests, or run commands that modify
  the environment unless explicitly requested.
- By default, limit work to inspecting files, editing code, and performing
  static analysis.
- If execution would be useful for validation, describe the proposed command
  but do not run it.
- An explicit request to execute something must come from the current prompt;
  do not infer permission from previous prompts or general project
  instructions.

## Reproducible analysis

- Make inputs, outputs, parameters, and random seeds explicit; preserve scripts
  or notebooks needed to reproduce derived results.
- Record organism, genome/transcript annotation release, transcript and
  isoform identifiers, strand, sequence source, and the exact region/window
  definition relative to AUG. Handle alternative start sites and transcript
  selection explicitly.
- For future structure analyses, consider ViennaRNA/RNAfold, MFE, RNA
  accessibility, and base-pairing probabilities. Record software/version,
  settings (including temperature and constraints), and sequence boundaries;
  these are candidate methods, not existing functionality.
- For comparisons across transcript regions or species, define matched regions,
  controls, orthology/alignment criteria, and handling of missing data. Treat
  predicted structure as a model, and distinguish association from causation.

## Decisions and reporting

- Document consequential scientific and technical choices close to the
  analysis (README or a short methods/decision note): rationale, alternatives
  considered, versions, and expected effect on interpretation.
- Report sample selection, exclusions, parameters, and limitations with each
  result so analyses can be reviewed and repeated.
- README files under `work/` must be short: include a brief introduction, optional Downloads section and a section with commands needed to execute the workflow (+ one or two optional lines explaining each command). Create subheader/section for each script starting with Step 1, Step 2 etc.
