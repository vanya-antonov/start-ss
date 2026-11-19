# start-ss

## Project

Early-stage bioinformatics research on whether RNA secondary structure upstream
of or around translation start codons (AUG) influences translation initiation in
mammalian transcripts, especially human genes. Questions include conservation,
local structure and initiation, and associations with translation efficiency.

## Repository structure

- Standalone, dated research experiments belong under `work/`; keep each
  experiment's scripts, README, raw data, and results together in its own
  subdirectory. Do not move experiment-specific files into shared pipeline
  directories unless the task explicitly calls for it.
- README files under `work/` must be concise and include the following sections:
  - **Introduction** — a brief description of the analysis and its objective.
  - **Downloads** *(optional)* — all downloads required to execute the workflow scripts.
  - **Workflow** — commands required to execute the workflow, with one or two optional explanatory lines for each command.
    - Create a separate subsection for each script, numbered sequentially: `Step 1`, `Step 2`, etc.  
- Research documentation is organized by dated folders under `docs/2026/`.

## Development practices

- Keep Python readable and modular; use clear names, small functions, and
  comments for non-obvious scientific choices. Follow existing project style.
- Add dependencies to the relevant environment file only when needed and
  document why.
- Do not commit raw or restricted datasets, credentials, or machine-specific
  paths. Record data source, release/version, and access terms.

### Python scripts under `work/`

- Use `argparse` for command-line argument parsing.
- Each script must accept at least one required input argument.
- When run without arguments, the script should print a short description of
  its purpose and usage, including the required arguments, rather than a
  lengthy traceback.
- Scripts under `work/` must be run from the repository root using their
  relative paths.

Example:

```sh
python work/2026/1111.gabarap_orthologs.ncbi/02_some_script.py SOME_INPUT
```

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
