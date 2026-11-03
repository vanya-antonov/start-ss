# start-ss

## Project

Early-stage bioinformatics research on whether RNA secondary structure upstream
of or around translation start codons (AUG) influences translation initiation in
mammalian transcripts, especially human genes. Questions include conservation,
local structure and initiation, and associations with translation efficiency.

## Repository today

- `README.md`: project title and Conda setup instructions.
- `environment.yml`: Python 3.14 and Django 6.1.2.
- `django/`: generated Django project scaffold (`manage.py`, settings, URLs,
  ASGI/WSGI); only the admin URL is configured.
- No analysis pipeline, application-specific code, notebooks, datasets, or
  tests are present yet. Do not assume these have been implemented.

## Development practices

- Keep Python readable and modular; use clear names, small functions, and
  comments for non-obvious scientific choices. Follow existing project style.
- Keep research analysis separate from Django configuration and interfaces.
  Add dependencies to `environment.yml` only when needed and document why.
- Do not commit raw or restricted datasets, credentials, or machine-specific
  paths. Record data source, release/version, and access terms.

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
