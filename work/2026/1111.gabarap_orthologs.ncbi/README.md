# NCBI GABARAP ortholog workflow

## Objective

Resolve a supplied human gene symbol to exactly one human NCBI Gene ID, retrieve
NCBI Datasets' complete calculated ortholog set for that reference, and retain
ortholog reports whose NCBI taxonomy IDs match the nonhuman species in
`scripts/mammals/constants.py`. Orthology is taken only from NCBI's assignment;
this workflow does not infer orthology from sequence similarity.

## NCBI API endpoints and retrieval

The workflow uses the NCBI Datasets v2 REST API at
`https://api.ncbi.nlm.nih.gov/datasets/v2`:

1. `GET /gene/symbol/{symbol}/taxon/9606/dataset_report` resolves the symbol in
   *Homo sapiens*. The response is retained in `data/`. Candidate reports are
   filtered to an exact case-insensitive symbol match; zero or multiple exact
   matches stop the analysis.
2. `GET /gene/id/{gene_id}/orthologs` retrieves the ortholog set for the
   resolved GeneID. No taxon restriction is supplied, so the complete set made
   available by NCBI is requested. This response is retained in `data/`.

The script selects records using the NCBI `tax_id`, preserving every returned
ortholog record in each configured target species. HTTP requests use a 60-second
timeout, retries with exponential backoff, and at most five requests per second
(the default NCBI Datasets rate limit). Existing raw response files are reused;
delete the relevant JSON file inside `data/` to refresh that response. Raw
responses may include records outside the project species list; the processed
table includes only configured targets.

NCBI's API documentation describes the API and its OpenAPI specification at
[NCBI Datasets v2 REST API](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/api/rest-api/).
The ortholog-set semantics are documented in [Download an ortholog data
package](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/genes/download-ortholog-data-package/)
and [Retrieve ortholog data and metadata](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/tutorials/download-ortholog-dataset/).

## Species included

All nonhuman keys in `SPECIES_DICT` are included. Species display names come
from that shared dictionary; the associated NCBI taxonomy IDs are configured
once in the experiment script because the shared constants currently contain
no taxonomy IDs.

| Key | Species | NCBI taxon ID |
| --- | --- | ---: |
| mouse | Mus musculus | 10090 |
| whale | Orcinus orca | 9733 |
| cow | Bos taurus | 9913 |
| sheep | Ovis aries | 9940 |
| panda | Ailuropoda melanoleuca | 9646 |
| manatee | Trichechus manatus latirostris | 127582 |
| elephant | Loxodonta africana | 9785 |
| dog | Canis lupus familiaris | 9615 |
| aardvark | Orycteropus afer afer | 9818 |
| horse | Equus caballus | 9796 |
| rabbit | Oryctolagus cuniculus | 9986 |
| marmoset | Callithrix jacchus | 9483 |

The reference is *Homo sapiens*, taxon ID 9606, and is not a target.

## Usage and reproducibility

From the repository root, run:

```sh
python work/2026/1111.gabarap_orthologs.ncbi/01_download_orthologs_from_ncbi.py GABARAP
```

The input symbol determines the output basename. For another gene, replace
`GABARAP` with its human symbol. The script uses only the Python standard
library and the repository's shared mammalian species constants. It records the
retrieval date in its console output; rerun the command to reproduce the
request using the current NCBI records. The raw JSON snapshots in `data/` make
the API inputs to the processed table inspectable and are reused on subsequent
runs.

## Output

The processed table is written to
`results/{SYMBOL}_ncbi_orthologs.tsv`, one row per returned ortholog. There are
no placeholder rows for target species lacking an NCBI ortholog record.

| Column | Meaning |
| --- | --- |
| `human_gene_symbol` | Input human symbol |
| `human_gene_id` | Resolved human NCBI GeneID |
| `species` | Species name from the shared project constants |
| `taxon_id` | Configured NCBI taxonomy ID used to assign the target species |
| `ortholog_gene_id` | Target NCBI GeneID |
| `ortholog_gene_symbol` | Target gene symbol, when present |
| `ortholog_description` | Target gene description, when present |
| `ortholog_type` | Orthology type, when the endpoint provides it |
| `refseq_accession` | Available RefSeq RNA accession, when present |

## Limitations

NCBI ortholog sets are computed assignments, not experimental proof of
functional equivalence. Coverage depends on the species and gene annotations
available to NCBI and may omit genes without an assignment. Gene symbols and
descriptions may differ across species or be absent. A symbol can be ambiguous;
the workflow deliberately stops unless exactly one exact human symbol match is
reported. RefSeq accessions are optional and may depend on endpoint annotation
details. The taxonomy ID map is specific to the species keys currently in the
shared constants and must be updated if those keys change. API content can
change over time; raw JSON captures the response used for each output.

## Retrieval date

This workflow definition was prepared on **2026-11-11**. The exact data
retrieval date is printed by each script run; record that date alongside any
reported result.
