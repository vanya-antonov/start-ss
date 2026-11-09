# Planned mammalian analysis

This document describes a planned cross-species analysis of RNA secondary
structures around annotated translation start codons. Canonical transcript
annotation retrieval is implemented; downstream orthology, sequence, and
structure analysis remains planned.

## Planned species and assemblies

The planned species are listed in the project order. Assemblies are UCSC
assembly names.

| No. | Species ID | Common name | Scientific name | UCSC assembly |
| --- | --- | --- | --- | --- |
| 1 | human | Human | *Homo sapiens* | hg38 |
| 2 | mouse | Mouse | *Mus musculus* | mm10 |
| 3 | whale | Killer whale | *Orcinus orca* | orcOrc1 |
| 4 | cow | Cow | *Bos taurus* | bosTau8 |
| 5 | sheep | Sheep | *Ovis aries* | oviAri3 |
| 6 | panda | Panda | *Ailuropoda melanoleuca* | ailMel1 |
| 7 | manatee | Manatee | *Trichechus manatus latirostris* | triMan1 |
| 8 | elephant | Elephant | *Loxodonta africana* | loxAfr3 |
| 9 | dog | Dog | *Canis lupus familiaris* | canFam3 |
| 10 | aardvark | Aardvark | *Orycteropus afer afer* | oryAfe1 |
| 11 | horse | Horse | *Equus caballus* | equCab2 |
| 12 | rabbit | Rabbit | *Oryctolagus cuniculus* | oryCun2 |
| 13 | marmoset | Marmoset | *Callithrix jacchus* | calJac3 |

## Planned workflow

### Download canonical transcript annotations

The workflow uses the Ensembl REST API at `/info/species` to match each
configured common and scientific name to an Ensembl production name, and
`/info/data` to determine the current release. Download the matching GTF
annotation manually before running the script. The script checks for that file
under `data/mammals/`, parses it in streaming mode, and never downloads,
deletes, or modifies the GTF. UCSC assembly IDs are not treated as Ensembl
production names. Only transcripts tagged `Ensembl_canonical` with both gene
and transcript biotypes `protein_coding` are included. CDS length is the sum
of the transcript's CDS feature lengths. Genome FASTA and GenBank files are not
downloaded.

### Step 1. Download the Ensembl GTF file

For the species you plan to process, use its Ensembl REST `/info/species`
record to get `name` (production name), `display_name`, and `assembly`; use
`/info/data` to get the current release. Set the shell variables below from
that record. For `species_name`, replace spaces in `display_name` with
underscores. Repeat this step for each species when processing `__ALL__`.

```bash
species_name="<display_name_with_underscores>"
production_name="<ensembl_production_name>"
assembly="<ensembl_assembly>"
release="<current_ensembl_release>"

wget -P data/mammals \
  "https://ftp.ensembl.org/pub/release-${release}/gtf/${production_name}/${species_name}.${assembly}.${release}.gtf.gz"
```

`wget -P` retains the GTF's Ensembl filename in `data/mammals/`; that is the
filename the script expects. The downloaded `*.gtf.gz` files are ignored by
Git.

### Step 2. Parse canonical transcript annotations

Pass a species ID from `SPECIES_ORDER` to process one species, or pass
`__ALL__` after downloading all 13 GTF files. Run from the repository root or
another working directory with the `start-ss-scripts` environment active:

```bash
python scripts/mammals/01_download_canonical_transcripts.py __ALL__
```

For example, process human annotations only with
`python scripts/mammals/01_download_canonical_transcripts.py human`.

It writes one tab-separated file per species to
`data/mammals/01_canonical_transcripts_{species_id}.tsv`, using the species IDs
in `SPECIES_ORDER`. Each file contains gene ID/name and biotype, transcript ID
and biotype, protein ID, and CDS length. REST, FTP, missing-species, malformed
GTF, and missing-canonical-annotation failures are logged and cause a non-zero
exit status. A species output is written only after its GTF has been completely
downloaded and validated.

### Planned downstream analysis

The proposed analysis will identify orthologous genes across these species,
retrieve annotated transcripts, and extract transcript-oriented 50-nt windows
around each annotated start codon: 23 nt upstream, the ATG start codon, and
24 nt downstream. The resulting windows will be used to analyze conservation
of predicted RNA secondary structures.

Implementation will need to document orthology criteria, transcript and
isoform selection, annotation and sequence sources, missing-data handling, and
the structure prediction software, version, and settings. Predicted structure
is a model; conservation or association in these data does not establish a
causal effect on translation initiation.

Additional mammalian scripts and inputs will live under `scripts/mammals/` and
`data/mammals/`. Derived results may be stored under `results/mammals/`.
