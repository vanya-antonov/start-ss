# Planned mammalian analysis

This document describes a planned cross-species analysis of RNA secondary
structures around annotated translation start codons. Canonical transcript
annotation retrieval is implemented; downstream orthology, sequence, and
structure analysis remains planned.

## Planned species and assemblies

The planned species are listed in the project order. Assemblies are UCSC
assembly names.

| Common name | Scientific name | UCSC assembly |
| --- | --- | --- |
| Human | *Homo sapiens* | hg38 |
| Mouse | *Mus musculus* | mm10 |
| Killer whale | *Orcinus orca* | orcOrc1 |
| Cow | *Bos taurus* | bosTau8 |
| Sheep | *Ovis aries* | oviAri3 |
| Panda | *Ailuropoda melanoleuca* | ailMel1 |
| Manatee | *Trichechus manatus latirostris* | triMan1 |
| Elephant | *Loxodonta africana* | loxAfr3 |
| Dog | *Canis lupus familiaris* | canFam3 |
| Aardvark | *Orycteropus afer afer* | oryAfe1 |
| Horse | *Equus caballus* | equCab2 |
| Rabbit | *Oryctolagus cuniculus* | oryCun2 |
| Marmoset | *Callithrix jacchus* | calJac3 |

## Planned workflow

### Download canonical transcript annotations

The downloader queries Ensembl BioMart for each species' canonical
(`transcript_is_canonical = 1`) transcript from protein-coding genes. It
resolves datasets by scientific name using the live BioMart dataset registry;
UCSC assembly IDs are not used as Ensembl dataset names. It does not retrieve
genome FASTA, GTF, or GenBank files.

Pass a species ID from `SPECIES_ORDER` to download one species, or pass
`__ALL__` to download all 13. Run from the repository root or another working
directory with the `start-ss-scripts` environment active:

```bash
python scripts/mammals/01_download_canonical_transcripts.py __ALL__
```

For example, download human annotations only with
`python scripts/mammals/01_download_canonical_transcripts.py human`.

It writes one tab-separated file per species to
`data/mammals/01_canonical_transcripts_{species_id}.tsv`, using the species IDs
in `SPECIES_ORDER`. Each file contains gene ID/name and biotype, transcript ID
and biotype, and protein ID and CDS length when BioMart provides them. Missing
optional fields are blank. API or dataset failures are logged, cause a
non-zero exit status, and do not replace that species' output with a partial
file.

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
