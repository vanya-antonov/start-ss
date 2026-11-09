# Planned mammalian analysis

This document describes a planned cross-species analysis of RNA secondary
structures around annotated translation start codons. The workflow is not yet
implemented: no mammalian scripts or data are present.

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

When implemented, mammalian scripts and inputs will live under
`scripts/mammals/` and `data/mammals/`. Derived results may be stored under
`results/mammals/`. No mammalian scripts or data have been implemented yet.
