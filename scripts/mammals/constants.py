"""Shared constants for the planned mammalian analysis pipeline."""

SPECIES_DICT = {
    "human": {
        "sort_order": 0,
        "title": "Human",
        "species": "Homo sapiens",
        "ncbi_taxon_id": 9606,
        "ucsc_assembly": "hg38",
    },
    "mouse": {
        "sort_order": 1,
        "title": "Mouse",
        "species": "Mus musculus",
        "ncbi_taxon_id": 10090,
        "ucsc_assembly": "mm10",
    },
    "whale": {
        "sort_order": 2,
        "title": "Killer whale",
        "species": "Orcinus orca",
        "ncbi_taxon_id": 9733,
        "ucsc_assembly": "orcOrc1",
    },
    "cow": {
        "sort_order": 3,
        "title": "Cow",
        "species": "Bos taurus",
        "ncbi_taxon_id": 9913,
        "ucsc_assembly": "bosTau8",
    },
    "sheep": {
        "sort_order": 4,
        "title": "Sheep",
        "species": "Ovis aries",
        "ncbi_taxon_id": 9940,
        "ucsc_assembly": "oviAri3",
    },
    "panda": {
        "sort_order": 5,
        "title": "Panda",
        "species": "Ailuropoda melanoleuca",
        "ncbi_taxon_id": 9646,
        "ucsc_assembly": "ailMel1",
    },
    "manatee": {
        "sort_order": 6,
        "title": "Manatee",
        "species": "Trichechus manatus latirostris",
        "ncbi_taxon_id": 127582,
        "ucsc_assembly": "triMan1",
    },
    "elephant": {
        "sort_order": 7,
        "title": "Elephant",
        "species": "Loxodonta africana",
        "ncbi_taxon_id": 9785,
        "ucsc_assembly": "loxAfr3",
    },
    "dog": {
        "sort_order": 8,
        "title": "Dog",
        "species": "Canis lupus familiaris",
        "ncbi_taxon_id": 9615,
        "ucsc_assembly": "canFam3",
    },
    "aardvark": {
        "sort_order": 9,
        "title": "Aardvark",
        "species": "Orycteropus afer afer",
        "ncbi_taxon_id": 1230840,
        "ucsc_assembly": "oryAfe1",
    },
    "horse": {
        "sort_order": 10,
        "title": "Horse",
        "species": "Equus caballus",
        "ncbi_taxon_id": 9796,
        "ucsc_assembly": "equCab2",
    },
    "rabbit": {
        "sort_order": 11,
        "title": "Rabbit",
        "species": "Oryctolagus cuniculus",
        "ncbi_taxon_id": 9986,
        "ucsc_assembly": "oryCun2",
    },
    "marmoset": {
        "sort_order": 12,
        "title": "Marmoset",
        "species": "Callithrix jacchus",
        "ncbi_taxon_id": 9483,
        "ucsc_assembly": "calJac3",
    },
}

SPECIES_ORDER = [
    key for key, _ in sorted(SPECIES_DICT.items(), key=lambda item: item[1]["sort_order"])
]
