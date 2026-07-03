# Snakemake_UShER_pango_workflow

**Reconstructing clade-founder genomes for Pango lineages**

This repository contains a reproducible Snakemake workflow for extracting SARS-CoV-2 Pango lineage subtrees from a global mutation-annotated tree and reconstructing the founder genome sequence for each lineage.

The workflow is designed as a step-by-step pipeline for users who may not be familiar with Snakemake. Each stage is documented with what it does, why it is needed, and what files it produces.

## Pipeline overview

This repository includes the following main workflow files:

* `Snakefile_pango_extract`

## What this pipeline is for

SARS-CoV-2 evolves as it spreads. Each transmission event can introduce mutations, and genomes sampled around the world form a large phylogenetic tree. Public databases distribute this tree in a compact format called a **mutation-annotated tree** or **MAT**. A MAT stores a phylogeny of millions of genomes, with mutations placed on the branches where they arose.

Related genomes are grouped into named lineages. Two commonly used naming systems are:

* **Nextstrain clades**, which are broad labels such as `20B` or `24E`
* **Pango lineages**, which are more detailed labels such as `AD.2` or `KP.3.1.1`

A single Nextstrain clade can contain hundreds of Pango lineages. This distinction is important because the pipeline performs lineage-specific analysis at the Pango level.

## What the pipeline does

The pipeline performs two main tasks for each selected Pango lineage:

1. **Extracts the relevant part of the global SARS-CoV-2 tree**

   The workflow identifies the subtree corresponding to a chosen Pango lineage and extracts the associated mutations and sample paths.

2. **Reconstructs the founder genome for that lineage**

   The workflow reconstructs a FASTA sequence representing the inferred founder genome of the selected Pango lineage. It also produces the list of substitutions that define the lineage founder.

The Pango lineages to analyze are listed in a configuration file, making it easy to run the workflow on multiple lineages in a reproducible way.

## Reproducibility

The analysis is fully reproducible using the provided Conda environment.

First, install Conda if it is not already available on your system. Then create the environment using:

```bash
conda env create -f environment.yml
```

This will create a Conda environment named:

```bash
SARS2-mut-spectrum
```

Activate the environment with:

```bash
conda activate SARS2-mut-spectrum
```

After activating the environment, the Snakemake workflows can be run using the provided Snakefiles and configuration files.


## Checking if the environment is loaded

```bash
matUtils --version
```


## Repository structure

```text
.
├── environment.yml
├── Snakefile_pango_extract
├── Snakefile_pango_founders
├── scripts/
│   └── 5b_pango_founder_fasta.py
└── config/
```
## Running the pipline 

Please check and update the *Config* file for choosing the clade. 


## Expected outputs

For each Pango lineage, the workflow produces lineage-specific outputs, including:

* the extracted subtree or lineage-specific tree information
* mutation paths for samples in the lineage
* the inferred founder genome in FASTA format
* the substitutions defining the lineage founder

## Intended users

This workflow is intended for researchers working with SARS-CoV-2 phylogenetics, Pango lineages, mutation-annotated trees, and lineage-specific founder genome reconstruction.
