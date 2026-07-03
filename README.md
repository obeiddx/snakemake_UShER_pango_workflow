# snakemake_UShER_pango_workflow
Reconstructing clade-founder genomes for Pango lineages


A step-by-step walkthrough of the Snakemake pipeline: what each step does, why it is there, and what it
produces. Written for readers who have never used Snakemake.
Pipelines: Snakefile_pango_extract and Snakefile_pango_founders · Script: scripts/
5b_pango_founder_fasta.py

What this pipeline is for ?

SARS-CoV-2 evolves as it spreads: each transmission can introduce mutations, and the genomes sampled
worldwide form a giant family tree. Public databases distribute this tree in a compact form called a mutation-
annotated tree (MAT): a phylogeny of millions of genomes with the mutations placed on the branches
where they arose.

Related genomes are grouped into named lineages. Two naming systems are in common use: Nextstrain
clades (coarse labels such as 20B or 24E ) and Pango lineages (fine-grained labels such as AD.2 or KP.
3.1.1 ). One Nextstrain clade can contain hundreds of Pango lineages. Section 6 shows exactly why that
difference matters for our analysis.

The two jobs of this pipeline are:
1. Extract the part of the global tree belonging to a chosen Pango lineage, and pull out its mutations and
2. sample paths.define it. Reconstruct the founder genome (a FASTA sequence) for that lineage, plus the list of substitutions that
Both are done per Pango lineage, and the lineages of interest are simply listed in a configuration file.
