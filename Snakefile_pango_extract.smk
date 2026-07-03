"""Standalone, Pango-keyed MAT extraction pipeline.

Keeps the tree/reference infrastructure of the main Snakefile verbatim
(get_mat_tree, get_ref_fasta, get_ref_gtf, edit_ref_gtf_R) and the per-clade
extraction chain, but:

  * selects clades by PANGO lineage (config: pango_clades), not Nextstrain;
  * extracts with `matUtils extract -c <pango>` (native clade extraction);
  * names every output by the Pango lineage;
  * stops at the extraction artefacts + founder FASTA (no counting / fitness /
    report / K=2 rules).

Per Pango lineage it produces, under results_{mat}/pango_extract/ :
    {pango}.pb                        extracted MAT subtree
    {pango}_tree.nwk                  Newick tree
    {pango}_mutations.tsv             translated (AA) mutations
    {pango}_sample_paths.tsv          sample paths (nt, incl. noncoding)
    {pango}_nt_mutations.csv          nucleotide mutations
    {pango}.fa                        founder FASTA (indels excluded)
    {pango}_ref_to_founder_muts.csv   ref -> founder mutations

Run (from the repo root, so scripts/ resolves):
    snakemake -s Snakefile_pango_extract --configfile config_pango_extract.yaml --cores 4

Requires on PATH (as the main pipeline does): matUtils (UShER).
"""

import os

configfile: "config_pango_extract.yaml"

MAT = config["current_mat"]
SUBDIR = config.get("pango_extract_dir", "pango_extract")
PANGO_CLADES = config.get("pango_clades", [])
if not PANGO_CLADES:
    raise WorkflowError(
        "No Pango lineages requested. Add a `pango_clades:` list to the config."
    )

# Pango names are alphanumerics + dots (no underscores), so {pango} does not
# greedily consume the `_tree.nwk` / `_mutations.tsv` / ... suffixes.
wildcard_constraints:
    pango=r"[A-Za-z0-9.]+",
    mat=r"[A-Za-z0-9_]+",


def out(fname):
    return os.path.join("results_{mat}", SUBDIR, fname)


rule all:
    input:
        expand(out("{pango}.pb"), mat=MAT, pango=PANGO_CLADES),
        expand(out("{pango}_tree.nwk"), mat=MAT, pango=PANGO_CLADES),
        expand(out("{pango}_mutations.tsv"), mat=MAT, pango=PANGO_CLADES),
        expand(out("{pango}_sample_paths.tsv"), mat=MAT, pango=PANGO_CLADES),
        expand(out("{pango}_nt_mutations.csv"), mat=MAT, pango=PANGO_CLADES),
        expand(out("{pango}.fa"), mat=MAT, pango=PANGO_CLADES),
        expand(out("{pango}_ref_to_founder_muts.csv"), mat=MAT, pango=PANGO_CLADES),


# =========================================================================
# Tree / reference infrastructure  (unchanged from the main Snakefile)
# =========================================================================
rule get_mat_tree:
    """Get the pre-built mutation-annotated tree."""
    params:
        url=lambda wc: config["mat_trees"][wc.mat],
    output:
        mat="results_{mat}/mat/mat_tree.pb.gz",
    shell:
        "curl {params.url} > {output.mat}"


rule get_ref_fasta:
    """Get the reference FASTA."""
    params:
        url=config["ref_fasta"],
    output:
        ref_fasta="results_{mat}/ref/ref.fa",
    shell:
        "wget -O - {params.url} | gunzip -c > {output.ref_fasta}"


rule get_ref_gtf:
    """Get the reference GTF."""
    params:
        url=config["ref_gtf"],
    output:
        ref_gtf="results_{mat}/ref/original_ref.gtf",
    shell:
        "wget -O - {params.url} | gunzip -c > {output.ref_gtf}"


rule edit_ref_gtf_R:
    input:
        gtf=rules.get_ref_gtf.output.ref_gtf,
    output:
        gtf="results_{mat}/ref/edited_ref_R.gtf",
    params:
        edits=config["add_to_ref_gtf"],
    run:
        import json, tempfile
        tf = tempfile.NamedTemporaryFile(delete=False, suffix=".json")
        tf.write(json.dumps(params.edits).encode())
        tf.close()
        shell("Rscript scripts/1_ref_gtf_minimal.R {input.gtf} {output.gtf} {tf.name}")


rule clade_founder_jsons:
    """Get JSONs with clade founders (indels not included)."""
    params:
        **config["clade_founders"],
    output:
        neher_json="results_{mat}/clade_founders_no_indels/clade_founders_neher.json",
        roemer_json="results_{mat}/clade_founders_no_indels/clade_founders_roemer.json",
    shell:
        """
        curl -L {params.neher_json} > {output.neher_json}
        curl -L {params.roemer_json} > {output.roemer_json}
        """


# =========================================================================
# Per-Pango extraction  (selection by Pango; outputs named by Pango)
# =========================================================================
rule mat_clade_by_pango:
    """Extract MAT for one Pango lineage via native clade extraction.

    Mirrors the main pipeline's `matUtils extract ... -O -o` invocation, but
    selects with `-c <pango>` (the tree's Pango annotation) instead of a
    sample list. Extracts the monophyletic clade (subtree incl. descendants).
    """
    input:
        mat=rules.get_mat_tree.output.mat,
    output:
        mat="results_{mat}/pango_extract/{pango}.pb",
    shell:
        r"""
        matUtils extract -i {input.mat} -c "{wildcards.pango}" -O -o {output.mat}
        """


rule translate_mat:
    """Translate mutations on the extracted MAT (AA mutations)."""
    input:
        mat=rules.mat_clade_by_pango.output.mat,
        ref_fasta=rules.get_ref_fasta.output.ref_fasta,
        ref_gtf=rules.edit_ref_gtf_R.output.gtf,
    output:
        tsv="results_{mat}/pango_extract/{pango}_mutations.tsv",
    shell:
        """
        matUtils summary \
            -i {input.mat} \
            -g {input.ref_gtf} \
            -f {input.ref_fasta} \
            -t {output.tsv}
        """


rule mat_sample_path:
    """Get sample paths on the extracted MAT (nt mutations incl. noncoding)."""
    input:
        mat=rules.mat_clade_by_pango.output.mat,
    output:
        tsv="results_{mat}/pango_extract/{pango}_sample_paths.tsv",
    shell:
        "matUtils extract -i {input.mat} --sample-paths {output.tsv}"


rule sample_path_to_nt_mutations:
    """Get all nucleotide mutations from sample paths."""
    input:
        tsv=rules.mat_sample_path.output.tsv,
    output:
        csv="results_{mat}/pango_extract/{pango}_nt_mutations.csv",
    script:
        "scripts/4_sample_path_to_nt_mutations.py"


rule extract_tree_per_clade:
    """Newick tree for the extracted Pango clade."""
    input:
        mat=rules.mat_clade_by_pango.output.mat,
    output:
        tree="results_{mat}/pango_extract/{pango}_tree.nwk",
    shell:
        "matUtils extract -i {input.mat} -t {output.tree}"


# =========================================================================
# Founder FASTA (Pango-keyed, indels excluded)  -- reuses scripts/5b_...
# =========================================================================
rule pango_founder_fasta_and_muts:
    input:
        roemer_json=rules.clade_founder_jsons.output.roemer_json,
        ref_fasta=rules.get_ref_fasta.output.ref_fasta,
    params:
        strict_ref_check=config.get("founder_strict_ref_check", True),
    output:
        fasta="results_{mat}/pango_extract/{pango}.fa",
        muts="results_{mat}/pango_extract/{pango}_ref_to_founder_muts.csv",
    script:
        "scripts/5b_pango_founder_fasta.py"
