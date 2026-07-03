"""Build a clade-founder FASTA + ref->founder mutation table for a single
Pango lineage, keyed *directly* off the Roemer pango-consensus summary JSON.

Indels (nucDeletions) are NOT applied to the sequence (founder "no indels"),
matching the convention of the main pipeline's clade_founder_fasta_and_muts.

Substitutions in the Roemer JSON are relative to Wuhan-Hu-1 (wuhCor1 /
MN908947.3), 1-based, format e.g. 'C241T' = ref C at site 241 -> T.

Written as a Snakemake `script:`; core logic is in build_founder() so it can
be unit-tested independently of Snakemake.
"""

import csv
import json
import re

SUB_RE = re.compile(r"^([ACGT])(\d+)([ACGT])$")


def read_ref_fasta(path):
    """Read a single-record reference FASTA -> (header, uppercase sequence)."""
    header = None
    seq_parts = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    raise ValueError(f"{path}: expected a single FASTA record")
                header = line[1:].strip()
            else:
                seq_parts.append(line.strip())
    if header is None:
        raise ValueError(f"{path}: no FASTA header found")
    return header, "".join(seq_parts).upper()


def write_fasta(path, name, seq, wrap=70):
    with open(path, "w") as f:
        f.write(f">{name}\n")
        for i in range(0, len(seq), wrap):
            f.write(seq[i : i + wrap] + "\n")


def build_founder(ref_fasta, roemer_json, pango, out_fasta, out_muts,
                  strict=True):
    """Apply Pango nucSubstitutions to the reference and write FASTA + muts.

    strict=True raises if any substitution's stated ref base disagrees with the
    reference FASTA (catches indexing / ref-version mismatches). Returns a dict
    summary for logging/testing.
    """
    _, ref = read_ref_fasta(ref_fasta)
    with open(roemer_json) as f:
        summary = json.load(f)

    if pango not in summary:
        # surface near matches to make typos obvious
        prefix = pango.split(".")[0]
        near = sorted(k for k in summary if k.startswith(prefix))[:20]
        raise KeyError(
            f"Pango lineage {pango!r} not in Roemer JSON "
            f"({len(summary)} lineages). Nearby keys: {near}"
        )

    entry = summary[pango]
    seq = list(ref)
    muts = []
    mismatches = []

    for s in entry.get("nucSubstitutions", []):
        if not s:
            continue
        m = SUB_RE.match(s)
        if not m:
            raise ValueError(f"{pango}: unparseable substitution {s!r}")
        ref_nt, pos, alt_nt = m.group(1), int(m.group(2)), m.group(3)
        if pos < 1 or pos > len(seq):
            raise ValueError(f"{pango}: site {pos} outside reference length {len(seq)}")
        obs = seq[pos - 1]
        if obs != ref_nt:
            mismatches.append((pos, ref_nt, obs))
        seq[pos - 1] = alt_nt
        muts.append((pos, ref_nt, alt_nt, s))

    if mismatches and strict:
        head = ", ".join(f"site {p}: JSON={r} ref={o}" for p, r, o in mismatches[:10])
        raise ValueError(
            f"{pango}: {len(mismatches)} substitution(s) disagree with reference "
            f"base (possible ref-version / indexing issue): {head}"
        )

    founder = "".join(seq)
    write_fasta(out_fasta, pango, founder)

    with open(out_muts, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["pango", "nextstrain_clade", "nt_site", "ref_nt",
                    "founder_nt", "nt_mutation"])
        nsclade = entry.get("nextstrainClade", "")
        for pos, ref_nt, alt_nt, s in sorted(muts):
            w.writerow([pango, nsclade, pos, ref_nt, alt_nt, s])

    return {
        "pango": pango,
        "nextstrain_clade": entry.get("nextstrainClade", ""),
        "n_substitutions": len(muts),
        "n_deletions_excluded": len([d for d in entry.get("nucDeletions", []) if d]),
        "n_ref_mismatches": len(mismatches),
        "seq_len": len(founder),
    }


if "snakemake" in globals():  # running under Snakemake
    summary = build_founder(
        ref_fasta=snakemake.input.ref_fasta,
        roemer_json=snakemake.input.roemer_json,
        pango=snakemake.wildcards.pango,
        out_fasta=snakemake.output.fasta,
        out_muts=snakemake.output.muts,
        strict=snakemake.params.get("strict_ref_check", True),
    )
    print(summary)
