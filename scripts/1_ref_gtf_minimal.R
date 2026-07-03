#!/usr/bin/env Rscript
# Minimal GTF editor (commented)
# Reads a GTF, inserts manual gene records, writes edited GTF.
# Works under Snakemake (as script:) or as a CLI script:
# Rscript edit_ref_gtf_minimal.R <in.gtf> <out.gtf> '{"GENE":[start,end],...}'

suppressPackageStartupMessages({
  # data.table for fast table IO & manipulation
  # jsonlite to parse JSON string for edits
  library(data.table)
  library(jsonlite)
})

# ---- Input handling: support Snakemake or command-line ----
# If run inside Snakemake via `script: "scripts/edit_ref_gtf_minimal.R"`,
# snakemake object is available and we read inputs from snakemake@...
if (exists("snakemake")) {
  in_gtf  <- snakemake@input[[1]]        # input GTF path provided by Snakemake
  out_gtf <- snakemake@output[[1]]       # output GTF path Snakemake expects
  edits_p <- snakemake@params[["edits"]] # edits param (could be R list or JSON string)
} else {
  # Otherwise get CLI args: input_gft, output_gtf, edits_json
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 3) stop("Usage: Rscript edit_ref_gtf_minimal.R <in.gtf> <out.gtf> <edits_json>")
  in_gtf  <- args[1]  # e.g. "results_mat/ref/original_ref.gtf"
  out_gtf <- args[2]  # e.g. "results_mat/ref/edited_ref.gtf"
  edits_p <- args[3]  # JSON string like '{"GENE1":[100,200],"GENE2":[500,800]}'
}

# ---- Parse edits param into an R named list: gene -> c(start,end) ----
parse_edits <- function(x) {
  # Accepts:
  #  - a named R list (e.g., passed from Snakemake), or
  #  - a single JSON string (CLI) that maps gene -> [start,end].
  if (is.null(x)) return(list())
  if (is.list(x)) return(lapply(x, function(v) as.integer(v[1:2])))
  if (is.character(x) && length(x) == 1) {
    parsed <- fromJSON(x, simplifyVector = TRUE) # jsonlite::fromJSON
    return(lapply(parsed, function(v) as.integer(v[1:2])))
  }
  stop("Unsupported edits format")
}
edits <- parse_edits(edits_p) # now edits is a named list: edits[["GENE"]] -> integer vector length 2

# ---- Read GTF into a data.table ----
# GTF has 9 tab-separated columns with no header; we name them explicitly here
cols <- c("seqname","source","feature","start","end","score","strand","frame","attribute")
# fread: very fast reader; quote="" avoids stripping quotes inside attributes
gtf <- fread(in_gtf, sep = "\t", header = FALSE, col.names = cols, quote = "")

# ensure numeric columns for later numeric comparisons
gtf[, `:=`(start = as.integer(start), end = as.integer(end))]

# ---- Basic consistency checks (keeps behavior identical to the longer version) ----
stopifnot(length(unique(gtf$seqname)) == 1L)  # ensure only one reference contig/seq in GTF
stopifnot(length(unique(gtf$source)) == 1L)   # ensure single source field (same formatting expectation)
trans_starts <- gtf[feature == "transcript", start]
# require transcript starts sorted (the original Python asserted this)
if (!identical(trans_starts, sort(trans_starts))) stop("transcript starts not sorted")

# store single seqname and source for creating new records
seqname <- unique(gtf$seqname)[1]
src     <- unique(gtf$source)[1]

# ---- Helper: build the 5 records (transcript, exon, CDS, start_codon, stop_codon) for a gene ----
build_rows <- function(gene, s, e) {
  s <- as.integer(s); e <- as.integer(e)
  if (s >= e) stop("start must be < end for ", gene)
  # attributes string for different features (keeps same attribute formatting as original)
  attrs <- function(feat) {
    if (feat == "transcript")
      sprintf('gene_id "%s"; transcript_id "%s";', gene, gene)
    else if (feat == "exon")
      sprintf('gene_id "%s"; transcript_id "%s"; exon_number "1"; exon_id "%s";', gene, gene, gene)
    else
      sprintf('gene_id "%s"; transcript_id "%s";', gene, gene)
  }
  # create rows as lists (keeps columns consistent with original script)
  rows <- list(
    list(seqname, src, "transcript", s, e, ".", "+", ".", attrs("transcript")),
    list(seqname, src, "exon",       s, e, ".", "+", ".", attrs("exon")),
    # CDS end uses e-3 (common convention to exclude stop codon in CDS), but never < s
    list(seqname, src, "CDS",        s, max(s, e - 3L), ".", "+", 0L, attrs("CDS")),
    list(seqname, src, "start_codon", s, s + 2L, ".", "+", 0L, attrs("start_codon")),
    list(seqname, src, "stop_codon",  e - 2L, e, ".", "+", 0L, attrs("stop_codon"))
  )
  dt <- rbindlist(lapply(rows, function(r) as.list(r)))
  setnames(dt, cols)
  dt[, `:=`(start = as.integer(start), end = as.integer(end))]
  dt
}

# ---- Insert edits in deterministic order (sorted by start) ----
if (length(edits) > 0) {
  # convert edits named list into a small table and sort by start position
  ed_df <- data.table(gene = names(edits),
                      start = as.integer(sapply(edits, `[`, 1)),
                      end   = as.integer(sapply(edits, `[`, 2)))
  setorder(ed_df, start)  # ensure deterministic insertion order
  
  # iterate over edits and insert GTF blocks at the right point
  for (i in seq_len(nrow(ed_df))) {
    g <- ed_df$gene[i]; s <- ed_df$start[i]; e <- ed_df$end[i]
    # find indices of transcript features to determine insertion index
    tx_idx <- which(gtf$feature == "transcript")
    # insert before the first transcript that starts AFTER s; otherwise append at end
    insert_at <- if (any(gtf$start[tx_idx] > s)) min(tx_idx[which(gtf$start[tx_idx] > s)]) else nrow(gtf) + 1L
    add <- build_rows(g, s, e)  # new records for this gene
    # three cases for rbind: insert at beginning, middle, or end
    if (insert_at == 1L) {
      gtf <- rbindlist(list(add, gtf), use.names = TRUE, fill = TRUE)
    } else if (insert_at > nrow(gtf)) {
      gtf <- rbindlist(list(gtf, add), use.names = TRUE, fill = TRUE)
    } else {
      gtf <- rbindlist(list(gtf[1:(insert_at - 1)], add, gtf[insert_at:.N]), use.names = TRUE, fill = TRUE)
    }
  }
}

# ---- Write edited GTF: tab-separated, no header, no quotes (match original format) ----
fwrite(gtf, out_gtf, sep = "\t", col.names = FALSE, quote = FALSE)
