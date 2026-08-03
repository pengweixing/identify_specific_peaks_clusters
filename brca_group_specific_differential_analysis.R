#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
})

option_list <- list(
  matrix = "data_m.Pancancer_peak.specific.mean_replicates.txt",
  groups = "sample_cancer_groups.tsv",
  outdir = "brca_subtype_vs_all_rest_deseq2",
  fdr = 0.05,
  logfc = 1,
  min_cpm = 1
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) %% 2 != 0) {
  stop("Arguments must be key-value pairs, e.g. --fdr 0.05", call. = FALSE)
}
if (length(args) > 0) {
  keys <- sub("^--", "", args[seq(1, length(args), by = 2)])
  vals <- args[seq(2, length(args), by = 2)]
  for (i in seq_along(keys)) {
    key <- keys[[i]]
    if (!key %in% names(option_list)) {
      stop("Unknown argument: --", key, call. = FALSE)
    }
    current <- option_list[[key]]
    option_list[[key]] <- if (is.numeric(current)) as.numeric(vals[[i]]) else vals[[i]]
  }
}

matrix_file <- option_list$matrix
group_file <- option_list$groups
outdir <- option_list$outdir
fdr_cutoff <- option_list$fdr
logfc_cutoff <- option_list$logfc
min_cpm <- option_list$min_cpm

target_groups <- c("Basal", "Her2", "LumA", "LumB")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "target_vs_rest"), showWarnings = FALSE, recursive = TRUE)

message("Reading group table: ", group_file)
sample_groups <- fread(group_file)
required_group_cols <- c("sample_name", "cancer_group")
if (!all(required_group_cols %in% names(sample_groups))) {
  stop("sample_cancer_groups.tsv must contain: sample_name, cancer_group", call. = FALSE)
}

sample_groups <- sample_groups[!cancer_group %in% c("BRCA_NA", "Normal")]

message("Reading matrix: ", matrix_file)
mat_dt <- fread(matrix_file)
if (!"peak_id" %in% names(mat_dt)) {
  setnames(mat_dt, 1, "peak_id")
}

sample_cols <- intersect(sample_groups$sample_name, names(mat_dt))
missing_samples <- setdiff(sample_groups$sample_name, names(mat_dt))
if (length(missing_samples) > 0) {
  warning("Dropping samples not found in matrix: ", paste(missing_samples, collapse = ", "))
}
sample_groups <- sample_groups[match(sample_cols, sample_name)]

if (length(sample_cols) < 2) {
  stop("Not enough usable samples found in the matrix.", call. = FALSE)
}

group_counts <- table(sample_groups$cancer_group)
target_counts <- group_counts[target_groups]
target_counts[is.na(target_counts)] <- 0
if (any(target_counts == 0)) {
  stop(
    "All four BRCA groups must have at least one sample. Counts: ",
    paste(target_groups, as.integer(target_counts), sep = "=", collapse = ", "),
    call. = FALSE
  )
}

message("Samples used by group:")
message(paste(names(group_counts), as.integer(group_counts), sep = "=", collapse = ", "))
fwrite(sample_groups, file.path(outdir, "samples_used.tsv"), sep = "\t")

counts <- as.matrix(mat_dt[, ..sample_cols])
storage.mode(counts) <- "numeric"
counts <- round(counts)
counts[counts < 0] <- 0
storage.mode(counts) <- "integer"
rownames(counts) <- mat_dt$peak_id

lib_size <- colSums(counts)
cpm <- t(t(counts) / lib_size * 1e6)
min_samples <- min(as.integer(target_counts))
keep <- rowSums(cpm >= min_cpm) >= min_samples

message("Regions before CPM filter: ", nrow(counts))
message("Regions after CPM filter: ", sum(keep))
counts <- counts[keep, , drop = FALSE]
cpm <- cpm[keep, , drop = FALSE]
peak_ids <- rownames(counts)

sample_info <- data.frame(
  row.names = sample_groups$sample_name,
  cancer_group = sample_groups$cancer_group
)

size_dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = sample_info,
  design = ~ 1
)
size_dds <- estimateSizeFactors(size_dds)
norm_counts <- counts(size_dds, normalized = TRUE)
log_norm_counts <- log2(norm_counts + 1)

target_means <- sapply(target_groups, function(group) {
  rowMeans(log_norm_counts[, sample_groups$cancer_group == group, drop = FALSE])
})
target_means_dt <- data.table(peak_id = peak_ids, target_means)
setnames(target_means_dt, target_groups, paste0("mean_log2normcount_", target_groups))
fwrite(
  target_means_dt,
  file.path(outdir, "brca_target_group_mean_log2normcount.tsv"),
  sep = "\t"
)

all_specific <- list()
summary_rows <- list()

for (target in target_groups) {
  message("Running DESeq2: ", target, " vs all rest")

  condition <- ifelse(sample_groups$cancer_group == target, "target", "rest")
  col_data <- data.frame(
    row.names = sample_groups$sample_name,
    condition = factor(condition, levels = c("rest", "target"))
  )

  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = col_data,
    design = ~ condition
  )
  sizeFactors(dds) <- sizeFactors(size_dds)
  dds <- DESeq(dds, quiet = TRUE)

  res <- results(
    dds,
    contrast = c("condition", "target", "rest"),
    alpha = fdr_cutoff
  )
  res_dt <- as.data.table(as.data.frame(res), keep.rownames = "peak_id")
  setnames(
    res_dt,
    c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"),
    c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
  )
  setorder(res_dt, padj, -log2FoldChange, na.last = TRUE)
  fwrite(
    res_dt,
    file.path(outdir, "target_vs_rest", paste0(target, "_vs_all_rest.tsv")),
    sep = "\t"
  )

  rest_mean <- rowMeans(log_norm_counts[, sample_groups$cancer_group != target, drop = FALSE])
  target_mean <- target_means_dt[[paste0("mean_log2normcount_", target)]]

  specific <- data.table(
    peak_id = peak_ids,
    target_group = target,
    log2FoldChange_vs_all_rest = res_dt$log2FoldChange[match(peak_ids, res_dt$peak_id)],
    padj_vs_all_rest = res_dt$padj[match(peak_ids, res_dt$peak_id)],
    mean_log2normcount_all_rest = rest_mean
  )
  mean_cols <- paste0("mean_log2normcount_", target_groups)
  specific <- cbind(specific, target_means_dt[, ..mean_cols])
  specific <- specific[
    !is.na(padj_vs_all_rest) &
      log2FoldChange_vs_all_rest >= logfc_cutoff &
      padj_vs_all_rest <= fdr_cutoff &
      target_mean > mean_log2normcount_all_rest
  ]
  setorder(specific, -log2FoldChange_vs_all_rest, padj_vs_all_rest)

  fwrite(
    specific,
    file.path(outdir, paste0(target, "_specific_regions.tsv")),
    sep = "\t"
  )
  all_specific[[target]] <- specific
  summary_rows[[target]] <- data.table(
    target_group = target,
    target_samples = as.integer(group_counts[[target]]),
    rest_samples = length(sample_cols) - as.integer(group_counts[[target]]),
    specific_regions = nrow(specific)
  )
}

summary_dt <- rbindlist(summary_rows)
fwrite(summary_dt, file.path(outdir, "summary.tsv"), sep = "\t")

combined <- rbindlist(all_specific, use.names = TRUE, fill = TRUE)
fwrite(combined, file.path(outdir, "all_group_specific_regions.tsv"), sep = "\t")

bed <- copy(combined[, .(peak_id, target_group)])
if (nrow(bed) > 0) {
  bed[, c("chrom", "coords") := tstrsplit(peak_id, ":", fixed = TRUE)]
  bed[, c("start", "end") := tstrsplit(coords, "-", fixed = TRUE)]
  bed[, coords := NULL]
  bed <- bed[, .(chrom, start, end, peak_id, target_group)]
}
fwrite(
  bed,
  file.path(outdir, "all_group_specific_regions.bed"),
  sep = "\t",
  col.names = FALSE
)

message("Done. Results written to: ", outdir)
message("Cutoffs: padj <= ", fdr_cutoff, ", log2FoldChange >= ", logfc_cutoff, ", min CPM >= ", min_cpm)
message("Note: replicate-mean counts were rounded to integers for DESeq2.")
