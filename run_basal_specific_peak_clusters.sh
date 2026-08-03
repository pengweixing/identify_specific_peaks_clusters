#!/usr/bin/env bash
set -euo pipefail

# Run from any directory; all relative paths are resolved from this script's folder.
workflow_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
cd "$workflow_dir"

raw_matrix="$workflow_dir/data_m.txt"
specific_matrix="$workflow_dir/data_m.Pancancer_peak.specific.txt"
merged_matrix="$workflow_dir/data_m.Pancancer_peak.specific.mean_replicates.txt"
sample_groups="$workflow_dir/sample_cancer_groups.tsv"
gene_annotation="$workflow_dir/hg38.refGene.gtf"
rscript="/home/pengwei/software/miniforge3/envs/py312/bin/Rscript"

if [[ ! -f "$raw_matrix" ]]; then
  echo "Missing input: $raw_matrix" >&2
  echo "Decompress it first with: gzip -dk data_m.txt.gz" >&2
  exit 1
fi

if [[ ! -f "$gene_annotation" ]]; then
  echo "Missing input: $gene_annotation" >&2
  echo "Decompress it first with: gzip -dk hg38.refGene.gtf.gz" >&2
  exit 1
fi

echo "[1/8] Filtering the raw matrix with the updated pan-cancer peak BED"
python3 filter_data_m_by_pancancer_peaks.py \
  --bed Pancancer_peak.specific.bed \
  --matrix "$raw_matrix" \
  --output "$specific_matrix"

echo "[2/8] Merging technical replicates"
python3 merge_replicates_by_mean.py \
  --input "$specific_matrix" \
  --output "$merged_matrix"

echo "[3/8] Calculating cancer-type CPM means and variances"
python3 summarize_cpm_by_cancer_type.py \
  --matrix "$merged_matrix" \
  --sample-groups "$sample_groups" \
  --mean-output data_m.Pancancer_peak.specific.mean_replicates.mean_cpm_by_cancer_type.tsv \
  --variance-output data_m.Pancancer_peak.specific.mean_replicates.var_cpm_by_cancer_type.tsv

echo "[4/8] Identifying BRCA subtype-specific peaks with DESeq2"
"$rscript" brca_group_specific_differential_analysis.R \
  --matrix "$merged_matrix" \
  --groups "$sample_groups" \
  --outdir brca_subtype_vs_all_rest_deseq2 \
  --fdr 0.05 \
  --logfc 1 \
  --min_cpm 1

echo "[5/8] Extracting all Basal-specific regions as BED"
python3 extract_specific_regions_bed.py \
  --specific-tsv brca_subtype_vs_all_rest_deseq2/Basal_specific_regions.tsv \
  --output brca_subtype_vs_all_rest_deseq2/Basal_specific_regions.bed

echo "[6/8] Clustering all Basal-specific peaks"
python3 cluster_peaks_by_distance.py \
  brca_subtype_vs_all_rest_deseq2/Basal_specific_regions.bed \
  brca_subtype_vs_all_rest_deseq2/Basal_specific_regions.distance_clusters.bed \
  --max-gap 5000

echo "[7/8] Ranking Basal-specific peak clusters"
python3 rank_basal_specific_clusters.py \
  --input brca_subtype_vs_all_rest_deseq2/Basal_specific_regions.distance_clusters.bed \
  --output brca_subtype_vs_all_rest_deseq2/Basal_specific_regions.distance_clusters.ranked.bed \
  --score-weight 0.40 \
  --log2fc-weight 0.15 \
  --fdr-weight 0.15 \
  --mean-cpm-weight 0.25 \
  --var-cpm-weight 0.05 \
  --variance-preference low

echo "[8/8] Plotting the top 20 clusters with 50 kb flanks"
mkdir -p Basal_specific_cluster_igv_plots_50kb
find Basal_specific_cluster_igv_plots_50kb -maxdepth 1 -type f -name '*.png' -delete
python3 plot_basal_specific_cluster_igv.py \
  --cluster-file brca_subtype_vs_all_rest_deseq2/Basal_specific_regions.distance_clusters.ranked.bed \
  --top-n 20 \
  --outdir Basal_specific_cluster_igv_plots_50kb \
  --flank 50000

echo "Workflow complete: $workflow_dir/Basal_specific_cluster_igv_plots_50kb"
