# Basal-specific peak clusters (One Example)

This directory contains the complete workflow for identifying cancer-specific regions from the pan-cancer peak set, clustering all
specific regions, ranking the clusters, and drawing the top IGV-like
figures.

## Download this repository

This repository uses Git LFS for `data_m.txt.gz`. Install Git and Git LFS first,
then clone the repository and download the LFS object:

```bash
git lfs install
git clone https://github.com/pengweixing/identify_specific_peaks_clusters.git
cd identify_specific_peaks_clusters
git lfs pull
```

After `git lfs pull`, `data_m.txt.gz` should be approximately 408 MB. If it is
only a small text pointer, Git LFS was not installed or the LFS download did not
finish; install Git LFS and run `git lfs pull` again.

## Run from beginning to end

The full count matrix and hg38 gene annotation are stored as gzip archives.
Decompress both before the first workflow run while keeping the compressed
archives:

```bash
gzip -dk data_m.txt.gz
gzip -dk hg38.refGene.gtf.gz
```

This creates `data_m.txt` and `hg38.refGene.gtf`, which are required by the
workflow. Then run:

```bash
./run_basal_specific_peak_clusters.sh
```

The script follows one fixed path:

1. Filter `data_m.txt` to the regions in the updated
   `Pancancer_peak.specific.bed`.
2. Merge technical replicate columns by sample barcode and arithmetic mean.
3. Calculate cancer-type CPM means and variances for cluster scoring.
4. Run DESeq2 for Basal, Her2, LumA, and LumB versus all other BRCA subtypes
   (`FDR <= 0.05`, `log2FC >= 1`, and `min CPM >= 1`).
5. Convert every row of `Basal_specific_regions.tsv` to BED. 
6. Join consecutive Basal-specific peaks separated by at most 5 kb and retain
   clusters containing at least three peaks.
7. Rank clusters using structural score, mean log2 fold change, FDR, mean CPM,
   and CPM stability.
8. Plot the top 20 clusters with 50 kb flanks.


## Main inputs and scripts

- `data_m.txt.gz`: compressed full TCGA ATAC count matrix. Run
  `gzip -dk data_m.txt.gz` to create the required `data_m.txt` input.
- `hg38.refGene.gtf.gz`: compressed hg38 RefGene annotation. Run
  `gzip -dk hg38.refGene.gtf.gz` to create the plotting input
  `hg38.refGene.gtf`.
- `Pancancer_peak.specific.bed`: input peak set by excluding all healthy peaks.
- `sample_cancer_groups.tsv`: sample-to-cancer/subtype assignments.
- `filter_data_m_by_pancancer_peaks.py`: extract specific-region matrix rows.
- `merge_replicates_by_mean.py`: merge technical replicates.
- `summarize_cpm_by_cancer_type.py`: calculate mean and variance CPM matrices.
- `brca_group_specific_differential_analysis.R`: BRCA subtype-versus-rest DESeq2.
- `extract_specific_regions_bed.py`: convert differential results to a standard
  three-column BED without occurrence metadata or filtering.
- `cluster_peaks_by_distance.py`: cluster nearby Basal-specific peaks.
- `rank_basal_specific_clusters.py`: rank the resulting clusters.
- `plot_basal_specific_cluster_igv.py`: generate IGV-like plots, including its
  own bigWig-reading and gene-annotation helper functions.

- `healthy.peak.bed`: derived from DOI: 10.1016/j.cell.2021.10.024

## External bigWig inputs (not included)

The `BRCA_subtype_random5_bigwigs/` and `selected_10/` directories are too large
to upload and are therefore **not included with this project**. The bigWig files
in both directories were downloaded from the NCI Genomic Data Commons TCGA
ATAC-seq publication page:

<https://gdc.cancer.gov/about-data/publications/ATACseq-AWG>

The plotting script expects the following directories beside
`run_basal_specific_peak_clusters.sh`:

```text
identify_specific_peaks_clusters/
├── BRCA_subtype_random5_bigwigs/
│   ├── Basal_1.bw
│   ├── Basal_2.bw
│   ├── ...
│   ├── Her2_1.bw
│   ├── LumA_1.bw
│   └── LumB_1.bw
└── selected_10/
    ├── ACC_1.bw
    ├── ACC_2.bw
    ├── ...
    ├── BLCA_1.bw
    ├── ...
    └── UCEC_10.bw
```

For `selected_10/`, ten samples were randomly selected for each cancer type
after download. Files were then renamed to `<cancer_type>_<number>.bw`, for
example `ACC_1.bw` through `ACC_10.bw`. The BRCA subtype tracks were similarly
renamed using their subtype labels, such as `Basal_1.bw`, `Her2_1.bw`,
`LumA_1.bw`, and `LumB_1.bw`.

Download, select, and rename these external bigWigs before running the plotting
step. The rest of the workflow can run without them, but step 8 will stop if
either directory is absent or contains no `.bw` files. The gene annotation is
provided separately as the compressed `hg38.refGene.gtf.gz` archive described
above.

## Main outputs

- `data_m.Pancancer_peak.specific.txt`
- `data_m.Pancancer_peak.specific.mean_replicates.txt`
- `data_m.Pancancer_peak.specific.mean_replicates.mean_cpm_by_cancer_type.tsv`
- `data_m.Pancancer_peak.specific.mean_replicates.var_cpm_by_cancer_type.tsv`
- `brca_subtype_vs_all_rest_deseq2/Basal_specific_regions.tsv`
- `brca_subtype_vs_all_rest_deseq2/Basal_specific_regions.bed`
- `brca_subtype_vs_all_rest_deseq2/Basal_specific_regions.distance_clusters.bed`
- `brca_subtype_vs_all_rest_deseq2/Basal_specific_regions.distance_clusters.ranked.bed`
- `Basal_specific_cluster_igv_plots_50kb/`

## Cluster ranking

Each metric is converted to a percentile rank. The weights are 0.40 structural
cluster score, 0.15 mean log2 fold change, 0.15 FDR significance, 0.25 mean CPM,
and 0.05 CPM stability; lower CPM variance is preferred.

BED intervals use zero-based, half-open coordinates. Each Basal bigWig track is
locally auto-scaled, while non-Basal tracks share a scale based on the average
local Basal maximum.
