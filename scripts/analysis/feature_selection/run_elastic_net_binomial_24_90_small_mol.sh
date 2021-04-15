#!/bin/bash

/allen/aics/gene-editing/RNA_seq/scRNAseq_SeeligCollaboration/2019_analysis/feature_selection/scripts/boelastic_net_binomial_24_90_small_mol.R \
    --input_exp "/allen/aics/gene-editing/RNA_seq/scRNAseq_SeeligCollaboration/2019_analysis/feature_selection/train_cm_entrez_only_normalized_counts_D24_D90.RData" \
    --metadata "/allen/aics/gene-editing/RNA_seq/scRNAseq_SeeligCollaboration/2019_analysis/feature_selection/train_cm_entrez_only_metadata_D24_D90.csv" \
    --lambda_seq "/allen/aics/gene-editing/RNA_seq/scRNAseq_SeeligCollaboration/2019_analysis/feature_selection/all_lambda_seq_log_D24_D90.RData" \
    --alpha 0.5 \
    --num_bootstrap 1000 \
    --replace FALSE \
    --subsample_fraction 0.8 \
    --out_dir "/allen/aics/gene-editing/RNA_seq/scRNAseq_SeeligCollaboration/2019_analysis/feature_selection/results/" \
    --analysis_name "D24_D90_small_mol"

