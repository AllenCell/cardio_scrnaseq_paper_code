#!/bin/bash

/allen/aics/gene-editing/RNA_seq/scRNAseq_SeeligCollaboration/2019_analysis/feature_selection/scripts/boelastic_net_binomial_early_late_small_mol.R \
    --input_exp "/allen/aics/gene-editing/RNA_seq/scRNAseq_SeeligCollaboration/2019_analysis/feature_selection/train_cm_entrez_only_normalized_counts.RData" \
    --metadata "/allen/aics/gene-editing/RNA_seq/scRNAseq_SeeligCollaboration/2019_analysis/feature_selection/train_cm_entrez_only_metadata.csv" \
    --lambda_seq "/allen/aics/gene-editing/RNA_seq/scRNAseq_SeeligCollaboration/2019_analysis/feature_selection/all_lambda_seq_log.RData" \
    --alpha 0.5 \
    --num_bootstrap 1000 \
    --replace FALSE \
    --subsample_fraction 0.8 \
    --out_dir "/allen/aics/gene-editing/RNA_seq/scRNAseq_SeeligCollaboration/2019_analysis/feature_selection/results/" \
    --analysis_name "early_late_small_mol"

