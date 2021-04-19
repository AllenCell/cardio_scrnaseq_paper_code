#!/usr/bin/env python


import fire
from pathlib import Path, PurePath
import quilt3


def distribute_feature_selection(
    pkg_dest="aics/wtc11_hipsc_cardiomyocyte_scrnaseq_d0_to_d90",
    readme="README.md",
    s3_bucket="s3://allencell",
    edit=False,
):

    # either edit package if it exists or create new
    if edit:
        p = quilt3.Package.browse(pkg_dest, registry=s3_bucket)
    else:
        p = quilt3.Package()

    feat_selection_dir = "/allen/aics/gene-editing/RNA_seq/scRNAseq_SeeligCollaboration/2019_analysis/feature_selection/"
    plot_csvs_dir = "/allen/aics/gene-editing/Manuscripts/Transcriptomics_2019/scrnaseq_supplement/paper_plots/feature_selection/"

    # package components
    # training and holdout data
    training_d12_24_14_26 = PurePath(f"{feat_selection_dir}/train_cm_entrez_only_normalized_counts.RData")

    p.set(f"feature_selection/training_holdout/{training_d12_24_14_26.name}", training_d12_24_14_26)

    training_d12_24_14_26_cell_metadata = PurePath(f"{feat_selection_dir}/train_cm_entrez_only_metadata.csv")
    p.set(f"feature_selection/training_holdout/{training_d12_24_14_26_cell_metadata.name}", training_d12_24_14_26_cell_metadata)

    training_d24_90 = PurePath(f"{feat_selection_dir}/train_cm_entrez_only_normalized_counts_D24_D90.RData")
    p.set(f"feature_selection/training_holdout/{training_d24_90.name}", training_d24_90)

    training_d24_90_cell_metadata = PurePath(f"{feat_selection_dir}/train_cm_entrez_only_metadata_D24_D90.csv")
    p.set(f"feature_selection/training_holdout/{training_d24_90_cell_metadata.name}", training_d24_90_cell_metadata)

    holdout_d12_24_14_26 = PurePath(f"{feat_selection_dir}/holdout_cm_entrez_only_normalized_counts.RData")
    p.set(f"feature_selection/training_holdout/{holdout_d12_24_14_26.name}", holdout_d12_24_14_26)

    holdout_d12_24_14_26_cell_metadata = PurePath(f"{feat_selection_dir}/holdout_cm_entrez_only_metadata.csv")
    p.set(f"feature_selection/training_holdout/{holdout_d12_24_14_26_cell_metadata.name}", holdout_d12_24_14_26_cell_metadata)

    holdout_d24_90 = PurePath(f"{feat_selection_dir}/holdout_cm_entrez_only_normalized_counts_D24_D90.RData")
    p.set(f"feature_selection/training_holdout/{holdout_d24_90.name}", holdout_d24_90)

    holdout_d24_90_cell_metadata = PurePath(f"{feat_selection_dir}/holdout_cm_entrez_only_metadata_D24_D90.csv")
    p.set(f"feature_selection/training_holdout/{holdout_d24_90_cell_metadata.name}", holdout_d24_90_cell_metadata)

    # lambda sequences
    lambda_seq_d12_24_14_26 = PurePath(f"{feat_selection_dir}/all_lambda_seq_log.RData")
    p.set(f"feature_selection/lambda_sequence/{lambda_seq_d12_24_14_26.name}", lambda_seq_d12_24_14_26)

    lambda_seq_d24_90 = PurePath(f"{feat_selection_dir}/all_lambda_seq_log_D24_D90.RData")
    p.set(f"feature_selection/lambda_sequence/{lambda_seq_d24_90.name}", lambda_seq_d24_90)

    # results
    protocol1_d12v24 = PurePath(f"{feat_selection_dir}/results/early_late_small_mol.RData")
    protocol1_d24v90 = PurePath(f"{feat_selection_dir}/results/D24_D90_small_mol.RData")
    protocol2_d14v26 = PurePath(f"{feat_selection_dir}/results/early_late_cytokine.RData")

    p.set(f"feature_selection/bootstrap_results/{protocol1_d12v24.name}", protocol1_d12v24)
    p.set(f"feature_selection/bootstrap_results/{protocol1_d24v90.name}", protocol1_d24v90)
    p.set(f"feature_selection/bootstrap_results/{protocol2_d14v26.name}", protocol2_d14v26)

    # prediction accuracy and other csv files used for figure making

    # selected gene groups
    selected_groups_d12v24 = PurePath(f"{plot_csvs_dir}/selected_gene_groups.csv")
    p.set("feature_selection/figure_csvs/selected_gene_groups_D12_D24.csv", selected_groups_d12v24)

    selected_groups_d24v90 = PurePath(f"{plot_csvs_dir}/selected_gene_groups_D24_D90.csv")
    p.set("feature_selection/figure_csvs/selected_gene_groups_D24_D90.csv", selected_groups_d24v90)

    # D12 vs D24 selected and random gene accuracy of day prediction in holdout cells
    accuracy_selected_genes_d12v24 = PurePath(f"{plot_csvs_dir}/lambda_sequence_stats_glmnet_ALL_20200123.csv")
    p.set("feature_selection/figure_csvs/selected_gene_accuracy_D12_D24.csv", accuracy_selected_genes_d12v24)

    accuracy_random_genes_d12v24 = PurePath(f"{plot_csvs_dir}/random_gene_stats_glmnet_20200123.csv")
    p.set("feature_selection/figure_csvs/random_gene_accuracy_D12_D24.csv", accuracy_random_genes_d12v24)

    # D24 vs. D90 selected and random gene accuracy of day prediction in holdout cells
    accuracy_selected_genes_d24v90 = PurePath(f"{plot_csvs_dir}/lambda_sequence_stats_glmnet_D24_D90_ALL_20210201.csv")
    p.set("feature_selection/figure_csvs/selected_gene_accuracy_D24_90.csv", accuracy_selected_genes_d24v90)

    accuracy_random_genes_d24v90 = PurePath(f"{plot_csvs_dir}/random_gene_stats_glmnet_D24_D90_ALL_20210201.csv")
    p.set("feature_selection/figure_csvs/random_gene_accuracy_D24_90.csv", accuracy_random_genes_d24v90)

    # D12 vs D24 gene selection frequency in bootstraps
    bootstrap_selection_frequency_d12v24 = PurePath(f"{plot_csvs_dir}/feature_selection_frequency.RData")
    p.set("feature_selection/figure_csvs/bootstrap_selection_frequency_D12_D24.csv", bootstrap_selection_frequency_d12v24)


    # D24 vs D90 gene selection frequency in bootstraps
    bootstrap_selection_frequency_d24v90 = PurePath(f"{plot_csvs_dir}/feature_selection_frequency_d24_d90.RData")
    p.set("feature_selection/figure_csvs/bootstrap_selection_frequency_D24_D90.csv", bootstrap_selection_frequency_d24v90)


    print(p)
    # p.push(pkg_dest, s3_bucket, message="feature selection")


if __name__ == "__main__":
    fire.Fire(distribute_feature_selection)
