#!/usr/bin/env python


import fire
from pathlib import Path, PurePath
import quilt3


def distribute_reproducibility(
    pkg_dest="aics/wtc11_hipsc_cardiomyocyte_scrnaseq_d0_to_d90",
    readme="README.md",
    s3_bucket="s3://allencell",
    edit=True,
):

    # either edit package if it exists or create new
    if edit:
        p = quilt3.Package.browse(pkg_dest, registry=s3_bucket)
    else:
        p = quilt3.Package()


    # package components
    # R object
    seurat_obj = PurePath("/allen/aics/gene-editing/Manuscripts/Transcriptomics_2019/scrnaseq_supplement/paper_plots/protocol_experiment_cell_line_comparison/both_protocols_seq_exp1_paper_figures.RData")
    p.set("reproducibility_D12_D14_D24_D26/seurat_object/reproducibility_D12_D14_D24_D26_paper_figures.RData", seurat_obj)

    # pairwise differential expression csvs
    de_dir = Path("/allen/aics/gene-editing/Manuscripts/Transcriptomics_2019/scrnaseq_supplement/paper_plots/protocol_experiment_cell_line_comparison/de_analysis_0.4/")

    for file_path in de_dir.glob("cluster_de_genes/*/*.csv"):
        de_file_name = PurePath(file_path).name

        p.set(f"reproducibility_D12_D14_D24_D26/cluster_differential_expression/pairwise_comparisons/{de_file_name}", file_path)

    print(p)
    p.push(pkg_dest, s3_bucket, message="reproducibility analysis")


if __name__ == "__main__":
    fire.Fire(distribute_reproducibility)
