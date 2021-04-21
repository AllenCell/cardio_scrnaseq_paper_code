#!/usr/bin/env python


import fire
from pathlib import Path, PurePath
import quilt3


def distribute_protocol1(
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
    seurat_obj = PurePath("/allen/aics/gene-editing/Manuscripts/Transcriptomics_2019/scrnaseq_supplement/paper_plots/D0_D12_D24_D90_paper_figures.RData")
    p.set("protocol1_D0_D12_D24_D90/seurat_object/D0_D12_D24_D90_paper_figures.RData", seurat_obj)

    # pairwise differential expression csvs
    de_dir = Path("/allen/aics/gene-editing/Manuscripts/Transcriptomics_2019/scrnaseq_supplement/paper_plots/small_molecule_cardio_timepoints_indepth/de_analysis_calcnorm/")

    for file_path in de_dir.glob("cluster_de_genes/*/*.csv"):
        de_file_name = PurePath(file_path).name

        p.set(f"protocol1_D0_D12_D24_D90/cluster_differential_expression/pairwise_comparisons/{de_file_name}", file_path)


    # go term csvs
    for x in de_dir.glob("*_go_bp*csv"):
        go_file_name = PurePath(x).name
        p.set(f"protocol1_D0_D12_D24_D90/cluster_differential_expression/enriched_go_terms/{go_file_name}", x)

    print(p)
    p.push(pkg_dest, s3_bucket, message="protocol1 analysis")


if __name__ == "__main__":
    fire.Fire(distribute_protocol1)
