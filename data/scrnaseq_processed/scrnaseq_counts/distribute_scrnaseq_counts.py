#!/usr/bin/env python


import fire
from pathlib import Path, PurePath
import quilt3


def distribute_counts(
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

    # package components
    counts_matrix = "/allen/aics/gene-editing/Manuscripts/Transcriptomics_2019/scrnaseq_supplement/data_package/raw_counts.mtx"
    genes = "/allen/aics/gene-editing/Manuscripts/Transcriptomics_2019/scrnaseq_supplement/data_package/genes.csv"
    cells = "/allen/aics/gene-editing/Manuscripts/Transcriptomics_2019/scrnaseq_supplement/data_package/cells.csv"
    cell_metadata = "/allen/aics/gene-editing/Manuscripts/Transcriptomics_2019/scrnaseq_supplement/data_package/cell_metadata.csv"

    p.set("scrnaseq_counts/raw_counts.mtx", counts_matrix)
    p.set("scrnaseq_counts/genes.csv", genes)
    p.set("scrnaseq_counts/cells.csv", cells)
    p.set("scrnaseq_counts/cell_metadata.csv", cell_metadata)

    print(p)
    # p.push(pkg_dest, s3_bucket, message="counts")


if __name__ == "__main__":
    fire.Fire(distribute_counts)
