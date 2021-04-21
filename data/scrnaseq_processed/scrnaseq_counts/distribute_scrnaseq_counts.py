#!/usr/bin/env python


import fire
from pathlib import Path, PurePath
import quilt3


def distribute_counts(
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
    counts_matrix = PurePath("/allen/aics/gene-editing/Manuscripts/Transcriptomics_2019/scrnaseq_supplement/data_package/raw_counts.mtx")
    genes = PurePath("/allen/aics/gene-editing/Manuscripts/Transcriptomics_2019/scrnaseq_supplement/data_package/genes.csv")
    cells = PurePath("/allen/aics/gene-editing/Manuscripts/Transcriptomics_2019/scrnaseq_supplement/data_package/cells.csv")
    cell_metadata = PurePath("/allen/aics/gene-editing/Manuscripts/Transcriptomics_2019/scrnaseq_supplement/data_package/cell_metadata.csv")

    p.set(f"scrnaseq_counts/{counts_matrix.name}", counts_matrix)
    p.set(f"scrnaseq_counts/{genes.name}", genes)
    p.set(f"scrnaseq_counts/{cells.name}", cells)
    p.set(f"scrnaseq_counts/{cell_metadata.name}", cell_metadata)

    print(p)
    p.push(pkg_dest, s3_bucket, message="counts")


if __name__ == "__main__":
    fire.Fire(distribute_counts)
