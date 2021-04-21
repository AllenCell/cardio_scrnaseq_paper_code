#!/usr/bin/env python


import fire
from pathlib import Path, PurePath
import pandas as pd
import quilt3


def distribute_scrnaseq_raw(
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


    # raw data manifest
    raw_csv = "/allen/aics/gene-editing/Manuscripts/Transcriptomics_2019/scrnaseq_supplement/data_package/scrnaseq_quilt/scrnaseq_raw/scrnaseq_data_raw.csv"

    raw_df = pd.read_csv(raw_csv)
    fastq_files = [PurePath(f) for f in raw_df.fastq_files]

    # fastq files
    for f in fastq_files:
        # print(f)
        p.set(f"scrnaseq_raw/fastq_files/{f.name}", f)

    # read assignment files
    read_assignment_files = [PurePath(r) for r in raw_df.read_assignment_files.unique()]
    for r in read_assignment_files:
        # print(r)
        p.set(f"scrnaseq_raw/read_assignment/{r.name}", r)

    # file manifest
    manifest = PurePath("/allen/aics/gene-editing/Manuscripts/Transcriptomics_2019/scrnaseq_supplement/data_package/scrnaseq_quilt/scrnaseq_raw/scrnaseq_data_raw.csv")
    p.set(f"scrnaseq_raw/manifest.csv", manifest)

    print(p)
    p.push(pkg_dest, s3_bucket, message="scrnaseq raw data")


if __name__ == "__main__":
    fire.Fire(distribute_scrnaseq_raw)
