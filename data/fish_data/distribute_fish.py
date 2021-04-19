#!/usr/bin/env python


import fire
from pathlib import Path, PurePath
import pandas as pd
import quilt3


def distribute_fish(
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



    # cell profiler files
    cellprofiler_features = PurePath("/allen/aics/gene-editing/FISH/2019/chaos/data/cp_20200303/merged_features/cp_features.csv")
    cellprofiler_counts = PurePath("/allen/aics/gene-editing/FISH/2019/chaos/data/cp_20200303/merged_features/image_object_counts.csv")
    p.set("rnafish/cellprofiler_output/cp_feature.csv", cellprofiler_features)
    p.set("rnafish/cellprofiler_output/image_object_counts.csv", cellprofiler_counts)

    # image manifest
    image_manifest = "/allen/aics/gene-editing/FISH/2019/chaos/data/cp_testing_zeiss_nonstructure/zeiss_image_set/fov_metadata/nonstructure_fov_manifest_for_quilt.csv"

    raw_image_df = pd.read_csv(image_manifest)
    raw_fovs = [PurePath(f) for f in raw_image_df.fov_path]

    # 3D czi files
    for f in raw_fovs:
        p.set(f"rnafish/raw_images/{f.name}", f)

    # 2D segmentations 
    segmented_tiffs = [PurePath(s) for s in raw_image_df.merged_2D_fov_tiff_path]
    for s in segmented_tiffs:
        p.set(f"rnafish/segmented_2d_tiffs/{s.name}", s)

    # file manifest
    manifest = PurePath(image_manifest)
    p.set(f"rnafish/manifest.csv", manifest)

    print(p)
    # p.push(pkg_dest, s3_bucket, message="rnafish data")


if __name__ == "__main__":
    fire.Fire(distribute_fish)
