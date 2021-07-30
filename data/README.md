# Data set for _A comprehensive analysis of gene expression changes in a high replicate and open-source dataset of differentiating hiPSC-derived cardiomyocytes_

This data package contains input data for analyses in the manuscript [_A comprehensive analysis of gene expression changes in a high replicate and open-source dataset of differentiating hiPSC-derived cardiomyocytes_](https://www.biorxiv.org/content/10.1101/2021.04.22.441027v1).

**Overview:** This dataset contains single cell RNA-sequencing data from in vitro hiPSC derived cardiomyocytes
and includes both raw (fastq) and processed data (count matrix).

## Cardiomyocyte differentiation and sample collection 
hiPSCs were differentiated into cardiomyocytes using two protocols: small molecule (Lian et al. 2012) and 
cytokine (Palpant et al. 2015). Samples were collected for sequencing at 4 time points: 1) D0 time point
before the initiation of differentiation, 2) D12/D14 after differentiation was initiated, 3) D24/D26, and 4) D90.

## Library preparation and sequencing
Single cell libraries were prepared using the Split-Seq method. (see manuscript for more details)

## Contents 
- `scrnaseq_raw`: scRNA-seq raw fastq and read assignment files
- `scrnaseq_counts`: scRNA-seq cell by gene count matrix
- `protocol1_D0_D12_D24_D90`: analysis of D0 stem cells and hiPSCs differentiated into cardiomyocytes using Protocol 1 (see manuscript); samples collected at D12, D24, and D90 after initiation of differentiation
- `reproducibility_D12_D14_D24_D26`: analysis of hiPSCs differentiated into cardiomyocytes using either Protocol 1 (D12 and D24) or Protocol 2 (D14 and D26) (see manuscript)
- `feature_selection`: inputs and outputs of bootstrapped sparse regression analysis to identify and rank genes that best predict cell time point
- `rnafish`: raw and segmented RNA-FISH images and transcript quantification

## Access
The data are programmatically accessible via `quilt`.

### Bulk download
To download the entire data set, install the `quilt` python package using
```bash
pip install quilt
```
and then
```python
import quilt3
b = quilt3.Bucket("s3://allencell")
b.fetch("aics/wtc11_hipsc_cardiomyocyte_scrnaseq_d0_to_d90/", "./")
```

### Download specific files or data sets
To download only certain individual files, navigate the web ui here to the specific file you are interested in, and use the `DOWNLOAD FILE` button in the upper right of the page.

To download specific folders/directories of data, similarly use the web ui to find the directory you want, and check the `<> CODE` tab at the top of the page for the python code that downloads that specific subset of data.

### Programmatic access
To access the data via the python quilt API, install `quilt` via `pip`, and then load the package with:

```python
pkg = quilt3.Package.browse(
    "aics/wtc11_hipsc_cardiomyocyte_scrnaseq_d0_to_d90",
    "s3://allencell",
)
```
Instructions for interacting with quilt packages in Python can be found [here](https://docs.quiltdata.com/walkthrough/getting-data-from-a-package).

## Citation
```
@article {Grancharova2021.04.22.441027,
    author = {Grancharova, Tanya and Gerbin, Kaytlyn A and Rosenberg, Alexander B and Roco, Charles M and Arakaki, Joy and DeLizzo, Colette and Dinh, Stephanie Q and Donovan-Maiye, Rory and Hirano, Matthew and Nelson, Angelique and Tang, Joyce and Theriot, Julie A and Yan, Calysta and Menon, Vilas and Palecek, Sean P and Seelig, Georg and Gunawardane, Ruwanthi N},
    title = {A comprehensive analysis of gene expression changes in a high replicate and open-source dataset of differentiating hiPSC-derived cardiomyocytes},
    elocation-id = {2021.04.22.441027},
    year = {2021},
    doi = {10.1101/2021.04.22.441027},
    publisher = {Cold Spring Harbor Laboratory},
    URL = {https://www.biorxiv.org/content/early/2021/04/23/2021.04.22.441027},
    eprint = {https://www.biorxiv.org/content/early/2021/04/23/2021.04.22.441027.full.pdf},
    journal = {bioRxiv}
}
```

## License
For questions on licensing please refer to https://www.allencell.org/terms-of-use.html.

## Contact
Allen Institute for Cell Science E-mail: cells@alleninstitute.org

## Feedback
Feedback on benefits and issues you discovered while using this data package is greatly appreciated. [Feedback Form](https://forms.gle/GUBC3zU5kuA8wyS17)

