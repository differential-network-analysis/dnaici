# A NOTE ON THE INPUT DEMO DATA USED IN DNAICI

## FOREWORD

To make the demo dataset more lightweight, only multi-omics data from chromosomes **18** and **19** in both untreated MCF7 cells and tamoxifen-resistant MCF7TR cells are included. For the data containing 23 chromosomes from MCF7 cells with untreated and one hour of E2 treated, as well as tamoxifen-resistant MCF7TR cells, we recommend that you download and use the [FULL DATASET](https://drive.google.com/file/d/1YbdZ7y5bRNqbP_4hVt6rcZM2Om1PoA-b/view). 

The file names for data from untreated MCF7 cells contain `t0` or `untreated`, while the data from one hour of E2 treated MCF7 cells and tamoxifen-resistant MCF7TR cells have file names containing `t1` and `tamr`, respectively.

ATTENTION: For full dataset, downloading and some steps of calculation might be time consuming, but you will get the same results as in our paper.

## hg19 data

**File**: `/hg19/hg19.chrom.sizes.clear.sorted`

**Format**:

```
...
chr18	78077248
chr19	59128983
...
```

**Description**: The file tells the names and sizes of each chromosome. It was downloaded [HERE](https://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/).

**File**: `/hg19/hg19-blacklist.bed`

**Format**:

```
...
chr18	96416	97552	Satellite_repeat	1000	.
chr18	105658	112233	Satellite_repeat	1000	.
...
```

**Description**: The black list file exclude a set of regions in the human genome that have anomalous, unstructured, high signal/read counts in next gen sequencing experiments independent of cell line and type of experiment. It was downloaded [HERE](https://www.encodeproject.org/annotations/ENCSR636HFF/).

**Usage**: Chromosome size file and black list file is required by module `DNAICI.preprocess_hic_homer2bed()` to make bed file of window bin based on Bedtools from Ref [[1]](https://academic.oup.com/bioinformatics/article/26/6/841/244688).

> **Note**: More recent hg38 block list can be found [HERE](https://www.nature.com/articles/s41598-019-45839-z).


## Hi-C interaction

**Example file**: `/hic_data/hicup_processed/tamr/chr18.tags.tsv`

**Format**:

```
...
chr18	27971	1	0.5	44	chr13	25318728	0	38
chr18	33411	0	0.5	35	chr18	20874150	0	44
...
```

**Description**: Hi-C raw data in MCF7 cells and tamoxifen-resistant MCF7TR cells. Each row of it represents a inter- or intra- chromosomal interaction. The data were obtained from Ref [[2]](https://www.nature.com/articles/s41467-019-09320-9).

**Usage**: Hi-C raw data are the input to the module `DNAICI.preprocess_hic_tag2homer()`. Software HOMER from Ref [[3]](https://www.cell.com/molecular-cell/pdf/S1097-2765(10)00366-7.pdf) is needed to filter and reorganize significant intra-chromosomal Hi-C interactions.

> **Note**: Folders 'untreated' and 'tamr' include Hi-C data in untreated MCF7 cells and tamoxifen-resistant MCF7TR cells, respectively. Due to size limitation, chr18.tags.tsv and chr19.tags.tsv are compressed and uploaded. Users need to unzip these two .gz files to calculate.


## Gene expresion

**Example file**: `/expression_data/mcf7_tamr_geneExp.bed`

**Format**:

```
...
chr18  58038295  58040008  -1.11102465870699  MC4R  4160  -  ENSG00000166603  0.0425870323229458
chr18  59000663  59223012  -1.410928086530052  CDH20  28316  +  ENSG00000101542  0.0178982516039949
...
```

**Description**: Microarray expression profiles of genes in untreated MCF7 cells and tamoxifen-resistant MCF7TR cells. Each row tells the gene name and its zscore value in a certain region. MicroRNAs of genes in zero and one hour of E2 treated MCF7 cells were from Ref [[4]](https://www.sciencedirect.com/science/article/pii/S0002944010600090), while RNA-Seq data in tamoxifen-resistant MCF7TR cells were obtained from Ref [[5]](https://www.nature.com/articles/s41556-020-0514-z).

**Usage**: Gene expresion are the input to the module `DNAICI.preprocess_omics_map2hic()` when preprocessing multi-omics datasets.

> **Note**: There is another file in this folder: `/expression_data/rna_count.txt`. It recorded the count value of RNA-seq from two repeat observations of MCF7 cells and MCF7TR cells. This file is the input to the dnaici.diegs_enrichment() function, which is used to select differential expressed genes.


## Nucleosome density

**Example file**: `/nucleosome_density_data/mcf7_untreated_DNas_200b.bed`

**Format**:

```
...
chr18	11800	12000	1.266	7.700865	7.57475
chr18	12000	12200	1.266	7.700865	7.57475
...
```

**Description**: DNase-seq and ATAC-Seq of nucleosome density in MCF7 and MCF7TR cells. Each row tells the zscore value in a certain region. DNase-seq of nucleosome density in zero and one hour of E2 treated MCF7 cells were from Ref [[2]](https://www.nature.com/articles/s41467-019-09320-9), while ATAC-Seq in tamoxifen-resistant MCF7TR cells were obtained from Ref [[5]](https://www.nature.com/articles/s41556-020-0514-z).

**Usage**: Nucleosome density are the input to the module `DNAICI.preprocess_omics_map2hic()` when preprocessing multi-omics datasets.


## Histone marker

**Example file**: `/histone_data/mcf7_tamr_ctcf_200b.bed`

**Format**:

```
...
chr18	12200	12400	-0.273	0.0
chr18	12400	12600	-0.273	0.0
...
```

**Description**: ChIP-Seq experiments of markers for enhancer (H3K27ac and H3K4me1), promoter (H3K4me3), repressor (H3K27me3 and H3K9me3), and insulator (CTCF) in MCF7 and MCF7TR cells. Each row tells the zscore value in a certain region. ChIP-Seq of histone markers in zero and one hour of E2 treated MCF7 cells were from Ref [[2]](https://www.nature.com/articles/s41467-019-09320-9), while ChIP-Seq in tamoxifen-resistant MCF7TR cells were obtained from Ref [[5]](https://www.nature.com/articles/s41556-020-0514-z).

**Usage**: Histone markers are the input to the module `DNAICI.preprocess_omics_map2hic()` when preprocessing multi-omics datasets.


## References

[1] Quinlan AR et al: BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics 2010, 26(6):841-842.

[2] Zhou Y et al: Temporal dynamic reorganization of 3D chromatin architecture in hormone-induced breast cancer and endocrine resistance. Nat Commun 2019, 10(1):1522.

[3] Heinz S et al: Simple Combinations of Lineage-Determining Transcription Factors Prime cis-Regulatory Elements Required for Macrophage and B Cell Identities. Molecular Cell 2010, 38(4):576-589.

[4] Cicatiello L et al: Estrogen Receptor alpha Controls a Gene Network in Luminal-Like Breast Cancer Cells Comprising Multiple Transcription Factors and MicroRNAs. American Journal of Pathology 2010, 176(5):2113-2130.

[5] Bi M et al: Enhancer reprogramming driven by high-order assemblies of transcription factors promotes phenotypic plasticity and breast cancer endocrine resistance. Nat Cell Biol 2020, 22(6):701-715.







