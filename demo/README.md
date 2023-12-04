# A NOTE ON THE INPUT DEMO DATA USED IN DNAICI

## FOREWORD

To make the demo dataset more lightweight, only multi-omics data from chromosomes 18 and 19 in both untreated MCF7 cells and tamoxifen-resistant MCF7TR cells are included. For the data containing 23 chromosomes from MCF7 cells with untreated and one hour of E2 treated, as well as tamoxifen-resistant MCF7TR cells, we recommend that you download and use the [FULL DATASET](https://drive.google.com/file/d/1YbdZ7y5bRNqbP_4hVt6rcZM2Om1PoA-b/view). 

ATTENTION: downloading and some steps of calculation might be time consuming, but you will get the same results as in our paper.

## hg19 data

**File**: `/hg19/hg19.chrom.sizes.clear.sorted`

**Format**:

```
chr1	249250621
chr2	243199373
...
```

**Description**: The file tells the names and sizes of each chromosome. It was downloaded [HERE](https://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/).

**File**: `/hg19/hg19-blacklist.bed`

**Format**:

```
chr1	564449	570371	High_Mappability_island	1000	.
chr1	724136	727043	Satellite_repeat	1000	.
...
```

**Description**: The black list file exclude a set of regions in the human genome that have anomalous, unstructured, high signal/read counts in next gen sequencing experiments independent of cell line and type of experiment. It was downloaded [HERE](https://www.encodeproject.org/annotations/ENCSR636HFF/).

**Usage**: Chromosome size file and black list file is required by module `DNAICI.preprocess_hic_homer2bed()` to make bed file of window bin based on Bedtools from Ref [[1]](https://academic.oup.com/bioinformatics/article/26/6/841/244688).

> **Note**: You can use more recent [hg38](https://www.nature.com/articles/s41598-019-45839-z) block list.


## Hi-C interaction

**Example file**: `/hic_data/hicup_processed/tamr/chr18.tags.tsv`

**Format**:

```
chr18	27971	1	0.5	44	chr13	25318728	0	38
chr18	33411	0	0.5	35	chr18	20874150	0	44
...
```

**Description**: Hi-C raw data in MCF7 cells and tamoxifen-resistant MCF7TR cells. Each row of it represents a inter- or intra- chromosomal interaction. The data were obtained from Ref [[2]](https://www.nature.com/articles/s41467-019-09320-9).

**Usage**: Hi-C raw data are the input to the module `DNAICI.preprocess_hic_tag2homer()`. Software HOMER from Ref [[3]](https://www.cell.com/molecular-cell/pdf/S1097-2765(10)00366-7.pdf) is needed to filter and reorganize significant intra-chromosomal Hi-C interactions.

> **Note**: folders 'untreated' and 'tamr' include Hi-C data in untreated MCF7 cells and tamoxifen-resistant MCF7TR cells.


## Gene expresion

**File**: `/demo/expression_data/mcf7_t0_geneExp.bed`

**Format**:

```
chrom	starts	ends	mcf7_0h_log2zscore	new_name	new_id	strand	ensemble_id	MCF-7_0h
chr1	65419	71585	-0.6984308889605118	OR4F5	79501	+	ENSG00000186092	290.678579
chr1	367640	368634	-0.8728869158218728	OR4F29	729759	+	ENSG00000284733	244.645389
...
```

**Description**: Microarray expression profiles of genes in zero and one hour of E2 treated MCF7 cells, as well as RNA-Seq data in tamoxifen-resistant MCF7TR cells. Each row tells the gene name and its zscore value in a certain region.

**Usage**: Gene expresion are the input to the module `DNAICI.preprocess_omics_map2hic()` when preprocessing multi-omics datasets.

**Note**: MicroRNAs of genes in zero and one hour of E2 treated MCF7 cells were from Ref [4](https://www.sciencedirect.com/science/article/pii/S0002944010600090), while RNA-Seq data in tamoxifen-resistant MCF7TR cells were obtained from Ref [5](https://www.nature.com/articles/s41556-020-0514-z).


## Nucleosome density

**File**: `/demo/nucleosome_density_data/mcf7_t0_DNas_200b.bed`

**Format**:

```
chr1	10400	10600	-0.891	0.0	0
chr1	13200	13400	-0.891	0.0	0
...
```

**Description**: DNase-seq and ATAC-Seq of nucleosome density in MCF7 and MCF7TR cells. Each row tells the zscore value in a certain region.

**Usage**: Nucleosome density are the input to the module `DNAICI.preprocess_omics_map2hic()` when preprocessing multi-omics datasets.

**Note**: DNase-seq of nucleosome density in zero and one hour of E2 treated MCF7 cells were from Ref [2](https://www.nature.com/articles/s41467-019-09320-9), while ATAC-Seq in tamoxifen-resistant MCF7TR cells were obtained from Ref [5](https://www.nature.com/articles/s41556-020-0514-z).


## Histone marker

**File**: `/demo/histone_data/mcf7_t0_ctcf_200b.bed`

**Format**:

```
chrom	starts	ends	TIME0_ctcf_log2zscore	TIME0_ctcf
chr1	521400	521600	1.410	5.43354
chr1	521600	521800	1.410	5.43354
...
```

**Description**: ChIP-Seq experiments of markers for enhancer (H3K27ac and H3K4me1), promoter (H3K4me3), repressor (H3K27me3 and H3K9me3), and insulator (CTCF) in MCF7 and MCF7TR cells. Each row tells the zscore value in a certain region.

**Usage**: Histone markers are the input to the module `DNAICI.preprocess_omics_map2hic()` when preprocessing multi-omics datasets.

**Note**: ChIP-Seq of histone markers in zero and one hour of E2 treated MCF7 cells were from Ref [[2]](https://www.nature.com/articles/s41467-019-09320-9), while ChIP-Seq in tamoxifen-resistant MCF7TR cells were obtained from Ref [[5]](https://www.nature.com/articles/s41556-020-0514-z).


## References

[1] Quinlan AR et al: BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics 2010, 26(6):841-842.

[2] Zhou Y et al: Temporal dynamic reorganization of 3D chromatin architecture in hormone-induced breast cancer and endocrine resistance. Nat Commun 2019, 10(1):1522.

[3] Heinz S et al: Simple Combinations of Lineage-Determining Transcription Factors Prime cis-Regulatory Elements Required for Macrophage and B Cell Identities. Molecular Cell 2010, 38(4):576-589.

[4] Cicatiello L et al: Estrogen Receptor alpha Controls a Gene Network in Luminal-Like Breast Cancer Cells Comprising Multiple Transcription Factors and MicroRNAs. American Journal of Pathology 2010, 176(5):2113-2130.

[5] Bi M et al: Enhancer reprogramming driven by high-order assemblies of transcription factors promotes phenotypic plasticity and breast cancer endocrine resistance. Nat Cell Biol 2020, 22(6):701-715.







