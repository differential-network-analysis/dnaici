# README for DNAICI USERS

## DNAICI: Differential Network Analysis in Intra-chromosomal Community Interaction

## FOREWORD

The README file provides instructions for analyzing [demo data](https://github.com/differential-network-analysis/dnaici/tree/master/demo). The default values for parameters (e.g. `cohort = 'untreated'`, `chromosome = ['chr18','chr19']`) are designed for demo data. For more details about demo data, please refer to [README.md](https://github.com/differential-network-analysis/dnaici/blob/master/demo/README.md).

The DNAICI pipeline consists of three main mudules
* `preprocess`: preprocess_hic, preprocess_omics
* `analysis`: cluster, enrichment, network, dieg
* `parameter`: estimate

The module `preprocess` contains two parts for the preprocessing of raw Hi-C data and other omics data. Four parts in the `analysis` module are used for integrated analysis of multi-omics data. The `parameter` module are used to estimate the optimal values of two tuning parameters. All the above modules depend on two support modules: `tools` and `bin`.

To obtain consistent results with those in our paper, we recommend that you download and use the [full data](https://drive.google.com/file/d/1YbdZ7y5bRNqbP_4hVt6rcZM2Om1PoA-b/view?usp=drive_link). Necessary adjustments to the parameters need to be made for analysis. Full data analysis is time consuming, please be careful and patient. 

## STEP 0. Load DNAICI package and initialize global parameters

HOMER and Bedtools are required by data preprocessing. Java environment and ModularityOptimizer.jar are required by clustering algorithm.

```python
from dnaici import dnaici
dna = dnaici.DNAICI(in_data_folder, out_data_folder, cohort = 'untreated', chromosome = ['chr18', 'chr19'], resolution = 500000)
```

> **Parameters:**
>
> ***in_data_folder***: Input data directory where Hi-C, gene expression, nucleosome density, histone marker, and chromosome regions hg19 are deposited.
>
> ***out_data_folder***: Output data directory where processed Hi-C, gene expression, nucleosome density, histone marker and other analysis results are deposited.
> 
> ***cohort***: The experimental condition for input data. The default is `'untreated'`, representing untreated MCF7 cell line. It can be changed to `'tamr'` for demo data, or `'t0'` and `'t1'` for full data.
>
> ***chromosome***: List of chromosomes you want to investigate. The default is `['chr18', 'chr19']`. If you want to get the complete results for 23 chromosomes, please let the input be `'whole_genome'`.
>
> ***resolution***: The resolution is used to divide the chromosome into equally sized window bin for analysis. The default is `500000`, representing 500 kilobase pair.

**Note**: Currently, DNAICI cannot be downloaded using `pip install dnaici`.


## STEP 1.1 Preprocess Hi-C data

The section is designed for data preprocessing of raw omics data. Firstly, applying HOMER to raw Hi-C data (HICUP) to filter and reorganize significant intra-chromosomal interactions. 

```python
dna.preprocess_hic_tag2homer(super_resolution, p_value = 0.1, zscore = 1.0)
```
> **Parameters:**
>
> ***super_resolution***: The range of signal averaging. In the manual of HOMER, it is recommended that super resolution be the same as resolution.
>
> ***p_value***: Cutoff for intra-chromosomal interactions. The default is `0.1`.
> 
> ***zscore***: Cutoff for intra-chromosomal interactions. The default is `1.0`.

The input data are located in `'/in_data_folder/hic_data/hicup_processed/cohort'` while the output data can be found in `'/out_data_folder/hic_data/hic_interaction_homer'`.

Then, BEDTools is used to convert HOMER exported significant interactions to bed format files and Hi-C adjacency matrices are build with predefined resolution. 

```python
dna.preprocess_hic_homer2bed(genome_version = 'hg19', fig_dpi = 300)
```
> **Parameters:**
>
> ***genome_version***: Chromosome regions such as hg19, hg38. The default is `'hg19'`.
>
> ***fig_dpi***: Figure resolution in dots per inch. The default is `300`.

The input data are previously selected significant intra-chromosomal Hi-C interactions from `dna.preprocess_hic_tag2homer()`, which are deposited in `'/out_data_folder/hic_data/hic_interaction_homer'`. The output data including heatmap of interactions, are stored in `'/out_data_folder/hic_data/hic_interaction_bed'`.

## STEP 1.2 Preprocess multi-omics data

For other omics data, they are first organized into bed format files.

```python
dna.preprocess_omics_map2hic(multi_omics)
```
> **Parameters:**
>
> ***multi_omics***: Only accept `'gene expression'`, `'nucleosome density'`, and `'histone marker'` currently. Other omics data can also be added if they have been appropriately prepared.

The input data come from `'/in_data_folder/multi_omics'` while the output data can be found in `'/out_data_folder/multi_omics/resolution/out_data'`.

Then, multi-omics datasets are mapped to the Hi-C adjacency matrices and heatmaps are drawn based on the genomic feature matrices.

```python
dna.preprocess_omics_heatmap(multi_omics, type_of_calculation, fig_dpi = 300)
```
> **Parameters:**
>
> ***multi_omics***: Only accept `'gene expression'`, `'nucleosome density'`, and `'histone marker'` currently. Other omics data can also be added if they have been appropriately prepared.
> 
> ***type_of_calculation***: Two options for the value of multi-omics data in the crossed window. It should be setted as `'mean'` or `'max'`.
> 
> ***fig_dpi***: Figure resolution in dots per inch. The default is `300`.

With the output of the previous function as input, the zscore matrix and heatmap are generated to `'/out_data_folder/multi_omics/resolution/out_plot'`


## STEP 2.1 Identify intra-chromosomal communities

After preparing the Hi-C interaction matrices, clustering algorithm (ModularityOptimizer.jar) is used to identify network communities.

```python
dna.cluster_for_hic(modularity_function = 1, resolution_parameter = 1.0, optimization_algorithm = 3, n_random_starts = 10, n_iterations = 100, random_seed = 0, print_output = 1)
```
> **Parameters:**
>
> ***modularity_function***: The modularity function for clustering method where `1 = standard` and `2 = alternative`. The default is `1`. 
> 
> ***resolution_parameter***: The resolution of clustering method. The default is `1.0`. 
> 
> ***optimization_algorithm***: The selection of optimization algorithm where `1 = original Louvain algorithm`, `2 = Louvain algorithm with multilevel refinement`, and `3 = SLM algorithm`. The default is `3`. 
>
> ***n_random_starts***: The number of random starts for clustering. The default is `10`.
> 
> ***n_iterations***: The number of iterations per random start. The default is `100`.
> 
> ***random_seed***: The seed of the random number generator. The default is `0`.
> 
> ***print_output***: Whether or not to print output to the console (`0 = no`; `1 = yes`). The default is `1`.

Then, cluter results like nodes and communities are organized into .tsv format, while networks are transformed into .json format. All these results can be found in `'/out_data_folder/hic_data/resolution/hic_community_data'`.

```python
dna.cluster_to_json()
```
Next, valid communities should be identified in all clusters. Hi-C interactions and other genomic features within each community need to be analyzed.

```python
dna.cluster_community_structure(minCluster_size = 20, fig_dpi = 300)
```
> **Parameters:**
>
> ***minCluster_size***: The minimum size (or number of interactions) of valid communities. It is recommeded to use module `estimate_community_size` to identify the optimal value at a specific window resolution. The default is `20`. Value smaller than `3` is not accepted.
> 
> ***fig_dpi***: Figure resolution in dots per inch. The default is `300`.

Through this, a summary of interactions and modularity score of network clustering and valid communities are exported to a table located in `'/out_data_folder/hic_data/resolution/hic_community_figures'`. 

Furthermore, a super network is constructed according to the clustering results, where the node size and the edge width indicate the community size and the number of community interactions.

```python
dna.cluster_super_network(minCluster_size = 20, fig_dpi = 300)
```
> **Parameters:**
>
> ***minCluster_size***: The minimum size (or number of interactions) of valid communities. It should be consistent with the `minCluster_size` in `cluster_community_structure`. The default is `20`. Value smaller than `3` is not accepted.
> 
> ***fig_dpi***: Figure resolution in dots per inch. The default is `300`. 

The exported super networks can be found in `'/out_data_folder/hic_data/resolution/hic_community_figures'`.

## STEP 2.2. Enrichment of genomic feature 

After dividing the chromosome into several valid structures, random permutation tests of enrichment are performed by comparing genomic features within and outside each community.

```python
dna.enrichment_permutation(permutation = 100, pval_cutoff = 0.01, fig_dpi = 300)
```
> **Parameters:**
>
> ***permutation***: Resampling times of random permutation test of genomic enrichment. The default is `100`.
> 
> ***pval_cutoff***: Cutoff for features in a community is significantly higher/lower than randomly selected samples. The default is `0.01`.
> 
> ***fig_dpi***: Figure resolution in dots per inch. The default is `300`. 

Then, draw both tval and pval of genomic feautre heatmaps from `enrichment_permutation`.

```python
dna.enrichment_heatmap(permutation = 100, fig_dpi = 300)
```
> **Parameters:**
>
> ***permutation***: Number of sampling when doing permutation test in `enrichment_permutation`. The default is `100`.
> 
> ***fig_dpi***: Figure resolution in dots per inch. The default is `300`. 

The permutation results are stored in `'/out_data_folder/hic_data/resolution/hic_community_figures'`.

## STEP 3 Differential network analysis
## STEP 3.1. Differential interacting nodes

Find centrality of nodes in Hi-C networks

**Required**: None

**Optional**: permutation (Default = 100)

```python
dna.network_centrality()
```

**Note**: The following analysis requires completing the above steps for two cohorts!

Screen for significant nodes between networks from two cohorts

**Required**: cohort1 and cohort2 with above steps completed are required by network comparison. Chromosome is required.

**Optional**: pval_cutoff (Default = 0.05), permutation (Default = 100), fig_dpi (Default = 300)

```python
dna.network_sigNodes(cohort1, cohort2, chromosome)
```

Compare significant nodes with the information from feature data 

**Required**: cohort1 and cohort2 with above steps completed are required by network comparison. Chromosome is required.

**Optional**: pval_cutoff (Default = 0.05)

```python
dna.network_comparison(cohort1, cohort2, chromosome)
```

## STEP 3.2. Enrichment and network analysis

**Note**: Before this step, users should go to DAVID website for gene enrichment analysis base on the DIEGs if they want to get figures of enrichment analysis.

Gene enrichment anaysis based on identified DIEGs from DIGs 

**Required**: cohort1 and cohort2 with above steps completed are required by network comparison

**Optional**: method (Default = 'relativeRatio'), pval_cutoff (Default = 0.05), fig_dpi (Default = 300)

```python
dna.diegs_enrichment(cohort1, cohort2)
```

Construct subnetworks based on selected DIEGs

**Required**: cohort1 and cohort2 with above steps completed are required by network comparison. Chromosome is required.

**Optional**: pval_cutoff (Default = 0.05)

```python
dna.diegs_subnetwork(cohort1, cohort2, chromosome)
```


## STEP S. Estimate tuning parameters

**Note**: If you want to explore the full input data including 23 chromosomes, the computing will be quite time consuming, be careful and patient :)

**Required**: Chromosome, cohort1 and cohort2 are required. cal_type is required and only 0, 1, 2, 3 are accepted:
    0: comparison between different resolution, 
    1: comparison between different super resolution with resolution equal to 50kb,
    2: comparison between different super resolution with resolution equal to 100kb,
    3: comparison between different super resolution with resolution equal to 500kb.

**Optional**: fig_dpi (Default = 300)

```python
dna.estimate_resolution(chromosome, cal_type, cohort1, cohort2)
```

**Note**: The following analysis requires the first two functions in STEP 2.1.Cluster Hi-C interactions to communities finished (`dna.cluster_for_hic()` and `dna.cluster_to_json()`).  It is recommended to get the optimal cutoff by using all 23 chromosomes

**Required**: cohort and chromosome are required.

**Optional**: cutoff4Proportion (float = 0.02), fig_dpi (Default = 300)

```python
dna.estimate_community_size(cohort, chromosome)
```

***

More details about parameters, functions, and modules can be found in the main function [dnaici.py](https://github.com/differential-network-analysis/dnaici/blob/master/dnaici/dnaici/dnaici.py). For supplementary files, a [HOMEPAGE](https://dnaici.github.io/dnaici/) about the paper is under development...



