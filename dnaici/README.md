# README for DNAICI USERS

## DNAICI: Differential Network Analysis in Intra-chromosomal Community Interaction

## FOREWORD

The README file provides instructions for analyzing [demo data](https://github.com/differential-network-analysis/dnaici/tree/master/demo). The default values for parameters (e.g. `cohort = 'untreated'`, `chromosome = ['chr18','chr19']`) are designed for demo data. For more details about demo data, please refer to [README.md](https://github.com/differential-network-analysis/dnaici/blob/master/demo/README.md).

The DNAICI pipeline consists of three main mudules
* `preprocess`: preprocess_hic, preprocess_omics
* `analysis`: cluster, enrichment, network, dieg
* `parameter`: estimate

The module `preprocess` contains two parts for the preprocessing of raw Hi-C data and other omics data. Four parts in the `analysis` module are used for integrated analysis of multi-omics data. The `parameter` module are used to estimate the optimal values of two tuning parameters. All the above modules depend on two support modules: `tools` and `bin`.

To obtain consistent results with those in our paper, we recommend that you download and use the [full data](https://drive.google.com/file/d/1YbdZ7y5bRNqbP_4hVt6rcZM2Om1PoA-b/view?usp=drive_link). Necessary adjustment to the parameters need to be made for analyzing the full data. 

## STEP 0. Load DNAICI package and initialize global parameters

HOMER and Bedtools are required by data preprocessing. Java environment and ModularityOptimizer.jar are required by clustering algorithm.

`from dnaici import dnaici`

`dna = dnaici.DNAICI(in_data_folder, out_data_folder)`

**Required**
    
    in_data_folder: the directory where Hi-C, gene expression, nucleosome density, histone marker, and chromosome regions hg19 are deposited.
    
    out_data_folder: the directory where you want to export files.
    
**Optional** 
    
    cohort: Default = 'untreated'. 
    
    chromosome: Default = ['chr18','chr19'].
    
    resolution: Default = 500000.

>**Note**: If you want to get the complete results for 23 chromosomes, please let the input of chromosome be 'whole_genome'. Whole genome analysis is time consuming, please be careful and patient.

## STEP 1. Preprocess multi_omics data
## STEP 1.1. Preprocess Hi-C data

**Required**: 'super resolution' is required by preprocess_hic_tag2homer

**Optional**: p_value (Default = 0.1), zscore (Default = 1.0)

`dna.preprocess_hic_tag2homer(super_resolution)`

**Required**: None

**Optional**: genome_version (Default = 'hg19'), fig_dpi (Default = 300)

`dna.preprocess_hic_homer2bed()`


## STEP 1.2 Preprocess multi-omics data

**Required**: 'multi_omics' is required by preprocess_omics_map2hic. Only accept 'gene expression', 'nucleosome density', and 'histone marker'

**Optional**: None

`dna.preprocess_omics_map2hic(multi_omics)`

**Required**: 'multi_omics', 'type_of_calculation' is required by preprocess_omics_heatmap. 'multi_omics' should be 'gene expression', 'nucleosome density', or 'histone marker'. 'type_of_calculation' should be 'mean' or 'max'

**Optional**: color_start (Default = -2), color_end (Default = 2), bar_start (Default = 0), bar_end (Default = 0), fig_dpi (Default = 300)

`dna.preprocess_omics_heatmap(multi_omics, type_of_calculation)`


## STEP 2 Identify intra-chromosomal communities
## STEP 2.1. Cluster Hi-C interactions to communities

Identify communities using Hi-C data  

**Required**: None

**Optional**: modularity_function (Default = 1), resolution_parameter (Default = 1.0), optimization_algorithm (Default = 3), n_random_starts (Default = 10), n_iterations (Default = 100), random_seed (Default = 0), print_output (Default = 1)

`dna.cluster_for_hic()`

Transform cluter results into json format

**Required**: None

**Optional**: None

`dna.cluster_to_json()`

Export the heatmaps for clustered Hi-C interactions and genomic features in each community

**Required**: None

**Optional**: minCluster_size (Default = 20), fig_dpi (Default = 300)

`dna.cluster_community_structure()`

**Note**: minCluster_size smaller than 3 is not accepted.

Make a super network based on cluster labels of nodes

**Required**: None

**Optional**: minCluster_size (Default = 20), fig_dpi (Default = 300)

`dna.cluster_super_network()`


## STEP 2.2. Enrichment of genomic feature 

Random permutation test of genomic features within each community

**Required**: None

**Optional**: permutation (Default = 100), pval_cutoff (Default = 0.01), fig_dpi (Default = 300)

`dna.enrichment_permutation()`

Heatmaps of genomic feature enrichment

**Required**: None

**Optional**: permutation (Default = 100), fig_dpi (Default = 300)

`dna.enrichment_heatmap()`


## STEP 3 Differential network analysis
## STEP 3.1. Differential interacting nodes

Find centrality of nodes in Hi-C networks

**Required**: None

**Optional**: permutation (Default = 100)

`dna.network_centrality()`

**Note**: The following analysis requires completing the above steps for two cohorts!

Screen for significant nodes between networks from two cohorts

**Required**: cohort1 and cohort2 with above steps completed are required by network comparison. Chromosome is required.

**Optional**: pval_cutoff (Default = 0.05), permutation (Default = 100), fig_dpi (Default = 300)

`dna.network_sigNodes(cohort1, cohort2, chromosome)`

Compare significant nodes with the information from feature data 

**Required**: cohort1 and cohort2 with above steps completed are required by network comparison. Chromosome is required.

**Optional**: pval_cutoff (Default = 0.05)

`dna.network_comparison(cohort1, cohort2, chromosome)`


## STEP 3.2. Enrichment and network analysis

**Note**: Before this step, users should go to DAVID website for gene enrichment analysis base on the DIEGs if they want to get figures of enrichment analysis.

Gene enrichment anaysis based on identified DIEGs from DIGs 

**Required**: cohort1 and cohort2 with above steps completed are required by network comparison

**Optional**: method (Default = 'relativeRatio'), pval_cutoff (Default = 0.05), fig_dpi (Default = 300)

`dna.diegs_enrichment(cohort1, cohort2)`

Construct subnetworks based on selected DIEGs

**Required**: cohort1 and cohort2 with above steps completed are required by network comparison. Chromosome is required.

**Optional**: pval_cutoff (Default = 0.05)

`dna.diegs_subnetwork(cohort1, cohort2, chromosome)`


## STEP S. Estimate tuning parameters

**Note**: If you want to explore the full input data including 23 chromosomes, the computing will be quite time consuming, be careful and patient :)

**Required**: Chromosome, cohort1 and cohort2 are required. cal_type is required and only 0, 1, 2, 3 are accepted:
    0: comparison between different resolution, 
    1: comparison between different super resolution with resolution equal to 50kb,
    2: comparison between different super resolution with resolution equal to 100kb,
    3: comparison between different super resolution with resolution equal to 500kb.

**Optional**: fig_dpi (Default = 300)

`dna.estimate_resolution(chromosome, cal_type, cohort1, cohort2)`

**Note**: The following analysis requires the first two functions in STEP 2.1.Cluster Hi-C interactions to communities finished (`dna.cluster_for_hic()` and `dna.cluster_to_json()`).  It is recommended to get the optimal cutoff by using all 23 chromosomes

**Required**: cohort and chromosome are required.

**Optional**: cutoff4Proportion (float = 0.02), fig_dpi (Default = 300)

`dna.estimate_community_size(cohort, chromosome)`

***

More details about parameters, functions, and modules can be found in the main function [dnaici.py](https://github.com/differential-network-analysis/dnaici/blob/master/dnaici/dnaici/dnaici.py). For supplementary files, a [HOMEPAGE](https://dnaici.github.io/dnaici/) about the paper is under development...



