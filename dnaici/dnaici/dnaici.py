import os
import sys
import numpy as np
import pandas as pd
import networkx as nx
import scipy.stats as stats
from .preprocess import Preprocess_hic_tag2homer
from .preprocess import Preprocess_hic_homer2bed
from .preprocess import Preprocess_omics_map2hic
from .preprocess import Preprocess_omics_heatmap
from .analysis import Cluster_for_hic
from .analysis import Cluster_to_json
from .analysis import Cluster_community_structure
from .analysis import Cluster_super_network
from .analysis import Enrichment_permutation
from .analysis import Enrichment_heatmap
from .analysis import Network_centrality
from .analysis import Network_sigNodes
from .analysis import Network_comparison
from .analysis import DIEGs_enrichment
from .analysis import DIEGs_subnetwork
from .parameter import Estimate_resolution
from .parameter import Estimate_community_size


class DNAICI():
    '''Differential Network Analysis in Intra-chromosomal Community Interaction'''
    def __init__(self, 
                 in_data_folder: str, 
                 out_data_folder: str,
                 cohort: str = 'untreated', 
                 chromosome: list = ['chr18','chr19'],
                 resolution: int = 500000
                 ):
        '''
        Initialize DNAICI with user-defined parameters
        ----------
        Parameters
        ----------
        in_data_folder : str
            DESCRIPTION. Input data directory, including Hi-C, gene expression, nucleosome density, 
            histone marker, and chromsome regions hg19 data in this folder.
        out_data_folder : str
            DESCRIPTION. output data directory including processed Hi-C, gene expression, nucleosome
            density, histone marker and computational results in it.
        cohort : str, optional
            DESCRIPTION. The experimental condition for input data. The default is 'untreated', 
            representing untreated MCF7 cell line.
        chromosome : str, optional
            DESCRIPTION. List of chromosomes you want to investigate. The default is ['chr18','chr19'], if 
            you want to examine all chromosomes, use 'whole_genome'.
        resolution : int, optional
            DESCRIPTION. The resolution is used to divide the chromosome we investigated into 
            equally sized window bin. The default is 500000, representing 500 kilobase pair.

        '''
        self.in_data_folder = in_data_folder
        self.cohort = cohort
        self.chromosome = chromosome
        self.resolution = resolution
        self.out_data_folder = out_data_folder
    
    
    def preprocess_hic_tag2homer(self,
                                 super_resolution: int,
                                 p_value: float = 0.1,
                                 zscore: float = 1.0
                                 ):
        '''
        -----------
        Description
        -----------
        Applying HOMER to HICUP data to filter and reorganize significant intra-chromosomal Hi-C 
        interactions. HICUP data should be stored in 'in_data_folder/hic_data'. Super resolution, 
        p-value, and z-score are parameters required by HOMER for selecting significant 
        interactions. Output data are stored in 'out_data_folder/hic_data/hic_interaction_homer'.
        ----------
        Parameters
        ----------
        super_resolution : int
            DESCRIPTION. The range of signal averaging, recommended to be the same as resolution.
        p_value : float, optional
            DESCRIPTION. Paremeter for filtering intra-chromosomal interactions. The default is 0.1.
        zscore : float, optional
            DESCRIPTION. Paremeter for filtering intra-chromosomal interactions. The default is 1.0.
        
        '''
        in_data_folder = self.in_data_folder
        cohort = self.cohort
        chromosome = self.chromosome
        resolution = self.resolution
        out_data_folder = self.out_data_folder
        # preprocessing Hi-C tags using homer
        Preprocess_hic_tag2homer.main(in_data_folder, cohort, chromosome, resolution, 
                                      super_resolution, p_value, zscore, out_data_folder)
    
    
    def preprocess_hic_homer2bed(self,
                                 genome_version: str = 'hg19',
                                 fig_dpi: float = 300
                                 ):
        '''
        -----------
        Description
        -----------
        Using BEDTools to convert HOMER exported significant interactions to bed format files. Input
        data are previously generated significant intra-chromosomal Hi-C interactions, which are 
        deposited in 'out_data_folder/hic_data/hic_interaction_homer'. Output data including heatmap
        of interactions, are stored in 'out_data_folder/hic_data/hic_interaction_bed'. Chromosome 
        region file are used as reference to make window bin bed file.
        ----------
        Parameters
        ----------
        genome_version : str, optional
            DESCRIPTION. Chromosome regions such as hg19 or hg38. The default is 'hg19'.
        fig_dpi : float, optional
            DESCRIPTION. Figure resolution in dots per inch. The default is 300.
        
        '''
        in_data_folder = self.in_data_folder
        cohort = self.cohort
        chromosome = self.chromosome
        resolution = self.resolution
        out_data_folder = self.out_data_folder
        # preprocessing Hi-C data to bed format
        Preprocess_hic_homer2bed.main(in_data_folder, cohort, chromosome, resolution, 
                                      out_data_folder, genome_version, fig_dpi)
        
        
    def preprocess_omics_map2hic(self,
                                 multi_omics: str
                                 ):
        '''
        -----------
        Description
        -----------
        Make multi-omics matrices basaed on Hi-C interactions. Input multi-omics data should be 
        stored in 'in_data_folder/'. The type of multi-omics data needs to be specified in input 
        parameter. For example, if multi_omics = 'gene expression', then the RNA-Seq should be 
        stored in 'in_data_folder/gene expression'. Output are stored in 
        'out_data_folder/multi-omics/out_data'.
        ----------
        Parameters
        ----------
        multi_omics : str
            DESCRIPTION. Different omics data from the same breast cancer cell system, including
            'gene expression', 'nucleosome density', 'histone marker' data.

        '''
        in_data_folder = self.in_data_folder
        cohort = self.cohort
        chromosome = self.chromosome
        resolution = self.resolution
        out_data_folder = self.out_data_folder
        # preprocessing multi-omics data by mapping to hic data
        Preprocess_omics_map2hic.main(in_data_folder, cohort, chromosome, resolution, multi_omics,
                                      out_data_folder)

    
    def preprocess_omics_heatmap(self,
                                 multi_omics: str,
                                 type_of_calculation: str,
                                 fig_dpi: float = 300
                                 ):
        '''
        -----------
        Description
        -----------
        After map multi-omics data to Hi-C data, use multi-omics data to make heatmap based on Hi-C 
        interaction matrix. The edge weight and node weight are computered based on methods 
        published in https://pubmed.ncbi.nlm.nih.gov/36220609/
        ----------
        Parameters
        ----------
        multi_omics : str
            DESCRIPTION. Different omics data from the same breast cancer cell system, including
            'gene expression', 'nucleosome density', 'histone marker' data.
        type_of_calculation : str. 
            DESCRIPTION. it should be setted as 'mean' or 'max'.
        fig_dpi : float, optional
            DESCRIPTION. Figure resolution in dots per inch. The default is 300.

        '''
        in_data_folder = self.out_data_folder
        cohort = self.cohort
        chromosome = self.chromosome
        resolution = self.resolution
        out_data_folder = self.out_data_folder
        # plot heatmaps of multi-omics data
        Preprocess_omics_heatmap.main(in_data_folder, cohort, chromosome, resolution, multi_omics,
                                      type_of_calculation, out_data_folder, fig_dpi)
    
    
    def cluster_for_hic(self,
                       modularity_function: int = 1,
                       resolution_parameter: float = 1.0,
                       optimization_algorithm: int = 3,
                       n_random_starts: int = 10,
                       n_iterations: int = 100,
                       random_seed: int = 0,
                       print_output: int = 1
                       ):
        '''
        -----------
        Description
        -----------
        Network clustering by using java -jar ModularityOptimizer.jar. The community for each intra-
        chromosomal interactions are obtained and exported. Export a table summarizing the number of
        edges and the score of communities for each chromosome.
        ----------
        Parameters
        ----------
        modularity_function : int, optional
            DESCRIPTION. The modularity function for clustering method where 1 = standard and
            2 = alternative. The default is 1. 
        resolution_parameter : float, optional
            DESCRIPTION. The resolution of clustering method. The default is 1.0. 
        optimization_algorithm : int, optional
            DESCRIPTION. The selection of optimization algorithm where 1 = original Louvain 
            algorithm, 2 = Louvain algorithm with multilevel refinement, and 3 = SLM algorithm. 
            The default is 3. 
        n_random_starts : int, optional
            DESCRIPTION. The number of random starts for clustering. The default is 10.
        n_iterations : int, optional
            DESCRIPTION. The number of iterations per random start. The default is 100.
        random_seed : int, optional
            DESCRIPTION. The seed of the random number generator. The default is 0.
        print_output : int, optional
            DESCRIPTION. Whether or not to print output to the console (0 = no; 1 = yes). The 
            default is 1.
            
        '''
        in_data_folder = self.out_data_folder
        cohort = self.cohort
        chromosome = self.chromosome
        resolution = self.resolution
        out_data_folder = self.out_data_folder
        # do network clustering in intra-chromosomal interaction matrices
        Cluster_for_hic.main(in_data_folder, cohort, chromosome, resolution, out_data_folder,
                            modularity_function, resolution_parameter, optimization_algorithm,
                            n_random_starts, n_iterations, random_seed, print_output)
    
    
    def cluster_to_json(self):
        '''
        -----------
        Description
        -----------
        Organize nodes and communities into .tsv format. Export the network to json format.
        ----------
        Parameters
        ----------
        No parameters are required.
        
        '''
        in_data_folder = self.out_data_folder
        cohort = self.cohort
        chromosome = self.chromosome
        resolution = self.resolution
        out_data_folder = self.out_data_folder
        # make clustering results to json format
        Cluster_to_json.main(in_data_folder, cohort, chromosome, resolution, out_data_folder)


    def cluster_community_structure(self,
                                    minCluster_size: int = 20,
                                    fig_dpi: float = 300
                                    ):
        '''
        -----------
        Description
        -----------
        Make and export 1) heatmap for classified Hi-C interactions in valid communities 2) heatmaps 
        for multi-omics features in each of clustered Hi-C interaction communicity. A summary of 
        interactions and modularity score of network clustering and valid communities filtering are 
        exported to a table. 
        ----------
        Parameters
        ----------
        minCluster_size : int, optional
            DESCRIPTION. The minimum size (or number of interactions) of valid communities. It is 
            recommeded to use modele 'estimate_community_size' to identify the optimal value at a 
            specific window resolution. The default is 20. Value smaller than 3 is not accepted.
        fig_dpi : float, optional
            DESCRIPTION. Figure resolution in dots per inch. The default is 300.
            
        '''
        in_data_folder = self.out_data_folder
        cohort = self.cohort
        chromosome = self.chromosome
        resolution = self.resolution
        out_data_folder = self.out_data_folder
        # plot heatmaps for valid communities
        Cluster_community_structure.main(in_data_folder, cohort, chromosome, resolution, 
                                         out_data_folder, minCluster_size, fig_dpi)


    def cluster_super_network(self,
                              minCluster_size: int = 20,
                              fig_dpi: float = 300
                              ):
        '''
        -----------
        Description
        -----------
        Make and export a super network of each chromosome. The node size and the edge width 
        indicate the community size and the number of community interactions, respectively. 
        ----------
        Parameters
        ----------
        minCluster_size : int, optional
            DESCRIPTION. The minimum size (or number of interactions) of valid communities. It is 
            recommeded to use modele 'estimate_community_size' to identify the optimal value at a 
            specific window resolution. The default is 20. Value smaller than 3 is not accepted. The
            value should be consistent with the parameter in 'cluster_community_structure'.
        fig_dpi : float, optional
            DESCRIPTION. Figure resolution in dots per inch. The default is 300.

        '''
        in_data_folder = self.out_data_folder
        cohort = self.cohort
        chromosome = self.chromosome
        resolution = self.resolution
        out_data_folder = self.out_data_folder
        # plot super networks for valid communities     
        Cluster_super_network.main(in_data_folder, cohort, chromosome, resolution, out_data_folder,
                                   minCluster_size, fig_dpi)


    def enrichment_permutation(self,
                               permutation: int = 100,
                               pval_cutoff: float = 0.01,
                               fig_dpi: float = 300
                               ):
        '''
        -----------
        Description
        -----------
        Perform random permutation tests of enrichment by comparing genomic features within and 
        outside the community.
        ----------
        Parameters
        ----------
        permutation : int, optional
            DESCRIPTION. Resampling times of random permutation test of genomic enrichment. The 
            default is 100.
        pval_cutoff : float, optional
            DESCRIPTION. Cutoff for features in a community is significantly higher/lower than 
            randomly selected samples. The default is 0.01.
        fig_dpi : float, optional
            DESCRIPTION. Figure resolution in dots per inch. The default is 300.

        '''
        in_data_folder = self.out_data_folder
        cohort = self.cohort
        chromosome = self.chromosome
        resolution = self.resolution
        out_data_folder = self.out_data_folder
        # permutation test for genomic enrichment 
        Enrichment_permutation.main(in_data_folder, cohort, chromosome, resolution, out_data_folder,
                                    permutation, pval_cutoff, fig_dpi)


    def enrichment_heatmap(self,
                           permutation: int = 100,
                           fig_dpi: float = 300
                           ):
        '''
        -----------
        Description
        -----------
        Draw both tval and pval of genomic feautre heatmaps from permutation test in network 
        clustering.
        ----------
        Parameters
        ----------
        permutation : int, optional
            DESCRIPTION. Number of sampling when doing permutation test in enrichment_permutation.
            The parameter is used for finding the correct input file. The default is 100.
        fig_dpi : float, optional
            DESCRIPTION. Figure resolution in dots per inch. The default is 300.

        '''
        in_data_folder = self.out_data_folder
        cohort = self.cohort
        chromosome = self.chromosome
        resolution = self.resolution
        out_data_folder = self.out_data_folder
        # draw permutation results for genomic enrichment
        Enrichment_heatmap.main(in_data_folder, cohort, chromosome, resolution, out_data_folder,
                                permutation, fig_dpi)


    def network_centrality(self,
                           permutation: int = 100
                           ):
        '''
        -----------
        Description
        -----------
        Calculate centrality of nodes in Hi-C networks
        ----------
        Parameters
        ----------
        permutation : int, optional
            DESCRIPTION. Number of random permurtations for assigning a p-val to each node's 
            centrality. The default is 100.

        '''
        in_data_folder = self.out_data_folder
        cohort = self.cohort
        chromosome = self.chromosome
        resolution = self.resolution
        out_data_folder = self.out_data_folder
        # calculate the node centrality of hic network
        Network_centrality.main(in_data_folder, cohort, chromosome, resolution, out_data_folder, 
                                permutation)


    def network_sigNodes(self,
                         cohort1: str,
                         cohort2: str,
                         chromosome: list,
                         pval_cutoff: float = 0.05,
                         permutation: int = 100,
                         fig_dpi: float = 300
                         ):
        '''
        -----------
        Description
        -----------
        Compare analyzed data and identify significant nodes in each chromosome between two groups 
        such as untreated and tamr. Both node centrality and genomic features are considered to .
        ----------
        Parameters
        ----------
        cohort1 : str
            DESCRIPTION. dataset 1.
        cohort2 : str
            DESCRIPTION. dataset 2.
        chromosome : list
            DESCRIPTION. The chromosome list users want to investigate. 'whole_genome' is 
            recommended. 
        pval_cutoff : float, optional
            DESCRIPTION. P-value for screen for significant nodes. The default is 0.05.
        permutation : int, optional
            DESCRIPTION. Number of sampling when doing permutation test in enrichment_permutation.
            The parameter is used for finding the correct input file. The default is 100.
        fig_dpi : float, optional
            DESCRIPTION. Figure resolution in dots per inch. The default is 300.

        '''
        in_data_folder = self.out_data_folder
        resolution = self.resolution
        out_data_folder = self.out_data_folder
        # compare and identify significant nodes
        Network_sigNodes.main(in_data_folder, chromosome, resolution, out_data_folder, cohort1,
                              cohort2, pval_cutoff, permutation, fig_dpi)
        

    def network_comparison(self,
                           cohort1: str,
                           cohort2: str,
                           chromosome: list,
                           pval_cutoff: float = 0.05
                           ):
        '''
        -----------
        Description
        -----------
        For the selected significantly changed nodes betwenn cohort 1 and cohort 2, link feature 
        data to each node and summarize both commons and differences.
        ----------
        Parameters
        ----------
        cohort1 : str
            DESCRIPTION. dataset 1.
        cohort2 : str
            DESCRIPTION. dataset 2.
        chromosome : list
            DESCRIPTION. The chromosome list users want to investigate. 'whole_genome' is 
            recommended. 
        pval_cutoff : float, optional
            DESCRIPTION. P-value for screen for significant nodes in network_sigNodes. The 
            parameter is used for finding the correct input file. The default is 0.05.

        '''
        in_data_folder = self.in_data_folder
        resolution = self.resolution
        out_data_folder = self.out_data_folder
        # link feature data to selected nodes
        Network_comparison.main(in_data_folder, chromosome, resolution, out_data_folder, cohort1, 
                                cohort2, pval_cutoff)


    def diegs_enrichment(self,
                         cohort1: str,
                         cohort2: str,
                         method: str = 'relativeRatio',
                         pval_cutoff: float = 0.05,
                         fig_dpi: float = 300
                         ):
        '''
        -----------
        Description
        -----------
        Find differentially interacting and expressed genes (DIEGs) from the selected significant 
        nodes. Gene enrichment analysis should be performed based on these DIEGs.
        ----------
        Parameters
        ----------
        cohort1 : str
            DESCRIPTION. dataset 1.
        cohort2 : str
            DESCRIPTION. dataset 2.
        method : str, optional
            DESCRIPTION. Methods (relativeRatio: relative ratio; foldChange: fold change) for
            selecting differentially expressed genes. The default is 'relativeRatio'.
        pval_cutoff : float, optional
            DESCRIPTION. P-value for screen for significant nodes in network_sigNodes. The 
            parameter is used for finding the correct input file. The default is 0.05.
        fig_dpi : float, optional
            DESCRIPTION. Figure resolution in dots per inch. The default is 300.

        '''
        in_data_folder = self.in_data_folder
        resolution = self.resolution
        out_data_folder = self.out_data_folder
        # find DIEGs from DIGs
        DIEGs_enrichment.main(in_data_folder, resolution, out_data_folder, cohort1, cohort2, method,
                              pval_cutoff, fig_dpi)


    def diegs_subnetwork(self,
                         cohort1: str,
                         cohort2: str,
                         chromosome: list,
                         pval_cutoff: float = 0.05
                         ):
        '''
        -----------
        Description
        -----------
        Construct subnetwork based on selected DIEGs. Genomic information (histone markers of 
        enhancer/repressor , gene expression) and gene names within each node are exported.       
        ----------
        Parameters
        ----------
        cohort1 : str
            DESCRIPTION. dataset 1.
        cohort2 : str
            DESCRIPTION. dataset 2.
        chromosome : list
            DESCRIPTION. The chromosome users want to investigate. 'whole_genome' is recommended. 
        pval_cutoff : float, optional
            DESCRIPTION. P-value for screen for significant nodes in network_sigNodes. The 
            parameter is used for finding the correct input file. The default is 0.05.

        '''
        in_data_folder = self.in_data_folder
        resolution = self.resolution
        out_data_folder = self.out_data_folder
        # export subnetwork from DIEGs
        DIEGs_subnetwork.main(in_data_folder, chromosome, resolution, out_data_folder, cohort1, 
                              cohort2, pval_cutoff)


    def estimate_resolution(self,
                            chromosome: list,
                            cal_type: int,
                            cohort1: str,
                            cohort2: str,
                            fig_dpi: float = 300
                            ):
        '''
        -----------
        Description
        -----------
        Determine the optimal resolution and super resolution by comparing the number of significant
        interactions and the distance of each interaction.
        ----------
        Parameters
        ----------
        chromosome: list
            DESCRIPTION. The chromosome you want to investigate, i.e. ['chr1', 'chr2', ...]. If you 
            want to check all the chromosome, please use 'whole_genome'. ATTENTION: whole genome 
            estimation is time consuming.
        cal_type : int
            DESCRIPTION. The type of calculation, where 
            0: comparison between different resolution, 
            1: comparison between different super resolution with resolution equal to 50kb,
            2: comparison between different super resolution with resolution equal to 100kb,
            3: comparison between different super resolution with resolution equal to 500kb.
        cohort1 : str
            DESCRIPTION. dataset 1.
        cohort2 : str
            DESCRIPTION. dataset 2.
        fig_dpi : float, optional
            DESCRIPTION. Figure resolution in dots per inch. The default is 300.

        '''
        
        in_data_folder = self.in_data_folder
        out_data_folder = self.out_data_folder
        # compare the hic information with different resolution and super resolution
        Estimate_resolution.main(in_data_folder, out_data_folder, chromosome, cal_type, cohort1, 
                                 cohort2, fig_dpi)


    def estimate_community_size(self,
                                cohort: str,
                                chromosome: list,
                                cutoff4Proportion: float = 0.02,
                                fig_dpi: float = 300
                                ):
        '''
        -----------
        Description
        -----------
        Determine the optimal resolution and super resolution by comparing the number of significant
        interactions and the distance of each interaction.
        ----------
        Parameters
        ----------
        cohort : str
            DESCRIPTION. dataset you want to investigate.
        chromosome: list
            DESCRIPTION. The chromosome you want to investigate, i.e. ['chr1', 'chr2', ...]. If you 
            want to check all the chromosome, please use 'whole_genome'. ATTENTION: whole genome 
            estimation is time consuming.
            time consuming.
        cutoff4Proportion : float, optional
            DESCRIPTION. The proportion is used to determine the minimal size of valid community.
            The default is 0.02.
        fig_dpi : float, optional
            DESCRIPTION. The default is 300.

        '''
        
        in_data_folder = self.out_data_folder
        resolution = self.resolution
        out_data_folder = self.out_data_folder
        
        Estimate_community_size.main(in_data_folder, cohort, chromosome, resolution, 
                                     out_data_folder, cutoff4Proportion, fig_dpi)


def main():
    DNAICI()


if __name__== '__main__':
    DNAICI()













