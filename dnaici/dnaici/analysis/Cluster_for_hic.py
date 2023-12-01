import numpy as np
import networkx as nx
import os
import pandas as pd
from dnaici.tools.data_functions import find_cutoff4percentage
import subprocess
import scipy.stats


def find_global_cutoff4matrix(chrom_strs, file_path, percentage, file_postprefix_string):
    
    '''
    find global cutoff value for all input adj-matrix based on selected the lowest percentage
    of data will be removed from all matrix
    input : chrom_strs, a list of chromosome strings cuch as chr1, chr2, ....
            file_path, input file folder location
            percentage, a percenage that used to filtering matrix such as percentage=0.1 then the lowest 10% of data will be removed
            file_postprefix_string, a string for input file name postprefix string such as '_zscore.matrix'
    return: gloabl_precentage_cutoff 
            a numerical value which can be used to as a cutoff value to remove the lowest percentage of data from all matrix.
    '''
    record_nonzeros_in_all_matrix = []
    print('Read:')
    for chrom in chrom_strs:
        hic_f = os.path.join(file_path, chrom + file_postprefix_string)
        print(hic_f)
        a_matrix = np.loadtxt(hic_f)
        A1 = a_matrix.copy()
        #here we assume it is undirect graph up/low diagnal elements are the same
        A=np.triu(A1).copy()
        #find all nonzeros in the matrix
        record_nonzeros_in_all_matrix.append(A[A!=0].copy())
    
    if (not np.isnan(percentage)) and (percentage != 0):
        all_values = np.sort(np.concatenate(record_nonzeros_in_all_matrix))
        len_of_values = all_values.shape[0]
        global_percentage_cutoff = all_values[int(len_of_values*percentage):][0]
    else:
        global_percentage_cutoff = np.nan
        all_values = np.nan
        len_of_values = np.nan 
    
    return global_percentage_cutoff, all_values, len_of_values


def find_significant_interactions(all_values, len_of_values, pval_cutoff=0.3):
    '''
    Convert Z-scores to P-values then select a globable cutoff value for selecting significant interactions based
    on predefined pval-cutoff value such as 0.3 in default
    '''
    #single tailed p-value
    #fit zscores to a normal distribution
    pval2zscores = scipy.stats.norm.sf(abs(all_values))
    #only consider right side zscores as significant interactions
    filter_conditions = (pval2zscores<pval_cutoff) & (all_values>0)
    global_pval_cutoff = min(all_values[filter_conditions])
    percentage_significant_interactions = len(all_values[all_values>=global_pval_cutoff])/len_of_values
    print(' %4.2f percentage of interactions with P-value < %4.2f '% ( percentage_significant_interactions*100, pval_cutoff )  )
    
    return global_pval_cutoff, pval_cutoff,percentage_significant_interactions
  

def export_edges4matrix(chrom_strs,file_path, isFind_global_cutoff, global_percentage_cutoff,\
                        withWeight,edges_rm_str,file_postprefix_string, percentage, out_path):
    
    '''
    Make input file for network clustering
    Input: chrom_strs, a list of chrom string such as chr1, chr2, ...
           file_path, input matrix file path
           isFind_global_cutoff, True or False for finding global cutoff based on lowest percentages
           global_percentagge_cutoff, a cutoff value will be used to filter lowest percentage of interactions
           withWeight, True or False, is export network with edge weight or not
           edges_rm_str, a string for edges remove such as _edges_rmGlob or _edges_rm
           file_postprefix_string, a string for file name postprefix such as _zscore.matrix
           percentage, percentage of the lowest interactions will be removed from matrix
           out_path, output file path
    Return:
           out_data_file, new adj-matirx after filtering,
           out_edge_file, network edges file
           record_edge_df, a dataframe for all edges
    '''
    record_num_edges_in_chrom = {} 
    record_adj_matrix4chrom = {}
    record_G4chrom = {}
    for chrom_str in chrom_strs:
        #load hic matrix
        hic_f = os.path.join(file_path, chrom_str + file_postprefix_string)
        #print(hic_f)
        a_matrix = np.loadtxt(hic_f)
        A1 = a_matrix.copy()
        #only conside updiagonl element in matrix
        A = np.triu(A1).copy()
        if not isFind_global_cutoff:
            #find cutoff value for the local lowest percentage of values in A
            cutoff_value4edges = find_cutoff4percentage(A, percentage)
        else:
            cutoff_value4edges = global_percentage_cutoff
        #find all edges in a matrix
        all_edges = len(A[A!=0]) 
        #filtering weak interactions based on lowest percentage
        if not np.isnan(cutoff_value4edges):
            print('Filteirng adj-matrix by weakest interactions')
            A[A<cutoff_value4edges] = 0
        #edges after filterring
        edges_after_filterring = len(A[A!=0])
        record_num_edges_in_chrom[chrom_str] = [all_edges, edges_after_filterring, cutoff_value4edges]
        #build a network based on filtered matrix
        G = nx.from_numpy_array(A)
        if True:
            max_weight = 20000
            #convert edge weight to positive values
            for u,v,d in G.edges(data = True):
                new_w = np.exp(d['weight'])
                if new_w > max_weight:
                    d['weight'] = max_weight
                else:
                    d['weight'] = new_w
        # — let us store the degree centralities for each nodes for a graph in a dictionary
        #deg_cen = {}
        #deg_cen = nx.degree_centrality(G)
        #export edges of each hic_interacton matrix after filtering lowest interactions
        if withWeight:
            out_file = chrom_str + edges_rm_str + str(percentage) + '_percent_4zscore.matrix'
            out_data_file = os.path.join(out_path, out_file)
            nx.write_weighted_edgelist(G, out_data_file, delimiter='\t')
        else:
            out_file = chrom_str + edges_rm_str + str(percentage) + 'percent_noWeight_4zscore.matrix'
            out_data_file = os.path.join(out_path, out_file)
            nx.write_edgelist(G,out_data_file, delimiter='\t', data=False)
        #export all edeges passed filtering condition in a file
        #here both up/dow diagnol elements of adj-matrix are exported 
        record_adj_matrix4chrom[chrom_str] = A.copy()
        record_G4chrom[chrom_str] = G.copy()
        print('Output 1 (%s):\n' %chrom_str, out_data_file, '\n')
    #export recorded edges before and aftering filtering
    #here both up and down diagnol elements are counted
    record_edge_df = pd.DataFrame.from_dict(record_num_edges_in_chrom).T.copy()
    #replace nan as zero
    record_edge_df.fillna(0, inplace=True)
    #record_edge_df=record_edge_df.astype(int)
    record_edge_df.columns = ['all_edges', 'edges_after_filtering', 'cutoff_value4edge']
    record_edge_df.all_edges = record_edge_df.all_edges.astype(int)
    record_edge_df.edges_after_filtering = record_edge_df.edges_after_filtering.astype(int)
    record_edge_df['percentage'] = record_edge_df.edges_after_filtering/record_edge_df.all_edges
    #out_edge_file='num_of'+ edges_rm_str+str(percentage) + 'percentage.edges'
    #out_edge_file=os.path.join(out_path,out_edge_file)
    ##record_edge_df.to_csv(out_edge_file, sep='\t')
    #print(out_edge_file)
    return out_data_file, record_edge_df, record_adj_matrix4chrom, record_G4chrom


def do_network_clustering4matrix(chrom_str, input_file,out_cluster_path,record_edge_df, modularity_function,
                                 resolution_parameter, optimization_algorithm, n_random_starts, n_iterations,random_seed, print_output):
    '''
    do network clustering and record modularity score for each cluster
    input: chrom_str, is chrom string such as chr1
           input_file, is input file name for edges in a adj-matrix
           out_cluster_path, is output file path
           record_edge_df, is a dataframe for storing of number of interactions, modularity scores of the chromosome
           other parameters are used by program (http://www.ludowaltman.nl/slm/) for performing network clustering
    for example, 
		modularity_function     Modularity function (1 = standard; 2 = alternative)
		resolution_parameter    Value of the resolution parameter
		optimization_algorithm  Algorithm for modularity optimization (1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM algorithm)
		n_random_starts Number of random starts
		n_iterations    Number of iterations per random start
		random_seed     Seed of the random number generator
		print_output    Whether or not to print output to the console (0 = no; 1 = yes)
    '''     
    output_file = os.path.join( out_cluster_path, os.path.basename(input_file).replace('_4zscore.matrix','_4communities.txt') )
    
    bin_folder = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../', 'bin'))
    cmd = 'java -jar ' + bin_folder + '/ModularityOptimizer.jar ' + input_file + ' ' \
        + output_file \
        + ' ' + str(modularity_function) + ' ' + str(resolution_parameter) + ' '+ str(optimization_algorithm) + ' ' \
        + str(n_random_starts) + ' ' + str(n_iterations) + ' ' + str(random_seed) + ' ' +str(print_output)
    
    print('Command:\n', cmd, '\n')
    print('Output 2 (%s):\n' %chrom_str, output_file, '\n')
    #results=os.system(cmd)
    #if results !=0:
        #print('Error in ', cmd)
    #export modularity score
    P = subprocess.Popen(cmd, stdout=subprocess.PIPE,stderr=subprocess.PIPE, shell=True)
    out, err = P.communicate()
    if len(err.decode('utf-8')) == 0:
        for li in out.decode('utf-8').split('\n'):
            if 'Maximum modularity in' in li:
                modularity_score = li.split(':')[-1].strip()
                record_edge_df.loc[chrom_str,'Mod_score'] = float(modularity_score)
    else:
        print(err.decode('utf-8'))
    
    return record_edge_df

########## START ###########
def main(in_data_folder,
         cohort,
         chromosome,
         resolution,
         out_data_folder,
         modularity_function = 1,       # 1 = standard; 2 = alternative
         resolution_parameter = 1.0,    # value of the resolution parameter
         optimization_algorithm = 3,    # 1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM algorithm
         n_random_starts = 10,          # number of random starts
         n_iterations = 100,            # number of iterations per random start
         random_seed = 0,               # seed of the random number generator
         print_output = 1):             # whether or not to print output to the console (0 = no; 1 = yes)
         
    percentage = 0
    #export edges with weight
    withWeight = False
    isFind_global_cutoff = False
    isFind_significant_interactions = False
    pval2zscore = 0.5   #default is 0.3
    #input hic adj-matrix and the corresponding window bin regions which is genrated in folder mcf7_hic_jbw
    #based on hic-interactions mapped to a predefined genomic window bin such as 500kb 
    file_postprefix_string = '_meanZscore.matrix'
    #path to hic adj-matrix
    bin_str = str(int(resolution/1000)) + 'kb'
    file_path = in_data_folder + '/hic_data/' + bin_str + '/hic_interaction_bed/' + cohort
    
    if chromosome == 'whole_genome':
        chrom_strs = ['chr'+str(i) for i in range(1,24)]
    else: 
        chrom_strs = chromosome
    
    #output path for data or figures will be generated by current analysis
    out_path = out_data_folder + '/hic_data/' + bin_str + '/hic_community_data/' + cohort
    if not os.path.exists(out_path):
        os.makedirs(out_path)
        print('Creat output folder:\n', out_path, '\n')
    #parameters for filtering weak edges in adj-matrix 
    if isFind_global_cutoff:
        edges_rm_str = '_edges_rmGlob_'
    else:
        edges_rm_str = '_edges_rm_'    
    #find a global cutoff value for rmeove the lowest percentages of interactions from genome-wide data
    if isFind_global_cutoff:
        #find a cutoff value for the weakest interactions
        global_percentage_cutoff, all_values, len_of_all_values = find_global_cutoff4matrix(chrom_strs, file_path, percentage, file_postprefix_string)
        if isFind_significant_interactions:
            #find a cutoff value for significant interactions where default pval_cutoff=0.3
            global_pval_cutoff, pval_cutoff, percentage_significant_interactions = find_significant_interactions(all_values, len_of_all_values, pval_cutoff=pval2zscore)
            global_percentage_cutoff = global_pval_cutoff
            percentage = "{:.2f}".format(1 - percentage_significant_interactions)
            print('Select (%5.2f percentags) significant interactions (P< %5.2f) for further study' % (percentage_significant_interactions*100,  pval_cutoff ))
        else:
            percentage = "{:.2f}".format(percentage)
            percentage_significant_interactions = 1- float(percentage)
    else:
        global_percentage_cutoff = np.nan
        
    #export edges passed filtering in a text file for being used by network clustering latter
    out_data_file, record_edge_df, record_adj_matrix4chrom, record_G4chrom = export_edges4matrix(chrom_strs,file_path, 
                                                                                                 isFind_global_cutoff, 
                                                                                                 global_percentage_cutoff,
                                                                                                 withWeight, edges_rm_str, 
                                                                                                 file_postprefix_string,
                                                                                                 percentage, out_path)    
    #out
    out_edge_file = 'num_of' + edges_rm_str + str(percentage) + 'percentage.edges'
    out_edge_file = os.path.join(out_path, out_edge_file)
    #run network clustering
    #network clustering based on filtered matrix
    for chrom_str in chrom_strs:
        if withWeight:
            input_file = os.path.join(out_path, chrom_str + edges_rm_str + str(percentage) + 'percent_4zscore.matrix')
        else:
            input_file = os.path.join(out_path, chrom_str + edges_rm_str + str(percentage) + 'percent_noWeight_4zscore.matrix')
        
        out_cluster_path = out_path #os.path.join('out_network_clusters',time_str,chrom_str)
        
        record_edge_df = do_network_clustering4matrix(chrom_str, input_file, out_cluster_path, record_edge_df, modularity_function,
                                                      resolution_parameter, optimization_algorithm, n_random_starts, n_iterations,random_seed, print_output) 
        
        #export recorded scoress
        print('Output 3 (%s):\n' %chrom_str, out_edge_file, '\n')
        record_edge_df.to_csv(out_edge_file, sep='\t')







