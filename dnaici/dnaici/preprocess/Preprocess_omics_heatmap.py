import numpy as np
import pandas as pd
from dnaici.tools.plot_tools import draw_plot
import matplotlib as mlp
import matplotlib.pyplot as plt
import os
mlp.use('agg')


def calculate_mean_of_pairwise_interaction(hic_matrix, omics_df):
    #assign mean value of multi-omics data to Hi-C interaction matrix
    new_matrix = np.zeros(hic_matrix.shape)
    for i in range(0, new_matrix.shape[0]):
        for j in range(i, new_matrix.shape[1]):
            if hic_matrix[i,j] != 0:
                #here df region index assume start from 1, but matrix index start from 0 
                tmp_region_i = i+1
                tmp_region_j = j+1
                tmp_val = (omics_df[omics_df.region==tmp_region_i].mean_value.to_numpy() + omics_df[omics_df.region==tmp_region_j].mean_value.to_numpy())[0]/2
                new_matrix[i,j] = tmp_val
                new_matrix[j,i] = tmp_val

    matrix2 = np.triu(new_matrix).copy()
    return matrix2.copy()


def calculate_weight_of_pairwise_interaction(hic_matrix, omics_df):
    #calculate a weight of two nodes' weights
    new_matrix = np.zeros(hic_matrix.shape)
    for i in range(0, new_matrix.shape[0]):
        for j in range(i, new_matrix.shape[1]):
            if hic_matrix[i,j] != 0:
                #here df region index assume start from 1, but matrix index start from 0
                tmp_region_i = i+1
                tmp_region_j = j+1
                tmp_val=max([0, (omics_df[omics_df.region==tmp_region_i].mean_value.to_numpy() + omics_df[omics_df.region==tmp_region_j].mean_value.to_numpy())[0]])
                new_matrix[i,j] = tmp_val
                new_matrix[j,i] = tmp_val
    
    matrix2 = np.triu(new_matrix).copy()
    return matrix2.copy()


def main(in_folder,
         cohort,
         chromosome,
         resolution,
         multi_omics,
         type_of_calculation,
         out_folder,
         fig_dpi
         ):
    
    bin_str = str(int(resolution/1000)) + 'kb'
    
    if multi_omics == 'gene expression':
        in_strs = 'expression_data'
        out_strs = ['_geneExp']
    elif multi_omics == 'nucleosome density':
        in_strs = 'nucleosome_density_data'
        out_strs = ['_DNas']
    elif multi_omics == 'histone marker':
        in_strs = 'histone_data'
        out_strs = ['_ctcf', '_h3k27ac', '_h3k27me3', '_h3k4me1', '_h3k4me3', '_h3k9me3']
    else:
        out_strs = []
        print('No such data type found :( Please set parameter \'multi_omics\' as \'gene expression\', \'nucleosome density\', or \'histone marker\'.')
        
    if chromosome == 'whole_genome':
        chrom_strs = ['chr'+str(i) for i in range(1,24)]
    else: 
        chrom_strs = chromosome
    
    for chrom_str in chrom_strs:
        for out_str in out_strs:
            data_type = chrom_str + out_str
            out_path = out_folder + '/' + in_strs + '/' + bin_str + '/out_plots/' + cohort
            if not os.path.exists(out_path):
                os.makedirs(out_path)            
            #read omics data
            in_path = in_folder + '/' + in_strs + '/' + bin_str + '/out_data/' + cohort
            in_file = os.path.join(in_path, chrom_str + '_' + bin_str + '_regions' + out_str + '_array.bed')
            print('Read multi-omics data: ', in_file)
            omics_df = pd.read_csv(in_file, sep='\t')
            #read Hi-C adj-matrix
            hic_path = in_folder + '/hic_data/' + bin_str + '/hic_interaction_bed/' + cohort
            hic_file = os.path.join(hic_path, chrom_str + '_meanZscore.matrix')
            print('Read Hi-C data: ', hic_file)
            hic_matrix = np.loadtxt(hic_file)
            
            if type_of_calculation == 'mean':
                #assign mean value to Hi-C interaction matrix
                matrix2 = calculate_mean_of_pairwise_interaction(hic_matrix, omics_df)
                cell_type4plot = data_type + '_mean'
                out_f = os.path.join(out_path, cell_type4plot+'_zscore.matrix')
            elif type_of_calculation == 'max':
                matrix2 = calculate_weight_of_pairwise_interaction(hic_matrix, omics_df)
                cell_type4plot = data_type + '_maxWeight'
                out_f = os.path.join(out_path, cell_type4plot+'_zscore.matrix')
            else:
                print('Not a valid input for \'type_of_calculation\' :( Please set parameter \'type_of_calculation\' as \'mean\' or \' max\'.')
                break
            #prepare for plotting of heatmap for the matrix
            matrix_start = 0
            matrix_end = matrix2.shape[1]
            
            #draw figure
            draw_plot(matrix2, cell_type4plot, matrix_start, matrix_end, -2, 2, 0, 0, out_path, fig_dpi=fig_dpi)
            
            #export omics matrix
            print('Export data: ',out_f)
            np.savetxt(out_f, matrix2, delimiter='\t', fmt='%10.5f')



