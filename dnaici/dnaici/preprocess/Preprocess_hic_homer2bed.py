import pandas as pd
import os
import numpy as np
import matplotlib as mlp
mlp.use('agg')
import matplotlib.pyplot as plt
import matplotlib.lines as lines
from PIL import Image
from dnaici.preprocess import Preprocess_window_bin_bed
from dnaici.tools.plot_tools import draw_plot2


def read_homerInteraction2df(f1, out_data_path, isExport=True):
    '''
    input: f1 file name of HOMER results, output_path output path of exported files
    return:
         left_df, for left node of pair of interactions
         right_df, for right node of pair of interactions
         out_left_file, file name of left node df
         out_right_file, file name of right node df
    Read exported significant interactions from HOMER and convert them to a pair of interactoins
    to left and right node dataframe or export these dataframe if it is required.
    '''
    #read significant interactoins
    df0 = pd.read_csv(f1, sep='\t')
    df0['idx'] = df0.index.to_list()

    #only keep intra-chromosome interactions
    df = df0[df0['chr(1)']==df0['chr(2)']].copy()

    #split left and right nodes to two dataframe
    left_df = df[['chr(1)', 'start(1)', 'end(1)','Z-score', 'LogP' ,'InteractionID']].copy()
    right_df = df[['chr(2)', 'start(2)', 'end(2)','Z-score', 'LogP' ,'InteractionID']].copy()

    left_df.sort_values(by='start(1)', inplace=True)
    right_df.sort_values(by='start(2)', inplace=True)

    out_left_file = os.path.join(out_data_path, os.path.basename(f1).replace('.txt', '_leftPos.bed.gz'))
    out_right_file = os.path.join(out_data_path, os.path.basename(f1).replace('.txt', '_rightPos.bed.gz'))
  
    #replace negative position by 1
    left_df.loc[left_df['start(1)']<0, 'start(1)'] = 1
    right_df.loc[right_df['start(2)']<0, 'start(2)'] = 1

    #export pair of edges to left and right position bed files
    if isExport:
        print(out_left_file)
        print(out_right_file)
        left_df.to_csv(out_left_file, sep='\t', index=False, header=None, compression='gzip')
        right_df.to_csv(out_right_file, sep='\t', index=False, header=None, compression ='gzip')
    return left_df.copy(), right_df.copy(), out_left_file, out_right_file


def map_left_right_node2defined_window_bins(in_bin_region_file, in_left_chrom_file, in_right_chrom_file, bin_str):
    '''
    map left and right chrom position of nodes to a predifined window bins such as 500kb window size
    input 500kb window bin regions in all chromosomes
    all files in bed format
    input:
        in_bin_region_file, a bed format of window bin size file
        in_left_chrom_file, a bed format of left node chromosome position
        in_right_chrom_file, a bed format of right node chromosome position
        bin_str: a postprefix of bin size string such as 500KB_bin , which will be used for file name of export files
    '''
    #in_bin_region_file='mcf7_hic_omer/chrom_region_500KB_window_bin.bed.gz'
    #in_left_chrom_file= out_left_file
    #in_right_chrom_file= out_right_file

    #set output file for mathced regions in 500kb window bins
    out_left_chrom_file = in_left_chrom_file.replace('.bed.gz', '_' + bin_str + '.bed')
    out_right_chrom_file = in_right_chrom_file.replace('.bed.gz', '_' + bin_str + '.bed')

    #use bed tools to assigne regino to window bins
    print('Mapp left and right nodes positions to window bin regions:\n  ', in_bin_region_file)
    cmd = 'bedtools intersect -a ' + in_bin_region_file + ' -b ' + in_left_chrom_file + ' -wao >' + out_left_chrom_file
    os.system(cmd)
    print(out_left_chrom_file)
    #
    cmd = 'bedtools intersect -a ' + in_bin_region_file + ' -b ' + in_right_chrom_file + ' -wao >' + out_right_chrom_file
    os.system(cmd)
    print(out_right_chrom_file)
    return out_left_chrom_file, out_right_chrom_file


def build_adj_matrix4mapped_interactions(left_region_file, right_region_file, chrom_str, weightGt=4, isFilterReads_count=False):
    '''
    Buld adj-matrix or count matrix for mapped interactions in predifined window bins
    input: left_region_file and right_region_file are two nodes' position files after bedtools intersect between significant interactions and predifined 
    window bins
            chrom_str, is the selected chromosome for current matrix such as chr1
    output: mean_adj_matrix, mean adj-matrix for zscores/read counts in pair-wise interactions
            a_matrix , total sum of zscores/read counts in adj-matrix
            count_matrix, count of number of regions/edges in a bin
            left_df and right_df are left and right nodes position in predfined window bins for a pair of significant interactions
            selected_left_df2, selected_right_df2 are left and right nodes postions in predfinied window bins after removing empty bins        
    '''
    #output datafram column names
    columns=['chrom','bin_start','bin_end','bin_id','numeric_chrom','chrom_region',
             'region_start','region_end','reads_count','pval','edges_id','overlaps']
    print(chrom_str)
    #left_region_file=out_left_chrom_file
    #right_region_file=out_right_chrom_file

    #read matched regions in 500kb window bins for both left and right nodes
    left_df=pd.read_csv(left_region_file,sep='\t',header=None)
    right_df=pd.read_csv(right_region_file,sep='\t',header=None)
    left_df.columns=columns
    right_df.columns=columns
    
    left_df.reads_count=left_df.reads_count.replace('.','0')
    right_df.reads_count=right_df.reads_count.replace('.','0')
    #filter df with a weight cutoff at reads_count and
    #find window bins of one chrom
    #edges_id=left_df.edges_id.unique()
    #chrom_str='chr1'
    if chrom_str=='chr23':
        chrom_str='chrX'

    if isFilterReads_count:
        print('Reads count < '+  str(weightGt) + ' are removed!' )
        selected_left_df=left_df[(left_df.chrom==chrom_str) & (left_df.reads_count.astype(float)>weightGt) ].copy()
        selected_right_df=right_df[ (right_df.chrom==chrom_str) & (right_df.reads_count.astype(float)>weightGt ) ].copy()
    else:
        selected_left_df=left_df[ left_df.chrom==chrom_str  ].copy()
        selected_right_df=right_df[ right_df.chrom==chrom_str  ].copy()
        
    #edges ids shall available in both nodes
    selected_edges_ids_right=selected_right_df.edges_id.unique()
    selected_edges_ids_left=selected_left_df.edges_id.unique()
    selected_edges_ids= list( set(selected_edges_ids_right) & set(selected_edges_ids_left))

    #selected_left_df=selected_left_df[selected_left_df.edges_id.apply(lambda x: x in set(selected_edges_ids))].copy()
    #selected_right_df=selected_right_df[selected_right_df.edges_id.apply(lambda x: x in set(selected_edges_ids))].copy()
  
    #remove unmapped window bins
    sel_edges_ids= [si for si in selected_edges_ids if si != '.']
    sel_edges_ids= np.array(sel_edges_ids,dtype=str)

    #create an initial adj-matrix for the chromosome based on length of window bins in a chrosome
    right_len_of_bins=selected_right_df.bin_id.max()
    left_len_of_bins=selected_left_df.bin_id.max()
    len_of_bins=max(right_len_of_bins, left_len_of_bins)
    print(len_of_bins)
    a_matrix=np.zeros((len_of_bins, len_of_bins))
    count_matrix=np.zeros((len_of_bins, len_of_bins))

    #filter bins without match to regions
    selected_left_df2=selected_left_df[selected_left_df.edges_id!='.'].copy()
    selected_right_df2=selected_right_df[selected_right_df.edges_id!='.'].copy()
    selected_right_df2.reads_count=selected_right_df2.reads_count.astype(float)

    #for each interaction/edge to find regions in bins
    for ei in sel_edges_ids:
        lf=selected_left_df2[selected_left_df2.edges_id==ei].copy()
        rt=selected_right_df2[selected_right_df2.edges_id==ei].copy()
        #in canse , one interaction is mapped to multiple window bins
        for idx,row_lf in lf.iterrows():
            for idx2,row_rt in rt.iterrows():
                #add total read counts/zscores to matched window bins, and assume both left and right nodes reads count are the same
                #here -1 is removed if input region index started from 0 such as jan2023 ??
                a_matrix[row_lf.bin_id-1,row_rt.bin_id-1] += row_rt.reads_count
                a_matrix[row_rt.bin_id-1,row_lf.bin_id-1] += row_rt.reads_count
                #count the number of regions mapped to the window bin
                count_matrix[row_lf.bin_id-1,row_rt.bin_id-1] +=1
                count_matrix[row_rt.bin_id-1,row_lf.bin_id-1] +=1

    #compute the mean of read counts/zscores in a window bin
    mean_adj_matrix=a_matrix/(count_matrix+0.0001)
    return mean_adj_matrix.copy(), a_matrix.copy(), count_matrix.copy(),left_df.copy(),right_df.copy(), selected_left_df2.copy(), selected_right_df2.copy()


def export_data_matrix_and_genomw_bin_file(cell, out_data_path,mean_adj_matrix, count_matrix, left_df, bin_str):
    '''
    Export the mean adj-matrix and count matirx to files
    and also export the window bin position file for the chromsome in bed format
    Input: cell, is chromosome name such as chr1
           count_matirx and mean_adj_matirx, are matrix for exporting
           out_data_path, is output file path
           bin_str , is the bin size string for output file such as 500KB_bin
           left_df, a complete left node dataframe which will be used to obtain all window bins in a chromosome
    '''
    out_matrix_file=os.path.join(out_data_path, cell +'_meanZscore.matrix')
    out_count_file=os.path.join(out_data_path, cell+ '_totalCount.matrix')
    print('Export matrix to file: ')
    print(out_matrix_file)
    print(out_count_file)
    np.savetxt(out_matrix_file, mean_adj_matrix, delimiter='\t',fmt='%10.5f')
    np.savetxt(out_count_file, count_matrix, delimiter='\t', fmt='%10.5f')

    #export genome window bin file for this chrom
    #find all genome window bin for this chromosome
    all_bins_df=left_df[['chrom','bin_start','bin_end','bin_id','numeric_chrom']].drop_duplicates().copy()
    if cell=='chr23':
        cell2='chrX'
    else:
        cell2=cell
    chrom_bins_df=all_bins_df[all_bins_df.chrom==cell2].copy()
    chrom_bins_df.sort_values(by=['numeric_chrom','bin_start'],inplace=True)
    out_bin_file=os.path.join(out_data_path, cell+ '_' + bin_str+ '_regions.tsv')
    print('Export genome window bin positino file :')
    print(out_bin_file)
    chrom_bins_df.to_csv(out_bin_file,sep='\t',index=False, header=None)
    return out_matrix_file, out_count_file, out_bin_file


def main(in_data_folder,
         cohort,
         chromosome,
         resolution,
         out_data_folder,
         genome_version,
         fig_dpi):
    #map interactions to window bins
    bin_str = str(int(resolution/1000)) + 'kb'
    #in_bin_region_file='genome/hg19/hg19_XY.'+bin_str+'.windows_bin.bed'
    in_bin_region_file = os.path.join(in_data_folder, genome_version, genome_version + '_XY.' + bin_str + '.windows_bin.bed')
    #make window bin bed file
    if not os.path.exists(in_bin_region_file):
        Preprocess_window_bin_bed.main(in_data_folder, resolution)
    #input Homer exported significant interactions
    homer_data_path = out_data_folder + '/hic_data/' + bin_str + '/hic_interaction_homer/' + cohort
    print(homer_data_path)
    #set output path
    out_data_path = out_data_folder + '/hic_data/' + bin_str + '/hic_interaction_bed/' + cohort
    if not os.path.exists(out_data_path):
        os.makedirs(out_data_path)
    #build adj-matrix
    if chromosome == 'whole_genome':
        chrom_strs = ['chr'+str(i) for i in range(1,24)]
    else: 
        chrom_strs = chromosome
    for chrom_str in chrom_strs:
        f1 = os.path.join(homer_data_path, cohort + '_significantInteractions_norm'+ bin_str + '_' + chrom_str + '.txt')
        #input data
        print(f1)
        left_df, right_df, out_left_file, out_right_file = read_homerInteraction2df(f1, out_data_path)
        #map interactions to window bins
        in_left_chrom_file = out_left_file
        in_right_chrom_file = out_right_file
        out_left_chrom_file, out_right_chrom_file =  map_left_right_node2defined_window_bins(in_bin_region_file, in_left_chrom_file, in_right_chrom_file, bin_str)
    
        #build adj-matrix
        left_region_file=out_left_chrom_file
        right_region_file=out_right_chrom_file
        mean_adj_matrix, a_matrix, count_matrix, left_df, right_df, selected_left_df2, selected_right_df2= build_adj_matrix4mapped_interactions(left_region_file, right_region_file, chrom_str)
    
        #draw heatmap
        matrix= mean_adj_matrix
        matrix2=matrix.copy()
        clip_value=6
        matrix2[matrix2>clip_value]=clip_value
        matrix2[matrix2<-clip_value]=-clip_value
        cell=chrom_str
        matrix_start=0
        matrix_end=matrix.shape[0]
        color_start=-clip_value
        color_end=clip_value
        bar_start=0
        bar_end=0
        matrix2 = np.triu(matrix2)
        draw_plot2(matrix2, cell, matrix_start, matrix_end, color_start, color_end, bar_start, bar_end, out_data_path, fig_dpi=fig_dpi)
        #export data matrix for mean adj-matrix and count-matrix
        out_matrix_file, out_count_file, out_bin_file= export_data_matrix_and_genomw_bin_file(cell, out_data_path, mean_adj_matrix, count_matrix, left_df, bin_str)
    return  out_matrix_file, out_count_file, out_bin_file  
















