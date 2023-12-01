import pandas as pd
import os
import glob
from dnaici.analysis.Enrichment_permutation import plot_heatmap4tval


def main(in_data_folder, 
         cohort, 
         chromosome, 
         resolution, 
         out_data_folder,
         permutation,
         fig_dpi):
        
    bin_str = str(int(resolution/1000)) + 'kb'

    if chromosome == 'whole_genome':
        chrom_strs = ['chr'+str(i) for i in range(1,24)]
    else: 
        chrom_strs = chromosome
    
    in_data_path = in_data_folder + '/hic_data/' + bin_str + '/hic_community_figures/' + cohort
    cluster_str = 'rm_0'
    
    for ci in chrom_strs:
        in_file = os.path.join(in_data_path, ci, ci+'_'+str(permutation)+'_loopPgt*'+cluster_str+'p*_tval.tsv')
        print(in_file, '\n')
        network_clusters_tval_file = glob.glob(in_file)
        
        print('Read files: ', '\n')
        print(network_clusters_tval_file, '\n')
        print('Export: ', '\n')
        for fi in network_clusters_tval_file:
            tmp_df = pd.read_csv(fi, sep='\t', index_col=0)	        
            #plot heat map
            colormap = 'BWY'
            tval_str = 'tval'
            fig_out_path = os.path.split(fi)[0]
            fig_out_name = os.path.basename(fi).replace('.tsv','').replace('_in_clusters','')
            fig_out_name0 = '_'.join(fig_out_name.split('_')[0:-1])
            fig_out_name1 = fig_out_name0 + '_' + tval_str
            color_clip_value = 20
            out_fig_file = plot_heatmap4tval(tmp_df, color_clip_value, colormap, fig_out_name1, fig_out_path, fig_dpi, tval_str=tval_str)
            
            tval_str = 'pval'
            fig_out_name1 = fig_out_name0+'_' + tval_str
            color_clip_value = 6
            out_fig_file = plot_heatmap4tval(tmp_df, color_clip_value, colormap, fig_out_name1, fig_out_path, fig_dpi, tval_str=tval_str)

    