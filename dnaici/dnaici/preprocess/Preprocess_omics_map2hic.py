import pandas as pd
import os
import numpy as np
import matplotlib as mlp
import matplotlib.pyplot as plt
import glob


def main(in_data_folder,
         cohort,
         chromosome,
         resolution,
         multi_omics,
         out_data_folder
         ):
        
    bin_str = str(int(resolution/1000)) + 'kb'
    
    if multi_omics == 'gene expression':
        in_strs = 'expression_data'
        out_strs = ['_geneExp']
    elif multi_omics == 'nucleosome density':
        in_strs = 'nucleosome_density_data'
        out_strs = ['_nucleoDens']
    elif multi_omics == 'histone marker':
        in_strs = 'histone_data'
        out_strs = ['_ctcf','_h3k27ac','_h3k27me3','_h3k4me1','_h3k4me3','_h3k9me3']
    else:
        out_strs = []
        print('No such data type found :( Please set parameter \'multi_omics\' as \'gene expression\', \'nucleosome density\', or \'histone marker\'')
        
    if chromosome == 'whole_genome':
        chrom_strs = ['chr'+str(i) for i in range(1,24)]
    else: 
        chrom_strs = chromosome
    
    for chrom_str in chrom_strs:
        #read region and multi-omics files
        for out_str in out_strs:
            out_str2 = out_str + '_array.bed'
            #read chromsome position of window bins
            hic_pos_f1 = out_data_folder + '/hic_data/' + bin_str + '/hic_interaction_bed/' + cohort + '/' + chrom_str + '_' + bin_str + '_regions.tsv'  
            #read mlti-omics data in 200b resolution
            
            pathname = in_data_folder + '/' + in_strs + '/*' + cohort + '*.bed'
            omics_f2 = glob.glob(pathname)[0]

            #top N percenage of singal will be takend the mean for pair-wise interactions
            topN_percent = 1.0
            out_path = out_data_folder + '/' + in_strs + '/' + bin_str + '/out_data/' + cohort
            if not os.path.exists(out_path):
                os.makedirs(out_path)
            #find multi-omics data for chrom regions
            out_file_postprefix = out_str + '.tsv'
            out_file = os.path.join(out_path, os.path.basename(hic_pos_f1).replace('.tsv', out_file_postprefix))
            print(out_file)
            #map to selected window bin regions 
            cmd = 'bedtools intersect -a ' + hic_pos_f1 + ' -b ' + omics_f2 + ' -wao >' + out_file
            result = os.system(cmd)
            if result != 0:
                print('Error in % ', cmd)
                exit(1)
            #read overlapping multi-omics in regions
            omics_df = pd.read_csv(out_file, sep='\s+', header=None, dtype='object')
            #read Hi-C position
            new_df = pd.DataFrame(columns=[0, 1, 2, 3, 4, 5, 6])
            region_df = pd.read_csv(hic_pos_f1, sep='\t', header=None)
            #for each Hi-C region find its omics data
            record_df_list=[]
            for index, rows in region_df.iterrows():
                indx = (omics_df[0] == rows[0]) & (omics_df[1].astype(int) == rows[1]) & (omics_df[2].astype(int) == rows[2])
                tmp_df = omics_df[indx].copy()
                tmp_df.replace('.', 0, inplace=True)
                topN = int(np.ceil(tmp_df[8].shape[0] * topN_percent))
                topN_array = np.sort(tmp_df[8].astype(float))[-topN:]
                rows[5] = topN_array.mean()
                rows[6] = tmp_df[8].astype(float).to_list()
                record_df_list.append(rows.to_frame().T.copy())

            new_df = pd.concat(record_df_list).copy()
            new_df.columns = ['chrom', 'starts', 'ends', 'region', 'chr', 'mean_value', 'value']

            #export regions DNas/geneExp/histnone data in selected window bin size
            out_file2 = out_file.replace(out_file_postprefix, out_str2)
            print(out_file2)
            new_df.to_csv(out_file2, sep='\t', index=False)










