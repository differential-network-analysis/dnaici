import os
import glob
import pandas as pd
import numpy as np
from ast import literal_eval
from dnaici.preprocess import Preprocess_window_bin_bed

def importData(file_DIEGs, file_DIGs, in_folder):
    
    DIEGs = pd.read_csv(file_DIEGs)
    print('Read file:', file_DIEGs, '\n')
    DIGs = pd.read_csv(file_DIGs)
    print('Read file:', file_DIGs, '\n')
    
    DIEGs_in_chr = DIGs[DIGs.gene.isin(DIEGs.gene)]
    DIEGs_in_chr = DIEGs_in_chr.reset_index(drop=True)
    
    loc = pd.DataFrame(columns=['starts_x','ends_x'])
    for ind in DIEGs_in_chr.index:
        gene, chrom = DIEGs_in_chr['gene'][ind], DIEGs_in_chr['chrom'][ind]
        if chrom == 'chrX':
            chrom = 'chr23'
            
        file = in_folder + '/%s/common_genes.tsv' %chrom
        common_genes = pd.read_csv(file, sep='\t')
        print('Read file:', file, '\n')
        
        tmp = common_genes[common_genes.gene==gene][['starts_x','ends_x']]
        loc = pd.concat([loc,tmp], ignore_index=True)
    
    DIEGs_in_chr[['start','end']] = loc[['starts_x','ends_x']]
    DIEGs_in_chr = DIEGs_in_chr.replace(to_replace='chrX', value='chr23')
    
    return DIEGs_in_chr


def findNodes(DIEGs_in_chr, in_data_folder, out_data_folder, resolution, cohort1, cohort2, pval_cutoff):
        
    file_bin = in_data_folder + '/hg19/hg19_XY.'+str(int(resolution/1000))+'kb.windows_bin.BlackListFiltered.bed'

    if not os.path.exists(file_bin):
        Preprocess_window_bin_bed.main(in_data_folder, resolution)
    
    bin_df = pd.read_csv(file_bin, sep='\t', names=['chrom','start','end','node','chrom_int'])
    print('Read file:', file_bin, '\n')
    bin_df = bin_df.replace(to_replace='chrX', value='chr23')
    
    file_sigNodes = out_data_folder + '/%s_vs_%s_pval_ls'%(cohort1, cohort2)+str(pval_cutoff)+'_selected_sigNodes.tsv'
    sigNodes_df = pd.read_csv(file_sigNodes, sep='\t')
    print('Read file:', file_sigNodes, '\n')
    sigNodes_df = pd.DataFrame(sigNodes_df['id'].str.split('_', expand=True).values, columns=['chrom', 'nodes'])
    sigNodes_df['nodes'] = sigNodes_df['nodes'].astype(int)
    sigNodes_df['nodes'] = sigNodes_df['nodes']+1
    
    node = []
    for ind in DIEGs_in_chr.index:
        chrom, start, end = DIEGs_in_chr['chrom'][ind], DIEGs_in_chr['start'][ind], DIEGs_in_chr['end'][ind]
        sub_bin_df = bin_df[bin_df.chrom.isin([chrom])]
        if len(sub_bin_df[sub_bin_df['start']-start <= 0]) > 0 and len(sub_bin_df[sub_bin_df['end']-end >= 0]) > 0:
            min_node = sub_bin_df[sub_bin_df['start']-start <= 0]['node'].values[-1]
            max_node = sub_bin_df[sub_bin_df['end']-end >= 0]['node'].values[0]
            #convert bin defined node index to network defined network
            n = np.arange(min_node, max_node+1, 1)
            
            sub_sigNodes_df = sigNodes_df[sigNodes_df.chrom.isin([chrom])]
            node.append(np.intersect1d(n, sub_sigNodes_df.nodes))
        else:
            node.append([])
    
    DIEGs_with_nodes = DIEGs_in_chr.copy()
    DIEGs_with_nodes['nodes'] = node
    
    return DIEGs_with_nodes


def countNodes(DIEGs_with_nodes, chrom, label_file, bin_str, cohort1, cohort2):
    
    sub_DIEGs_with_nodes = DIEGs_with_nodes[DIEGs_with_nodes.chrom.isin([chrom])]
    
    nodes = []
    for ind in sub_DIEGs_with_nodes.index:
        if len(sub_DIEGs_with_nodes['nodes'][ind])>0:
            for n in sub_DIEGs_with_nodes['nodes'][ind]:
                if n not in nodes:
                    nodes.append(n)
    
    nodes_for_DIEGs = pd.DataFrame({'nodes':nodes})
    
    label_file_1 = label_file + '/hic_data/'+ bin_str + '/hic_community_data/' + cohort1 + '/%s_edges_rm_0percent_noWeight_4communities_class.tsv'%chrom
    label_file_2 = label_file + '/hic_data/'+ bin_str + '/hic_community_data/' + cohort2 + '/%s_edges_rm_0percent_noWeight_4communities_class.tsv'%chrom
    
    label_df_1 = pd.read_csv(label_file_1, sep='\t')   
    label_df_1['nodes'] = label_df_1['nodes']+1
    label_df_1['labels'] = label_df_1['labels']+1
    label_df_2 = pd.read_csv(label_file_2, sep='\t')   
    label_df_2['nodes'] = label_df_2['nodes']+1
    label_df_2['labels'] = label_df_2['labels']+1
    
    nodes_with_label_1 = label_df_1[label_df_1.nodes.isin(nodes_for_DIEGs['nodes'])]
    nodes_with_label_2 = label_df_2[label_df_2.nodes.isin(nodes_for_DIEGs['nodes'])]
    
    return nodes_for_DIEGs, nodes_with_label_1, nodes_with_label_2


def geneInNode(DIEGs_with_nodes, nodes_for_DIEGs):
    
    genes = []
    for i in nodes_for_DIEGs.nodes:
        g = []
        for j in range(len(DIEGs_with_nodes)):
            if i in DIEGs_with_nodes.nodes[j]:
                g.append(DIEGs_with_nodes.gene[j])
        genes.append(g)    
    
    return genes


def networkMatrix(nodes_for_DIEGs, nodes_pairs):
    
    n = len(nodes_for_DIEGs)
    adj_mat = pd.DataFrame(columns=['node1','node2','edge'])
    
    for i in range(n-1):
        for j in range(i+1,n):
            n1 = nodes_for_DIEGs.nodes[i]
            n2 = nodes_for_DIEGs.nodes[j]
            if str(n1)+':'+str(n2) in nodes_pairs or str(n2)+':'+str(n1) in nodes_pairs:
                edge = 1
            else:
                edge = 0
            tmp = pd.DataFrame(data={'node1':[n1],'node2':[n2],'edge':[edge]})
            adj_mat = pd.concat([adj_mat,tmp], ignore_index=True)
                
    return adj_mat


def countFeature(out_data_folder, bin_str, cohort, chrom, feature_str1, feature_str2, nodes_with_label, genes):
    
    feature_df1_DIEGs = pd.DataFrame({'feature1':[[]]*len(nodes_with_label)})
    feature_df2_DIEGs = pd.DataFrame({'feature2':[[]]*len(nodes_with_label)})
    # enhancer
    for feature in feature_str1:
        feature_file =  out_data_folder + '/histone_data/' + bin_str + '/out_data/' + cohort + '/%s_%s_regions_'%(chrom, bin_str) + feature + '_array.bed'
        print('Read file:', feature_file, '\n')
        feature_df_f1 = pd.read_csv(feature_file, sep='\t', converters={'value': literal_eval})
        tmp_f1 = pd.DataFrame({'feature1':feature_df_f1[feature_df_f1.region.isin(nodes_with_label.nodes)]['value'].values})
        feature_df1_DIEGs = feature_df1_DIEGs + tmp_f1
    # repressor
    for feature in feature_str2:
        feature_file =  out_data_folder + '/histone_data/' + bin_str + '/out_data/' + cohort + '/%s_%s_regions_'%(chrom, bin_str) + feature + '_array.bed'
        print('Read file:', feature_file, '\n')
        feature_df_f2 = pd.read_csv(feature_file, sep='\t', converters={'value': literal_eval})
        tmp_f2 = pd.DataFrame({'feature2':feature_df_f2[feature_df_f2.region.isin(nodes_with_label.nodes)]['value'].values})
        feature_df2_DIEGs = feature_df2_DIEGs + tmp_f2
    # gene_expression
    feature_file =  out_data_folder + '/expression_data/' + bin_str + '/out_data/' + cohort + '/%s_%s_regions_'%(chrom, bin_str) + 'geneExp_array.bed'
    feature_df_f3 = pd.read_csv(feature_file, sep='\t', converters={'value': literal_eval})
    feature_df3_DIEGs = pd.DataFrame({'feature3':feature_df_f3[feature_df_f3.region.isin(nodes_with_label.nodes)]['value'].values})
    
    feature_df_DIEGs = pd.DataFrame()
    
    feature_df_DIEGs['label'] = nodes_with_label.labels.values
    
    feature_df_DIEGs['genes'] = genes
    
    feature_df_DIEGs['mean_enhancer'] = feature_df1_DIEGs.feature1.apply(np.mean)
    feature_df_DIEGs['max_enhancer'] = feature_df1_DIEGs.feature1.apply(np.max)
    
    feature_df_DIEGs['mean_repressor'] = feature_df2_DIEGs.feature2.apply(np.mean)
    feature_df_DIEGs['max_repressor'] = feature_df2_DIEGs.feature2.apply(np.max)
    
    feature_df_DIEGs['mean_expression'] = feature_df3_DIEGs.feature3.apply(np.mean)
    feature_df_DIEGs['max_expression'] = feature_df3_DIEGs.feature3.apply(np.max)
    
    feature_df_DIEGs.index=nodes_with_label.nodes
    
    return feature_df_DIEGs


def main(in_data_folder, 
         chromosome,
         resolution, 
         out_data_folder,
         cohort1,
         cohort2,
         pval_cutoff):
    
    bin_str = str(int(resolution/1000)) + 'kb'
    out_folder = out_data_folder + '/differential_network_analysis/' + bin_str
    
    if chromosome == 'whole_genome':
        chrom_strs = set(['chr'+str(i) for i in range(1,24)])
    else: 
        chrom_strs = set(chromosome)
    
    if glob.glob(os.path.join(out_folder, '*_DIEGs_'+'%s_resolution.tsv'%bin_str)) == []:
        print('No DIEGs found in the input chromosome :( Please move to another chromosome')
    else:
        file_DIEGs = glob.glob(os.path.join(out_folder, '*_DIEGs_'+'%s_resolution.tsv'%bin_str))[0]
        file_DIGs = glob.glob(os.path.join(out_folder, '*_DIGs_'+'%s_resolution.tsv'%bin_str))[0]
    
        DIEGs_in_chr = importData(file_DIEGs, file_DIGs, out_folder)
        DIEGs_with_nodes = findNodes(DIEGs_in_chr, in_data_folder, out_folder, resolution, cohort1, cohort2, pval_cutoff)
        
        chroms = set(DIEGs_with_nodes['chrom'])
        co_chroms = chroms & chrom_strs
        
        if len(co_chroms) == 0:
            print('No DIEGs found in the input chromosome :( Please move to another chromosome')
        else:
            print('Start analysis...', '\n')
            writer1 = pd.ExcelWriter(out_folder + '/network_matrix_%s.xlsx'%cohort1)
            writer2 = pd.ExcelWriter(out_folder + '/network_matrix_%s.xlsx'%cohort2)
            writer3 = pd.ExcelWriter(out_folder + '/feature_in_nodes_%s.xlsx'%cohort1)
            writer4 = pd.ExcelWriter(out_folder + '/feature_in_nodes_%s.xlsx'%cohort2)
            
            for chrom in chroms:
                
                print('investigating ' + chrom + ' ...', '\n')
                
                #label_file = 'mcf7_hic/mcf7_tcc/'
                nodes_for_DIEGs, nodes_with_label_1, nodes_with_label_2 = countNodes(DIEGs_with_nodes, chrom, out_data_folder, bin_str, cohort1, cohort2)
                genes = geneInNode(DIEGs_with_nodes, nodes_for_DIEGs)
                
                if len(nodes_for_DIEGs) > 0:
                    # output common nodes
                    out_file = out_folder + '/' + chrom + '/common_nodes_for_DIEGs.tsv'
                    print('Export file:', out_file, '\n')
                    nodes_for_DIEGs.to_csv(out_file, sep='\t', index=False)
                    # output edges between nodes for constructing networks
                    ## cohort 1
                    file_1 = out_folder + '/' + chrom + '/%s_only_nodes_pairs.tsv'%cohort1
                    print('Read file:', file_1, '\n')
                    nodes_pairs_1 = pd.read_csv(file_1)['%s_only_nodes_pairs'%cohort1].values
                    
                    adj_mat_1 = networkMatrix(nodes_for_DIEGs, nodes_pairs_1)
                    out_file1 = out_folder + '/' + chrom + '/network_matrix_%s.tsv'%cohort1
                    print('Export file:', out_file1, '\n')
                    adj_mat_1.to_csv(out_file1, sep='\t', index=True)
                    adj_mat_1.to_excel(writer1, sheet_name=chrom, index=False)
                    
                    ## cohort 2
                    file_2 = out_folder + '/' + chrom + '/%s_only_nodes_pairs.tsv'%cohort2
                    print('Read file:', file_2, '\n')
                    nodes_pairs_2 = pd.read_csv(file_2)['%s_only_nodes_pairs'%cohort2].values
                    
                    adj_mat_2 = networkMatrix(nodes_for_DIEGs, nodes_pairs_2)
                    out_file2 = out_folder + '/' + chrom + '/network_matrix_%s.tsv'%cohort2
                    print('Export file:', out_file2, '\n')
                    adj_mat_2.to_csv(out_file2, sep='\t', index=True)
                    adj_mat_2.to_excel(writer2, sheet_name=chrom, index=False)
            
                    # output genomic features (mean z-score) of nodes for coloring
                    feature_str1 = ['h3k4me1', 'h3k4me3', 'h3k27ac']
                    feature_str2 = ['h3k27me3', 'h3k9me3']
                    ## cohort 1
                    feature_1 = countFeature(out_data_folder, bin_str, cohort1, chrom, feature_str1, feature_str2, nodes_with_label_1, genes)
                    out_file1 = out_folder + '/' + chrom + '/feature_in_nodes_%s.tsv'%cohort1
                    print('Export file:', out_file1, '\n')
                    feature_1.to_csv(out_file1, sep='\t', index=True)
                    feature_1.to_excel(writer3, sheet_name=chrom)
                    ## cohort 2
                    feature_2 = countFeature(out_data_folder, bin_str, cohort2, chrom, feature_str1, feature_str2, nodes_with_label_2, genes)
                    out_file2 = out_folder + '/' + chrom + '/feature_in_nodes_%s.tsv'%cohort2
                    print('Export file:', out_file2, '\n')
                    feature_2.to_csv(out_file2, sep='\t', index=True)
                    feature_2.to_excel(writer4, sheet_name=chrom)
            # excel file finished
            writer1.close()
            writer2.close()
            writer3.close()
            writer4.close()

    print('Analysis finished! We recommend Cytoscape for drawing the network figures')


