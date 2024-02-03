import os
import glob
import pandas as pd
import numpy as np
import multiprocessing as mp
from dnaici.tools import data_functions
from dnaici.tools.plot_tools import draw_plot2


def add_TvalSign2log10Pval(in_df,tval_str,cluster_str='',is_transpose=True,min_p=10**(-6)):
    ''' in_df is a dataframe contains both tval and pval for each feature
        tval_str is 'tval' for return tvalues, if it is 'pval' then return tval signed abs(log10(pval))
        cluster_str is additional string that need to be added in the search 
    ''' 
    if is_transpose:
        new_df=in_df.iloc[:,in_df.columns.str.contains(cluster_str+tval_str)].T.copy()
    else:
        new_df=in_df.iloc[:,in_df.columns.str.contains(cluster_str+tval_str)].copy()
    if tval_str=='pval':
        #min_p=10**(-6)
        if is_transpose:
            tval_df=in_df.iloc[:,in_df.columns.str.contains(cluster_str+'tval')].T.copy()
        else:
            tval_df=in_df.iloc[:,in_df.columns.str.contains(cluster_str+'tval')].copy()
        tval_sign_df=tval_df.apply(lambda x: x/np.abs(x)).copy()
        pval_df=new_df.apply(lambda x:np.abs( np.log10(x+min_p))).copy()
        if not pval_df.empty:
            #convert pval to abs(log10(pval)) and take a negative or positive sign from Tvals
            out_df=pd.DataFrame(data=pval_df.values*tval_sign_df.values,columns=pval_df.columns,index=pval_df.index).copy()
        else:
            out_df=new_df.copy()
    else:
        out_df=new_df.copy()
        
    return out_df.copy()


def plot_heatmap4tval(in_df,color_clip_value,colormap, fig_out_name, fig_out_path,fig_dpi=100,tval_str='tval'):
    ''' Plot a heatmap for exported Tvalues from permutation test
        where columns of in_df with tval are T-values
        color_clip_value is the maximum value allowed in heatmap
        colormap is a color map for heat map such as GWR, or BWY
        tval_str is used to select plot either tval or pval of a dataframe
    '''
    ##in_file='mcf7_hic/mcf7_tcc/t1/out/out_network_cluster/../out_figs_maxWeight/t1/chr11/chr11_1000_loopPgt_0.01_in_clusters_edges_rmGlob_0percent_noWeight_tval.tsv'
    ##in_df=pd.read_csv(in_file,sep='\t',index_col=0)
    #new_df=in_df.iloc[:,in_df.columns.str.contains(tval_str)].T.copy()
    #if tval_str=='pval':
    #   min_p=0.000001
    #   tval_df=in_df.iloc[:,in_df.columns.str.contains('tval')].T.copy()
    #   tval_sign_df=tval_df.apply(lambda x: x/np.abs(x)).copy()
    #   pval_df=new_df.apply(lambda x:np.abs( np.log10(x+min_p))).copy()
    #   #convert pval to abs(log10(pval)) and take a negative or positive sign from Tvals
    #   new_df=pd.DataFrame(data=pval_df.values*tval_sign_df.values,columns=pval_df.columns,index=pval_df.index).copy()
    new_df=add_TvalSign2log10Pval(in_df,tval_str)
    #convert matrix to numeric 
    matrix=new_df.to_numpy().astype(float)
    clip_value=float(color_clip_value)
    matrix[matrix>clip_value]=clip_value
    matrix[matrix<-clip_value]=-clip_value
    matrix_start=0
    matrix_end=matrix.shape[0]
    color_start=-clip_value
    color_end=clip_value
    cell=fig_out_name 
    start=matrix_start
    end=matrix_end
    bar_start=0
    bar_end=0
    line_location=[]
    xticks=new_df.columns.str.replace('_','').to_list()
    yticks=new_df.index.str.replace('_cn','C').str.replace('_'+tval_str,'').to_list()
    out_fig=draw_plot2(matrix, cell, start, end, color_start, color_end, bar_start, bar_end,fig_out_path, colors=colormap,plotXtick=xticks,plotLines=line_location,plotYtick=yticks,fig_dpi=fig_dpi)
    
    return out_fig


def find_geneInfo4cluster(in_cluster_data_files, in_chrom_gene_file, id_sep_string):
    ''' find gene information for network clusters.
        Input: chrom-sepcific network clustering results files, chrom-specific gene expression file, 
            separate string for edge/interaction ID such as : for node1:node2
        output:
            a dictionary of dataframes for clusters of a network with gene name added in each cluster
            a dictionary of dataframes for genes and their expressions in a cluster
    '''
    #in_cluster_data_files=glob.glob(os.path.join(cluster_path0,'*'+percent_str+'*_data.tsv'))
    cluster_label=['_cn1_','_cn2_','_cn3_','_cn4_','_cn5_','_cn6_','_cn7_','_cn8_','_cn9_','_cn10_',
                   '_cn11_','_cn12_','_cn13_','_cn14_','_cn15_','_cn16_','_cn17_','_cn18_']
    #----- START ----
    #loop in a number of network clusters from chrom ci
    all_dfs={}
    ci_id2file={}
    for fi in in_cluster_data_files:
        print(fi, '\n')
        tmp_df=pd.read_csv(fi,sep='\t',index_col=0)
        tmp_df['id']=tmp_df.index
        #find cluster label
        ci_id=[ci for ci in cluster_label if ci in fi]
        if len(ci_id)==1:
            all_dfs[ci_id[0]]=tmp_df.copy()
            ci_id2file[ci_id[0]]=fi
        else:
            print(ci_id, '\n')
            input('Error in file name Stop ')
    #in each network cluster find its genes and do enrichment test
    #input gene info/expression file
    #in_chrom_gene_file=glob.glob(os.path.join('mcf7_expression/bead_array_exp',data_folder_str,time_str, 'out_data',ci+'_*_geneExp.tsv'))
    if len(in_chrom_gene_file)==1:
        gene_df=pd.read_csv(in_chrom_gene_file[0],sep='\t',header=None)
        gene_df.columns=['chrom','chrom_start','chrom_end','chrom_region','numerical_chrom','chrom_gene','gene_start','gene_end','expression','gene_name',
                         'gene_id','strand','ensemble','reads','gene_length']
        gene_df['gene_expression']=gene_df.gene_name.astype(str) +':'+gene_df.expression.astype(str)
    else:
        print(in_chrom_gene_file, '\n')
        input('Error in gene file')
    
    print('Export: \n')
    #loop in network clusters to get gene info 
    sep_str=id_sep_string #':'
    new_all_dfs={}
    new_uqGenes_in_cluster={}
    for ki in all_dfs.keys():
        print(ki, '\n')
        tmp_df=all_dfs[ki]
        tmp_gene_list=[]
        for idx,row in tmp_df.iterrows():
            row_id=[int(i) for i in row.id.split(sep_str)]
            tmp_gene=gene_df[gene_df.chrom_region.apply(lambda x : x in set(row_id))].gene_expression.to_list()
            tmp_gene_list.append([row.id, tmp_gene])
        
        tmp_gene_df=pd.DataFrame(data=tmp_gene_list)
        tmp_gene_df.columns=['id','gene_expression']
        new_all_dfs[ki]=pd.merge(tmp_df,tmp_gene_df,on='id',how='left').copy()
        #remove empty cell
        new_all_dfs[ki].gene_expression=new_all_dfs[ki].gene_expression.apply(lambda x: [i for i in x if i!='.'+sep_str+'.'])
        out_file=ci_id2file[ki].replace('_data.tsv','_data_info.tsv')
        #export new data gene matrix
        new_all_dfs[ki].to_csv(out_file,sep='\t',index=False)
        print(out_file, '\n')
        #export uniq gene name in this cluster
        out_gene_file=ci_id2file[ki].replace('_data.tsv','_uniq_gene.tsv')
        print(out_gene_file, '\n')
        tmp_genes= new_all_dfs[ki].gene_expression.to_list()
        t=[]
        for i in tmp_genes:
            t=t+i
        tmp_uq_gene=list(set(t))
        tmp_out_df=pd.DataFrame(data=tmp_uq_gene)
        tmp_out_df2=tmp_out_df[0].str.split(sep_str,expand=True)
        tmp_out_df2=tmp_out_df2[tmp_out_df2[0]!='.'].copy()
        tmp_out_df2[1]=tmp_out_df2[1].astype(float)
        tmp_out_df2.sort_values(by=1,ascending=False,inplace=True)
        tmp_out_df2.to_csv(out_gene_file,sep='\t',header=None,index=False)  
        new_uqGenes_in_cluster[ki]=tmp_out_df2.copy() 
    
    return new_all_dfs.copy() , new_uqGenes_in_cluster.copy()


def row_calculation(x1,x2,tmp_vect,isMAX=True):
    ''' caculate mean or a weight between two values in tmp_vect
        x1 , x2 are index (start from 1) of values in tmp_vect
    '''
    #here the calculation shall be the same as calculate_weight_of_pairwise_interaction()
    # and  calculate_mean_of_pairwise_interaction() in  make_hic_map4weight.py 
    tmp_list=np.array([tmp_vect[int(x1)-1],tmp_vect[int(x2)-1]]).sum()
    if isMAX:
        tmp_out=np.max([0, tmp_list])
    else:
        tmp_out=tmp_list/2
    
    return tmp_out


def find_genome_wide_data4features(in_data_folder, bin_str, chrom_strs, time_str, feature_str):
    ''' find genome wide data (max or mean) of values in an edge for selected 8 features
        here, we assume pair-wise interactions for all available nodes
        Input, chrom information
        Output, dictionary of all features in dataframe , and dictionary of all pair-wise nodes features in dataframe
    '''
    #read all data files without interaction matrix, here all features are available even if there is not an interaction in Hi-C adj-matrix
    record_data_dict={}
    #feature_str=['_geneExp_','_DNas_','_h3k4me1_','_h3k4me3_','_h3k27ac_','_h3k27me3_','_h3k9me3_','_ctcf_']
    #initialize dictionary
    record_features_dict={}
    for fs in feature_str:
        record_features_dict[fs]=[]
    #loop in chroms
    for chrom_str in chrom_strs:
        #read full data for each feature
        record_data_dict[chrom_str] = data_functions.read_data_files_without_interaction(in_data_folder, bin_str, chrom_str, time_str)
        tmp_keys = record_data_dict[chrom_str].keys()
        #loop in each feature
        for tk in tmp_keys:
            #find corresponding feature key
            selected_feature = [fi for fi in feature_str if fi in tk]
            if len(selected_feature) == 1:
                #make a pair-wise interaction matrix for selected feature vector
                region_df = record_data_dict[chrom_str][tk].copy()
                #make a matrix index for all updiagnol pair-wise interaction positions
                #assume the fourth columns is the bin window reginos for each chromosome
                len_of_bins = region_df.iloc[-1,3]+1
                product_lists = []
                for i in range(1,len_of_bins):
                    for j in range(i,len_of_bins):
                        product_lists.append((region_df.chr[i-1],i,j))
                #make a pair-wise position index with three columns, chrom, bin_i, bin_j
                region_pair_df = pd.DataFrame(data=product_lists).copy()
                #compute a max or mean values in a pair-wise interaction
                #and assume the feature values are stored in mean_DNas column
                tmp_vect = region_df.mean_value.to_numpy()
                region_pair_df[3] = region_pair_df.apply(lambda x: row_calculation(x[1],x[2],tmp_vect),axis=1)
                region_pair_df[4] = region_pair_df.apply(lambda x: row_calculation(x[1],x[2],tmp_vect,isMAX=False),axis=1)
                region_pair_df.columns = ['bin_chrom','bin_i','bin_j','bin_max','bin_mean']
                record_features_dict[selected_feature[0]] += [region_pair_df.copy()]
            else:
                print(selected_feature, '\n')
                input('Error ?? ')
                pass
        
    out_record_features_dict = {}
    for ki in record_features_dict.keys():
        out_record_features_dict[ki] = pd.concat(record_features_dict[ki].copy())
    
    return record_data_dict, out_record_features_dict


def remove_str4cluster_label(x):
    ''' here we assume x =_cn1_ , _cn2_, .. that contains _cn#_ then
        return an index by a predifined dictionary
    '''
    cluster_index={'_cn1_':1, '_cn2_':2, '_cn3_':3, '_cn4_':4, '_cn5_':5, '_cn6_':6, '_cn7_':7,'_cn8_':8, '_cn9_':9, '_cn10_':10,
                   '_cn11_':11, '_cn12_':12, '_cn13_':13, '_cn14_':14, '_cn15_':15, '_cn16_':16, '_cn17_':17, '_cn18_':18}
    #tmp_x=x.replace('_cn','').replace('_','')
    if x in cluster_index.keys():
        tmp_x=cluster_index[x]
    else:
        tmp_x=x
    #print(tmp_x)
    return tmp_x


def do_permutation_test4features_in_cluster(ci, in_cluster_data_files,out_record_features_dict, chrom_cluster_dict,genome_feature_str,feature_str, loop4test, pval_cutoff):
    ''' Perform permutation test for features in a cluster against randomly selected samples from genome-wide data
        Input:
            ci, is the chromsome
            in_cluster_data_files, a list of network clustering results files
        chrom_cluster_dict, is generated by function find_geneInfo4cluster(), where a dataframe of all features in each network cluster are stored in a dictionary
        ome_feature_str, is a list of all genome features in each network cluster such as ['_geneExp_','_nucleoDens_','_h3k4me1_','_h3k4me3_','_h3k27ac_','_h3k27me3_','_h3k9me3_','_ctcf_']
        ture_str, is a string for indicating a computtion of values between two nodes in an edge such as maxWeight, or mean
        p4test, is the number of permutations will be done in order to estimate an expected p-values for ttest enrichment of a feature between the cluster and randomly sampled data.
        l_cutoff, is a p-value cutoff value to count the number of ttest with pval greater than this pval_cuttoff,
        Output:
            a dataframe contains all expected pvals and tvals from the permutation tests for features in network clusters. 
            Here the smaller the expected pval indicates the more significant the feature enrichment in a cluster when compared against randomly drawed samples. 	 
    '''
    #loop4test=1000
    #pval_cutoff=0.01
    record_ki_dfs=[]
    ttest_column_str=['tval','pval']
    #sort network cluster number
    keys2chrom_cluster=list(chrom_cluster_dict.keys())
    #print(keys2chrom_cluster)
    keys2chrom_cluster.sort(key=remove_str4cluster_label)
    for ki in keys2chrom_cluster:
        print(ki, '\n')
        tmp_df=chrom_cluster_dict[ki].copy()
        tmp_df_ids=set(tmp_df.id.to_list())
        record_fi_tval_pval={}
        num_of_processes=len(genome_feature_str)
        pool = mp.Pool(processes=num_of_processes)
        #parallel loop in features
        tmp_record_fi_tval_pval = pool.map(data_functions.call_permutation_test4a_feature,[(ci,genome_feature_str[loop],tmp_df.copy(),tmp_df_ids,out_record_features_dict.copy(),loop4test,pval_cutoff,feature_str,ttest_column_str) for loop in range(0,num_of_processes)],1)
        pool.close()      
        #combine all features to one dictionary
        for ti in tmp_record_fi_tval_pval:
            record_fi_tval_pval |= ti
        #convert dictionary to dataframe
        record_fi_df=pd.DataFrame(data=record_fi_tval_pval,dtype=float).copy()
        record_fi_df=record_fi_df.T
        record_fi_df.iloc[:,0]=record_fi_df.iloc[:,0].map('{:4.2f}'.format)
        record_fi_df.iloc[:,1]=record_fi_df.iloc[:,1].map('{:5.6f}'.format)
        record_fi_df.columns=[ki+ttest_column_str[0], ki+ttest_column_str[1]]
        record_ki_dfs.append(record_fi_df.T.copy())
    #combine the dataframes of all clusters
    all_ki_df=pd.concat(record_ki_dfs,axis=0).T.copy()
    #export results for the chromosome
    out_path0,out_file0= os.path.split(in_cluster_data_files[0])
    out_file=os.path.join(out_path0, ci+'_'+str(loop4test) +'_loopPgt_'+str(pval_cutoff)+ '_in_clusters_' + '_'.join(out_file0.split('_')[3:]))
    out_file=out_file.replace('_data.tsv','_tval.tsv')
    all_ki_df.to_csv(out_file,sep='\t')
    print(out_file, '\n')
    #plot heat map
    color_clip_value=20
    colormap='BWY'
    fig_dpi=300
    fig_out_path=os.path.split(out_file)[0]
    fig_out_name=os.path.basename(out_file).replace('.tsv','').replace('_in_clusters','')
    out_fig_file= plot_heatmap4tval(all_ki_df,color_clip_value,colormap,fig_out_name, fig_out_path,fig_dpi=fig_dpi)
    
    return all_ki_df.copy(),out_file,out_fig_file


def main(in_data_folder, 
         cohort, 
         chromosome, 
         resolution, 
         out_data_folder,
         permutation,
         pval_cutoff,
         fig_dpi):
    #file parameters
    bin_str = str(int(resolution/1000)) + 'kb'
    
    if chromosome == 'whole_genome':
        chrom_strs = ['chr'+str(i) for i in range(1,24)]
    else: 
        chrom_strs = chromosome
        
    percent_str = '0'
    feature_str ='maxWeight'
    #feature_str='mean'
    '''
    isFind_global_cutoff=False
    if isFind_global_cutoff:
        edges_rm_str = '_edges_rmGlob_'
    else:
        edges_rm_str = '_edges_rm_'

    in_path = in_data_folder + '/hic_data/' + bin_str

    hic_f_path = in_path + '/hic_interaction_bed/' + cohort
    
    #hic file name postprefix string
    '''
    hic_f_postprefix_file_name = '_meanZscore.matrix'
    
   
    chrom_region_postprefix_file_name = '_'+ bin_str +'_regions.tsv'
    #here only one of them to be True
    isExport_GW_features = True
    isREAD_GW_features = False
    #out path for combined features by genome-wide pair-wise interactions
    out_path1 = out_data_folder + '/hic_data/' + bin_str + '/hic_combined_features/' + cohort    
    if not os.path.exists(out_path1):
        os.makedirs(out_path1)
        print('Create output folder for combined features: \n', out_path1, '\n')
    
    genome_feature_str = ['_geneExp_','_nucleoDens_','_ctcf_','_h3k4me1_','_h3k4me3_','_h3k27ac_','_h3k27me3_','_h3k9me3_']
    #read all data files without interaction matrix, here all features are available even if there is not an interaction in Hi-C adj-matrix
    
    if False:
        record_matrix={}
        record_region_pair_df={}
        record_files={}
        for chrom_str in chrom_strs:
            print(chrom_str, '\n')
            hic_f_path0 = in_data_folder + '/hic_data/' + bin_str + '/hic_interaction_bed/' + cohort
            hic_f = os.path.join(hic_f_path0, chrom_str + hic_f_postprefix_file_name)
            chrom_region_file = os.path.join(hic_f_path0, chrom_str + chrom_region_postprefix_file_name)
            print(hic_f, '\n')
            print(chrom_region_file, '\n')
            all_matrix, new_region_pair_df, all_files = data_functions.read_all_data_files(out_data_folder, bin_str, feature_str, chrom_str, cohort, hic_f, chrom_region_file)
            record_matrix[chrom_str]=all_matrix.copy()
            record_region_pair_df[chrom_str]=new_region_pair_df.copy()
            record_files[chrom_str]=all_files.copy()
        input('Click for continute ...')


    if isExport_GW_features:
        print('Compute and export genome-wide features: \n')
        # This pair-wise nodes interaction for features in genome-wise interactions are the same for all types called signifncat interactions if
        # they are in the same window bin size such as 500kb bins and unders the same condition such as t0. Thus, it only needs to be generated 
        # once for all the same window bins in the same condition
        record_data_dict = {}
        record_data_dict, out_record_features_dict = find_genome_wide_data4features(in_data_folder, bin_str, chrom_strs, cohort, genome_feature_str)
        #export combined features for pair-wise nodes to files
        for ki in out_record_features_dict.keys():
            out_file0 = ki + cohort + '_nodes_pair' + chrom_region_postprefix_file_name
            out_file = os.path.join(out_path1, out_file0)
            out_record_features_dict[ki]['id'] = out_record_features_dict[ki].bin_i.astype(str) + ':' + out_record_features_dict[ki].bin_j.astype(str)
            out_record_features_dict[ki].to_csv(out_file, sep='\t', index=False)
            print(out_file, '\n')
    
    if isREAD_GW_features:
        print('Read genome-wide features: \n')
        out_record_features_dict = {}
        for fi in genome_feature_str:
            in_file0 = fi + cohort + '_nodes_pair' + chrom_region_postprefix_file_name
            in_file = os.path.join(out_path1, in_file0)
            print(in_file, '\n')
            out_record_features_dict[fi] = pd.read_csv(in_file, sep='\t') 
            out_record_features_dict[fi]['id'] = out_record_features_dict[fi].bin_i.astype(str)+ ':'+ out_record_features_dict[fi].bin_j.astype(str)
    #find gene information for bins
    #loop in all chromosomes
    record_out_df = {}
    record_out_file = {}
    for ci in chrom_strs:
        #ci=chroms_str[0]
        print(ci, '\n')
        #results of network clustering path
        #here cluster data file has to be generated by "read_clustered_network_homer.py" first!
        cluster_path = in_data_folder + '/hic_data/' + bin_str + '/hic_community_figures/' + cohort + '/' + ci
        in_cluster_data_files = glob.glob(os.path.join(cluster_path, '*'+percent_str+'*_data.tsv'))
        chrom_path = in_data_folder + '/expression_data/' + bin_str + '/out_data/' + cohort
        in_chrom_gene_file = glob.glob(os.path.join(chrom_path, ci+'_*_geneExp.tsv'))
        id_sep_string = ':'
        print('Read gene information in clusters: \n')
        print(in_chrom_gene_file, '\n')
        chrom_cluster_dict, chrom_cluster_uqGene_dict = find_geneInfo4cluster(in_cluster_data_files,in_chrom_gene_file ,id_sep_string )
        #print(chrom_cluster_dict)
        #print(chrom_cluster_uqGene_dict)
        #loop in each cluster do permutation test for each featurs
        if len(in_cluster_data_files) > 0:
            all_ki_df, out_file, out_fig_file = do_permutation_test4features_in_cluster(ci, in_cluster_data_files, out_record_features_dict, chrom_cluster_dict, genome_feature_str, feature_str, permutation, pval_cutoff)
            record_out_df[ci] = all_ki_df.copy()
            record_out_file[ci] = [out_file,out_fig_file]
        else:
            print('Network clustering results are not found for ', ci)
            record_out_df[ci] = []
            record_out_file[ci] = ''
            #input('Click to confinue ...')

