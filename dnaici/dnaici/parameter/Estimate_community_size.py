import numpy as np
import pandas as pd
import scipy.stats
import matplotlib.pyplot as plt
import os


def importData(path, time, chrom):
    
    file_clsuter = os.path.join(path, '%s_edges_rm_0percent_noWeight_4communities_class.tsv' % chrom)
    print('Read:', file_clsuter, '\n')
    cluster = pd.read_csv(file_clsuter, sep='\t')
    
    file_edge = os.path.join(path, '%s_edges_rm_0percent_noWeight_4zscore.matrix' % chrom)
    print('Read:', file_edge, '\n')
    edge = pd.read_csv(file_edge, sep='\t', header=None)
    
    return cluster, edge


def countEdges(cluster, edge):
    
    clusters_labels = cluster.labels.unique()
    
    num_of_edges_in_cluster_dict = {}

    for i in clusters_labels:
        num_of_edges_in_cluster_dict[i]=0
    
    for idx, row in edge.iterrows():
        label2node_row0 = cluster.labels[cluster.nodes==row[0]].to_list()[0]
        label2node_row1 = cluster.labels[cluster.nodes==row[1]].to_list()[0]
        if label2node_row0 == label2node_row1:
            num_of_edges_in_cluster_dict[label2node_row0] += 1
    
    num_of_edges_in_cluster = []
    for e in num_of_edges_in_cluster_dict.values():
        if e > 0:
            num_of_edges_in_cluster.append(e)
    
    return num_of_edges_in_cluster


def combineEdges(path, time, chrom_strs):
    
    num_of_all_edges = []
    
    for chrom in chrom_strs:
        print('Diving into %s now :) ...' %chrom)
        cluster, edge = importData(path, time, chrom)
        num_of_edges_in_cluster = countEdges(cluster, edge)
        num_of_all_edges += num_of_edges_in_cluster
    
    return num_of_all_edges


def normalization(num_of_all_edges):
    
    normalized_all_edges = [x/np.mean(num_of_all_edges) for x in num_of_all_edges]
    
    return normalized_all_edges


def getCutoff(cutoff4Proportion, num_of_all_edges, normalized_all_edges):
    
    arr = np.asarray(normalized_all_edges)
    
    ind = (np.abs(arr - cutoff4Proportion)).argmin()
    
    return num_of_all_edges[ind]


def disPlot(normalized_all_edges, num_of_all_edges, cutoff4Proportion, cutoff, time_str, out_path, fig_dpi, chrom_len):
    
    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(6, 3))
    n_bins = 20
    
    axs[0].hist(normalized_all_edges, n_bins, density=True, histtype='bar', color='cadetblue', alpha=0.8)
    axs[0].axvline(x = cutoff4Proportion, color = 'k', ls = '--')
    axs[0].text(cutoff4Proportion+0.4, 0.9, 'p=%s'%cutoff4Proportion)
    axs[0].legend(frameon=False)
    kde = scipy.stats.gaussian_kde(normalized_all_edges)
    x = np.linspace(0, max(normalized_all_edges), 200)
    axs[0].plot(x, kde(x), color='gray', ls='--')
    axs[0].set_title('Distribution of $p_{rj}$')
    axs[0].set_xlabel('Proportion')
    
    axs[0].set_ylabel('Frequency')
    
    axs[1].hist(num_of_all_edges, n_bins, histtype='bar', color='salmon', alpha=0.8)
    axs[1].axvline(x = cutoff, color = 'k', ls = '--')
    axs[1].text(cutoff+10, chrom_len, 'cutoff=%d'%cutoff)
    axs[1].legend(frameon=False)
    axs[1].set_title('$e_{rj}$, time=%s'%time_str)
    #axs[1].sharey(axs[0])
    #axs[1].set_xlim([0, 400])
    axs[1].set_xlabel('Edges')
    
    for axs in fig.get_axes():
        axs.grid(False)
    fig.tight_layout()
            
    plt.savefig(out_path + 'cutoff_cohort_%s_proportion_%s.jpg'%(time_str,cutoff4Proportion), dpi=fig_dpi)
    plt.close(fig)


def main(in_data_folder,
         cohort,
         chromosome,
         resolution,
         out_data_folder,
         cutoff4Proportion,
         fig_dpi):
    
    # using the proportion to find the cutoff for given resolution
    bin_str = str(int(resolution/1000)) + 'kb'
    
    if chromosome == 'whole_genome':
        chrom_strs = ['chr'+str(i) for i in range(1,24)]
    else: 
        chrom_strs = chromosome
    
    chrom_len = len(chrom_strs)
    
    in_path = in_data_folder + '/hic_data/' + bin_str + '/hic_community_data/' + cohort
    num_of_all_edges = combineEdges(in_path, cohort, chrom_strs)
    num_of_all_edges.sort()
    normalized_all_edges = normalization(num_of_all_edges)
    
    cutoff = getCutoff(cutoff4Proportion, num_of_all_edges, normalized_all_edges)
    print('The cutoff for %dkb is %d under dataset %s' %(int(resolution/1000), cutoff, cohort), '\n')
    
    print('Start drawing the distribution of edges', '\n')
    
    out_path = out_data_folder + '/hic_data/estimate_community_size/'
    if not os.path.exists(out_path):
        os.makedirs(out_path)
        print('Create:', out_path, '\n')
    
    disPlot(normalized_all_edges, num_of_all_edges, cutoff4Proportion, cutoff, cohort, out_path, fig_dpi, chrom_len)    
    
    
    
    
    
    
    
    
    
    
    
    