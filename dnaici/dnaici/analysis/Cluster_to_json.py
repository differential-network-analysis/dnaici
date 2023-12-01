import pandas as pd
import numpy as np
import networkx as nx
import matplotlib as mlp
import json
import os


def nodes2cluster_label_df(in_path, cluster_f):
    '''
    Input:
        in_path, input file path
        cluster_f, file name of clustering results from network clustering ,  where
        cluster label of all nodes are shown in one column
    Return:
        cluster_df, a dataframe with two columns nodes, labels
        num_of_nodes, number of nodes in the cluster
    '''
    #read clustered clustering label for each node
    cluster_c=np.loadtxt(os.path.join(in_path,cluster_f) )
    num_of_nodes=cluster_c.shape[0]
    
    #make a dataframe for nodes and cluster labels
    nodes=[i for i in range(0,num_of_nodes)]
    cluster_df=pd.DataFrame()
    cluster_df['nodes']=nodes
    cluster_df['labels']=cluster_c
    
    return cluster_df.copy(), num_of_nodes


def color_code4clusters(cluster_df, max_clusters=2000):
    '''
    make a color code for cluster labels,
    Input: 
         cluster_df, is a dataframe with two columns, nodes, labels
         max_clusters, is the maximum number of clusters allowed in color map
    return:
         colors, a list of colors for all nodes
    '''
    #make color for clusters
    color_code={0: 'red', 
                1: 'green',
                2: 'Lime',
                3: 'blue',
                4: 'aliceblue',
                5: 'aquamarine',
                6: 'antiqueWhite',
                7:  'Aqua',
                8:  'Coral',
                9:  'Brown',
                10:  'BlueViolet'}
    #fill the rest of cluster as black color
    for i in range(11,max_clusters):
        color_code[i]='black'
    
    #return colors for all nodes
    colors=[]
    for idx,row in cluster_df.iterrows():
        colors.append(color_code[row.labels.astype(int)])
        
    return colors


def adj_matrix4network_view(num_of_nodes, in_df, out_path, withWeight, colors, nodes_size, edges_color, chrom_str, edges_rm_str, percent_str, cluster_df):
    ''' Make json file and node cluster label file for Cytoscape view.
    Input:
        num_of_nodes, number of nodes in a graph
        in_df, a dataframe contains all pair-wise edges in a graph, where with two columns, left node, right node of an edge
        out_path, a file path contains all network clustering results such as out_network_clusters
        withWeight, True or False, to export weight or no-weight edges for graph view
        colors,  a list color code for cluster labels which is exported by color_code4clusters
        nodes_size, a list or single value for node size
        edges_color, a list or single value for edge color
        chrom_str, string of a chromosome such as chr1
        edge_rm_str, string of edege filtering condition such as  _edges_rmGlob_ or _edges_rm_
        percent_str, string of percentage of lowest interactions will be removed such as 0.1 
    Return:
        G, a graph from networkx basad on adj-matrix, A
        A, an adj-matrix of graph G
        out_network_file, output file name for networks in json format which can be viewed in Cytoscape
        out_node_f, output file of clustering labels of all nodes which can be used in Cytoscape for coloring nodes in different clusters.
    '''
    #make adjancy matrix for the edges of original interactoins , then convert to Graph
    adj_matrix=np.zeros((num_of_nodes,num_of_nodes))
    if withWeight:
        for idx,row in in_df.iterrows():
            adj_matrix[row[0].astype(int),row[1].astype(int)]=row[2]
    else:
        for idx,row in in_df.iterrows():
            adj_matrix[row[0].astype(int),row[1].astype(int)]=1
    #build a network by the a-matrix
    A=adj_matrix.copy()
    G = nx.from_numpy_array(A)
    #graph options
    options = {'node_color': colors,
               'node_size': nodes_size,
               'width': 0.1 ,
               'with_labels': False,
               'font_weight': 'bold',
               'edge_color': edges_color}
    
    nx.draw_spring(G, **options)
    #export networks and nodes colors for cytoscape
    #export clustered node labels
    if withWeight:
        out_node_f=chrom_str + edges_rm_str + percent_str + 'percent_4communities_class.tsv'
    else:
        out_node_f=chrom_str +edges_rm_str + percent_str + 'percent_noWeight_4communities_class.tsv'
    
    out_node_f = os.path.join(out_path, out_node_f)
    print('Output 1 (%s):\n' %chrom_str, out_node_f, '\n')
    cluster_df.to_csv(out_node_f, sep='\t', index=False)
    #export graph data to cytoscape format
    t=nx.cytoscape_data(G) 
    #change int to str for nodes names
    for ti in t['elements']['edges']:
        ti['data']['source']=str(ti['data']['source'])
        ti['data']['target']=str( ti['data']['target'])
    #export to cytoscape json format
    if withWeight:
        out_network_file = chrom_str + edges_rm_str + percent_str + 'percent_network.json'
    else:
        out_network_file = chrom_str + edges_rm_str + percent_str + 'percent_noWeight_network.json'
    
    out_network_file = os.path.join(out_path, out_network_file)
    print('Output 2 (%s):\n' %chrom_str, out_network_file, '\n')
    with open(out_network_file, 'w', encoding='utf-8') as f:
        json.dump(t, f, ensure_ascii=False, indent=4)
    
    return G.copy(), A.copy(), out_network_file, out_node_f


def main(in_data_folder,
         cohort,
         chromosome,
         resolution,
         out_data_folder):
        
    isFind_global_cutoff = False
    percent_str = '0'
    #whether exported edges with weight
    withWeight = False
    
    if isFind_global_cutoff:
        edges_rm_str = '_edges_rmGlob_'
    else:
        edges_rm_str = '_edges_rm_'

    if chromosome == 'whole_genome':
        chrom_strs = ['chr'+str(i) for i in range(1,24)]
    else: 
        chrom_strs = chromosome
    
    bin_str = str(int(resolution/1000)) + 'kb'
    
    in_path = in_data_folder + '/hic_data/' + bin_str + '/hic_community_data/' + cohort
    out_path = out_data_folder + '/hic_data/' + bin_str + '/hic_community_data/' + cohort
    
    for chrom_str in chrom_strs:
        #load network cluststering results for the nodes
        if withWeight:
            cluster_f = chrom_str + edges_rm_str + percent_str + 'percent_4communities.txt'
        else:
            cluster_f = chrom_str + edges_rm_str + percent_str + 'percent_noWeight_4communities.txt'
        f = os.path.join(in_path, cluster_f)
        print('Read 1 (%s):\n' %chrom_str, f, '\n')
        #assignn nodes and its cluster label to a dataframe
        cluster_df, num_of_nodes = nodes2cluster_label_df(in_path, cluster_f)
        
        #generate color code for all cluster labels
        colors = color_code4clusters(cluster_df)
        
        #load hi-c interaction matrix
        if withWeight:
            f = os.path.join(in_path, chrom_str + edges_rm_str + percent_str + 'percent_4zscore.matrix')
        else:
            f = os.path.join(in_path, chrom_str + edges_rm_str + percent_str + 'percent_noWeight_4zscore.matrix')
        print('Read: 2 (%s):\n' %chrom_str, f, '\n')
        #read all edges in a graph which was used to do network clustering and exported by function export_edges4matrix
        in_df = pd.read_csv(f, sep='\t', header=None)
        
        #set up graph nodes and edges property
        nodes_size = 1000
        edges_color = 'blue'
        #export graph and the corresponding nodes cluster label to json files which can be viewed by Cytoscape later
        G, adj_matrix, out_network_file, out_node_file = adj_matrix4network_view(num_of_nodes, in_df, out_path, withWeight, 
                                    colors, nodes_size, edges_color, chrom_str, edges_rm_str, percent_str, cluster_df)

