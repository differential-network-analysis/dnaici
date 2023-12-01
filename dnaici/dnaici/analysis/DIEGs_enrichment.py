import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def importData(all_file, selected_chroms):
    
    gene_list_all = []
    
    for chrom in selected_chroms:
        
        filename = os.path.join(all_file, chrom, 'common_genes.tsv')
        common_genes_df = pd.read_csv(filename, sep = '\t')
        gene_list = list(common_genes_df.gene)
        gene_list_all += gene_list
        
    gene_set = set(gene_list_all)
    
    return gene_set


def foldChange(gene_expression_file):
    
    gene_expression = pd.read_csv(gene_expression_file, header=0, sep='\t')
    gene_expression['MCF7'] = gene_expression.MCF7_RNA_rep1+gene_expression.MCF7_RNA_rep2+10e-1
    gene_expression['MCF7TR'] = gene_expression.MCF7TR_RNA_rep1+gene_expression.MCF7TR_RNA_rep2+10e-1
    
    expression_df = pd.DataFrame(data={'foldChange':np.divide(gene_expression['MCF7TR'].values,gene_expression['MCF7'].values), 'gene':gene_expression.GeneSymbol})
    
    return expression_df


def relativeRatio(gene_expression_file):
    
    gene_expression = pd.read_csv(gene_expression_file, header=0, sep='\t')
    gene_expression['MCF7'] = gene_expression.MCF7_RNA_rep1+gene_expression.MCF7_RNA_rep2+10e-1
    gene_expression['MCF7TR'] = gene_expression.MCF7TR_RNA_rep1+gene_expression.MCF7TR_RNA_rep2+10e-1
    gene_expression['mean'] = (gene_expression['MCF7']+gene_expression['MCF7TR'])/2
    gene_expression['difference'] = gene_expression['MCF7TR']-gene_expression['MCF7']
    
    expression_df = pd.DataFrame(data={'relativeRatio':np.divide(gene_expression['difference'].values,gene_expression['mean'].values), 'gene':gene_expression.GeneSymbol})
    
    return expression_df


def selectDEGs(gene_expression_file, method):
    
    if method == 'foldChange':
        expression_df = foldChange(gene_expression_file)
        DEGs_df = expression_df.loc[(expression_df['foldChange'] >= 2 ) | (expression_df['foldChange'] <= 0.5)]
        
    elif method == 'relativeRatio':
        expression_df = relativeRatio(gene_expression_file)
        DEGs_df = expression_df.loc[(expression_df['relativeRatio'] >= 2/3 ) | (expression_df['relativeRatio'] <= -2/3)]     
    
    DEGs_df = DEGs_df.rename(columns={'GeneSymbol': 'gene'})
    
    return DEGs_df


def importDavid(file):
    
    david_results = pd.read_excel(file)
    david_results['-log10PValue'] = -np.log10(david_results['PValue'])
    david_results = david_results[['Category','Term','-log10PValue','Fold Enrichment']]
    
    david_results_GO_BP = david_results[david_results.Category.str.contains('GOTERM_BP')].head(16)
    david_results_GO_MF = david_results[david_results.Category.str.contains('GOTERM_MF')].head(15)
    david_results_KEGG = david_results[david_results.Category.str.contains('KEGG')].head(15)
    david_results_REACTOME = david_results[david_results.Category.str.contains('REACTOME')].head(16)
    david_results_REACTOME = david_results_REACTOME.loc[david_results_REACTOME['Term'] != 'TRP channels']
            
    return david_results_GO_BP, david_results_GO_MF, david_results_KEGG, david_results_REACTOME


def enrichmentPlot(out_folder, df, name, fig_dpi):
    
    df = df[::-1]
    
    plt.rcParams['font.family'] = 'Arial'
    fig, ax = plt.subplots()
    ax2 = ax.twiny()
    
    ax.barh(df['Term'], df['Fold Enrichment'], color='#E5C6DD', alpha=1)
    ax2.scatter(df['-log10PValue'], df['Term'], color='#53ADA3')
    
    ax.tick_params(axis='y', which='both', left=False, right=False, )
    
    ax.xaxis.set_label_position('top')
    ax.xaxis.tick_top()
    ax.set_xlabel('Fold Enrichment')
    
    ax2.xaxis.set_label_position('bottom')
    ax2.xaxis.tick_bottom()
    ax2.set_xlabel('-log10(P-value)')
    #ax.axvline(x = -np.log10(0.05), color = 'gray', ls = '--')
    out_david = out_folder + '/david_'+name+'.jpg'
    plt.savefig(out_david, bbox_inches='tight', dpi=fig_dpi)
    print('Export figure:', out_david, '\n')
    plt.close(fig)


def main(in_data_folder, 
         resolution, 
         out_data_folder,
         cohort1,
         cohort2,
         method,
         pval_cutoff,
         fig_dpi):
    
    bin_str = str(int(resolution/1000)) + 'kb'

    # use relative ratio/fold change to select DEGs among DIGs
    print('Use', method, 'method to select DEGs among DIGs', '\n')
    
    in_folder = in_data_folder + '/expression_data/'
    rna_count_file = in_folder + 'rna_count.txt'
    all_DEGs_df = selectDEGs(rna_count_file, method)
    print('Read file:', rna_count_file, '\n')
    
    out_folder = out_data_folder + '/differential_network_analysis/' + bin_str
    nodes_file = os.path.join(out_folder, '%s_vs_%s_pval_ls'%(cohort1,cohort2)+str(pval_cutoff)+'_selected_sigNodes.tsv')
    
    if not os.path.exists(nodes_file):
        print('No DIEGs found in the input chromosome :( Please move to another chromosome')
    else:
        selected_df = pd.read_csv(nodes_file, sep='\t')
        print('Read file:', nodes_file, '\n')
        
        selected_chroms = selected_df.id.str.split('_',expand=True)[0].unique()
        gene_set = importData(out_folder, selected_chroms)
        expression_bed = in_folder + 'mcf7_%s_geneExp.bed'%cohort1
        expression_df = pd.read_csv(expression_bed, sep='\t', header=0)
        print('Read file:', expression_bed, '\n')
        
        gene_df = expression_df[expression_df.new_name.isin(gene_set)]
        out_genes_chromosomes = gene_df[['chrom', 'new_name']].drop_duplicates()
        out_genes_chromosomes.rename(columns = {'new_name':'gene'}, inplace = True)
        out_genes_chromosomes_cp = out_genes_chromosomes.copy()
        out_genes_chromosomes_cp = out_genes_chromosomes_cp.set_index('chrom')
        out_genes_file = os.path.join(out_folder, str(len(out_genes_chromosomes_cp)) + '_DIGs_%s_resolution.tsv'%bin_str)
        out_genes_chromosomes_cp.to_csv(out_genes_file)
        print('Export file:', out_genes_file, '\n')
        
        #gene_list_select = pd.read_csv(out_genes_file, header=0, sep=',')
        select_DEGs_df = pd.merge(all_DEGs_df, out_genes_chromosomes, on='gene')
        select_DEGs_df = select_DEGs_df.set_index('gene')
        print('There are ' + str(len(select_DEGs_df)) + ' DEGs in the ' + str(len(out_genes_chromosomes)) + ' genes we selected', '\n')
        # export these genes to txt file for david analysis again
        out_DIEGs_file = os.path.join(out_folder, str(len(select_DEGs_df)) + '_DIEGs_%s_resolution.tsv'%bin_str)
        select_DEGs_df.to_csv(out_DIEGs_file)
        print('Export file:', out_DIEGs_file, '\n')
        
        # enrichment analysis based on david results for DEGs in the DIGs
        david_file = out_folder + '/david_for_'+str(len(select_DEGs_df))+'_genes.xlsx'
        if os.path.exists(david_file):
            david_results_GO_BP, david_results_GO_MF, david_results_KEGG, david_results_REACTOME = importDavid(david_file)
            enrichmentPlot(out_folder, david_results_GO_BP, 'GO_BP', fig_dpi)
            enrichmentPlot(out_folder, david_results_GO_MF, 'GO_MF', fig_dpi)
            enrichmentPlot(out_folder, david_results_KEGG, 'KEGG', fig_dpi)
            enrichmentPlot(out_folder, david_results_REACTOME, 'REACTOME', fig_dpi)
        else:
            print('No DAVID results found, please go to DAVID website for enrichment analysis :)', '\n')
    
    
    
    
    
    
    












