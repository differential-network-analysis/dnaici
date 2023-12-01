import os
import pandas as pd


def main(in_data_folder,
         resolution
         ):
    #input selected sorted chromosome size file
    chrom = 'hg19' 
    f1 = in_data_folder + '/' + chrom + '/' + chrom + '.chrom.sizes.clear.sorted'
    #use bedtools to make bin bed file based on selected window size
    window = int(resolution/1000)
    rmBlackList = True
    
    out_bed_file = in_data_folder + '/' + chrom + '/' + chrom + '_XY.' + str(window) + 'kb.windows.bed'
    cmd = 'bedtools makewindows -g ' + f1 +  \
          ' -w ' + str(resolution) + ' > ' + out_bed_file
    print(cmd)
    os.system(cmd)
    
    df1 = pd.read_csv(out_bed_file, sep='\t', header=None)
    
    df1[1] = df1[1] + 1
    df1.to_csv(out_bed_file, index=False, sep='\t', header=None)
    
    uq_chrm = df1.iloc[:, 0].unique()
    #for each chrom region assign index number and numeric chromsome number
    record_chroms = []
    for ui in uq_chrm:
        tmp_df = df1[df1.iloc[:,0] == ui].copy()
        #here controls index starts from 0 or 1 !!
        tmp_index = [i for i in range(1, tmp_df.shape[0]+1)]
        tmp_chrm = ui.replace('chr', '')
        if tmp_chrm =='X':
            tmp_chrm = 23
        elif tmp_chrm == 'Y':
            tmp_chrm = 24
        elif tmp_chrm == 'M':
            tmp_chrm = 25
        else:
            tmp_chrm = int(tmp_chrm)
        
        tmp_chrm = [tmp_chrm]*tmp_df.shape[0]   
        tmp_df.insert(3, 3, tmp_index)
        tmp_df.insert(4, 4, tmp_chrm)
        if ui not in ['chrY', 'chrM']:
            record_chroms.append(tmp_df.copy())
    
    all_df = pd.concat(record_chroms).copy()
    out_bed_file2 = out_bed_file.replace('.bed', '_bin.bed')
    all_df.to_csv(out_bed_file2, sep='\t', index=False, header=None)
    
    if rmBlackList:
        out_bed_file3 = out_bed_file2.replace('.bed', '.BlackListFiltered.bed')
        #remove black list from the bins
        cmd2 = 'bedtools intersect -v -a ' + out_bed_file2 + \
               ' -b '+ in_data_folder + '/' + chrom +'/'+ chrom+'-blacklist.bed '+ \
               ' > ' + out_bed_file3
    print(cmd2)
    os.system(cmd2)







