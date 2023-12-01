import os
import subprocess


def submit_job_in_subprocess(record_cmds):
    processes = []
    for cmd in record_cmds:
        processes.append(subprocess.Popen(cmd, shell=True))
    for p in processes:
        p.communicate()


def main(in_data_folder,
         cohort,
         chromosome,
         resolution,
         super_resolution,
         p_value,
         zscore,
         out_data_folder
         ):
    
    bin_str = str(int(resolution/1000)) + 'kb'
    in_data_path = in_data_folder + '/hic_data/hicup_processed/' + cohort
    out_data_path = out_data_folder + '/hic_data/' + bin_str + '/hic_interaction_homer/' + cohort
    
    if not os.path.exists(out_data_path):
        os.makedirs(out_data_path)
        print('Create: ', out_data_path)
    
    if chromosome == 'whole_genome':
        chrom_strs = ['chr'+str(i) for i in range(1,24)]
    else: 
        chrom_strs = chromosome
    
    for i in range(0, len(chrom_strs)):
        if chrom_strs[i] == 'chr23':
            chrom_strs[i] = 'chrX'
    
    record_cmds = []
    all_out_files = []
    
    for chrom_str in chrom_strs:
        out_file0 = 'mcf7_' + cohort + '_significantInteractions_norm' + bin_str + '_' + chrom_str + '.txt'
        out_file = os.path.join(out_data_path, out_file0)
        cmd = 'analyzeHiC ' + in_data_path + ' -res ' + str(resolution) + ' -superRes ' + str(super_resolution) + ' -norm -nolog -center -chr ' + chrom_str + ' -interactions ' + out_file + ' -nomatrix -pvalue ' + str(p_value) + ' -zscore ' + str(zscore)
        all_out_files.append(out_file)
        
        print(cmd)
        #os.system(cmd)
        record_cmds.append(cmd)
    
    submit_job_in_subprocess(record_cmds)
    
    out_lk_file = [i for i in all_out_files if 'chrX' in i]
    if len(out_lk_file) > 0: 
        print(out_lk_file)
        source_file = os.path.basename(out_lk_file[0])
        target_file = source_file.replace('chrX', 'chr23')
        cmd1 = 'ln -s ' + source_file + ' ' + target_file
        print(cmd1)
        cmd2 = 'mv -f ' + target_file + ' ' + out_data_path
        print(cmd2)
        os.system(cmd1)
        os.system(cmd2)





