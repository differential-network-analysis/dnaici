import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import shutil
from dnaici.preprocess import Preprocess_hic_tag2homer


def countDistence(file, resolution, time):
    
    f = open(file, 'r')
    lines = f.readlines()
    
    ind = lines[0].split('	').index('Distance')
    
    distance = []
    for x in lines:
        distance.append(x.split('	')[ind])
    f.close()
    
    #distance = [np.log(float(i)) for i in distance[1:]]
    distance = [float(i) for i in distance[1:]]
    
    info = resolution.split(',')
    bin_size = info[0]
    super_resolution = info[1]
    
    bin_size_kb = str(int(int(bin_size)/1000))
    super_resolution_kb = str(int(int(super_resolution)/1000))
    
    Distance = pd.DataFrame({'Distance': distance, 
                             'bin size(kb), super resolution(kb)': [bin_size_kb+','+super_resolution_kb]*len(distance), 
                             'time': [time]*len(distance)})
    return Distance


def combineChromosome(resolution, cohort, in_data_folder, out_data_folder, chrom_strs):
    
    info = resolution.split(',')
    bin_size = info[0]
    super_resolution = info[1]
    
    path = out_data_folder + '/hic_data/estimate_resolution/hic_interaction_homer_%skb_%skb'%(str(int(int(bin_size)/1000)), str(int(int(super_resolution)/1000)))
    if not os.path.exists(path):
        os.makedirs(path)
        print('Create:', path, '\n')
    if not os.path.exists(path + '/' + cohort):
        Preprocess_hic_tag2homer.main(in_data_folder, cohort, chrom_strs, int(bin_size),
                                      int(super_resolution), 0.1, 1.0, out_data_folder+'/temp')
    path_temp = out_data_folder+'/temp/hic_data/'+str(int(int(bin_size)/1000))+'kb/hic_interaction_homer/'+cohort
    if os.path.exists(path_temp) and not os.path.exists(path+'/'+cohort):
        shutil.move(path_temp, path)
        shutil.rmtree(out_data_folder+'/temp')
        
    comChr = pd.DataFrame()
    
    for chrom in chrom_strs:
              
        file = path + '/' + cohort + '/mcf7_' + cohort + '_significantInteractions_norm%skb_%s.txt' %(str(int(int(bin_size)/1000)), chrom)
        Distance = countDistence(file, resolution, cohort)
        
        comChr = pd.concat([comChr, Distance], ignore_index=True)
        
    return comChr


def combineResolutions(resolutions, cohort, in_data_folder, out_data_folder, chrom_strs):
    
    comRes = pd.DataFrame()
    
    for resolution in resolutions:
        
        comChr = combineChromosome(resolution, cohort, in_data_folder, out_data_folder, chrom_strs)
        comRes = pd.concat([comRes, comChr], ignore_index=True)
        
    return comRes


def combineTime(resolutions, cohort_list, in_data_folder, out_data_folder, chrom_strs):
    
    comTime = pd.DataFrame()
    
    for cohort in cohort_list:
        
        comRes = combineResolutions(resolutions, cohort, in_data_folder, out_data_folder, chrom_strs)
        comTime = pd.concat([comTime, comRes], ignore_index=True)
                
    return comTime


def countInteraction(file):
    
    f = open(file, 'r')
    Interaction = len(f.readlines()) - 1
    
    return Interaction


def combineInteraction(resolutions, cohort, out_data_folder, chrom_strs):
    
    comInteraction = []
    
    for resolution in resolutions:
        
        meanInteraction = []
        for chrom in chrom_strs:
        
            info = resolution.split(',')
            bin_size = info[0]
            super_resolution = info[1]
            
            path = out_data_folder + '/hic_data/estimate_resolution/hic_interaction_homer_%skb_%skb'%(str(int(int(bin_size)/1000)), str(int(int(super_resolution)/1000)))
            file = path + '/' + cohort + '/mcf7_' + cohort + '_significantInteractions_norm%skb_%s.txt' %(str(int(int(bin_size)/1000)), chrom)
            Interaction = countInteraction(file)
            
            meanInteraction.append(Interaction)
            
        comInteraction.append(np.mean(meanInteraction))
    
    return comInteraction


def violinPlot(comTime, comInteraction0, comInteraction1, resolutions, out_data_folder, cohort1, cohort2, fig_dpi):
    
    sns.set_theme(style="whitegrid")
    my_colors = {cohort1: '#A4D3EE', cohort2: '#F6F6BC'}
    fig, ax = plt.subplots()
    sns.boxplot(data=comTime, x="bin size(kb), super resolution(kb)", y="Distance", hue="time", palette=my_colors)
    ax.grid(False)
    
    ax2 = plt.twinx()
    sns.lineplot(data=comInteraction0, linestyle='dashed', color='#A4D3EE', linewidth=3, ax=ax2)
    sns.lineplot(data=comInteraction1, linestyle='dashed', color='#F6F6BC', linewidth=3, ax=ax2)
    ax2.grid(False)
    ax2.set_ylabel('Mean number of significant interactions')
    #ax2.set(ylim=(-min(comInteraction1)/10, max(comInteraction0)+2000))
    
    bin_size1 = resolutions[0].split(',')[0]
    bin_size2 = resolutions[-1].split(',')[0]
    
    out_path = out_data_folder + '/hic_data/estimate_resolution/'
    if not os.path.exists(out_path):
        os.makedirs(out_path)
    
    if bin_size1 == bin_size2:
        s = 'same_resolution_%skb'%str(int(int(bin_size1)/1000))
    else:
        s = 'different_resolution'
    
    fig.tight_layout()
    plt.savefig(out_path + 'chr_all_' + s + '.jpg', dpi=fig_dpi)
    plt.close(fig)


def main(in_data_folder, 
         out_data_folder,
         chromosome,
         cal_type,
         cohort1, 
         cohort2,
         fig_dpi
         ):
    
    cohort_list = [cohort1, cohort2]
    if chromosome == 'whole_genome':
        chrom_strs = ['chr'+str(i) for i in range(1,24)]
    else: 
        chrom_strs = chromosome

    if cal_type == 0:
    # 0 represents comparison between different resolution
    
        print('Let us first check how precisely we can observe chromosomes')
        print('The super resolution is set to be the same as the bin size as recommended by the manual of HOMER :)')
            
        resolutions = ['500000,500000', '100000,100000', '50000,50000', '20000,20000', '10000,10000']
        
        comTime = combineTime(resolutions, cohort_list, in_data_folder, out_data_folder, chrom_strs)
        comInteraction0 = combineInteraction(resolutions, cohort1, out_data_folder, chrom_strs)
        comInteraction1 = combineInteraction(resolutions, cohort2, out_data_folder, chrom_strs)
        violinPlot(comTime, comInteraction0, comInteraction1, resolutions, out_data_folder, cohort1, cohort2, fig_dpi)
    
    elif cal_type == 1:
    # 1 represents comparison between different super with resolution equal to 50kb

        print('Change super resolution might help')
        print('Forgot to say, the bin size is fixed at 50kb :)')
        
        resolutions = ['50000,50000', '50000,100000', '50000,200000', '50000,300000', '50000,500000']
        
        comTime = combineTime(resolutions, cohort_list, in_data_folder, out_data_folder, chrom_strs)
        comInteraction0 = combineInteraction(resolutions, cohort1, out_data_folder, chrom_strs)
        comInteraction1 = combineInteraction(resolutions, cohort2, out_data_folder, chrom_strs)
        violinPlot(comTime, comInteraction0, comInteraction1, resolutions, out_data_folder, cohort1, cohort2, fig_dpi)

    elif cal_type == 2:
    # 2 represents comparison between different super with resolution equal to 100kb
    
        print('Change super resolution might help')
        print('Forgot to say, the bin size is fixed at 1000kb :)')
        
        resolutions = ['100000,100000', '100000,200000', '100000,300000', '100000,400000', '100000,500000']
        
        comTime = combineTime(resolutions, cohort_list, in_data_folder, out_data_folder, chrom_strs)
        comInteraction0 = combineInteraction(resolutions, cohort1, out_data_folder, chrom_strs)
        comInteraction1 = combineInteraction(resolutions, cohort2, out_data_folder, chrom_strs)
        violinPlot(comTime, comInteraction0, comInteraction1, resolutions, out_data_folder, cohort1, cohort2, fig_dpi)
    
    elif cal_type == 3:
    # 3 represents comparison between different super with resolution equal to 500kb
    
        print('Change super resolution might help')
        print('Forgot to say, the bin size is fixed at 5000kb :)')
    
        resolutions = ['500000,100000', '500000,200000', '500000,300000', '500000,500000', '500000,1000000']
        
        comTime = combineTime(resolutions, cohort_list, in_data_folder, out_data_folder, chrom_strs)
        comInteraction0 = combineInteraction(resolutions, cohort1, out_data_folder, chrom_strs)
        comInteraction1 = combineInteraction(resolutions, cohort2, out_data_folder, chrom_strs)
        violinPlot(comTime, comInteraction0, comInteraction1, resolutions, out_data_folder, cohort1, cohort2, fig_dpi)
    
    else:
        print('Wrong cla_type, only 0, 1, 2, and 3 are accepted :)')
    





