from PIL import Image
import os
import matplotlib.lines as lines
import numpy as np
import matplotlib as mlp
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from itertools import product
from sklearn.cluster import KMeans
mlp.use('agg')

#combine images
def combine_multiple_images(files,num_times_width,num_times_height,out_fig):
  #initialize new image size
  img=Image.open(files[-1])
  width=img.width*num_times_width
  height=img.height*num_times_height
  final_img=Image.new('RGB',(width,height),'#FFFFFF')
  x_offset=0
  y_offset=0
  #loop from files to combine images where we assume all file names are sorted by chromosomes
  loop=0
  for fi in files:
    img=Image.open(fi)
    #change y_offset when it reaches end of width
    if loop>=num_times_width:
       y_offset += img.size[1]
       x_offset =0
       loop=0
    final_img.paste(img,(x_offset,y_offset))
    x_offset += img.size[0]
    loop +=1
  #export new image
  #resized_img=final_img.resize((16,8))
  final_img.save(out_fig,dpi=(300,300))
  return final_img,out_fig


def draw_plot(matrix, cell, start, end, color_start, color_end, bar_start, bar_end,out_path, fig_dpi=100):
    plot_data = np.triu(matrix[start:end, start:end])
    #map min and max vales to 5 levels for both negative and positive values
    cvals = [color_start, color_start*0.8,color_start*0.6,color_start*0.4, color_start*0.2, 0 , 0.2 * color_end , 0.4*color_end,0.6 * color_end, 0.8 * color_end, color_end]
    ##generate color map for Green-> White >Red
    colors=[(0,1,0), (1,1,1), (1,0,0)]
    cmap= mlp.colors.LinearSegmentedColormap.from_list("GWR",colors,N=len(cvals))
    fig, ax = plt.subplots()
    line = lines.Line2D([bar_start - start, bar_end - start], [bar_start - start + 5, bar_end - start + 5],
                        lw=4, color='black', axes=ax)
    im=ax.imshow(plot_data, cmap=cmap,vmin=color_start,vmax=color_end)
    ax.set_title(cell.replace('_', ' '))
    ax.add_line(line)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    plt.axis('off')
    fig.colorbar(im)
    plt.draw()
    out_file=os.path.join(out_path,cell+'.jpg')
    plt.savefig(out_file,dpi=fig_dpi)
    plt.close(fig)
    return out_file


def draw_plot2(matrix, cell, start, end, color_start, color_end, bar_start, bar_end,out_path, colors='GWR', plotXtick=[],plotLines=[], plotYtick=[],fig_dpi=100):
    '''matrix is a numpy 2-D matrix
       cell is title and figure name
       start, end are the start and end position for the input matrix
       color_start, color_end are the minimum and maximum values for the color code
       bar_star, bar_end are the bar position for the plot
       colors are selected color map either GWR (green, white, red) or WR (white , red) 
       plotXtick is a list of xtick labels that will be added on plot
       plotLines is a list of line postions that will be added on plot
    '''
    if (matrix.shape[0]==matrix.shape[1]) & (len(plotLines)>0) :
       plot_data = np.triu(matrix[start:end, start:end])
    else:
       plot_data=matrix.copy()
   
    #map min and max vales to 5 levels for both negative and positive values
    cvals1 = [color_start, color_start*0.8,color_start*0.6,color_start*0.4, color_start*0.2, 0 , 0.2 * color_end , 0.4*color_end,0.6 * color_end, 0.8 * color_end, color_end]
    ##generate color map for Green-> White >Red
    colors1=[(0,1,0), (1,1,1),(1,0,0)]
    cmap1=mlp.colors.LinearSegmentedColormap.from_list("GWR",colors1,N=len(cvals1))
   
    #generate color map for Blue -> white-> Yellow
    colors3=[(0, 0, 0.7), (1,1,1), (1,0.5,0)]
    cmap3=mlp.colors.LinearSegmentedColormap.from_list("BWY",colors3,N=len(cvals1))

    #assume color start=0 for white->Red color map 
    cvals2=  [ 0,  0.2 * color_end , 0.4*color_end,0.6 * color_end, 0.8 * color_end, color_end]
    ## coloe map for White -> Red
    colors2=[(1,1,1),(1,0,0)]
    cmap2= mlp.colors.LinearSegmentedColormap.from_list("WR",colors2, N=len(cvals2))

    #select color map
    if colors=='GWR':
       out_cmap=cmap1
    elif colors=='WR':
       out_cmap=cmap2
    elif colors=='BWY':
       out_cmap=cmap3

    fig, ax = plt.subplots()
    #fig, ax=plt.subplots(1, figsize=(6,5))
    if len(plotLines)==0:
       #one line plot
       line = lines.Line2D([bar_start - start, bar_end - start], [bar_start - start + 5, bar_end - start + 5],
                        lw=4, color='black', axes=ax)
    else:
       #plot multi lines
       plt.hlines(plotLines,xmin=ax.get_xticks()[0]-1,xmax=matrix.shape[1],lw=0.5,color='black')

    if plot_data.shape[0]==plot_data.shape[1]:
      #a diagonal matrix
      im=ax.imshow(plot_data, cmap=out_cmap,vmin=color_start,vmax=color_end)
    else:
      #non diagonal 2D matrix
      im=ax.imshow(plot_data, cmap=out_cmap,vmin=color_start,vmax=color_end,interpolation='nearest', aspect='auto')

    ax.set_title(cell.replace('_', ' ').replace(' in clusters edges ','').replace(' noWeight','').replace(' edges','')) 
    if len(plotLines)==0 :
       ax.add_line(line)
    #ax.xaxis.set_visible(False)
    #ax.yaxis.set_visible(False)

    #add xtick label if it is needed
    if len(plotXtick)>0:
       ticks=plotXtick
       tmp_xticks_position=ax.get_xticks()
       if len(tmp_xticks_position)>len(ticks):
           xticks_position=tmp_xticks_position[1:len(ticks)+1]
       elif len(tmp_xticks_position)<len(ticks):
           xticks_position=[i for i in range(0,9)]
       else:
           xticks_position=tmp_xticks_position
       ax.set_xticks(xticks_position)
       #print(ticks)
       #print(xticks_position)
       ax.set_xticklabels(ticks, rotation=20, fontsize=8)
       
       left_plotLines=np.array([0]+plotLines)
       right_plotLines=np.array(plotLines+[plot_data.shape[0]])
       min_delta_plotLines=int( (right_plotLines-left_plotLines).min()/2)
       ax.set_yticks( np.concatenate( ( np.array(plotLines)-min_delta_plotLines , np.array([matrix.shape[0]])-min_delta_plotLines) )  )
       ax.set_yticklabels( [i for i in range(1,len(plotLines)+2)],fontsize=8)
    else:
       plt.axis('off')
       ax.xaxis.set_visible(False)
 
    if len(plotYtick)>0:
       ticks=plotYtick
       tmp_yticks_position=ax.get_yticks()
       if len(tmp_yticks_position)>len(ticks):
          ytick_position=tmp_yticks_position[1:len(ticks)+1]
       elif len(tmp_yticks_position)<len(ticks): 
          yticks_position=[i for i in range(0,len(ticks))]
       else:
          yticks_position=tmp_yticks_position
       #print(ticks)
       #print(yticks_position)
       ax.set_yticks(yticks_position)
       ax.set_yticklabels(ticks,rotation=0, fontsize=8)
    else:
       ax.yaxis.set_visible(False)

    #show color bar and plot figure
    fig.colorbar(im)
    plt.draw()
    out_file=os.path.join(out_path,cell+'.jpg')
    print('Output: \n', out_file, '\n')
    plt.savefig(out_file,dpi=fig_dpi)
    plt.close(fig)
    return out_file


def plot_elbow(chrom_str, num_of_clusters,new_region_pair_df,out_path):
  #elbow curve for kmeans clustering
  print('Estimate the number of clusters ....')
  plt.clf()
  Nc = range(1, num_of_clusters)
  kmeans = [KMeans(n_clusters=i) for i in Nc]
  score = [kmeans[i].fit(new_region_pair_df).score(new_region_pair_df) for i in range(len(kmeans))]
  plt.plot(Nc,score)
  plt.vlines([5,10],plt.gca().get_yticks().min(), plt.gca().get_yticks().max(),colors='red')
  plt.xlabel('Number of Clusters')
  plt.ylabel('Score')
  plt.title('Elbow Curve')
  out_fig=os.path.join(out_path, chrom_str+'_Elbow_curve_kmeans.jpg')
  print(out_fig)
  plt.savefig(out_fig)
  return out_fig


def plot_clustered_heatmap(chrom_str,feature_str,num_of_clusters, sub_df, color_clip_value,color_type,out_path,edge_feature_str='',fig_dpi=100):
  '''
    Kmeans-clustring of an imput dataframe and show it in a color coded heatmap 
  '''
  #num_of_clusters=3
  kmeans = KMeans(n_clusters=num_of_clusters).fit(sub_df)

  sub_df['kmeans']=kmeans.labels_
  sorted_new_region_pair_df=sub_df.sort_values(by='kmeans').copy()
  new_data_df=sorted_new_region_pair_df.copy()

  #find line position for each cluster
  new_data_df['line_location']=[i for i in range(1,new_data_df.shape[0]+1)]
  line_location=[]
  for i in range(0,num_of_clusters-1):
     line_location.append(new_data_df[new_data_df.kmeans==i].line_location.max())

  new_data_df.drop(['kmeans','line_location'],axis=1,inplace=True)

  #show heatmap of clustered rows
  matrix=new_data_df.to_numpy()
  clip_value=color_clip_value
  matrix[matrix>clip_value]=clip_value
  matrix[matrix<-clip_value]=-clip_value
  matrix_start=0
  matrix_end=matrix.shape[0]
  color_start=-clip_value
  color_end=clip_value
  cell=chrom_str+'_C'+str(num_of_clusters)+ edge_feature_str
  start=0
  end=0
  bar_start=0
  bar_end=0
  #out_path='out_figs_'+ feature_str
  xticks=new_data_df.columns.str.replace('_'+feature_str+'_zscore','').str.replace(chrom_str.split('_')[0]+'_','').to_list()
  out_fig=draw_plot2(matrix, cell, start, end, color_start, color_end, bar_start, bar_end,out_path, colors='GWR',plotXtick=xticks,plotLines=line_location,fig_dpi=fig_dpi)
  return out_fig

#this function for plot purpose
def find_cluster_label_in_index(idx, record_clustered_sub_df):
    '''find cluster label for each pair-wise interaction based on network clustering of subgraphs, if not find then return -1
       idx, is the pair of nodes such as 1:13
       record_clustered_sub_df is a dictionary of network clusters where pair of nodes already groupped to their corresponding clusters.
    '''
    is_find=[]
    for ki in record_clustered_sub_df.keys():
        if idx in record_clustered_sub_df[ki].index:
           is_find.append(ki)
        else :
           is_find.append(-1)
    out_find=np.array(is_find)
    out_find2=out_find[out_find>=0].copy()
    if len( out_find2) >0:
        return out_find2[0]+1
    else:
        return -1


def plot_network_clusteringLabel_heatmap(chrom_str,num_of_clusters, a_matrix, new_region_pair_df , record_clustered_sub_df,out_path,edge_feature_str='',fig_dpi=100):
  ''' input a numpy matrix represents 2-D interaction heatmap a_matrix
      Then find network clustering labels for each pair-wise interactions (new_region_pair_df)
      based on not subgraphs (record_clustered_sub_df)
  '''
  #assume the first one is Hi-C matrix
  #then build a cluster matrix based on network clustering labels
  #a_matrix=all_matrix[0].copy()
  cluster_matrix=np.zeros(a_matrix.shape)
  new_region_pair_df['location']=new_region_pair_df.index.to_list()
  new_region_pair_df[['row_i','col_j']]=new_region_pair_df.location.str.split(':',expand=True).astype(int)
  #for each pair of interactions, find its relevant network clustering label
  for idx,row in new_region_pair_df.iterrows():
    if a_matrix[row.row_i-1,row.col_j-1] !=0:
       #do not show filtered inereactions in heatmap
       cluster_matrix[row.row_i-1,row.col_j-1]=find_cluster_label_in_index(idx,record_clustered_sub_df)

  # define a color map for all cluster labels, where -1 means missing cluster label, 0 means no interactions ,the rest numbers are the cluster label from 1 to 16 
  color_map1 = {-1: np.array([0,0,0]), #black 
             0: np.array([255,255,255]), #white
             1: np.array([0, 255, 0]), # green
             3: np.array([255, 0, 0]), # red
             2: np.array([255, 255,0]), #yellow
             4: np.array([153, 0, 153]), #purpe
             5: np.array([51, 51,255]), #dark blue
             6: np.array([255, 204, 204]), #pink
             7: np.array([0, 0, 255]), # blue 
             8: np.array([0, 255,255]), #light blue
             9: np.array([153, 51, 255]),
             10: np.array([160,160,160]), #gray
             11: np.array([153,255,204]), #light green
             12: np.array([153,0,76]),
             13: np.array([0 ,128, 255]),
             14: np.array([229,255,204]),
             15: np.array([153,153,0]),
             16: np.array([51, 0, 0 ]), # black	
             17: np.array([255, 0, 255]), # magenta
             18: np.array([192, 192, 192]), # silver
             19: np.array([240, 230, 140]), # khaki
             20: np.array([128, 0, 0]), # maroon
             21: np.array([128, 128, 0]), # olive
             22: np.array([165, 42, 42]), # brown
             23: np.array([240, 128, 128]), # light coral
             24: np.array([143, 188, 143]), # dark sea green
             25: np.array([221, 160, 221]), # plum
             26: np.array([210, 105, 30]), # chocolate
             27: np.array([188, 143, 143]), # rosy brown
             28: np.array([255, 255, 240]), # ivory
             29: np.array([95, 158, 160]), # cadet blue
             30: np.array([216, 191, 216]), # thistle
             31: np.array([255, 250, 250]) # snow	     
            }

  color_map2 = {-1: np.array([0,0,0]), #black 
             0: np.array([255,255,255]), #white
             1: np.array([0, 255, 0]), # green
             2: np.array([255, 0,0]) #red
            }
  if num_of_clusters>2:
     color_map=color_map1
  else:
     color_map=color_map2

  #normalize color map for figure legend purpose
  norm_color_map={}
  for ki in color_map.keys():
     norm_color_map[ki]= list(color_map[ki]/255.0)+[1]

  # make a 3d numpy array that has a color channel dimension  
  data=cluster_matrix.copy()
  data_3d = np.ndarray(shape=(data.shape[0], data.shape[1], 3), dtype=int)
  for i in range(0, data.shape[0]):
    for j in range(0, data.shape[1]):
        data_3d[i][j] = color_map[data[i][j]]

  # display the plot 
  fig, ax = plt.subplots(1,1)
  im=ax.imshow(data_3d)
  plt.axis('off')

  # add cluster numbers to the plot 
  # thanks to tmdavison answer here https://stackoverflow.com/a/40890587/7871710
  if False:
    for i in range(0, data.shape[0]):
      for j in range(0, data.shape[1]):
        c = data[j,i]
        ax.text(i, j, str(c), va='center', ha='center')

  #add figure legend
  values = np.array(np.unique(data.ravel()),int)
  # get the colors of the values, according to the 
  # colormap used by imshow
  # create a patch (proxy artist) for every color 
  patches = [ mpatches.Patch(color=norm_color_map[values[i]], label="C {l}".format(l=values[i]) ) for i in range(len(values)) ]
  # put those patched as legend-handles into the legend
  plt.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0. )
  out_fig=os.path.join(out_path, chrom_str+'_'+str(num_of_clusters)+'clusters'+edge_feature_str + '_heatmap.jpg')
  plt.savefig(out_fig,dpi=fig_dpi)
  plt.close(fig)
  print('Output: \n', out_fig, '\n')
  return out_fig, cluster_matrix.copy()






