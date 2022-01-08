import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import warnings
import os
import sys
import scipy.spatial as sp
import time

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning) 
from pandas.core.common import SettingWithCopyWarning
from shapely.geometry import MultiPoint, Point, Polygon, LineString,MultiPolygon
from shapely.ops import polygonize,unary_union
from scipy.spatial import Voronoi, voronoi_plot_2d,ConvexHull, convex_hull_plot_2d
from scipy.fft import fft, ifft
from PIL import Image
from matplotlib.path import Path
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from functions import *

n_images_9normal=7                 #in realtà il numero di immagini viste è n_images-1
n_images_3normal=64
n_images_3tumor=11
n_images_5tumor=13
n_images_7tumor=6
n_images_8tumor=6
n_images_9tumor=18
n_images_12normal=24
n_images_12tumor=2
n_images_13tumor=9
n_images_tot=160

E=[]
colors=[]
edges=[]
xs=[]
fig,ax=plt.subplots()
legend_elements =[Line2D([0], [0], marker='o', color='w', label='9-norm',markerfacecolor='darkviolet',markeredgecolor='red',markersize=5),
                    Line2D([0], [0], marker='o', color='w', label='9-tum',markerfacecolor='slategrey',markeredgecolor='blue', markersize=5),
                    Line2D([0], [0], marker='o', color='w', label='13-tum',markerfacecolor='slategrey',markeredgecolor='purple', markersize=5),
                    Line2D([0], [0], marker='o', color='w', label='3-norm',markerfacecolor='darkviolet',markeredgecolor='black', markersize=5),
                    Line2D([0], [0], marker='o', color='w', label='3-tum',markerfacecolor='slategrey',markeredgecolor='orange', markersize=5),
                    Line2D([0], [0], marker='o', color='w', label='5-tum',markerfacecolor='slategrey',markeredgecolor='green', markersize=5),
                    Line2D([0], [0], marker='o', color='w', label='7-tum',markerfacecolor='slategrey',markeredgecolor='royalblue', markersize=5),
                    Line2D([0], [0], marker='o', color='w', label='8-tum',markerfacecolor='slategrey',markeredgecolor='olive', markersize=5),
                    Line2D([0], [0], marker='o', color='w', label='12-norm',markerfacecolor='darkviolet',markeredgecolor='tomato', markersize=5),
                    Line2D([0], [0], marker='o', color='w', label='12-tum',markerfacecolor='slategrey',markeredgecolor='yellow', markersize=5)]

h=0
ax.legend(handles=legend_elements)

columns=['9-norm','9-tum','13-tum','3-norm','3-tum','5-tum','7-tum','8-tum','12-norm','12-tum']
indexis=np.arange(63)+1
indexis.tolist()
output_frame=pd.DataFrame(index=indexis,columns=columns)

#calcolo l'<<energia>> per tutte le immagini di ogni cartella

for i in range(1,n_images_9normal):

    xs.append(h)
    image_N=str(i)
    df=pd.read_csv(r'D:\matteo_citterio\risultati_segment\ROI_9_INT_Normal\Results_'+image_N+'.csv')            #apri il csv
    num_nuclei=len(df.axes[0])
    funct=create_nuclei_function(num_nuclei,df)	
    b=get_distances_and_quantities_vectors(funct,num_nuclei,'img '+image_N,simple=True)                                   #calcola matrice delle distanze

    e=compute_energy(b,num_nuclei,df,image_N)/num_nuclei
    output_frame['9-norm'].iloc[i-1]=e                                                                          #l'indice delle foto parte da 1
    E.append(e)
    
    colors.append('darkviolet')
    edges.append('red')
    h+=1

for i in range(1,n_images_9tumor):

    xs.append(h)
    image_N=str(i)
    df=pd.read_csv(r'D:\matteo_citterio\risultati_segment\ROI_9_INT_Tumor\Results_'+image_N+'.csv')            #apri il csv
    num_nuclei=len(df.axes[0])
    funct=create_nuclei_function(num_nuclei,df)	
    b=get_distances_and_quantities_vectors(funct,num_nuclei,'img '+image_N,simple=True)                                   #calcola matrice delle distanze
   
    e=compute_energy(b,num_nuclei,df,image_N)/num_nuclei
    output_frame['9-tum'].iloc[i-1]=e
    E.append(e)
    colors.append('slategrey')
    edges.append('blue')
    h+=1

for i in range(1,n_images_13tumor):

    xs.append(h)
    image_N=str(i)
    df=pd.read_csv(r'D:\matteo_citterio\risultati_segment\ROI_13_INT_Tumor\Results_'+image_N+'.csv')            #apri il csv
    num_nuclei=len(df.axes[0])
    funct=create_nuclei_function(num_nuclei,df)	
    b=get_distances_and_quantities_vectors(funct,num_nuclei,'img '+image_N,simple=True)                                   #calcola matrice delle distanze
    
    e=compute_energy(b,num_nuclei,df,image_N)/num_nuclei
    output_frame['13-tum'].iloc[i-1]=e                                                                          #l'indice delle foto parte da 1
    E.append(e)
    
    colors.append('slategrey')
    edges.append('purple')

    h+=1

for i in range(1,n_images_3normal):

    xs.append(h)
    image_N=str(i)
    df=pd.read_csv(r'D:\matteo_citterio\risultati_segment\ROI_3_INT_Normal\Results_'+image_N+'.csv')            #apri il csv
    num_nuclei=len(df.axes[0])
    funct=create_nuclei_function(num_nuclei,df)	
    b=get_distances_and_quantities_vectors(funct,num_nuclei,'img '+image_N,simple=True)                                   #calcola matrice delle distanze
  
    e=compute_energy(b,num_nuclei,df,image_N)/num_nuclei
    output_frame['3-norm'].iloc[i-1]=e
    E.append(e)
    colors.append('darkviolet')
    edges.append('black')

    h+=1

for i in range(1,n_images_3tumor):

    xs.append(h)
    image_N=str(i)
    df=pd.read_csv(r'D:\matteo_citterio\risultati_segment\ROI_3_INT_Tumor\Results_'+image_N+'.csv')            #apri il csv
    num_nuclei=len(df.axes[0])
    funct=create_nuclei_function(num_nuclei,df)	
    b=get_distances_and_quantities_vectors(funct,num_nuclei,'img '+image_N,simple=True)                                   #calcola matrice delle distanze

    e=compute_energy(b,num_nuclei,df,image_N)/num_nuclei
    output_frame['3-tum'].iloc[i-1]=e
    E.append(e)
    colors.append('slategrey')
    edges.append('orange')
    h+=1

for i in range(1,n_images_5tumor):

    xs.append(h)
    image_N=str(i)
    df=pd.read_csv(r'D:\matteo_citterio\risultati_segment\ROI_5_INT_Tumor\Results_'+image_N+'.csv')            #apri il csv
    num_nuclei=len(df.axes[0])
    funct=create_nuclei_function(num_nuclei,df)	
    b=get_distances_and_quantities_vectors(funct,num_nuclei,'img '+image_N,simple=True)                                   #calcola matrice delle distanze

    e=compute_energy(b,num_nuclei,df,image_N)/num_nuclei
    output_frame['5-tum'].iloc[i-1]=e
    E.append(e)
    colors.append('slategrey')
    edges.append('green')
    h+=1

for i in range(1,n_images_7tumor):

    xs.append(h)
    image_N=str(i)
    df=pd.read_csv(r'D:\matteo_citterio\risultati_segment\ROI_7_INT_Tumor\Results_'+image_N+'.csv')            #apri il csv
    num_nuclei=len(df.axes[0])
    funct=create_nuclei_function(num_nuclei,df)	
    b=get_distances_and_quantities_vectors(funct,num_nuclei,'img '+image_N,simple=True)                                   #calcola matrice delle distanze

    e=compute_energy(b,num_nuclei,df,image_N)/num_nuclei
    output_frame['7-tum'].iloc[i-1]=e
    E.append(e)
    colors.append('slategrey')
    edges.append('royalblue')
    h+=1

for i in range(1,n_images_8tumor):

    xs.append(h)
    image_N=str(i)
    df=pd.read_csv(r'D:\matteo_citterio\risultati_segment\ROI_8_INT_Tumor\Results_'+image_N+'.csv')            #apri il csv
    num_nuclei=len(df.axes[0])
    funct=create_nuclei_function(num_nuclei,df)	
    b=get_distances_and_quantities_vectors(funct,num_nuclei,'img '+image_N,simple=True)                                   #calcola matrice delle distanze

    e=compute_energy(b,num_nuclei,df,image_N)/num_nuclei
    output_frame['8-tum'].iloc[i-1]=e
    E.append(e)
    colors.append('slategrey')
    edges.append('olive')
    h+=1

for i in range(1,n_images_12tumor):

    xs.append(h)
    image_N=str(i)
    df=pd.read_csv(r'D:\matteo_citterio\risultati_segment\ROI_12_INT_Tumor\Results_'+image_N+'.csv')            #apri il csv
    num_nuclei=len(df.axes[0])
    funct=create_nuclei_function(num_nuclei,df)	
    b=get_distances_and_quantities_vectors(funct,num_nuclei,'img '+image_N,simple=True)                                   #calcola matrice delle distanze

    e=compute_energy(b,num_nuclei,df,image_N)/num_nuclei
    output_frame['12-tum'].iloc[i-1]=e
    E.append(e)
    colors.append('slategrey')
    edges.append('yellow')
    h+=1

for i in range(1,n_images_12normal):

    xs.append(h)
    image_N=str(i)
    df=pd.read_csv(r'D:\matteo_citterio\risultati_segment\ROI_12_INT_Normal\Results_'+image_N+'.csv')            #apri il csv
    num_nuclei=len(df.axes[0])
    funct=create_nuclei_function(num_nuclei,df)	
    b=get_distances_and_quantities_vectors(funct,num_nuclei,'img '+image_N,simple=True)                                   #calcola matrice delle distanze
    #j=get_nearest_neighbors(b,89,num_nuclei)
    #a=89-1

    #circle1 = plt.Circle((df['XM'].iloc[a],df['YM'].iloc[a]),get_distance(df['XM'].iloc[a],df['YM'].iloc[a],df['XM'].iloc[int(j[3]-1)],df['YM'].iloc[int(j[3]-1)]) , color='b',fill=False)
    #ax.add_patch(circle1)
    #ax.scatter(df['XM'].iloc[a],df['YM'].iloc[a],s=30,c='purple')
    #ax.scatter(df['XM'],df['YM'],s=1,c='k')
    #ax.scatter(df['XM'].iloc[j[:]-1],df['YM'].iloc[j[:]-1],s=30,c='r')

    e=compute_energy(b,num_nuclei,df,image_N)/num_nuclei
    output_frame['12-norm'].iloc[i-1]=e
    E.append(e)
    colors.append('darkviolet')
    edges.append('tomato')
    h+=1

len_E=len(E)
E=np.asarray(E)
indici=np.argsort(E)
print(output_frame)
output_frame.to_csv("nematic_energy_6neigh_anglenoFE.csv")

for i in range(0,len_E):

    color=colors[indici[i]]
    edge=edges[indici[i]]
    ax.scatter(i,E[indici[i]],c=color,s=15,edgecolors=edge)

plt.show()