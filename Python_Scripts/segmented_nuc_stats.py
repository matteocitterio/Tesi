import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import warnings
import os
import sys
import scipy.spatial as sp
import time
import statistics

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

colors=[]
columns=['9-norm','9-tum','13-tum','3-norm','3-tum','5-tum','7-tum','8-tum','12-norm','12-tum']
img_ns=[7,64,11,13,6,6,18,24,2,9]
img_ns_string=['9_INT_Normal','3_INT_Normal','3_INT_Tumor','5_INT_Tumor','7_INT_Tumor','8_INT_Tumor','9_INT_Tumor','12_INT_Normal','12_INT_Tumor','13_INT_Tumor']
indexis=np.arange(63)+1
indexis.tolist()
output_frame=pd.DataFrame(index=indexis,columns=columns)

for eta in range(0,len(img_ns)):

    for h in range(1,img_ns[eta]):

        image_N=str(h)
        df=pd.read_csv(r'D:\matteo_citterio\risultati_segment\ROI_'+img_ns_string[eta]+'\Results_'+image_N+'.csv')            #apri il csv
        num_nuclei=len(df.axes[0])

        if not os.path.exists(r'D:\matteo_citterio\nuke_stats\ROI_'+img_ns_string[eta]+'\img_'+image_N):                      #fai la cartella se non esiste
	
	        os.mkdir(r'D:\matteo_citterio\nuke_stats\ROI_'+img_ns_string[eta]+'\img_'+image_N)

        num_bins=int(num_nuclei/11)

        quantity_distr_list=['Area','Perim.','AR','Round','Solidity','Angle','Feret','Circ.']
        means_and_std=np.zeros((2,len(quantity_distr_list)))

        for i in range(0,len(quantity_distr_list)):

            quant=quantity_distr_list[i]
            print('folder: ',img_ns_string[eta],'img:',image_N,'/',img_ns[eta]-1,' quant:',quant)
            fig,ax=plt.subplots()                                                                                       #area
            bins=np.linspace(df[quant].min(),df[quant].max(),num_bins)
            ax.hist(df[quant],bins)
            ax.set_title(quant+'distr.')
            means_and_std[0,i]=df[quant].mean()
            means_and_std[1,i]=statistics.stdev(df[quant]) 
            fig.savefig(r'D:\matteo_citterio\nuke_stats\ROI_'+img_ns_string[eta]+'\img_'+image_N+'\\'+quant+'_hist.png')

            for j in range(0,len(quantity_distr_list)):

                print('plotting vs',j,'/',len(quantity_distr_list))

                if (i!=j):

                    fig,ax=plt.subplots()  
                    ax.plot(df[quant],df[quantity_distr_list[j]],'o', mfc='red')
                    ax.set_xlabel(quant)
                    ax.set_ylabel(quantity_distr_list[j])
                    ax.set_title(quant+'vs'+quantity_distr_list[j])
                    fig.savefig(r'D:\matteo_citterio\nuke_stats\ROI_'+img_ns_string[eta]+'\img_'+image_N+'\\'+quant+' vs '+quantity_distr_list[j]+'.png')


    # output_frame['9-norm'].iloc[i-1]=e                                                                          #l'indice delle foto parte da 1
    # E.append(e)
    
    # colors.append('darkviolet')
    # edges.append('red')
    # h+=1