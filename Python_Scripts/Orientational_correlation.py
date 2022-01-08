from numpy import ma
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import warnings
import os
import sys
from pandas.core.indexes import multi
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

colors=[]
columns=['9-norm','9-tum','13-tum','3-norm','3-tum','5-tum','7-tum','8-tum','12-norm','12-tum']
# img_ns=[7,64,11,13,6,6,18,24,2,9]
real_img_n=[64,11,13,6,6,7,18,24,2,9]

img_ns_string=['ROI_3_INT_Normal', 'ROI_3_INT_Tumor', 'ROI_5_INT_Tumor', 'ROI_7_INT_Tumor', 'ROI_8_INT_Tumor', 'ROI_9_INT_Normal', 'ROI_9_INT_Tumor', 'ROI_12_INT_Normal', 'ROI_12_INT_Tumor', 'ROI_13_INT_Tumor']
# indexis=np.arange(63)+1
# indexis.tolist()
# output_frame_means=pd.DataFrame(index=indexis,columns=np.arange(6)+1)
# output_frame_std=pd.DataFrame(index=indexis,columns=columns)

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

# ax.legend(handles=legend_elements)

columns=['3-norm','3-tum','5-tum','7-tum','8-tum','9-norm','9-tum','12-norm','12-tum','13-tum']
indexis=np.arange(63)+1
indexis.tolist()
output_frame=pd.DataFrame(index=indexis,columns=columns)
output_frame_distances=pd.DataFrame(index=indexis,columns=columns)

multi_df=pd.read_csv('multipliers.csv')
pixel_sizes_df=pd.read_csv('pixel_sizes.csv')


for eta in range(0,len(img_ns_string)):

    for i in range(0,real_img_n[eta]-1):

        multis=multi_df[img_ns_string[eta]].iloc[i]
        pixel=pixel_sizes_df[img_ns_string[eta]].iloc[i]
        coeff=multis*pixel
        print('folder:',img_ns_string[eta],'img',i,'/',real_img_n[eta])
        df=pd.read_csv(r'D:\\matteo_citterio\\risultati_segment\\'+img_ns_string[eta]+'\\Results_'+str(i+1)+'.csv')   
        dj=pd.read_csv('D:\\matteo_citterio\\plot_al_fly\\nematics\\'+img_ns_string[eta]+'\\'+str(i+1)+'_corr.csv')
        max_len=dj['max'].iloc[0]*1.5
        # e=compute_energy_Voronoi_NN(df)
        # print(e)
        e=compute_energy_Voronoi_NN_with_metrics(df,max_len,coeff)
        print(e)
        output_frame[columns[eta]].iloc[i]=e
        output_frame_distances[columns[eta]].iloc[i]=dj['max'].iloc[0]

output_frame.to_csv('Orientational_correlation.csv')
output_frame_distances.to_csv('max_gr.csv')

