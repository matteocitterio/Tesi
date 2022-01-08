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

colors=[]
columns=['9-norm','9-tum','13-tum','3-norm','3-tum','5-tum','7-tum','8-tum','12-norm','12-tum']
img_ns=[7,64,11,13,6,6,18,24,2,9]
real_img_n=[64,11,13,6,6,7,18,24,2,9]
#real_img_n=[18,24,2,9]
#real_img_n=[7, 18, 9, 64, 11, 13, 6, 6, 24, 2]
img_ns_string=['ROI_3_INT_Normal', 'ROI_3_INT_Tumor', 'ROI_5_INT_Tumor', 'ROI_7_INT_Tumor', 'ROI_8_INT_Tumor', 'ROI_9_INT_Normal', 'ROI_9_INT_Tumor', 'ROI_12_INT_Normal', 'ROI_12_INT_Tumor', 'ROI_13_INT_Tumor']
# indexis=np.arange(63)+1
# indexis.tolist()
# output_frame_means=pd.DataFrame(index=indexis,columns=np.arange(6)+1)
# output_frame_std=pd.DataFrame(index=indexis,columns=columns)
#img_ns_string=['ROI_9_INT_Tumor', 'ROI_12_INT_Normal', 'ROI_12_INT_Tumor', 'ROI_13_INT_Tumor']

quant='AR'
disegno=False
for eta in range(0,len(img_ns_string)): 

    nuc=[]

    for i in range(1,real_img_n[eta]):

        df=pd.read_csv(r'D:\matteo_citterio\risultati_segment\\'+img_ns_string[eta]+'\Results_'+str(i)+'.csv')
        nuc.append(len(df.axes[0]))

    inx_max=max(nuc)
    indexis=np.arange(inx_max)+1
    indexis.tolist()
    columns=np.arange(real_img_n[eta]-1)+1
    columns.tolist()
    output_frame_ar=pd.DataFrame(index=indexis,columns=columns)
    output_frame_area=pd.DataFrame(index=indexis,columns=columns)

    for h in range(1,real_img_n[eta]):                                                                                              #ciclo su tutte le immagini di un folder

        print('folder:', img_ns_string[eta], 'img:', h,'/',real_img_n[eta]-1)
        image_N=str(h)
        ars,area,bho,bho1=voronoi_and_quantities_per_folder(img_ns_string[eta],h,h,disegno)

        for i in range(0,len(ars)):

            output_frame_ar[h].iloc[i]=ars[i]

        for i in range(0,len(area)):
            output_frame_area[h].iloc[i]=area[i]

    print(output_frame_ar)
    output_frame_ar.to_csv("AR_per_image\AR_voronoi"+img_ns_string[eta]+".csv")
    output_frame_area.to_csv("Area_ratio_per_image\Area_ratio_voronoi"+img_ns_string[eta]+".csv")

