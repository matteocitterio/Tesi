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
import math

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning) 
from pandas.core.common import SettingWithCopyWarning
from shapely.geometry import MultiPoint, Point, Polygon, LineString,MultiPolygon
from shapely.ops import polygonize,unary_union
from scipy.spatial import Voronoi, voronoi_plot_2d,ConvexHull, convex_hull_plot_2d
from scipy.fft import fft, ifft
import scipy.stats as stats
from PIL import Image
from matplotlib.path import Path
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from functions import *

colors=[]
columns=['9-norm','9-tum','13-tum','3-norm','3-tum','5-tum','7-tum','8-tum','12-norm','12-tum']
img_ns=[7,64,11,13,6,6,18,24,2,9]
real_img_n=[7, 18, 9, 64, 11, 13, 6, 6, 24, 2]
img_ns_string=['9_INT_Normal','9_INT_Tumor','13_INT_Tumor','3_INT_Normal','3_INT_Tumor','5_INT_Tumor','7_INT_Tumor','8_INT_Tumor','12_INT_Normal','12_INT_Tumor']
quant='AR'
divisors=[60,200,170,900,140,100,44,26,135,28]

fig2,an=plt.subplots()
ar_ratio=[]
gdl=[]
ar_means=[]

colors=['green','black','black','green','black','black','black','black','green','black']

correlations=[]

print('AR')

for eta in range(0,len(img_ns_string)): 

    df=pd.read_csv(r'D:\matteo_citterio\Python Scripts\AR_per_image\AR_voronoiROI_'+img_ns_string[eta]+'.csv')
    dj=pd.read_csv(r'D:\matteo_citterio\Python Scripts\AR_per_image\AR_'+img_ns_string[eta]+'.csv')

    num_nuclei_voro=len(df.axes[0])
    num_nuclei_nuke=len(dj.axes[0])

    voro_data=[]
    nuclear_data=[]

    for i in range(1,real_img_n[eta]):                  #ciclo a colonna che vuol dire a immagine

        for j in range(0,num_nuclei_voro):

            if(not math.isnan(df[str(i)].iloc[j])):

                voro_data.append(df[str(i)].iloc[j])

        for j in range(0,num_nuclei_nuke):

            if(not math.isnan(dj[str(i)].iloc[j])):

                nuclear_data.append(dj[str(i)].iloc[j])

    tempo=[]
    for i in range(0,len(voro_data)):

        if(voro_data[i]<5):
            tempo.append(voro_data[i])

    voro_data=tempo

    tempo=[]
    for i in range(0,len(voro_data)):
        tempo.append(nuclear_data[i])

    nuclear_data=tempo

    correlations.append(stats.spearmanr(voro_data,nuclear_data))
    print(img_ns_string[eta],':',stats.spearmanr(voro_data,nuclear_data))

print('AREA')

for eta in range(0,len(img_ns_string)): 

    df=pd.read_csv(r'D:\matteo_citterio\Python Scripts\Area_ratio_per_image\Area_ratio_voronoiROI_'+img_ns_string[eta]+'.csv')
    dj=pd.read_csv(r'D:\matteo_citterio\Python Scripts\Area_ratio_per_image\Area_ratio_'+img_ns_string[eta]+'.csv')

    num_nuclei_voro=len(df.axes[0])
    num_nuclei_nuke=len(dj.axes[0])

    voro_data=[]
    nuclear_data=[]

    for i in range(1,real_img_n[eta]):                  #ciclo a colonna che vuol dire a immagine

        for j in range(0,num_nuclei_voro):

            if(not math.isnan(df[str(i)].iloc[j])):

                voro_data.append(df[str(i)].iloc[j])

        for j in range(0,num_nuclei_nuke):

            if(not math.isnan(dj[str(i)].iloc[j])):

                nuclear_data.append(dj[str(i)].iloc[j])

    tempo=[]
    for i in range(0,len(voro_data)):

        if(voro_data[i]<5):
            tempo.append(voro_data[i])

    voro_data=tempo

    tempo=[]
    for i in range(0,len(voro_data)):
        tempo.append(nuclear_data[i])

    nuclear_data=tempo

    correlations.append(stats.spearmanr(voro_data,nuclear_data))
    print(img_ns_string[eta],':',stats.spearmanr(voro_data,nuclear_data))