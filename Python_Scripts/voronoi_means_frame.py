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
outliers=[3000,1000,1250,600,1000,1500,1000,3500,1000,2600]


fig2,an=plt.subplots()
ar_ratio=[]
gdl=[]
ar_means=[]

colors=['green','black','black','green','black','black','black','black','green','black']
columns=['9-norm','9-tum','13-tum','3-norm','3-tum','5-tum','7-tum','8-tum','12-norm','12-tum']

indexis=np.arange(63)+1
indexis.tolist()
output_frame_means=pd.DataFrame(index=indexis,columns=columns)
output_frame_std=pd.DataFrame(index=indexis,columns=columns)

output_frame_Area_means=pd.DataFrame(index=indexis,columns=columns)
output_frame_Area_std=pd.DataFrame(index=indexis,columns=columns)

for eta in range(0,len(img_ns_string)): 

    print(img_ns_string[eta])

    df=pd.read_csv(r'D:\matteo_citterio\Python Scripts\Area_ratio_per_image\Area_ratio_voronoiROI_'+img_ns_string[eta]+'.csv')
    dj=pd.read_csv(r'D:\matteo_citterio\Python Scripts\AR_per_image\AR_voronoiROI_'+img_ns_string[eta]+'.csv')

    num_nuclei=len(df.axes[0])


    for i in range(1,real_img_n[eta]):                  #ciclo a colonna che vuol dire a immagine

        output_frame_means[columns[eta]].iloc[i-1]=dj[str(i)].mean()
        output_frame_std[columns[eta]].iloc[i-1]=dj[str(i)].std()

        output_frame_Area_means[columns[eta]].iloc[i-1]=df[str(i)].mean()
        output_frame_Area_std[columns[eta]].iloc[i-1]=df[str(i)].std()

output_frame_means.to_csv("AR_voronoi_mean.csv")
output_frame_std.to_csv("AR_voronoi_std.csv")
output_frame_Area_means.to_csv("Area_ratio_voronoi_mean.csv")
output_frame_Area_std.to_csv("Area_ratio_voronoi_std.csv")
