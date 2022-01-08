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
from matplotlib.patches import Ellipse
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import skimage.io
from skimage.measure import EllipseModel
from functions import *

df=pd.read_csv('D:\matteo_citterio\\risultati_segment\ROI_3_INT_Normal\Results_21.csv')
posizioni=pd.read_csv('D:\matteo_citterio\\risultati_segment\ROI_3_INT_Normal\contour20.csv')
multi_df=pd.read_csv('multipliers.csv')
pixel_sizes_df=pd.read_csv('pixel_sizes.csv')
multi=multi_df['ROI_3_INT_Normal'].iloc[21-1]
coeff_sistemazione=(multi)#*(pixel_sizes_df['ROI_13_INT_Tumor'].iloc[8-1])
fig,ax=plt.subplots()
num_nuclei=len(df.axes[0])

img=skimage.io.imread('D:\matteo_citterio\SCAN - Copia\ROI_3_INT_Normal\\21.tif')
img_star=Image.open(r'D:\matteo_citterio\SCAN - Copia\\ROI_3_INT_Normal\\21.tif')
a,b=img_star.size
diagonal=((a**2)+(b**2))**0.5
print(a,b)
ax.imshow(img)
rois=[]
rois_val=[]
nuclei=create_nuclei_array(num_nuclei,df)[0]*coeff_sistemazione			#crea un array con le pos dei cdm dei nuclei dalla df di stardist


labels=rois_detections(nuclei,0,df,posizioni,rois,rois_val)	#modifica rois e ritorna la struttura labels

vor=[]
rois_restricted=[]
membranes_restricted=[]

voronoi_on_single_roi(labels,vor,rois_restricted,membranes_restricted,rois)	#riempie le prime tre liste in input
volumes=voronoi_vertices_inside_membrane(vor,membranes_restricted,ax,False)	

df['XM']=df['XM']*multi
df['YM']=df['YM']*multi

temp_list_of_things=[]

for l in range(0,len(volumes)):
		
	for eta in range(0,len(volumes[l])):

		temp_list_of_things.append(volumes[l][eta])

ars=[]
nuke_ars=[]

for i in range(0,len(volumes[1])):          #le macchie

    for j in range(0,len(volumes[1][i])):   

        numero=len(volumes[1][i][j])
        points=np.zeros((numero,2))
        x,y=np.zeros((numero)),np.zeros((numero))
        for l in range(0,numero):

            x[l]=volumes[1][i][j][l]
            y[l]=volumes[2][i][j][l]
            points[l,0]=volumes[1][i][j][l]
            points[l,1]=volumes[2][i][j][l]

        ell = EllipseModel()
        ell.estimate(points)

        print(j,'/',len(volumes[1][i]))

        if(ell.params is not None):


            volumes[3][i][j]=np.around(volumes[3][i][j],decimals=2)
            volumes[4][i][j]=np.around(volumes[4][i][j],decimals=2)
            df['XM']=np.around(df['XM'],decimals=2)
            df['YM']=np.around(df['YM'],decimals=2)
            nuke_ars.append(df.loc[(df['XM']==volumes[3][i][j]) & (df['YM']==volumes[4][i][j])]['AR'])

            xc, yc, alpha, beta, theta = ell.params
            ars.append(max(alpha,beta)/min(alpha,beta))



temp1,temp2=[],[]
for i in range(0,len(ars)):
    if(ars[i]<4):
        temp1.append(ars[i])
        temp2.append(nuke_ars[i])

fig,ax=plt.subplots()
ax.plot(temp1,temp2,'o')
print(stats.spearmanr(temp1,temp2))


# print(len(volumes[0]))
# temp=[]
# temp_vert=np.zeros(())
# for i in range(0,len(volumes[0])):

#     for j in range(0,len(volumes[0][i])):

#         temp.append(volumes[0][i][j])

# fig,ax=plt.subplots()
# print(len(temp))
# ax.hist(temp,bins=int(len(temp)/10))


plt.show()