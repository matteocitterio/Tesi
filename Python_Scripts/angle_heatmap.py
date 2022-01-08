import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import warnings
import os
import sys
import scipy.spatial as sp
from scipy.spatial import distance
import skimage.io
import skimage.color

from pandas.core.common import SettingWithCopyWarning
from shapely.geometry import MultiPoint, Point, Polygon, LineString,MultiPolygon
from shapely.ops import polygonize,unary_union
from scipy.spatial import Voronoi, voronoi_plot_2d,ConvexHull, convex_hull_plot_2d
from PIL import Image
from matplotlib.path import Path
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl
from functions import *


img_folder='ROI_13_INT_Tumor'
img_folder_2='ROI_5_INT_Tumor'
img_N1=1
img_N2=11

multi=10/3
multi2=10/3

multi_df=pd.read_csv('multipliers.csv')
pixel_sizes_df=pd.read_csv('pixel_sizes.csv')
pixel=multi

img=skimage.io.imread('D:\matteo_citterio\SCAN - Copia\\'+img_folder+'\\'+str(img_N1)+'.tif')
df=pd.read_csv('D:\matteo_citterio\\risultati_segment\\'+img_folder+'\Results_'+str(img_N1)+'.csv')
gray_image = skimage.color.rgb2gray(img)

img2=skimage.io.imread('D:\matteo_citterio\SCAN - Copia\\'+img_folder_2+'\\'+str(img_N2)+'.tif')
df2=pd.read_csv('D:\matteo_citterio\\risultati_segment\\'+img_folder_2+'\Results_'+str(img_N2)+'.csv')
gray_image2 = skimage.color.rgb2gray(img2)

fig,ax=plt.subplots(1,2)
ax[0].imshow(gray_image,cmap='gray')
ax[0].set_title('Patient 13',fontsize=13)
ax[1].set_title('Patient 5',fontsize=13)
s=ax[0].scatter(df['XM']*multi,df['YM']*multi,c=df['Angle'],vmin=0,vmax=180,cmap='twilight',s=22)
ax[1].imshow(gray_image2,cmap='gray')
s1=ax[1].scatter(df2['XM']*multi2,df2['YM']*multi2,c=df2['Angle']*multi2,vmin=0,vmax=180,cmap='twilight',s=22)
cb=fig.colorbar(s)
cb.set_label('$\\theta$',fontsize=13)
ax[0].set_xlabel('[pixel]',fontsize=13)
ax[0].set_ylabel('[pixel]',fontsize=13)
ax[1].set_xlabel('[pixel]',fontsize=13)

fig,ax=plt.subplots(1,2)
vmax=max(max(df['AR']),max(df2['AR']))
ax[0].imshow(gray_image,cmap='gray')
ax[0].set_title('Patient 13')
ax[1].set_title('Patient 5')
s=ax[0].scatter(df['XM']*multi,df['YM']*multi,c=df['AR'],vmin=1,vmax=vmax,cmap='twilight',s=10)
ax[1].imshow(gray_image2,cmap='gray')
s1=ax[1].scatter(df2['XM']*multi2,df2['YM']*multi2,c=df2['AR'],vmin=1,vmax=vmax,cmap='twilight',s=10)
cb=fig.colorbar(s)
cb.set_label('$AR$')
ax[0].set_xlabel('[pixel]')
ax[0].set_ylabel('[pixel]')
ax[1].set_xlabel('[pixel]')
ax[1].set_ylabel('[pixel]')

plt.suptitle('AR heatmap')

plt.show()