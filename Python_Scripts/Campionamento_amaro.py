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
from matplotlib.path import Path
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

image_N='35'
fig,ax=plt.subplots()
img=Image.open(r'D:\matteo_citterio\SCAN - Copia\ROI_3_INT_Normal\\'+image_N+'.tif')
a,b=img.size
diagonal=((a**2)+(b**2))**0.5


df=pd.read_csv(r'D:\matteo_citterio\risultati_segment\ROI_3_INT_Normal\Results_'+image_N+'.csv')

num_nuclei=len(df.axes[0])
data_size=num_nuclei                     #dimensione dei dati generati
print(a,b)
num_bins=int(num_nuclei/13)

# alpha=np.random.random(data_size)*a        #uniforme reali
# beta=np.random.random(data_size)*b
# gamma=np.random.random(data_size)+1

# punti=np.zeros((data_size,3))
# h=0

# print('Creating points')
# for i in range(0,len(alpha)):

# 	punti[i][0]=alpha[i]
# 	punti[i][1]=beta[i]
# 	h+=1
# df['simula']=alpha
# df['simulb']=beta
funct=create_nuclei_function(num_nuclei,df)			    #crea un array con le pos dei cdm dei nuclei dalla df di stardist

df.to_csv(r'D:\matteo_citterio\risultati_segment\ROI_3_INT_Normal\Results.txt', header=None, index=None, sep=' ',columns=['XM','YM']  ,mode='a')
# df.to_csv(r'D:\matteo_citterio\risultati_segment\ROI_3_INT_Normal\simul.txt', header=None, index=None, sep=' ',columns=['simula','simulb']  ,mode='a')
# array=get_distances_and_quantities_vectors(funct,num_nuclei,'exp')	#serve per trovare il minimo

density_matlab()

dj=pd.read_csv(r'D:\matteo_citterio\Python Scripts\temp2.csv')
print(dj)
os.remove(r'D:\matteo_citterio\risultati_segment\ROI_3_INT_Normal\Results.txt')
os.remove(r'D:\matteo_citterio\Python Scripts\temp2.csv')
# density_check()
# dh=pd.read_csv(r'D:\matteo_citterio\Python Scripts\temp2.csv')
# print(dh)
# os.remove(r'D:\matteo_citterio\risultati_segment\ROI_3_INT_Normal\simul.txt')
# os.remove(r'D:\matteo_citterio\Python Scripts\temp2.csv')
# bins_density,density_profile=density_profile_function_amara(funct,num_nuclei,a,b,num_bins_density)
dj.columns =['g(R)', 'R','c']
bins_density=dj['R']
density_profile=dj['g(R)']
# dh.columns =['g(R)', 'R','counts']
# bins_density2=dh['R']
# density_profile2=dh['g(R)']
# counts2=dh['counts']
df1=pd.DataFrame(density_profile)
k1=20
df1=df1.rolling(k1+3).mean()

size=len(dj.axes[0])
x=np.zeros((size))
y=np.zeros((size))

for j in range(0,size):
    x[j]=bins_density[j]
    y[j]=df1.iloc[j]


shift=np.count_nonzero(np.isnan(y))
print(shift)
y = y[~np.isnan(y)]

maxi=np.argsort(y)
print(x[maxi[-1]+shift])


ax.axvline(x[maxi[-1]+shift],c='red')
ax.plot(bins_density,df1,'-',c='green',label='scan3_1')
fig.suptitle('density function')
ax.set_ylabel('g(r)')
ax.set_xlim(-1,size)
ax.set_xlabel('r [pixel]')	
ax.legend()
ax.axhline(1,c='k')
plt.tight_layout()
plt.show()
