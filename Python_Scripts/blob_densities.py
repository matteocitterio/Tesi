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

from pandas.core.common import SettingWithCopyWarning
from shapely.geometry import MultiPoint, Point, Polygon, LineString,MultiPolygon
from shapely.ops import polygonize,unary_union
from scipy.spatial import Voronoi, voronoi_plot_2d,ConvexHull, convex_hull_plot_2d
from PIL import Image
from matplotlib.path import Path
from functions import *

#per 12 norm: 10/3, per 13_tum:10/3 o 2

df=pd.read_csv('D:\matteo_citterio\\risultati_segment\ROI_13_INT_Tumor\Results_8.csv')
posizioni=pd.read_csv('D:\matteo_citterio\\risultati_segment\ROI_13_INT_Tumor\contour7.csv')
num_nuclei=len(df.axes[0])

img=skimage.io.imread('D:\matteo_citterio\SCAN - Copia\ROI_13_INT_Tumor\\8.tif')
img_star=Image.open(r'D:\matteo_citterio\SCAN - Copia\ROI_13_INT_Tumor\\8.tif')
a,b=img_star.size
diagonal=((a**2)+(b**2))**0.5
print(a,b)

rois=[]
rois_val=[]
nuclei=create_nuclei_array(num_nuclei,df)[0]			#crea un array con le pos dei cdm dei nuclei dalla df di stardist



random_size_number=4000                                #numero dei nuclei fake da generare

random_nuclei=np.zeros((random_size_number,2))

random_nuclei2=np.zeros((10000,2))

random_size_number2=4000

for i in range(0,random_size_number):                   #genero punti random

    temp=np.zeros((2))
    temp[0]=np.random.uniform(0,1)*a
    temp[1]=np.random.uniform(0,1)*b
    random_nuclei[i]=temp

for i in range(0,random_size_number2):                   #genero punti random

    temp=np.zeros((2))
    temp[0]=np.random.uniform(0,1)*a
    temp[1]=np.random.uniform(0,1)*b
    random_nuclei2[i]=temp


print(num_nuclei)                                       #nuclei reali

labels=rois_detections(random_nuclei,0,df,posizioni,rois,rois_val)	#modifica rois e ritorna la struttura labels: vede quali nuclei simulati sono nel contorno
total_points=[]
for i in range(len(rois)):                               #total_points ha dentro tutti i nuclei fake nel contorno
    for j in range(len(rois[i])):
        total_points.append(rois[i][j])

print('len:', len(total_points))
print(((len(total_points)**2)-len(total_points))/2)

rois2=[]
labels2=rois_detections(random_nuclei2,0,df,posizioni,rois2,rois_val)	#modifica rois e ritorna la struttura labels: vede quali nuclei simulati sono nel contorno
total_points2=[]
for i in range(len(rois2)):                               #total_points ha dentro tutti i nuclei fake nel contorno
    for j in range(len(rois2[i])):
        total_points2.append(rois2[i][j])

punti_random=np.zeros((len(total_points),3))             #per distances vector pnti fake
for i in range(0,len(total_points)):

	punti_random[i][0]=total_points[i][0]
	punti_random[i][1]=total_points[i][1]

punti_random2=np.zeros((len(total_points2),3))             #per distances vector pnti fake
for i in range(0,len(total_points2)):

	punti_random2[i][0]=total_points2[i][0]
	punti_random2[i][1]=total_points2[i][1]

punti_real=np.zeros((num_nuclei,3))                      #per distance vector punti real
for i in range(0,num_nuclei):
    punti_real[i][0]=nuclei[i][0]*3.33
    punti_real[i][1]=nuclei[i][1]*3.33
    punti_real[i][2]=df['Angle'].iloc[i]

num_bins=500
#distance_vector_simu=get_distances_and_quantities_vectors(punti_random,len(total_points),'simu',simple=True)
#distance_vector_real=get_distances_and_quantities_vectors(punti_real,num_nuclei,'real',simple=True)

distance_vector_simu=get_distances_vector(punti_random,len(total_points),'simu')
distance_vector_simu2=get_distances_vector(punti_random2,len(total_points2),'simu2')
distance_vector_real=get_distances_and_angles(punti_real,num_nuclei,'real')

step=diagonal/num_bins
#num_iterations=int((max(a,b)/(2*step)))
num_iterations=int(diagonal/step)
print('n bins',num_bins)
bins=np.linspace(0,diagonal,num_bins)

digitized_simu= np.digitize(distance_vector_simu[0,:], bins)					#ritorna indici che servono per binnare poi il mio array AL DATO BINNAGGIO
digitized_simu2= np.digitize(distance_vector_simu2[0,:], bins)
digitized_real= np.digitize(distance_vector_real[0,:], bins)

simu_bin_rays=np.zeros((int(num_iterations)))                         #il primo elemento sarà per definiziione nullo: numero di elementi ad una distanza zero (escluso se stesso)
simu_bin_rays2=np.zeros((int(num_iterations))) 
real_bin_rays=np.zeros((int(num_iterations)))
nem_energy=np.zeros((num_iterations))    
cum=0
cum2=0
cum_real=0
xs=np.zeros((int(num_iterations)))                           #serve per il plot

for i in range(1,int(num_iterations)):                       #perchè digitized ritorna il bin a cui appartiene, bins invece ha gli estremi

    print(i,'/',int(num_iterations))
    xs[i]=bins[i]
    simu_bin_rays[i]=cum+len(distance_vector_simu[0,:][digitized_simu == i])         #[i-1] semplicemente per come si contano gli elementi
    cum+=len(distance_vector_simu[0,:][digitized_simu == i])
    simu_bin_rays2[i]=cum2+len(distance_vector_simu2[0,:][digitized_simu2 == i])         #[i-1] semplicemente per come si contano gli elementi
    cum2+=len(distance_vector_simu2[0,:][digitized_simu2 == i])
    real_bin_rays[i]=cum_real+len(distance_vector_real[0,:][digitized_real == i])
    nem_energy[i]=np.mean(distance_vector_real[1,:][digitized_real == i])
    cum_real+=len(distance_vector_real[0,:][digitized_real == i])

norm_simu_1=((len(total_points)**2)-len(total_points))/2
norm_simu_2=((len(total_points2)**2)-len(total_points2))/2
norm_real=((num_nuclei**2)-num_nuclei)/2
for i in range(0,int(num_iterations)):
    simu_bin_rays[i]=simu_bin_rays[i]/norm_simu_1
    simu_bin_rays2[i]=simu_bin_rays2[i]/norm_simu_2
    real_bin_rays[i]=real_bin_rays[i]/norm_real

quoz=np.zeros((int(num_iterations)))
quoz2=np.zeros((int(num_iterations)))

for i in range(0,int(num_iterations)):

    if(simu_bin_rays[i]==0):
        quoz[i]=0
        quoz2[i]=0

    else:

        quoz[i]=real_bin_rays[i]/simu_bin_rays[i]
        quoz2[i]=simu_bin_rays2[i]/simu_bin_rays[i]

print('max real:',max(distance_vector_real[0,:]))
print('num nuke', num_nuclei)
print('dfmax',max(df['YM']))
print('dfmin',min(df['YM']))

fig,ax=plt.subplots(1,4)

ax[0].plot(xs,simu_bin_rays,'o',c='red',label='simulation')
ax[0].plot(xs,real_bin_rays,'o',c='green',label='real')
ax[0].set_xlabel('R')
ax[0].set_ylabel('av. #counts per point')
ax[0].legend()

ax[1].plot(xs,quoz,'-',c='blue')
ax[1].set_xlabel('R')
ax[1].set_ylabel('g(R)')
ax[1].axhline(1,ls='--',color='red')

ax[2].plot(xs,quoz2,'-',c='blue')
ax[2].set_title('control')
ax[2].axhline(1,ls='--',color='red')

ax[3].plot(xs,nem_energy,'-',c='purple',label='nem_energy')
ax[3].axhline(2/np.pi,ls='--',color='red')
#ax[3].plot(xs,simu_bin_rays2,'o',c='blue',label='simulation2')

ax[3].legend()

fig,ax=plt.subplots()

ax.imshow(img)
ax.scatter(punti_random[:,0],punti_random[:,1],c='k',s=5)
ax.scatter(punti_real[:,0],punti_real[:,1])

plt.show()