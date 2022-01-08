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
from functions import *

data_size=3000                     #dimensione dei dati generati
side_1=1936
side_2=1460
num_bins=1000

alpha=np.random.random(data_size)*side_1        #uniforme reali
beta=np.random.random(data_size)*side_2
gamma=np.random.random(data_size)+1

punti=np.zeros((data_size,3))
h=0
print('Creating points')
for i in range(0,len(alpha)):

	punti[i][0]=alpha[i]
	punti[i][1]=beta[i]
	h+=1

print('points created')

image_diagonal=((side_1**2)+(side_2**2))**0.5				#diagonale dell'immagine
print('getting distances')
a=get_distances_and_quantities_vectors(punti,data_size,'real',simple=True)	#matrice delle distanze e vettore con grandezze
step=image_diagonal/num_bins
num_iterations=int((side_2/(2*step)))
print('n bins',num_bins)
bins=np.linspace(0,image_diagonal,num_bins)

print('min',a[0,:].min())
print('mean',a[0,:].mean())
print('distanze totali:',len(a[0,:]))

digitized= np.digitize(a[0,:], bins)					#ritorna indici che servono per binnare poi il mio array AL DATO BINNAGGIO

bin_rays=np.zeros((num_iterations))                           #il primo elemento sarà per definiziione nullo: numero di elementi ad una distanza zero (escluso se stesso)
theor_bin_rays=np.zeros((num_iterations))
cum=0
xs=np.zeros((num_iterations))

area=side_1*side_2
density=data_size/area
per_poly=np.zeros((num_iterations))
for i in range(1,num_iterations):                           #perchè digitized ritorna il bin a cui appartiene, bins invece ha gli estremi
    per_poly[i]=(bins[i])

    print(i)
    xs[i]=bins[i]
    bin_rays[i]=cum+len(a[0,:][digitized == i])      #[i-1] semplicemente per come si contano gli elementi
    cum+=len(a[0,:][digitized == i])
    theor_bin_rays[i]=(1/area)*((np.pi*area*(bins[i]**2))-((4/3)*(side_1+side_2)*(bins[i]**3))+(1/2)*(1.165)*(bins[i]**4))




bin_rays=bin_rays/data_size
theo=theor_bin_rays*density*(0.5)


diff=(theo-bin_rays)



h=1.58
eta=1.19
rimasuglis=[]
hs=[]

for j in range(0,100):
    
    # print(h)
    hs.append(h)
    rimasuglio=0
    temp_h=np.zeros((num_iterations))

    for l in range(1,num_iterations):
        temp_h[l]=-(bins[l]**3)*h*data_size/(area*1000)
        temp_h[l]+=(bins[l]**2)*(side_1+side_2)*eta*data_size/(area*2*1000)

    for l in range(1,num_iterations):
        rimasuglio+=abs(temp_h[l]-diff[l])

    rimasuglis.append(rimasuglio)
    h+=0.001

print(bin_rays[-1],theo[-1])

# print(np.polyfit(per_poly,diff,deg=3))
pol=np.polyfit(per_poly,diff,deg=3)

temp_diff=np.zeros((num_iterations))

for k in range(1,num_iterations):

    temp_diff[k]=pol[0]*(bins[k]**3)+pol[1]*(bins[k]**2)+pol[2]*(bins[k]**1)+pol[3]*(bins[k]**0)

fig,ax=plt.subplots(1,3)

ax[0].plot(xs,bin_rays,'o',c='red',label='simulation')
ax[0].plot(xs,theo,'o',c='green',label='model')
ax[0].set_xlabel('R')
ax[0].set_ylabel('av. #counts per point')
ax[0].legend()
# ax[1].plot(xs,diff,'o',c='blue')

# theo=theo-temp_diff

quoz=bin_rays/theo
print(quoz[0])

# ax[1].plot(xs,temp_h,'o',c='purple')
ax[1].plot(xs,quoz,'o',c='green')
ax[2].plot(xs,diff,'o',c='purple')
plt.show()

# bin_rays=[len(a[0,:][digitized == i]) for i in range(1, num_bins+1)]			#conteggi a raggio

# dr=(image_diagonal)/num_bins
# print('dr',dr)
# numer_of=int((smallest_side/(dr*2))-1)
# print('nuovo numero calcoalto',numer_of)
# print('numero input',num_bins)
# print(smallest_side)

# bin_rays_not_modi=[]

# for j in range(0,numer_of):


# 	print(j,'/',numer_of)

# 	bin_rays_not_modi.append(bin_rays[j])
# 	bin_rays[j]=bin_rays[j]/analytic_function(side1,side2,(bins[j]+2*dr),1,data_size,dr=dr/2)
	
