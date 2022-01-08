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

fig,ax=plt.subplots()

medie=[]
diagonal=((1460**2)+(1936**2))**0.5
for i in range(0,36):	#insitu

	print(i)

	df=pd.read_csv('/home/matteo/Scrivania/tesi/codice/processed2/insitu/Results_'+str(i+1)+'.csv')
	num_nuclei=len(df.axes[0])

	funct=create_nuclei_function(num_nuclei,df)			    #crea un array con le pos dei cdm dei nuclei dalla df di stardist
	array=get_distances_and_quantities_vectors(funct,num_nuclei)	#serve per trovare il minimo
	sorted=np.sort(array[0,:])								#trovo il passo cellulare
	min=next(x for x in sorted if x!=0)
	k1=10										#moltiplicatore bin dentisty
	k2=2										#moltiplicatore bin cosi
	num_bins_density=int(diagonal/10)*k1		#lo step è 5 volte la minima distanzas
	num_bins_corr=int(diagonal/min)*k2			#lo step è 1/2 della minima distanza
	bins_density,density_profile=density_profile_function(funct,num_nuclei,1936,1460,num_bins_density)
	medie.append(density_profile)

array=np.zeros((len(medie[0])))
for i in range(0,len(medie[0])):

	medie_reduced=np.zeros((36))

	for j in range(0,36):

		medie_reduced[j]=medie[j][i]

	array[i]=medie_reduced.mean()
df=pd.DataFrame(density_profile)
df=df.rolling(k1+3).mean()
	
# ax.plot(bins_density,array,'o-',c='tab:green',ms=3,label='insitu')
ax.plot(bins_density,df,'-',c='green',label='insitu')

medie=[]

for i in range(0,9):	#mixed

	print(i)

	df=pd.read_csv('/home/matteo/Scrivania/tesi/codice/processed2/mixed/Results_'+str(i+1)+'.csv')
	num_nuclei=len(df.axes[0])
	funct=create_nuclei_function(num_nuclei,df)			    #crea un array con le pos dei cdm dei nuclei dalla df di stardist
	array=get_distances_and_quantities_vectors(funct,num_nuclei)	#serve per trovare il minimo
	sorted=np.sort(array[0,:])								#trovo il passo cellulare
	min=next(x for x in sorted if x!=0)
	k1=10										#moltiplicatore bin dentisty
	num_bins_density=int(diagonal/10)*k1		#lo step è 5 volte la minima distanzas
	bins_density,density_profile=density_profile_function(funct,num_nuclei,1936,1460,num_bins_density)
	medie.append(density_profile)
array=np.zeros((len(medie[0])))
for i in range(0,len(medie[0])):

	medie_reduced=np.zeros((9))

	for j in range(0,9):

		medie_reduced[j]=medie[j][i]

	array[i]=medie_reduced.mean()
df1=pd.DataFrame(density_profile)
df1=df1.rolling(k1+3).mean()
	
# ax.plot(bins_density,array,'o-',c='tab:blue',ms=3,label='mixed')
ax.plot(bins_density,df1,'-',c='blue',label='mixed')

medie=[]

for i in range(0,18):	#infiltrating

	print(i)

	df=pd.read_csv('/home/matteo/Scrivania/tesi/codice/processed2/infiltrating/Results_'+str(i+1)+'.csv')
	num_nuclei=len(df.axes[0])
	funct=create_nuclei_function(num_nuclei,df)			    #crea un array con le pos dei cdm dei nuclei dalla df di stardist
	array=get_distances_and_quantities_vectors(funct,num_nuclei)	#serve per trovare il minimo
	sorted=np.sort(array[0,:])								#trovo il passo cellulare
	min=next(x for x in sorted if x!=0)
	k1=10										#moltiplicatore bin dentisty
	num_bins_density=int(diagonal/10)*k1		#lo step è 5 volte la minima distanzas
	bins_density,density_profile=density_profile_function(funct,num_nuclei,1936,1460,num_bins_density)
	medie.append(density_profile)

array=np.zeros((len(medie[0])))
for i in range(0,len(medie[0])):

	medie_reduced=np.zeros((18))

	for j in range(0,18):

		medie_reduced[j]=medie[j][i]

	array[i]=medie_reduced.mean()
df2=pd.DataFrame(density_profile)
df2=df2.rolling(k1+3).mean()
	
# ax.plot(bins_density,array,'o-',c='tab:purple',ms=3,label='infiltrating')
ax.plot(bins_density,df2,'-',c='purple',label='infiltrating')

fig.suptitle('density function')
ax.set_ylabel('C(r)')
ax.set_xlabel('r [pixel]')	
ax.legend()
ax.axhline(1,c='k')
ax.set_xlim(-1,300)
plt.tight_layout()
plt.show()
