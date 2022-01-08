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
from matplotlib.animation import FuncAnimation
from matplotlib.animation import PillowWriter 
from matplotlib.animation import FFMpegWriter
from functions import *

fig,ax=plt.subplots(2,1,figsize=(16,9))

kolor=['tab:blue','g','r','c','m','y','k','green','orange','purple']
h=0
i=0
folder='ROI_3_INT_Normal'

def do_the_correlation(i,folder):

	kw='AR'
	image_N=i

	img=Image.open(r'D:\matteo_citterio\SCAN - Copia\\'+folder+'\\'+str(image_N)+'.tif')												#
	a,b=img.size
	diagonal=((a**2)+(b**2))**0.5

	ax[0].clear()
	ax[1].clear()

	print(i)
	df=pd.read_csv(r'D:\matteo_citterio\risultati_segment\\'+folder+'\Results_'+str(image_N)+'.csv')
	num_nuclei=len(df.axes[0])
	print(num_nuclei)
	funct=create_nuclei_function(num_nuclei,df)			    		#crea un array con le pos dei cdm dei nuclei dalla df di stardist
	# array=get_distances_and_quantities_vectors(funct,num_nuclei)	#serve per trovare il minimo

	# sorted=np.sort(array[0,:])								#trovo il passo cellulare
	# min=next(x for x in sorted if x!=0)

	# k1=10										#moltiplicatore bin dentisty
	# k2=2										#moltiplicatore bin cosi

	num_bins=int(num_nuclei/13)

	# num_bins_density=int(diagonal/min)*k1		#lo step è 5 volte la minima distanzas
	# num_bins_corr=int(diagonal/min)*k2			#lo step è 1/2 della minima distanza

	df.to_csv(r'D:\matteo_citterio\risultati_segment\ROI_3_INT_Normal\Results.txt', header=None, index=None, sep=' ',columns=['XM','YM']  ,mode='a')

	bins, bin_means, bin_rays, media, media_quad = pair_correlation_function_for_iterations(funct,num_nuclei,a,b,num_bins,kw)

	density_matlab()
	
	dj=pd.read_csv(r'D:\matteo_citterio\Python Scripts\temp2.csv')
	print(dj)
	os.remove(r'D:\matteo_citterio\risultati_segment\ROI_3_INT_Normal\Results.txt')
	os.remove(r'D:\matteo_citterio\Python Scripts\temp2.csv')

	dj.columns =['g(R)', 'R','c']
	bins_density=dj['R']
	density_profile=dj['g(R)']

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


	ax[0].axvline(x[maxi[-1]+shift],c='red',ls='--')
	ax[0].plot(bins_density,df1,'-',c='green',label='scan3_1')
	ax[0].set_title('density function')
	ax[0].set_ylabel('g(r)')
	ax[0].set_xlim(-1,size)
	ax[0].set_xlabel('r [pixel]')	
	ax[0].legend()
	ax[0].axhline(1,c='k')

	ax[1].plot(bins,(bin_means-(media**2))/(media_quad-(media**2)),'o-',ms=3,c=kolor[0],label='insitu'+str(i+1))
	ax[1].set_xlim(-1,size)
	ax[1].set_xlabel('r [pixel]')
	ax[1].set_ylabel('$(<f_i \cdot f_j>-<f>^2)$/$\sigma^2}$')
	ax[1].axvline(x[maxi[-1]+shift],c='red',ls='--')

	fig.suptitle(kw)

	plt.tight_layout()
	fig.savefig('D:\matteo_citterio\plot_al_fly\spatial autocorrelations\\3_INT_Normal\AR\\'+str(image_N)+'.png')
	# plt.show()

####################################################################################################################################################################################################################

	# right_xlim=300
	# understat=[]
	# understat_bins=[]
	# understat_brown=[]
	# understat_bins_brown=[]
	# temp=np.zeros((int(right_xlim/(diagonal/(int(diagonal/min)*5)))))
	# temp_dos=np.zeros((int(right_xlim/(diagonal/(int(diagonal/min)*5)))))
	# temp_tres=np.zeros((int(right_xlim/(diagonal/(int(diagonal/min)*5)))))

	# for j in range(0,int(right_xlim/(diagonal/(int(diagonal/min)*5)))):

	# 	temp[j]=bin_means[j]
	# 	temp_dos[j]=density_profile[j]
	# 	temp_tres[j]=bin_rays[j]

	# max1=np.nanmax(temp)
	# min1=np.nanmin(temp)
	# max1=(max1-(media**2))/(media_quad-(media**2))
	# min1=(min1-(media**2))/(media_quad-(media**2))
	# max2=np.nanmax(temp_dos)
	# max3=np.nanmax(temp_tres)

	# print(max1,min1)

	# for j in range(0,len(bins)):

	# 	if (bin_rays[j]<=30):
	# 		understat.append(bin_means[j])
	# 		understat_bins.append(bins[j])
	# 	elif(bin_rays[j]<=45):
	# 		understat_brown.append(bin_means[j])
	# 		understat_bins_brown.append(bins[j])

	# df=pd.DataFrame(density_profile)
	# df=df.rolling(k1+3).mean()

	# df[0]=df[0].fillna(0)
	# d1=df[0].iloc[0]-df[0].iloc[1]

	# print(d1)
	# j=1
	# while(d1>=0):

	# 	j+=1
	# 	d1=df[0].iloc[j]-df[0].iloc[j+1]		#massimo locale		

	# d2=0
	# while(d2>=0):

	# 	j+=1
	# 	d2=df[0].iloc[j]-df[0].iloc[j+1]
		
	# 	print(d2)

	# x_primi_vicini=bins[j]

	# print('dopo il ciclo',d2)
	# print('x vicino', x_primi_vicini)


	# ax[1][1].plot(bins,bin_rays,'o',c=kolor[0],label=str(i+1))
	# ax[1][1].set_xlabel('r [pixel]')
	# ax[1][1].set_ylabel('#counts')
	# ax[1][1].set_ylim(-10,max3)
	# ax[1][1].set_xlim(-10,right_xlim)
	# ax[1][1].axhline(30,c='green',label='30, mean stat threshold')
	# ax[1][0].axhline(0,c='green')
	# ax[1][0].axhline(1,c='blue')
	# ax[1][0].axvline(x_primi_vicini,c='red')
	# ax[0][0].plot(bins_density,density_profile,'o-',ms=3,c='green',label='data')
	# ax[0][0].plot(bins_density,df,'-',c='red',label='moving av.')
	# ax[0][0].axvline(x_primi_vicini,c='red')
	# ax[0][0].legend()
	# ax[0][0].set_xlabel('r [pixel]')
	# ax[0][0].set_ylabel('C(r)')
	# ax[0][0].axhline(1)
	# ax[0][0].set_xlim(-5,right_xlim)
	# ax[0][0].set_ylim(-0.3,max2+(max2/4))
	# ax[1][0].legend()
	# ax[0][1].scatter(funct[:,0],funct[:,1],s=1,c='k')
	# ax[0][1].set_xlabel('Y [pixel]')
	# ax[0][1].set_ylabel('X [pixel]')
	fig.suptitle(kw)
	# h+=1

for j in range(0,63):

	do_the_correlation(j+1,folder)

# plt.show()
# i=1
# ani = FuncAnimation(fig, do_the_correlation(i,folder), frames=200, interval=63, repeat=False)
# i+=1

# f = "D:\matteoinsitu.mp4" 
# writervideo = FFMpegWriter(fps=0.5) 
# ani.save(f, writer=writervideo) 