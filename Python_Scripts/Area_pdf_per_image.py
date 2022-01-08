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

from PIL import Image
import matplotlib as mpl
from matplotlib.path import Path
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from functions import *
import scipy.stats as stats

colors=[]
columns=['9-norm','9-tum','13-tum','3-norm','3-tum','5-tum','7-tum','8-tum','12-norm','12-tum']
img_ns=[7,64,11,13,6,6,18,24,2,9]
real_img_n=[7, 18, 9, 64, 11, 13, 6, 6, 24, 2]
img_ns_string=['9_INT_Normal','9_INT_Tumor','13_INT_Tumor','3_INT_Normal','3_INT_Tumor','5_INT_Tumor','7_INT_Tumor','8_INT_Tumor','12_INT_Normal','12_INT_Tumor']
patient_string=['9-Normal','9-Tumor','13-Tumor','3-Normal','3-Tumor','5-Tumor','7-Tumor','8-Tumor','12-Normal','12-Tumor']
img_ns_tumor_string=['9_INT_Tumor','13_INT_Tumor','3_INT_Tumor','5_INT_Tumor','7_INT_Tumor','8_INT_Tumor','12_INT_Tumor']
real_img_n_tumor=[ 18, 9,11, 13, 6, 6,2]
real_img_n_normal=[7,64,24]
img_ns_normal_string=['9_INT_Normal','3_INT_Normal','12_INT_Normal']
quant='AR'
divisors=[90,120,115,1000,180,130,25,25,80,30]#,15,145,15]
ar_ratio=[]
gdl=[]
ar_means=[]
ar_sigmas=[]

colors=['green','black','black','green','black','black','black','black','green','black']

fig,ax=plt.subplots()
fig3,bx=plt.subplots()
cm=truncate_colormap(plt.cm.get_cmap('RdPu'), 0.4, 1)
color_gradients = cm(np.arange(0,10)/10) 

for eta in range(0,len(img_ns_string)): 

	print(img_ns_string[eta])

	dj=pd.read_csv(r'D:\matteo_citterio\Python Scripts\AR_per_image\AR_voronoiROI_'+img_ns_string[eta]+'.csv')
	#dj=pd.read_csv(r'D:\matteo_citterio\Python Scripts\Area_ratio_per_image\Area_ratio_voronoiROI_'+img_ns_string[eta]+'.csv')

	num_nuclei=len(dj.axes[0])

	data=[]
	area_ratio=[]

	for i in range(1,real_img_n[eta]):                  #ciclo a colonna che vuol dire a immagine

		for j in range(0,num_nuclei):

			if(not math.isnan(dj[str(i)].iloc[j])):

				data.append(dj[str(i)].iloc[j])
				#area_ratio.append(dj[str(i)].iloc[j])

	print(len(data))

	data=np.asarray(data)
	data=(data)/(np.nanmean(data))

	# df = pd.DataFrame(data=data)
	# results, bin_edges=pd.qcut(df[0],q=10,retbins=True)
	# bin_counts=pd.qcut(df[0], q=10).value_counts()

	# print(pd.qcut(df[0], q=10))

	ar_means.append(np.average(data))
	ar_sigmas.append(np.std(data))

	y,x=np.histogram(data,int(len(data)/divisors[eta]),density=True)
	hist_dist = stats.rv_histogram(np.histogram(data,int(len(data)/divisors[eta]),density=True))

	xnew=np.zeros((x.size+1))
	step=x[1]-x[0]

	for j in range(1,x.size+1):

		xnew[j]=x[j-1]

	x=xnew

	params=stats.lognorm.fit(data, loc=0)
	params_gamma=stats.gamma.fit(data, loc=0)

	if (img_ns_string[eta]=='9_INT_Normal'):

		params2=stats.lognorm.fit(data)

	else:

		params2=stats.lognorm.fit(data, loc=0)

	# print(params)
	print(params)
	gdl.append(params2[0])

	pdf_fitted =stats.lognorm.pdf(x, *params)
	pdf_fitted2 =stats.lognorm.pdf(x, *params2)
	ax.plot(x,hist_dist.pdf(x),'-o',color=color_gradients[eta], label='Patient '+patient_string[eta],lw=3)
	#ax.plot(x, pdf_fitted, color='r',label='gamma fit',lw=3)
	#ax.plot(x, pdf_fitted2, color='k',label='chi2 fit',ls='--',lw=3)

	bx.plot(x,hist_dist.pdf(x),'-o', label='data PDF',lw=3)
	bx.plot(x, pdf_fitted, color='r',label='gamma fit',lw=3)
	bx.plot(x, pdf_fitted2, color='k',label='chi2 fit',ls='--',lw=3)
	bx.set_yscale('log')
	# bx.set_xscale('log')
	bx.legend()

	# ax.hist(data,int(len(data)/divisors[eta]),density=True,label='data')

	ax.set_xlabel('area/<area>')
	ax.set_ylabel('Normalized counts')
	ax.set_title('Voronoi area distributions')
	ax.set_xlim(-0.1,6)


	bx.set_xlabel('Area [$\mu m ^{2}$]')
	bx.set_ylabel('Normalized counts')
	bx.set_title('patient:'+img_ns_string[eta]+' Area distribution LOG')

	ax.legend()
	# fig.savefig('D:\matteo_citterio\plot_al_fly\distribuzioni\Distribuzioni a paziente\\AR\\patient_'+img_ns_string[eta]+'_AR_distribution.png')

	y_ar,x_ar=np.histogram(data,int(len(data)/divisors[eta]))
	sampled=stats.lognorm.rvs(*params,size=len(data))
	sampled2=stats.lognorm.rvs(*params2,size=len(data))
	y_sampled,x_sampled=np.histogram(sampled,int(len(data)/divisors[eta]))
	y_sampled2,x_sampled2=np.histogram(sampled2,int(len(data)/divisors[eta]))

	qs=np.linspace(0,1,11)
	qs=qs.tolist()
	bin_extremes=stats.lognorm.ppf(qs, *params)
	expectation_value=int(len(data)/10)
	counts=pd.cut(data,bin_extremes).value_counts()


	results_table = pd.DataFrame(zip(bin_extremes, counts),columns=['bins', 'counts'])
	print(results_table)
	print('The expectation:',expectation_value)

	temp=0
	for mu in range(0,len(counts)):

		temp+=(((counts[mu]-expectation_value)**(2))/expectation_value)

	print('chisquared_lognorm:',temp)

	print('\n')
	print('gamma pvalue:',stats.ks_2samp(y_ar,y_sampled)[1])
	print('with lognorm cdf',stats.ks_1samp(data,cdf=stats.lognorm.cdf,args=params)[1])
	print('with gamma cdf',stats.ks_1samp(data,cdf=stats.gamma.cdf,args=params_gamma)[1])
	print('\n')


legend_elements =[Line2D([0], [0], marker='o', color='w', label='Tumor',markerfacecolor='black',markeredgecolor='black',markersize=5),
					Line2D([0], [0], marker='o', color='w', label='Normal',markerfacecolor='green',markeredgecolor='green', markersize=5)]



# an.plot(gdl,ar_ratio,'o',color=colors)
plt.show()

