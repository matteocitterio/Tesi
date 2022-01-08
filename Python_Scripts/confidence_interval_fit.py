import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import warnings
import os
import sys
from scipy.integrate._ivp.radau import P
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
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl
from functions import *



N=10000                     #numero di iterazioni

# fig,ax=plt.subplots(1,2)

coeff=0.69389538
intercetta=-0.69473173
y_pred=[]

cmap = plt.cm.get_cmap('Blues',10)
cm=truncate_colormap(plt.cm.get_cmap('Blues'), 0.4, 1)
colors = cm(np.arange(0,10)/10) 

sigmas=[]
mean_ars=[]
variab_ars=[]
p_vals=[]
ks=[]
ranger=10
fixing=1
fig,ax=plt.subplots()
fig2,ax2=plt.subplots()

for j in range(0,ranger):

	s=0.2#0.001+0.005*j
	sigmas.append(0.4+0.05*j)
	ars=[]
	aspect_ratios=[]

	alpha=beta=gamma=1



	for i in range(0,N):

		n1=np.random.normal(0,2)
		n2=np.random.normal(0,2)
		n3=np.random.normal(0,2)

		magnitude=np.sqrt(n1**2+n2**2+n3**2)

		n1=n1/magnitude
		n2=n2/magnitude
		n3=n3/magnitude




		#print(fixing,'-',j,'/',ranger,':',i,'/',N)

		a=np.random.lognormal(0,s)
		b=np.random.lognormal(0,s)
		c=np.random.lognormal(0.4+0.05*j,s)
		
		ars.append(max(a,b,c)/min(a,b,c))


		eta=(n1**2)*((1/b**2)+(1/c**2))+(n2**2)*((1/a**2)+(1/c**2))+(n3**2)*((1/a**2)+(1/b**2))
		gamma=(n1**2)/((b**2)*(c**2))+(n2**2)/((a**2)*(c**2))+(n3**2)/((a**2)*(b**2))

		beta1=(eta/2)-np.sqrt(((eta/2)**2)-gamma)
		beta2=(eta/2)+np.sqrt(((eta/2)**2)-gamma)

		aspect_ratio=np.sqrt(beta2/beta1)
		aspect_ratios.append(aspect_ratio)

	num_bins=300
	bins=np.linspace(0,15,num_bins)

	
	aspect_ratios=np.asarray(aspect_ratios)

	bins_cs=np.linspace(1,6,150)
	y,x=np.histogram(aspect_ratios,bins_cs,density=True)

	#hist_dist_a = stats.rv_histogram(aspect_ratios,bins_cs,density=True)

	xnew=np.zeros((x.size+1))+1
	step=x[1]-x[0]
	for h in range(1,x.size+1):

		xnew[h]=x[h-1]

	x=xnew

	params=stats.gamma.fit(aspect_ratios,loc=1)
	pdf_fitted =stats.gamma.pdf(x, *params)
	ax.plot(x, pdf_fitted,c=colors[j],lw=3)

	ax.set_xlabel('$AR$',fontsize=13)
	ax.set_ylabel('Gamma-fitted PDF',fontsize=13)
	ax.set_xlim(0.9,4.1)
	#ax.set_title('Patient:'+img_ns_string[eta]+' AR distribution')
	# ax.set_title('Normalized Voronoi AR distributions')
	ax.set_title('Nuclei $AR$ Gamma fit',fontsize=13)


	aspect_ratios_mean=np.nanmean(aspect_ratios)
	variab_ars.append(np.nanstd(aspect_ratios))
	mean_ars.append(np.nanmean(aspect_ratios))
	aspect_ratio_normalized=(aspect_ratios-1)/(aspect_ratios_mean-1)

	bins_cs=np.linspace(0,6,150)
	y,x=np.histogram(aspect_ratio_normalized,bins_cs,density=True)

	#hist_dist_a = stats.rv_histogram(aspect_ratios,bins_cs,density=True)

	xnew=np.zeros((x.size+1))
	step=x[1]-x[0]
	for h in range(1,x.size+1):

		xnew[h]=x[h-1]

	x=xnew

	params=stats.gamma.fit(aspect_ratio_normalized,loc=0)
	pdf_fitted =stats.gamma.pdf(x, *params)
	ax2.plot(x, pdf_fitted,c=colors[j],lw=3)
	pdf_fitted =stats.gamma.pdf(x, *params)
	print(params)

	ax2.set_xlabel('$(AR-1)/(\overline{AR}-1)$',fontsize=13)
	ax2.set_ylabel('Gamma-fitted PDF',fontsize=13)
	ax2.set_xlim(-0.1,5)
	#ax.set_title('Patient:'+img_ns_string[eta]+' AR distribution')
	# ax.set_title('Normalized Voronoi AR distributions')
	ax2.set_title('Normalized nuclei $AR$ Gamma fit',fontsize=13)


	y_ar,x_ar=np.histogram(aspect_ratio_normalized,bins)
	params=stats.gamma.fit(aspect_ratio_normalized,loc=0)
	ks.append(params[0])
	sampled=stats.chi2.rvs(*params,size=N)
	y_sampled,x_sampled=np.histogram(sampled,bins)

	y_pred.append(coeff*np.nanmean(aspect_ratios)+intercetta)

	qs=np.linspace(0,1,11)
	qs=qs.tolist()
	bin_extremes=stats.gamma.ppf(qs, *params)
	expectation_value=int(len(aspect_ratio_normalized)/10)
	counts=pd.cut(aspect_ratio_normalized,bin_extremes).value_counts()

	results_table = pd.DataFrame(zip(bin_extremes, counts),columns=['bins', 'counts'])
	print(results_table)
	print('The expectation:',expectation_value)

	temp=0
	for mu in range(0,len(counts)):

		temp+=(((counts[mu]-expectation_value)**(2))/expectation_value)

	# print('chisquared_lognorm:',temp)

	# p_vals.append(stats.ks_2samp(y_ar,y_sampled)[1])
	p_vals.append(temp)

	ars=np.asarray(ars)
	ars_normalized=(ars-1)/(np.average(ars)-1)

	# ax[0].hist(aspect_ratio_normalized,120,density=True,label='s='+str(s),alpha=0.5)
	# ax[0].set_xlabel('(AR-1)/(<AR>-1)')
	# ax[0].set_ylabel('PDF')
	# ax[0].legend()
	# ax[0].set_title('Simulated 2d ellipse aspect ratio')
	# ax[1].hist(ars_normalized,160,density=True,label='s='+str(s))
	# ax[1].set_title('Input 3D aspect ratio')
	# ax[1].set_xlabel('(AR-1)/(<AR>-1)')

params=(1.86, 0, 0.55)
pdf_fitted =stats.gamma.pdf(x, *params)
ax2.plot(x,pdf_fitted,c='red',lw=3,label='k=1.86')
ax2.legend(fontsize=13)

fig,ax=plt.subplots()
ax.plot(sigmas,mean_ars,'o',mfc='tab:blue')
ax.set_xlabel('Input 3D lognorm $s-$parameter')
ax.set_ylabel('Output 2D $\overline{AR}$')
ax.set_title('Simulation of the $\overline{AR}$ range')
ax.axhline(1.3,ls='--',color='red')
ax.axhline(1.7,ls='--',color='red')




cm=truncate_colormap(plt.cm.get_cmap('Blues'), 0.4, 1)

fig,ax=plt.subplots()
s=ax.scatter(mean_ars,variab_ars,c=sigmas, vmin=min(sigmas),vmax=max(sigmas),cmap=cm)
ax.plot(mean_ars,y_pred,color='red',label='experimental data linear regression')
cb=fig.colorbar(s)
cb.set_label('3D lognormal mean ')
ax.set_xlabel('$\overline{AR}$')
ax.legend()
ax.set_ylabel('$\sigma_{AR}$')
ax.set_title('$\overline{AR}-\sigma_{AR}$ relation')

fig,ax=plt.subplots()
s=ax.scatter(mean_ars,ks,c=sigmas, vmin=min(sigmas),vmax=max(sigmas),cmap=cm)
cb=fig.colorbar(s)
cb.set_label('3D lognormal $s-$parameter ')
ax.set_xlabel('$\overline{AR}$')
ax.axvline(1.3,ls='--',color='red')
ax.axvline(1.7,ls='--',color='red')
ax.set_ylabel('2D fitted  $k-$parameter')
ax.set_title('Two-dimensional $k$-parameter')

fig,ax=plt.subplots()
s=ax.scatter(mean_ars,p_vals,c=sigmas, vmin=min(sigmas),vmax=max(sigmas),cmap=cm)
cb=fig.colorbar(s)
cb.set_label('3D lognormal $s-$parameter ')
ax.set_xlabel('$\overline{AR}$')
ax.set_ylabel('$\chi^{2}$')
ax.axhline(16.9,ls='--',color='red')
ax.set_title('Goodness of fit')

plt.show()