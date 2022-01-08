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

sigmas=[]
mean_ars=[]
variab_ars=[]
p_vals=[]
ks=[]
ranger=10
fixing=1.4


for j in range(0,ranger):

    s=0.001+0.05*j
    s1=0.001+0.005*j
    sigmas.append(s)
    ars=[]
    aspect_ratios=[]

    for i in range(0,N):

        print('due s',fixing,'-',j,'/',ranger,':',i,'/',N)

        a=stats.lognorm.rvs(s)
        b=stats.lognorm.rvs(s1)
        c=stats.lognorm.rvs(s)*fixing
        
        # a=a1/np.sqrt(b1*c1)
        # b=b1/np.sqrt(a1*c1)
        # c=c1/np.sqrt(a1*b1)

        # ass.append(a)
        # bs.append(b)
        # cs.append(c)
        ars.append(max(a,b,c)/min(a,b,c))

        n1=np.random.normal(0,2)
        n2=np.random.normal(0,2)
        n3=np.random.normal(0,2)

        magnitude=np.sqrt(n1**2+n2**2+n3**2)

        n1=n1/magnitude
        n2=n2/magnitude
        n3=n3/magnitude

        eta=(n1**2)*((1/b**2)+(1/c**2))+(n2**2)*((1/a**2)+(1/c**2))+(n3**2)*((1/a**2)+(1/b**2))
        gamma=(n1**2)/((b**2)*(c**2))+(n2**2)/((a**2)*(c**2))+(n3**2)/((a**2)*(b**2))

        beta1=(eta/2)-np.sqrt(((eta/2)**2)-gamma)
        beta2=(eta/2)+np.sqrt(((eta/2)**2)-gamma)

        aspect_ratio=np.sqrt(beta2/beta1)
        aspect_ratios.append(aspect_ratio)

    num_bins=300
    bins=np.linspace(0,15,num_bins)

    aspect_ratios=np.asarray(aspect_ratios)
    aspect_ratios_mean=np.nanmean(aspect_ratios)
    variab_ars.append(np.nanstd(aspect_ratios))
    mean_ars.append(np.nanmean(aspect_ratios))
    aspect_ratio_normalized=(aspect_ratios-1)/(aspect_ratios_mean-1)

    y_ar,x_ar=np.histogram(aspect_ratio_normalized,bins)
    params=stats.chi2.fit(aspect_ratio_normalized,floc=0)
    ks.append(params[0])
    sampled=stats.chi2.rvs(*params,size=N)
    y_sampled,x_sampled=np.histogram(sampled,bins)

    y_pred.append(coeff*np.nanmean(aspect_ratios)+intercetta)

    p_vals.append(stats.ks_2samp(y_ar,y_sampled)[1])

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

fig,ax=plt.subplots()
ax.plot(sigmas,mean_ars,'o',mfc='purple',mec='k')
ax.set_xlabel('3D lognorm $\sigma$')
ax.set_ylabel('<AR>')
ax.set_title('vol+fixed: '+str(fixing))

cm=plt.cm.get_cmap('RdYlBu')

fig,ax=plt.subplots()
s=ax.scatter(mean_ars,variab_ars,c=sigmas, vmin=min(sigmas),vmax=max(sigmas),cmap=mpl.cm.cool)
ax.scatter(mean_ars,y_pred,c='red',label='exp-data linearity')
cb=fig.colorbar(s)
cb.set_label('$\sigma_{lognorm}$ 3D ')
ax.set_xlabel('<AR>')
ax.legend()
ax.set_ylabel('$\sigma_{AR}$')
ax.set_title('Atia relation fixed val')
plt.suptitle('lognormal 3D input: fixed: '+str(fixing))

fig,ax=plt.subplots()
s=ax.scatter(mean_ars,ks,c=sigmas, vmin=min(sigmas),vmax=max(sigmas),cmap=mpl.cm.cool)
cb=fig.colorbar(s)
cb.set_label('$\sigma_{lognorm}$ 3D ')
ax.set_xlabel('<AR>')
ax.axvline(1.3,ls='--',color='red')
ax.axvline(1.8,ls='--',color='red')
ax.set_ylabel('fitted df 2D distribution-fixed')
ax.set_title('fitted df, fixed: '+str(fixing))

fig,ax=plt.subplots()
s=ax.scatter(mean_ars,p_vals,c=sigmas, vmin=min(sigmas),vmax=max(sigmas),cmap=mpl.cm.cool)
cb=fig.colorbar(s)
cb.set_label('$\sigma_{lognorm}$ 3D ')
ax.set_xlabel('<AR>')
ax.set_ylabel('fit p_value')
ax.axhline(0.05,ls='--',color='red')
ax.set_title('godness of fit, fixed: '+str(fixing))

plt.show()