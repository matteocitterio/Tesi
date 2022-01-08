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


sigmas=[]


p_vals=[]
ks=[]
ranger=100
diffs=[]
mean_ars_s=[]
variab_ars_s=[]
y_pred_s=[]


for l in range(0,40):

    fixing=l
    diff=[]
    mean_ars=[]
    variab_ars=[]
    y_pred=[]
    un=1
    du=un
    tri=un*(l/3)

    for j in range(0,ranger):

        s=0.01+0.01*j
        mu=un*s
        mu2=du*s
        mu3=tri*s
        sigmas.append(s)
        ars=[]
        aspect_ratios=[]

        for i in range(0,N):

            print(fixing,'-',j,'/',ranger,':',i,'/',N)

            a=np.random.lognormal(mean=mu,sigma=s)
            b=np.random.lognormal(mean=mu2,sigma=s)
            c=np.random.lognormal(mean=mu3,sigma=s)
        
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
        diff.append(abs((coeff*np.nanmean(aspect_ratios)+intercetta)-(np.nanstd(aspect_ratios))))

        p_vals.append(stats.ks_2samp(y_ar,y_sampled)[1])

        ars=np.asarray(ars)
        ars_normalized=(ars-1)/(np.average(ars)-1)

    diffs.append(diff)
    mean_ars_s.append(mean_ars)
    variab_ars_s.append(variab_ars)
    y_pred_s.append(y_pred)

fig,ax=plt.subplots()
for i in range(0,40):

    ax.plot(mean_ars_s[i],variab_ars_s[i],'o',label=str(i/3))
    ax.plot(mean_ars_s[i],y_pred_s[i],'o',color='red')
    
ax.set_xlabel('<AR>')
ax.set_ylabel('$\sigma_{AR}$')
ax.legend()

fig,ax=plt.subplots()

for i in range(0,40):

    ax.plot(mean_ars_s[i],diffs[i],'o',label=str(i/3))


ax.set_xlabel('<AR>')
ax.set_ylabel('diff')
ax.legend()

plt.show()