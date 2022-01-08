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
import scipy.stats as stats
from PIL import Image
from matplotlib.path import Path
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d import Axes3D
from functions import *

c=1                         #semiasse maggiore
a=0.5                       #semiasse minore
N=100000                     #numero di iterazioni

per_provax=np.zeros((N))
per_provaz=np.zeros((N))
per_provay=np.zeros((N))


data=[]


# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')

df=pd.read_csv(r'D:\matteo_citterio\Python Scripts\AR_per_image\AR_3_INT_Normal.csv')
num_nuclei=len(df.axes[0])
for i in range(1,64):

    for j in range(0,num_nuclei):

        if(not math.isnan(df[str(i)].iloc[j])):
            
            data.append(df[str(i)].iloc[j])

data=np.asarray(data)
diffs=[]
data_mean=np.average(data)
data_normalized=(data-1)/(data_mean-1)

bins=np.linspace(.999,4.5,200)
y,x=np.histogram(data_normalized,bins,density=True)
sigmas=[]
p_vals=[]

for l in range(0,1):

    
    
    for j in range(0,20):

        df1=2.4+0.1*j
        df2=df1
        sigmas.append(df1)

        ass=[]
        bs=[]
        cs=[]
        ars=[]
        aspect_ratios=[]

        for i in range(0,N):

            # alpha=np.random.uniform(0,2*np.pi)                  #3 angoli di eulero
            # u=np.uniform(0,1)
            # beta=np.arccos(1-2*u)
            # gamma=np.random.uniform(0,2*np.pi)

            print(l,':', j,':',i)

            a=1
            b=1+stats.chi2.rvs(df=df1)
            c=1+stats.chi2.rvs(df=df2)
            # a=1
            # b=1.2
            # c=2

            ass.append(a)
            bs.append(b)
            cs.append(c)
            ars.append(max(b,c)/min(b,c))

            n1=np.random.normal(0,2)
            n2=np.random.normal(0,2)
            n3=np.random.normal(0,2)

            magnitude=np.sqrt(n1**2+n2**2+n3**2)

            n1=n1/magnitude
            n2=n2/magnitude
            n3=n3/magnitude

            per_provax[i]=n1
            per_provay[i]=n2
            per_provaz[i]=n3

            eta=(n1**2)*((1/b**2)+(1/c**2))+(n2**2)*((1/a**2)+(1/c**2))+(n3**2)*((1/a**2)+(1/b**2))
            gamma=(n1**2)/((b**2)*(c**2))+(n2**2)/((a**2)*(c**2))+(n3**2)/((a**2)*(b**2))

            beta1=(eta/2)-np.sqrt(((eta/2)**2)-gamma)
            beta2=(eta/2)+np.sqrt(((eta/2)**2)-gamma)

            aspect_ratio=np.sqrt(beta2/beta1)
            aspect_ratios.append(aspect_ratio)

        aspect_ratios=np.asarray(aspect_ratios)
        aspect_ratios_mean=np.average(aspect_ratios)
        aspect_ratio_normalized=(aspect_ratios-1)/(aspect_ratios_mean-1)

        ys,xs=np.histogram(aspect_ratio_normalized,bins,density=True)
        diffs.append(np.cumsum(abs(y-ys))[-1])
        p_vals.append(stats.mannwhitneyu(data_normalized,aspect_ratio_normalized)[1])

    fig,ax=plt.subplots()
    ax.plot(sigmas,p_vals)
    ax.set_xlabel('gdl')
    ax.set_ylabel('p_val')
    

    fig,ax=plt.subplots()
    ax.plot(sigmas,diffs)
    ax.set_xlabel('gdl')
    ax.set_ylabel('scarto')
    plt.suptitle('gdl c:'+str(df2))

plt.show()

fig,ax=plt.subplots()
ax.plot(sigmas,diffs)
plt.show()
# ax.scatter(per_provax,per_provay,per_provaz)
# ax.set_title('Simulated Spherical versor')
# plt.show()

fig,ax=plt.subplots(1,3)
ax[0].hist(ass,100)
ax[0].set_title('a')
ax[1].hist(bs,100)
ax[1].set_title('b')
ax[2].hist(cs,100)
ax[2].set_title('c')
fig.suptitle('simulated input parameters')
# plt.show()

fig,ax=plt.subplots(1,2)

aspect_ratios=np.asarray(aspect_ratios)
aspect_ratios_mean=np.average(aspect_ratios)
aspect_ratio_normalized=(aspect_ratios-1)/(aspect_ratios_mean-1)

fig.suptitle('s1:'+str(sigma1)+' s2:'+str(sigma2))

ax[0].hist(data_normalized,120,density=True,label='real')
ax[0].hist(aspect_ratio_normalized,120,density=True,label='Simulation')
ax[0].legend()
ax[0].set_title('Simulated 2d ellipse aspect ratio')
ax[1].hist(ars,160,density=True)
ax[1].set_title('Input 3D aspect ratio')
# y,x=np.histogram(data,120,density=True)
# ys,xs=np.histogram(aspect_ratios,120,density=True)

plt.show()