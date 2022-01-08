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
N=5000                     #numero di iterazioni

per_provax=np.zeros((N))
per_provaz=np.zeros((N))
per_provay=np.zeros((N))


data=[]
ass=[]
bs=[]
cs=[]
ars=[]
aspect_ratios=[]

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')

df=pd.read_csv(r'D:\matteo_citterio\Python Scripts\AR_per_image\AR_13_INT_Tumor.csv')
num_nuclei=len(df.axes[0])
for i in range(1,8):

    for j in range(0,num_nuclei):

        if(not math.isnan(df[str(i)].iloc[j])):
            
            data.append(df[str(i)].iloc[j])

data=np.asarray(data)

data_mean=np.average(data)
data_normalized=(data-1)/(data_mean-1)


for i in range(0,N):



    print(i)

    sigma1=.48
    sigma2=.48

    df=3.1

    #a=1+stats.chi2.rvs(df)
    #b=1+stats.chi2.rvs(df)
    #c=1+stats.chi2.rvs(df)

    s=0.2

    # s1=np.random.normal(0,2)
    # s2=np.random.normal(0,2)
    # s3=np.random.normal(0,2)

    # magnitude_s=np.sqrt(s1**2+s2**2+s3**2)

    # s1=s1/magnitude_s
    # s2=s2/magnitude_s
    # s3=s3/magnitude_s

    # strain=stats.lognorm.rvs(s)

    # a=strain*s1
    # b=strain*s2
    # c=strain*s3

    # a=stats.lognorm.rvs(s)
    # b=stats.lognorm.rvs(s)
    # c=stats.lognorm.rvs(s)
    a=np.random.normal(loc=1,scale=s)
    b=np.random.normal(loc=1,scale=s)
    c=np.random.normal(loc=1,scale=s)


    ass.append(a)
    bs.append(b)
    cs.append(c)
    ars.append(max(a,b,c)/min(a,b,c))

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

fig,ax=plt.subplots(1,3)
ax[0].hist(ass,100)
ax[0].set_title('a')
ax[1].hist(bs,100)
ax[1].set_title('b')
ax[2].hist(cs,100)
ax[2].set_title('c')
fig.suptitle('simulated input parameters')

fig,ax=plt.subplots(1,2)

aspect_ratios=np.asarray(aspect_ratios)
print(min(aspect_ratios))
aspect_ratios_mean=np.nanmean(aspect_ratios)
print(aspect_ratios_mean)
aspect_ratio_normalized=(aspect_ratios-1)/(aspect_ratios_mean-1)

ars=np.asarray(ars)
plt.suptitle('lognormal 3D input df: '+str(s))

num_bins=300
bins=np.linspace(0,5,num_bins)

ax[0].hist(aspect_ratio_normalized,bins,density=True,label='Simulation',alpha=0.5)
y_ar,x_ar=np.histogram(aspect_ratio_normalized,bins)


ax[0].set_xlabel('(AR-1)/(<AR>-1)')
ax[0].set_ylabel('PDF')

ax[0].set_title('Simulated 2d ellipse aspect ratio')

ars_normalized=(ars-1)/(np.average(ars)-1)


ax[1].hist(ars_normalized,160,density=True)
ax[1].set_title('Input 3D aspect ratio')
ax[1].set_xlabel('(AR-1)/(<AR>-1)')

params=stats.chi2.fit(aspect_ratio_normalized,floc=0)
x=np.linspace(-0.01,5,300)
pdf_fitted=stats.chi2.pdf(x,*params)
sampled=stats.chi2.rvs(*params,size=N)
# ax[0].hist(sampled,bins,density=True,label='sampled',alpha=0.5)
y_sampled,x_sampled=np.histogram(sampled,bins)

ax[0].legend()

ax[0].plot(x,pdf_fitted,color='red')
print('params',params)
print('params_no_loc',params)
print('params',stats.ks_2samp(y_ar,y_sampled))

# print(stats.chisquare(f_obs=y_ar,f_exp=y_sampled,ddof=3))
# print(np.cumsum(y_ar)[-1],np.cumsum(pdf_fitted)[-1])
plt.show()