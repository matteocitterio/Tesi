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
from matplotlib.path import Path
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d import Axes3D
from functions import *
import scipy.stats as stats

N=100000                     #numero di iterazioni

per_provax=np.zeros((N))
per_provaz=np.zeros((N))
per_provay=np.zeros((N))


data=[]


#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')

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

bins=np.linspace(-0.1,4.5,150)
y,x=np.histogram(data_normalized,bins,density=True)
sigmas=[]

sigma1=.3+10*0.01
sigmas.append(sigma1)

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

    print(i)

    sigma2=.4

    #a=1
    #b=1+(np.random.normal(1,sigma1)**2)
    #c=1+(np.random.normal(1,sigma2)**2)
    a=stats.lognorm.rvs(0.06)
    b=stats.lognorm.rvs(0.06)
    c=stats.lognorm.rvs(0.06)

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


#ax.scatter(per_provax,per_provay,per_provaz,color='tab:blue')
#ax.set_title('Simulated $\hat{n}$ unit vector')
#ax.set_xlabel('x')
#ax.set_ylabel('y')
#ax.set_zlabel('z')
#plt.show()

#fig,ax=plt.subplots()
#ax.plot(sigmas,diffs)
#plt.show()


fig,ax=plt.subplots(1,3)
bins_ass=np.linspace(0.7,1.3,150)
bins_cs=np.linspace(0.7,1.3,150)
y,x=np.histogram(ass,bins_ass,density=True)
y,xs=np.histogram(cs,bins_cs,density=True)
hist_dist_a = stats.rv_histogram(np.histogram(ass,bins_ass,density=True))
hist_dist_b = stats.rv_histogram(np.histogram(bs,bins_ass,density=True))
hist_dist_c = stats.rv_histogram(np.histogram(cs,bins_cs,density=True))
xnew=np.zeros((x.size+1))+0.7
step=x[1]-x[0]
for j in range(1,x.size+1):

    xnew[j]=x[j-1]
x=xnew
xnew=np.zeros((xs.size+1))+0.7
step=xs[1]-xs[0]
for j in range(1,xs.size+1):

    xnew[j]=xs[j-1]
xs=xnew
ax[0].plot(x,hist_dist_a.pdf(x),'-o',color='tab:purple', ms=8,mew=1,lw=3)
ax[1].plot(x,hist_dist_b.pdf(x),'-o',color='tab:purple', ms=8,mew=1,lw=3)
ax[2].plot(xs,hist_dist_c.pdf(xs),'-o',color='tab:purple', ms=8,mew=1,lw=3)
ax[0].set_title('a',size=13)
ax[1].set_title('b',size=13)
ax[2].set_title('c',size=13)
ax[0].set_ylabel('Semiaxes PDF',size=13)

# aspect_ratios=np.asarray(aspect_ratios)
# aspect_ratios_mean=np.average(aspect_ratios)
# aspect_ratio_normalized=(aspect_ratios-1)/(aspect_ratios_mean-1)

# fig.suptitle('s1:'+str(sigma1)+' s2:'+str(sigma2))

ys,x=np.histogram(aspect_ratio_normalized,bins,density=True)
diffs.append(np.cumsum(abs(y-ys))[-1])
hist_dist = stats.rv_histogram(np.histogram(aspect_ratio_normalized,bins,density=True))
hist_dist_data=stats.rv_histogram(np.histogram(data_normalized,bins,density=True))
xnew=np.zeros((x.size+1))-.1
step=x[1]-x[0]

for j in range(1,x.size+1):

    xnew[j]=x[j-1]

x=xnew

params=stats.gamma.fit(aspect_ratio_normalized)
pdf_fitted=stats.gamma.pdf(x, *params)

fig,ax=plt.subplots(1,2)

# ax[0].plot(x,hist_dist_data.pdf(x),'-o',color='tab:blue', ms=8,mew=1, label='Data PDF',lw=3)

ax[0].plot(x,hist_dist.pdf(x),'-*',color='tab:orange', ms=8,mew=1, label='Simulation PDF',lw=3)
ax[0].plot(x,pdf_fitted,ls='-',lw=3,color='red', label='k-Gamma, k=2')
# ax[0].hist(data_normalized,100,color='tab:blue',density=True,label='Patient 3 - Normal tissue')
# ax[0].hist(aspect_ratio_normalized,100,color='red',alpha=0.5,density=True,label='Simulation')
ax[0].set_xlabel('$(AR-1)/(\overline{AR}-1)$',fontsize=13)
ax[0].set_ylabel('PDF',fontsize=13)
ax[0].legend(fontsize=13)
ax[0].set_title('2D distribution',fontsize=13)


# ax[1].hist(ars,160,density=True,color='tab:blue')
ax[1].set_xlabel('$AR$',fontsize=13)
ax[1].set_title('3D distribution',fontsize=13)
print(stats.mannwhitneyu(data_normalized,aspect_ratio_normalized))
# y,x=np.histogram(data,120,density=True)
# ys,xs=np.histogram(aspect_ratios,120,density=True)

bisn=np.linspace(0.95,1.4,150)
y,x=np.histogram(ars,bisn,density=True)
hist_dist_3d=stats.rv_histogram(np.histogram(ars,bisn,density=True))
xnew=np.zeros((x.size+1))+0.95
step=x[1]-x[0]
for j in range(1,x.size+1):

    xnew[j]=x[j-1]
x=xnew
ax[1].plot(x,hist_dist_3d.pdf(x),'-o',color='tab:purple', ms=8,mew=1, label='3D PDF',lw=3)
ax[1].legend(fontsize=13,loc='upper left')

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

print('chisquared_lognorm:',temp)

plt.show()