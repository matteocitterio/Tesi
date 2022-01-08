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
from functions import *

colors=[]
columns=['9-norm','9-tum','13-tum','3-norm','3-tum','5-tum','7-tum','8-tum','12-norm','12-tum']
img_ns=[7,64,11,13,6,6,18,24,2,9]
real_img_n=[7, 18, 9, 64, 11, 13, 6, 6, 24, 2]
img_ns_string=['9_INT_Normal','9_INT_Tumor','13_INT_Tumor','3_INT_Normal','3_INT_Tumor','5_INT_Tumor','7_INT_Tumor','8_INT_Tumor','12_INT_Normal','12_INT_Tumor']
img_ns_tumor_string=['9_INT_Tumor','13_INT_Tumor','3_INT_Tumor','5_INT_Tumor','7_INT_Tumor','8_INT_Tumor','12_INT_Tumor']
real_img_n_tumor=[ 18, 9, 11, 13, 6, 6, 2]
real_img_n_normal=[7,64,24]
img_ns_normal_string=['9_INT_Normal','3_INT_Normal','12_INT_Normal']
quant='AR'
divisors=[60,200,170,1000,140,100,44,26,135,28]

divisor=570

fig2,an=plt.subplots()
ar_ratio=[]
gdl=[]

colors=['green','black','black','green','black','black','black','black','green','black']

tumor=[]
area_ratio=[]

for eta in range(0,len(img_ns_tumor_string)): 

    print(img_ns_tumor_string[eta])

    df=pd.read_csv(r'D:\matteo_citterio\Python Scripts\AR_per_image\AR_'+img_ns_tumor_string[eta]+'.csv')

    num_nuclei=len(df.axes[0])

    for i in range(1,real_img_n_tumor[eta]):                  #ciclo a colonna che vuol dire a immagine

        for j in range(0,num_nuclei):

            if(not math.isnan(df[str(i)].iloc[j])):

                tumor.append(df[str(i)].iloc[j])

normal=[]

for eta in range(0,len(img_ns_normal_string)): 

    print(img_ns_normal_string[eta])

    df=pd.read_csv(r'D:\matteo_citterio\Python Scripts\AR_per_image\AR_'+img_ns_normal_string[eta]+'.csv')

    num_nuclei=len(df.axes[0])

    for i in range(1,real_img_n_normal[eta]):                  #ciclo a colonna che vuol dire a immagine

        for j in range(0,num_nuclei):

            if(not math.isnan(df[str(i)].iloc[j])):

                normal.append(df[str(i)].iloc[j])

print(len(tumor))
print(len(normal))

################################################################################################################################################

y,x=np.histogram(tumor,int(len(tumor)/divisor),density=True)

print(divisors[eta])

hist_dist = stats.rv_histogram(np.histogram(tumor,int(len(tumor)/divisor),density=True))

xnew=np.zeros((x.size+1))+1
step=x[1]-x[0]

for j in range(1,x.size+1):

    xnew[j]=x[j-1]

x=xnew

fig,ax=plt.subplots()
fig3,bx=plt.subplots()

params=stats.gamma.fit(tumor, loc=1)

if (img_ns_string[eta]=='9_INT_Normal'):

    params2=stats.chi2.fit(tumor)

else:

    params2=stats.chi2.fit(tumor, loc=1)

# print(params)
print(params2)
gdl.append(params2[0])

pdf_fitted =stats.gamma.pdf(x, *params)
pdf_fitted2 =stats.chi2.pdf(x, *params2)
ax.plot(x,hist_dist.pdf(x),'-o', label='data PDF',lw=3)
ax.plot(x, pdf_fitted, color='r',label='gamma fit',lw=3)
ax.plot(x, pdf_fitted2, color='k',label='chi2 fit',ls='--',lw=3)

bx.plot(x,hist_dist.pdf(x),'-o', label='data PDF',lw=3)
bx.plot(x, pdf_fitted, color='r',label='gamma fit',lw=3)
bx.plot(x, pdf_fitted2, color='k',label='chi2 fit',ls='--',lw=3)
bx.set_yscale('log')
    # bx.set_xscale('log')
bx.legend()

    # ax.hist(data,int(len(data)/divisors[eta]),density=True,label='data')

ax.set_xlabel('AR')
ax.set_ylabel('Normalized counts')
ax.set_title('Tumor AR distribution')

bx.set_xlabel('AR')
bx.set_ylabel('Normalized counts')
bx.set_title('Tumor AR distribution LOG')

ax.legend()
    # fig.savefig('D:\matteo_citterio\plot_al_fly\distribuzioni\Distribuzioni a paziente\\AR\\patient_'+img_ns_string[eta]+'_AR_distribution.png')

####################################################################################################################

divisor_normal=1300

y,x=np.histogram(normal,int(len(normal)/divisor_normal),density=True)

hist_dist = stats.rv_histogram(np.histogram(normal,int(len(normal)/divisor_normal),density=True))

xnew=np.zeros((x.size+1))+1
step=x[1]-x[0]

for j in range(1,x.size+1):

    xnew[j]=x[j-1]

x=xnew

fig,ax=plt.subplots()
fig3,bx=plt.subplots()

params=stats.gamma.fit(normal, loc=1)

if (img_ns_string[eta]=='9_INT_Normal'):

    params2=stats.chi2.fit(normal)

else:

    params2=stats.chi2.fit(normal, loc=1)

# print(params)
print(params)
gdl.append(params2[0])

pdf_fitted =stats.gamma.pdf(x, *params)
pdf_fitted2 =stats.chi2.pdf(x, *params2)
ax.plot(x,hist_dist.pdf(x),'-o', label='data PDF',lw=3)
ax.plot(x, pdf_fitted, color='r',label='gamma fit',lw=3)
ax.plot(x, pdf_fitted2, color='k',label='chi2 fit',ls='--',lw=3)

bx.plot(x,hist_dist.pdf(x),'-o', label='data PDF',lw=3)
bx.plot(x, pdf_fitted, color='r',label='gamma fit',lw=3)
bx.plot(x, pdf_fitted2, color='k',label='chi2 fit',ls='--',lw=3)
bx.set_yscale('log')
    # bx.set_xscale('log')
bx.legend()

    # ax.hist(data,int(len(data)/divisors[eta]),density=True,label='data')

ax.set_xlabel('AR')
ax.set_ylabel('Normalized counts')
ax.set_title('Normal AR distribution')

bx.set_xlabel('AR')
bx.set_ylabel('Normalized counts')
bx.set_title('Normal AR distribution LOG')

ax.legend()
    # fig.savefig('D:\matteo_citterio\plot_al_fly\distribuzioni\Distribuzioni a paziente\\AR\\patient_'+img_ns_string[eta]+'_AR_distribution.png')

tumor_array=np.asarray(tumor)
normal_array=np.asarray(normal)
turmor_normalized=(tumor_array-1)/(np.average(tumor_array)-1)
normal_normalized=(normal_array-1)/(np.average(normal_array)-1)

fig,ax=plt.subplots()
ax.hist(turmor_normalized,160,density=True,label='tumor')
ax.hist(normal_normalized,160,density=True,label='normal',alpha=0.5)
ax.set_xlabel('(AR-1)/(<AR>-1)')
ax.set_ylabel('PDF')
params=stats.chi2.fit(normal_normalized,loc=0)
params2=stats.chi2.fit(turmor_normalized,loc=0)

print('tum',params2[0])
print('norm',params[0])

ax.legend()

print(stats.mannwhitneyu(turmor_normalized,normal_normalized))
print(stats.ks_2samp(turmor_normalized,normal_normalized))


plt.show()

