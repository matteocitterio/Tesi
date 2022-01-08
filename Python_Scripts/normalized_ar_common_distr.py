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

from scipy.stats.stats import _attempt_exact_2kssamp
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
import scipy.stats as stats

colors=[]
columns=['9-norm','9-tum','13-tum','3-norm','3-tum','5-tum','7-tum','8-tum','12-norm','12-tum']
img_ns=[7,64,11,13,6,6,18,24,2,9]
real_img_n=[7, 18, 9, 64, 11, 13, 6, 6, 24, 2]
img_ns_string=['9_INT_Normal','9_INT_Tumor','13_INT_Tumor','3_INT_Normal','3_INT_Tumor','5_INT_Tumor','7_INT_Tumor','8_INT_Tumor','12_INT_Normal','12_INT_Tumor']
quant='AR'
divisors=[38,170,170,900,140,100,44,26,80,34]

ar_ratio=[]
gdl=[]
ar_means=[]

colors=['green','black','black','green','black','black','black','black','green','black']

fig,ax=plt.subplots()
fig2,ax2=plt.subplots()
fig3,ax3=plt.subplots()
cmap = plt.cm.get_cmap('Blues',10)
kolors = ['#0000FF','#195BCD','#73A3F6','#191980','#1B50BA','#135EF5']
ls=['-','--','-.',':']
outliers=[3000,1000,1250,600,1000,1500,1000,3500,1000,2600]
for eta in range(0,len(img_ns_string)): 

    print(img_ns_string[eta])

    df=pd.read_csv(r'D:\matteo_citterio\Python Scripts\Area_ratio_per_image\Area_ratio_'+img_ns_string[eta]+'.csv')
    dj=pd.read_csv(r'D:\matteo_citterio\Python Scripts\Area_ratio_per_image\Area_ratio_'+img_ns_string[eta]+'.csv')

    num_nuclei=len(df.axes[0])

    data=[]
    area_ratio=[]

    for i in range(1,real_img_n[eta]):                  #ciclo a colonna che vuol dire a immagine

        for j in range(0,num_nuclei):

            if(not math.isnan(df[str(i)].iloc[j])):

                data.append(df[str(i)].iloc[j])
                area_ratio.append(dj[str(i)].iloc[j])

    print(len(data))

    tempo=[]
    for i in range(0,len(data)):

        if(data[i]<outliers[eta]):
            tempo.append(data[i])

    data=tempo

    ar_means.append(np.average(data))
    # data=np.asarray(data)
    # data=(data-1)/(np.average(data)-1)

    y,x=np.histogram(data,int(len(data)/divisors[eta]),density=True)
    hist_dist = stats.rv_histogram(np.histogram(data,int(len(data)/divisors[eta]),density=True))

    xnew=np.zeros((x.size+1))
    step=x[1]-x[0]

    for j in range(1,x.size+1):

        xnew[j]=x[j-1]

    x=xnew

    params=stats.lognorm.fit(data, loc=0)

    if (img_ns_string[eta]=='9_INT_Normal'):

        params2=stats.chi2.fit(data)

    else:

        params2=stats.chi2.fit(data, loc=0)

    print(params)
    #print(params2[0])
    gdl.append(params2[0])

    pdf_fitted =stats.lognorm.pdf(x, *params)
    pdf_fitted2 =stats.chi2.pdf(x, *params2)

    if (eta==0):
        ax.plot(x,hist_dist.pdf(x),'*',color='#CC0000', ms=8,mew=1, label='Patient 9 Normal tissue',lw=3)
        # ax.plot(x, pdf_fitted,ls=ls[eta], color='red',label='Patient 9 Normal tissue ',lw=3)
    if (eta==1):
        ax.plot(x,hist_dist.pdf(x),'o',color='#606060', ms=8,mew=1, label='Patient 9 Tumor tissue',lw=3)
        # ax.plot(x, pdf_fitted,ls=ls[eta], color='red',label='Patient 9 Tumor tissue ',lw=3)
    if (eta==2):
        print('pass')
        #ax.plot(x,hist_dist.pdf(x),'o',color='#606060', ms=8,mew=1, label='Patient 13 Tumor tissue',lw=3)
        # ax.plot(x, pdf_fitted,ls=ls[eta], color='red',label='Patient 13 Tumor tissue ',lw=3) 
    if (eta==3):
        ax2.plot(x,hist_dist.pdf(x),'*',color='#CC0000', ms=8,mew=1, label='Patient 3 Normal tissue',lw=3)
        # ax.plot(x, pdf_fitted,ls=ls[eta], color='red',label='Patient 3 Normal tissue ',lw=3)  
    if (eta==4):
        ax2.plot(x,hist_dist.pdf(x),'o',color='#606060', ms=8,mew=1, label='Patient 3 Tumor tissue',lw=3)
        # ax.plot(x, pdf_fitted,ls=ls[eta], color='red',label='Patient 3 Normal tissue ',lw=3)  
    if(eta==8):    
        ax3.plot(x,hist_dist.pdf(x),'*',color='#CC0000', ms=8,mew=1, label='Patient 12 Normal tissue',lw=3)
    if(eta==9):    
        ax3.plot(x,hist_dist.pdf(x),'o',color='#606060', ms=8,mew=1, label='Patient 12 Tumor tissue',lw=3)
    #ax.plot(x, pdf_fitted2, color='k',label='chi2 fit',ls='--',lw=3)

    # ax.hist(data,int(len(data)/divisors[eta]),density=True,label='data')

    # ax.set_xlabel('$\left(AR-1\\right)/\left( \overline{AR}-1\\right)$')
    ax.set_xlabel('Area $[\mu m ^{2}]$')
    ax.set_ylabel('Normalized counts')
    ax.set_xlim(-1,130)
    ax2.set_xlabel('Area $[\mu m ^{2}]$')
    ax2.set_ylabel('Normalized counts')
    ax2.set_xlim(-1,100)
    ax3.set_xlabel('Area $[\mu m ^{2}]$')
    ax3.set_ylabel('Normalized counts')
    ax3.set_xlim(-1,100)
    #ax.set_title('Patient:'+img_ns_string[eta]+' AR distribution')
    # ax.set_title('Normalized Voronoi AR distributions')
    ax.set_title('Patient 9')
    ax2.set_title('Patient 3')
    ax3.set_title('Patient 12')



    ax.legend()
    ax2.legend()
    ax3.legend()
    # fig.savefig('D:\matteo_citterio\plot_al_fly\distribuzioni\Distribuzioni a paziente\\AR\\patient_'+img_ns_string[eta]+'_AR_distribution.png')

    y_ar,x_ar=np.histogram(data,int(len(data)/divisors[eta]))
    sampled=stats.chi2.rvs(*params,size=len(data))
    sampled2=stats.chi2.rvs(*params2,size=len(data))
    y_sampled,x_sampled=np.histogram(sampled,int(len(data)/divisors[eta]))
    y_sampled2,x_sampled2=np.histogram(sampled2,int(len(data)/divisors[eta]))


    


legend_elements =[Line2D([0], [0], marker='o', color='w', label='Tumor',markerfacecolor='black',markeredgecolor='black',markersize=5),
                    Line2D([0], [0], marker='o', color='w', label='Normal',markerfacecolor='green',markeredgecolor='green', markersize=5)]

h=0
# an.legend(handles=legend_elements)

# an.plot(gdl,ar_ratio,'o',color=colors)
plt.show()

