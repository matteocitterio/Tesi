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
patient_string=['9-Normal','9-Tumor','13-Tumor','3-Normal','3-Tumor','5-Tumor','7-Tumor','8-Tumor','12-Normal','12-Tumor']
divisors=[1,1,1,1,1,1,1,1,1,1]

ar_ratio=[]
gdl=[]
ar_means=[]

colors=['green','black','black','green','black','black','black','black','green','black']

fig,ax=plt.subplots()
fig2,ax2=plt.subplots()
cmap = plt.cm.get_cmap('Blues',10)
cm=truncate_colormap(plt.cm.get_cmap('Blues'), 0.4, 1)
color_gradients = cm(np.arange(0,10)/10) 
kolors = ['#0000FF','#195BCD','#73A3F6','#191980','#1B50BA','#135EF5']
ls=['-','--','-.',':']
number=[150,100,250,250,150,150,150,120,120,100]
outliers=[3000,1000,1250,600,1000,1500,1000,3500,1000,2600]
for eta in range(0,len(img_ns_string)): 

    print(img_ns_string[eta])

    df=pd.read_csv(r'D:\matteo_citterio\Python Scripts\AR_per_image\AR_'+img_ns_string[eta]+'.csv')
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

        if(data[i]<8):
            tempo.append(data[i])

    data=tempo

    ar_means.append(np.average(data))
    data=np.asarray(data)                   #normalizzazione
    data=(data-1)/(np.average(data)-1)

    y,x=np.histogram(data,int(len(data)/divisors[eta]),density=True)
    hist_dist = stats.rv_histogram(np.histogram(data,int(len(data)/divisors[eta]),density=True))

    xnew=np.zeros((x.size+1))
    step=x[1]-x[0]

    for j in range(1,x.size+1):

        xnew[j]=x[j-1]

    x=xnew

    params=stats.gamma.fit(data, loc=1)

    if (img_ns_string[eta]=='9_INT_Normal'):

        params2=stats.chi2.fit(data)

    else:

        params2=stats.chi2.fit(data, loc=1)

    print(params)
    #print(params2[0])
    gdl.append(params2[0])

    pdf_fitted =stats.gamma.pdf(x, *params)
    pdf_fitted2 =stats.chi2.pdf(x, *params2)

    ax.plot(x, pdf_fitted,color=color_gradients[eta], label='Patient '+patient_string[eta],lw=3)

    # ax.hist(data,int(len(data)/divisors[eta]),density=True,label='data')

    # ax.set_xlabel('$\left(AR-1\\right)/\left( \overline{AR}-1\\right)$')
    ax.set_xlabel('$AR$')
    ax.set_ylabel('Gamma-fitted PDF')
    ax.set_xlim(-0.1,4.6)
    #ax.set_title('Patient:'+img_ns_string[eta]+' AR distribution')
    # ax.set_title('Normalized Voronoi AR distributions')
    ax.set_title('Nuclei $AR$ Gamma fit')
    ax.legend()
    # fig.savefig('D:\matteo_citterio\plot_al_fly\distribuzioni\Distribuzioni a paziente\\AR\\patient_'+img_ns_string[eta]+'_AR_distribution.png')

    y_ar,x_ar=np.histogram(data,int(len(data)/divisors[eta]))
    sampled=stats.chi2.rvs(*params,size=len(data))
    sampled2=stats.chi2.rvs(*params2,size=len(data))
    y_sampled,x_sampled=np.histogram(sampled,int(len(data)/divisors[eta]))
    y_sampled2,x_sampled2=np.histogram(sampled2,int(len(data)/divisors[eta]))



    bins_cs=np.linspace(-0.1,4.5,number[eta])
    y,xs=np.histogram(data,bins_cs,density=True)
    hist_dist_a = stats.rv_histogram(np.histogram(data,bins_cs,density=True))

    ax2.plot(x,hist_dist_a.pdf(x),'o',color=color_gradients[eta],label='Patient '+patient_string[eta], ms=4,mew=1,lw=3)
    ax2.set_xlabel('$AR-1/\overline{AR}-1$')
    ax2.set_ylabel('PDF')
    ax2.set_xlim(-0.1,4.5)
    #ax.set_title('Patient:'+img_ns_string[eta]+' AR distribution')
    # ax.set_title('Normalized Voronoi AR distributions')
    ax2.set_title('Nuclei $AR$ distribution')
    ax2.legend()


legend_elements =[Line2D([0], [0], marker='o', color='w', label='Tumor',markerfacecolor='black',markeredgecolor='black',markersize=5),
                    Line2D([0], [0], marker='o', color='w', label='Normal',markerfacecolor='green',markeredgecolor='green', markersize=5)]

h=0
# an.legend(handles=legend_elements)

# an.plot(gdl,ar_ratio,'o',color=colors)
plt.show()

