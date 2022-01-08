import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import warnings
import os
import sys
import scipy.spatial as sp
import time

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning) 
from pandas.core.common import SettingWithCopyWarning
from shapely.geometry import MultiPoint, Point, Polygon, LineString,MultiPolygon
from shapely.ops import polygonize,unary_union
from scipy.spatial import Voronoi, voronoi_plot_2d,ConvexHull, convex_hull_plot_2d
from scipy.fft import fft, ifft
from PIL import Image
from matplotlib.path import Path
from functions import *

n_images_9normal=7                 #in realtà il numero di immagini viste è n_images-1
n_images_3normal=64
n_images_3tumor=11
n_images_5tumor=13
n_images_7tumor=6
n_images_8tumor=6
n_images_9tumor=18
n_images_12normal=24
n_images_12tumor=2
n_images_13tumor=9

columns=['9-norm','9-tum','13-tum','3-norm','3-tum','5-tum','7-tum','8-tum','12-norm','12-tum']
indexis=np.arange(63)+1
indexis.tolist()
output_frame=pd.DataFrame(index=indexis,columns=columns)

for i in range(1,n_images_9normal):

    print('9-norm: ',i)

    image_N=str(i)
    df=pd.read_csv(r'D:\matteo_citterio\risultati_segment\ROI_9_INT_Normal\Results_'+image_N+'.csv')

    num_nuclei=len(df.axes[0])                   #dimensione dei dati generati  
    funct=create_nuclei_function(num_nuclei,df)			    #crea un array con le pos dei cdm dei nuclei dalla df di stardist
    df.to_csv(r'D:\matteo_citterio\risultati_segment\ROI_3_INT_Normal\Results.txt', header=None, index=None, sep=' ',columns=['XM','YM']  ,mode='a')

    density_matlab()

    dj=pd.read_csv(r'D:\matteo_citterio\Python Scripts\temp2.csv')
    
    os.remove(r'D:\matteo_citterio\risultati_segment\ROI_3_INT_Normal\Results.txt')
    os.remove(r'D:\matteo_citterio\Python Scripts\temp2.csv')
    dj.columns =['g(R)', 'R','c']
    bins_density=dj['R']
    density_profile=dj['g(R)']

    df1=pd.DataFrame(density_profile)
    k1=20
    df1=df1.rolling(k1+3).mean()

    size=len(dj.axes[0])
    x=np.zeros((size))
    y=np.zeros((size))

    for j in range(0,size):
        x[j]=bins_density[j]
        y[j]=df1.iloc[j]

    shift=np.count_nonzero(np.isnan(y))

    y = y[~np.isnan(y)]
    maxi=np.argsort(y)

    output_frame['9-norm'].iloc[i-1]=x[maxi[-1]+shift]    

for i in range(1,n_images_9tumor):

    print('9-tumor: ',i)

    image_N=str(i)
    df=pd.read_csv(r'D:\matteo_citterio\risultati_segment\ROI_9_INT_Tumor\Results_'+image_N+'.csv')

    num_nuclei=len(df.axes[0])                   #dimensione dei dati generati  
    funct=create_nuclei_function(num_nuclei,df)			    #crea un array con le pos dei cdm dei nuclei dalla df di stardist
    df.to_csv(r'D:\matteo_citterio\risultati_segment\ROI_3_INT_Normal\Results.txt', header=None, index=None, sep=' ',columns=['XM','YM']  ,mode='a')

    density_matlab()

    dj=pd.read_csv(r'D:\matteo_citterio\Python Scripts\temp2.csv')
    os.remove(r'D:\matteo_citterio\risultati_segment\ROI_3_INT_Normal\Results.txt')
    os.remove(r'D:\matteo_citterio\Python Scripts\temp2.csv')
    dj.columns =['g(R)', 'R','c']
    bins_density=dj['R']
    density_profile=dj['g(R)']

    df1=pd.DataFrame(density_profile)
    k1=20
    df1=df1.rolling(k1+3).mean()

    size=len(dj.axes[0])
    x=np.zeros((size))
    y=np.zeros((size))

    for j in range(0,size):
        x[j]=bins_density[j]
        y[j]=df1.iloc[j]

    shift=np.count_nonzero(np.isnan(y))

    y = y[~np.isnan(y)]
    maxi=np.argsort(y)

    output_frame['9-tum'].iloc[i-1]=x[maxi[-1]+shift]    

for i in range(1,n_images_13tumor):

    print('13-tum: ',i)

    image_N=str(i)
    df=pd.read_csv(r'D:\matteo_citterio\risultati_segment\ROI_13_INT_Tumor\Results_'+image_N+'.csv')

    num_nuclei=len(df.axes[0])                   #dimensione dei dati generati  
    funct=create_nuclei_function(num_nuclei,df)			    #crea un array con le pos dei cdm dei nuclei dalla df di stardist
    df.to_csv(r'D:\matteo_citterio\risultati_segment\ROI_3_INT_Normal\Results.txt', header=None, index=None, sep=' ',columns=['XM','YM']  ,mode='a')

    density_matlab()

    dj=pd.read_csv(r'D:\matteo_citterio\Python Scripts\temp2.csv')
    os.remove(r'D:\matteo_citterio\risultati_segment\ROI_3_INT_Normal\Results.txt')
    os.remove(r'D:\matteo_citterio\Python Scripts\temp2.csv')
    dj.columns =['g(R)', 'R','c']
    bins_density=dj['R']
    density_profile=dj['g(R)']

    df1=pd.DataFrame(density_profile)
    k1=20
    df1=df1.rolling(k1+3).mean()

    size=len(dj.axes[0])
    x=np.zeros((size))
    y=np.zeros((size))

    for j in range(0,size):
        x[j]=bins_density[j]
        y[j]=df1.iloc[j]

    shift=np.count_nonzero(np.isnan(y))

    y = y[~np.isnan(y)]
    maxi=np.argsort(y)

    output_frame['13-tum'].iloc[i-1]=x[maxi[-1]+shift]    

for i in range(1,n_images_3normal):

    print('3-norm: ',i)

    image_N=str(i)
    df=pd.read_csv(r'D:\matteo_citterio\risultati_segment\ROI_3_INT_Normal\Results_'+image_N+'.csv')

    num_nuclei=len(df.axes[0])                   #dimensione dei dati generati  
    funct=create_nuclei_function(num_nuclei,df)			    #crea un array con le pos dei cdm dei nuclei dalla df di stardist
    df.to_csv(r'D:\matteo_citterio\risultati_segment\ROI_3_INT_Normal\Results.txt', header=None, index=None, sep=' ',columns=['XM','YM']  ,mode='a')

    density_matlab()

    dj=pd.read_csv(r'D:\matteo_citterio\Python Scripts\temp2.csv')
    os.remove(r'D:\matteo_citterio\risultati_segment\ROI_3_INT_Normal\Results.txt')
    os.remove(r'D:\matteo_citterio\Python Scripts\temp2.csv')
    dj.columns =['g(R)', 'R','c']
    bins_density=dj['R']
    density_profile=dj['g(R)']

    df1=pd.DataFrame(density_profile)
    k1=20
    df1=df1.rolling(k1+3).mean()

    size=len(dj.axes[0])
    x=np.zeros((size))
    y=np.zeros((size))

    for j in range(0,size):
        x[j]=bins_density[j]
        y[j]=df1.iloc[j]

    shift=np.count_nonzero(np.isnan(y))

    y = y[~np.isnan(y)]
    maxi=np.argsort(y)

    output_frame['3-norm'].iloc[i-1]=x[maxi[-1]+shift]    

for i in range(1,n_images_3tumor):

    print('3-tum: ',i)

    image_N=str(i)
    df=pd.read_csv(r'D:\matteo_citterio\risultati_segment\ROI_3_INT_Tumor\Results_'+image_N+'.csv')

    num_nuclei=len(df.axes[0])                   #dimensione dei dati generati  
    funct=create_nuclei_function(num_nuclei,df)			    #crea un array con le pos dei cdm dei nuclei dalla df di stardist
    df.to_csv(r'D:\matteo_citterio\risultati_segment\ROI_3_INT_Normal\Results.txt', header=None, index=None, sep=' ',columns=['XM','YM']  ,mode='a')

    density_matlab()

    dj=pd.read_csv(r'D:\matteo_citterio\Python Scripts\temp2.csv')
    os.remove(r'D:\matteo_citterio\risultati_segment\ROI_3_INT_Normal\Results.txt')
    os.remove(r'D:\matteo_citterio\Python Scripts\temp2.csv')
    dj.columns =['g(R)', 'R','c']
    bins_density=dj['R']
    density_profile=dj['g(R)']

    df1=pd.DataFrame(density_profile)
    k1=20
    df1=df1.rolling(k1+3).mean()

    size=len(dj.axes[0])
    x=np.zeros((size))
    y=np.zeros((size))

    for j in range(0,size):
        x[j]=bins_density[j]
        y[j]=df1.iloc[j]

    shift=np.count_nonzero(np.isnan(y))

    y = y[~np.isnan(y)]
    maxi=np.argsort(y)

    output_frame['3-tum'].iloc[i-1]=x[maxi[-1]+shift]    

for i in range(1,n_images_5tumor):

    print('5-tum: ',i)

    image_N=str(i)
    df=pd.read_csv(r'D:\matteo_citterio\risultati_segment\ROI_5_INT_Tumor\Results_'+image_N+'.csv')

    num_nuclei=len(df.axes[0])                   #dimensione dei dati generati  
    funct=create_nuclei_function(num_nuclei,df)			    #crea un array con le pos dei cdm dei nuclei dalla df di stardist
    df.to_csv(r'D:\matteo_citterio\risultati_segment\ROI_3_INT_Normal\Results.txt', header=None, index=None, sep=' ',columns=['XM','YM']  ,mode='a')

    density_matlab()

    dj=pd.read_csv(r'D:\matteo_citterio\Python Scripts\temp2.csv')
    os.remove(r'D:\matteo_citterio\risultati_segment\ROI_3_INT_Normal\Results.txt')
    os.remove(r'D:\matteo_citterio\Python Scripts\temp2.csv')
    dj.columns =['g(R)', 'R','c']
    bins_density=dj['R']
    density_profile=dj['g(R)']

    df1=pd.DataFrame(density_profile)
    k1=20
    df1=df1.rolling(k1+3).mean()

    size=len(dj.axes[0])
    x=np.zeros((size))
    y=np.zeros((size))

    for j in range(0,size):
        x[j]=bins_density[j]
        y[j]=df1.iloc[j]

    shift=np.count_nonzero(np.isnan(y))

    y = y[~np.isnan(y)]
    maxi=np.argsort(y)

    output_frame['5-tum'].iloc[i-1]=x[maxi[-1]+shift]    

for i in range(1,n_images_7tumor):

    print('7-tum: ',i)

    image_N=str(i)
    df=pd.read_csv(r'D:\matteo_citterio\risultati_segment\ROI_7_INT_Tumor\Results_'+image_N+'.csv')

    num_nuclei=len(df.axes[0])                   #dimensione dei dati generati  
    funct=create_nuclei_function(num_nuclei,df)			    #crea un array con le pos dei cdm dei nuclei dalla df di stardist
    df.to_csv(r'D:\matteo_citterio\risultati_segment\ROI_3_INT_Normal\Results.txt', header=None, index=None, sep=' ',columns=['XM','YM']  ,mode='a')

    density_matlab()

    dj=pd.read_csv(r'D:\matteo_citterio\Python Scripts\temp2.csv')
    os.remove(r'D:\matteo_citterio\risultati_segment\ROI_3_INT_Normal\Results.txt')
    os.remove(r'D:\matteo_citterio\Python Scripts\temp2.csv')
    dj.columns =['g(R)', 'R','c']
    bins_density=dj['R']
    density_profile=dj['g(R)']

    df1=pd.DataFrame(density_profile)
    k1=20
    df1=df1.rolling(k1+3).mean()

    size=len(dj.axes[0])
    x=np.zeros((size))
    y=np.zeros((size))

    for j in range(0,size):
        x[j]=bins_density[j]
        y[j]=df1.iloc[j]

    shift=np.count_nonzero(np.isnan(y))

    y = y[~np.isnan(y)]
    maxi=np.argsort(y)

    output_frame['7-tum'].iloc[i-1]=x[maxi[-1]+shift]    

for i in range(1,n_images_8tumor):

    print('8-tum: ',i)

    image_N=str(i)
    df=pd.read_csv(r'D:\matteo_citterio\risultati_segment\ROI_8_INT_Tumor\Results_'+image_N+'.csv')

    num_nuclei=len(df.axes[0])                   #dimensione dei dati generati  
    funct=create_nuclei_function(num_nuclei,df)			    #crea un array con le pos dei cdm dei nuclei dalla df di stardist
    df.to_csv(r'D:\matteo_citterio\risultati_segment\ROI_3_INT_Normal\Results.txt', header=None, index=None, sep=' ',columns=['XM','YM']  ,mode='a')

    density_matlab()

    dj=pd.read_csv(r'D:\matteo_citterio\Python Scripts\temp2.csv')
    os.remove(r'D:\matteo_citterio\risultati_segment\ROI_3_INT_Normal\Results.txt')
    os.remove(r'D:\matteo_citterio\Python Scripts\temp2.csv')
    dj.columns =['g(R)', 'R','c']
    bins_density=dj['R']
    density_profile=dj['g(R)']

    df1=pd.DataFrame(density_profile)
    k1=20
    df1=df1.rolling(k1+3).mean()

    size=len(dj.axes[0])
    x=np.zeros((size))
    y=np.zeros((size))

    for j in range(0,size):
        x[j]=bins_density[j]
        y[j]=df1.iloc[j]

    shift=np.count_nonzero(np.isnan(y))

    y = y[~np.isnan(y)]
    maxi=np.argsort(y)

    output_frame['8-tum'].iloc[i-1]=x[maxi[-1]+shift]    

for i in range(1,n_images_12normal):

    print('12-norm: ',i)

    image_N=str(i)
    df=pd.read_csv(r'D:\matteo_citterio\risultati_segment\ROI_12_INT_Normal\Results_'+image_N+'.csv')

    num_nuclei=len(df.axes[0])                   #dimensione dei dati generati  
    funct=create_nuclei_function(num_nuclei,df)			    #crea un array con le pos dei cdm dei nuclei dalla df di stardist
    df.to_csv(r'D:\matteo_citterio\risultati_segment\ROI_3_INT_Normal\Results.txt', header=None, index=None, sep=' ',columns=['XM','YM']  ,mode='a')

    density_matlab()

    dj=pd.read_csv(r'D:\matteo_citterio\Python Scripts\temp2.csv')
    os.remove(r'D:\matteo_citterio\risultati_segment\ROI_3_INT_Normal\Results.txt')
    os.remove(r'D:\matteo_citterio\Python Scripts\temp2.csv')
    dj.columns =['g(R)', 'R','c']
    bins_density=dj['R']
    density_profile=dj['g(R)']

    df1=pd.DataFrame(density_profile)
    k1=20
    df1=df1.rolling(k1+3).mean()

    size=len(dj.axes[0])
    x=np.zeros((size))
    y=np.zeros((size))

    for j in range(0,size):
        x[j]=bins_density[j]
        y[j]=df1.iloc[j]

    shift=np.count_nonzero(np.isnan(y))

    y = y[~np.isnan(y)]
    maxi=np.argsort(y)

    output_frame['12-norm'].iloc[i-1]=x[maxi[-1]+shift] 

for i in range(1,n_images_12tumor):

    print('12-tum: ',i)

    image_N=str(i)
    df=pd.read_csv(r'D:\matteo_citterio\risultati_segment\ROI_12_INT_Tumor\Results_'+image_N+'.csv')

    num_nuclei=len(df.axes[0])                   #dimensione dei dati generati  
    funct=create_nuclei_function(num_nuclei,df)			    #crea un array con le pos dei cdm dei nuclei dalla df di stardist
    df.to_csv(r'D:\matteo_citterio\risultati_segment\ROI_3_INT_Normal\Results.txt', header=None, index=None, sep=' ',columns=['XM','YM']  ,mode='a')

    density_matlab()

    dj=pd.read_csv(r'D:\matteo_citterio\Python Scripts\temp2.csv')
    os.remove(r'D:\matteo_citterio\risultati_segment\ROI_3_INT_Normal\Results.txt')
    os.remove(r'D:\matteo_citterio\Python Scripts\temp2.csv')
    dj.columns =['g(R)', 'R','c']
    bins_density=dj['R']
    density_profile=dj['g(R)']

    df1=pd.DataFrame(density_profile)
    k1=20
    df1=df1.rolling(k1+3).mean()

    size=len(dj.axes[0])
    x=np.zeros((size))
    y=np.zeros((size))

    for j in range(0,size):
        x[j]=bins_density[j]
        y[j]=df1.iloc[j]

    shift=np.count_nonzero(np.isnan(y))

    y = y[~np.isnan(y)]
    maxi=np.argsort(y)

    output_frame['12-tum'].iloc[i-1]=x[maxi[-1]+shift]       

output_frame.to_csv("max_ray_density.csv")