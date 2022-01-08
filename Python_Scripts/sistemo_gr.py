import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import warnings
import os
import sys
import scipy.spatial as sp
from scipy.spatial import distance
import skimage.io

from scipy import stats
from pandas.core.common import SettingWithCopyWarning
from shapely.geometry import MultiPoint, Point, Polygon, LineString,MultiPolygon
from shapely.ops import polygonize,unary_union
from scipy.spatial import Voronoi, voronoi_plot_2d,ConvexHull, convex_hull_plot_2d
from PIL import Image
from matplotlib.path import Path
from scipy.optimize import curve_fit
from functions import *


# for i in range(0,10):
#     g_r_and_nematic_AR('ROI_3_INT_Tumor',i+1,7000,10/3,False,700)
def func(x, a, b, c):
    return a * np.exp(-b * x) + c

scales=[]
fig,ax=plt.subplots()
xer=0

columns=['3_INT_Normal','3_INT_Tumor','5_INT']
indexis=np.arange(63)+1
indexis.tolist()
output_frame=pd.DataFrame(index=indexis,columns=columns)

		

# for i in range(0,63):

#     x=g_r_and_nematic_AR('ROI_3_INT_Normal',i+1,7000,1/30,False,1600,True)
# ############################################################################################################################################################################################################
# for i in range(0,10):

#     x=g_r_and_nematic_AR('ROI_3_INT_Tumor',i+1,7000,1/30,False,1500,True)

# ############################################################################################################################################################################################################

# for i in range(0,12):

#     x=g_r_and_nematic_AR('ROI_5_INT_Tumor',i+1,7000,1/30,False,1200,True)

# # ############################################################################################################################################################################################################

# for i in range(0,5):

#     x=g_r_and_nematic_AR('ROI_7_INT_Tumor',i+1,10000,1/30,False,1200,True)

############################################################################################################################################################################################################

for i in range(0,5):

    x=g_r_and_nematic_AR('ROI_8_INT_Tumor',i+1,10000,1/30,False,800,True)
    xs,ys=x[0],x[3]

############################################################################################################################################################################################################

for i in range(3,6):

    x=g_r_and_nematic_AR('ROI_9_INT_Normal',i+1,10000,1/30,False,800,True)
    xs,ys=x[0],x[3]

# for i in range(0,17):

#     x=g_r_and_nematic_AR('ROI_9_INT_Tumor',i+1,10000,1/30,False,1300,True)
#     xs,ys=x[0],x[3]

# for i in range(0,1):

#     x=g_r_and_nematic_AR('ROI_12_INT_Tumor',i+1,10000,1/30,False,1300,True)

# for i in range(0,23):

#     x=g_r_and_nematic_AR('ROI_12_INT_Normal',i+1,10000,1/30,False,1200,True)


for i in range(0,8):

    x=g_r_and_nematic_AR('ROI_13_INT_Tumor',i+1,10000,1/30,False,1100,True)

# output_frame.to_csv('D:\matteo_citterio\plot_al_fly\\nematics\\exponential_drop.csv')

# for i in range(0,12):

#     x=g_r_and_nematic_AR('ROI_5_INT_Tumor',i+1,100,1,1/30,False,10000,False)

#     xs,ys=x[0],x[3]
#     popt, pcov = curve_fit(func, xs, ys)
#     scales.append(1/popt[1])
#     ax.scatter(np.random.normal(1,0.05),1/popt[1],color='black')

# for i in range(0,5):

#     x=g_r_and_nematic_AR('ROI_7_INT_Tumor',i+1,100,1,1/30,False,10000,False)

#     xs,ys=x[0],x[3]
#     popt, pcov = curve_fit(func, xs, ys)
#     scales.append(1/popt[1])
#     ax.scatter(np.random.normal(1,0.05),1/popt[1],color='black')

# for i in range(0,5):

#     x=g_r_and_nematic_AR('ROI_8_INT_Tumor',i+1,100,1,1/30,False,10000,False)

#     xs,ys=x[0],x[3]
#     popt, pcov = curve_fit(func, xs, ys)
#     scales.append(1/popt[1])
#     ax.scatter(np.random.normal(1,0.05),1/popt[1],color='black')

# for i in range(0,6):

#     x=g_r_and_nematic_AR('ROI_9_INT_Normal',i+1,100,1,1/30,False,10000,False)

#     xs,ys=x[0],x[3]
#     popt, pcov = curve_fit(func, xs, ys)
#     scales.append(1/popt[1])
#     ax.scatter(np.random.normal(0,0.05),1/popt[1],color='green')

# for i in range(0,17):

#     x=g_r_and_nematic_AR('ROI_9_INT_Tumor',i+1,100,1,1/30,False,10000,False)

#     xs,ys=x[0],x[3]
#     popt, pcov = curve_fit(func, xs, ys)
#     scales.append(1/popt[1])
#     ax.scatter(np.random.normal(1,0.05),1/popt[1],color='black')

# for i in range(0,23):

#     x=g_r_and_nematic_AR('ROI_12_INT_Normal',i+1,100,1,1/30,False,10000,False)

#     xs,ys=x[0],x[3]
#     popt, pcov = curve_fit(func, xs, ys)
#     scales.append(1/popt[1])
#     ax.scatter(np.random.normal(0,0.05),1/popt[1],color='green')

# for i in range(0,1):

#     x=g_r_and_nematic_AR('ROI_12_INT_Tumor',i+1,100,1,1/30,False,10000,False)

#     xs,ys=x[0],x[3]
#     popt, pcov = curve_fit(func, xs, ys)
#     scales.append(1/popt[1])
#     ax.scatter(np.random.normal(1,0.05),1/popt[1],color='black')

# for i in range(0,8):

#     x=g_r_and_nematic_AR('ROI_13_INT_Tumor',i+1,100,1,1/30,False,10000,False)

#     xs,ys=x[0],x[3]
#     popt, pcov = curve_fit(func, xs, ys)
#     scales.append(1/popt[1])
#     ax.scatter(np.random.normal(1,0.05),1/popt[1],color='black')

ax.set_ylabel('scale')
ax.set_xticks([0,1])
ax.set_xticklabels(['Norm','Tum'])
ax.set_title('scale')
plt.show()