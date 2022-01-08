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
import skimage.color

from pandas.core.common import SettingWithCopyWarning
from shapely.geometry import MultiPoint, Point, Polygon, LineString,MultiPolygon
from shapely.ops import polygonize,unary_union
from scipy.spatial import Voronoi, voronoi_plot_2d,ConvexHull, convex_hull_plot_2d
from PIL import Image
from matplotlib.path import Path
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl
from functions import *


img_folder='ROI_13_INT_Tumor'



for i in range(2,3):

    multi=2

    print(i)
    img=skimage.io.imread('D:\matteo_citterio\SCAN - Copia\\'+img_folder+'\\'+str(i+1)+'.tif')
    df=pd.read_csv('D:\matteo_citterio\\risultati_segment\\'+img_folder+'\Results_'+str(i+1)+'.csv')


    fig,ax=plt.subplots()
    ax.imshow(img)
    ax.set_title(str(i+1))
    ax.scatter(df['XM']*multi,df['YM']*multi,c='red',s=10)



plt.show()