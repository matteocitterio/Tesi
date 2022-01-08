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
from matplotlib.patches import Ellipse
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import skimage.io
from skimage.measure import EllipseModel
from functions import *
img_ns_string=['ROI_9_INT_Normal','ROI_9_INT_Tumor','ROI_13_INT_Tumor','ROI_3_INT_Normal','ROI_3_INT_Tumor','ROI_5_INT_Tumor','ROI_7_INT_Tumor','ROI_8_INT_Tumor','ROI_12_INT_Normal','ROI_12_INT_Tumor']
kolumns=['9-norm','9-tum','13-tum','3-norm','3-tum','5-tum','7-tum','8-tum','12-norm','12-tum']
real_img_n=[64,11,13,6,6,7,18,24,2,9]
real_img_n=[7,18,9,64,11,13,6,6,24,2]



for eta in range(0,len(img_ns_string)):


    img_folder=img_ns_string[eta]
    header=kolumns[eta]
    img_N=real_img_n[eta]-1
    img_start=1
    disegno=False
    xar,ars_binned,xarea,area_binned,fluc_ar,fluc_area=voronoi_and_nuclei_per_folder(header,img_folder,img_start,img_N,disegno)

    indexis=np.arange(80)
    indexis.tolist()
    columns=['xar','ars_binned','xareas','areas_binned','fluc_ar','fluc_area']
    output_frame=pd.DataFrame(index=indexis,columns=columns)
    output_frame['xar']=xar
    output_frame['ars_binned']=ars_binned
    output_frame['xareas']=xarea
    output_frame['areas_binned']=area_binned
    output_frame['fluc_ar']=fluc_ar
    output_frame['fluc_area']=fluc_area

    output_frame.to_csv('D:\\matteo_citterio\\Python Scripts\\Voronoi_nuke_correlations\\'+img_ns_string[eta]+'80.csv')



plt.show()