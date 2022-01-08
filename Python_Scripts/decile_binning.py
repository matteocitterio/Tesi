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
import matplotlib as mpl
from matplotlib.path import Path
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from functions import *
import scipy.stats as stats
data=[]
area_ratio=[]
dj=pd.read_csv(r'D:\matteo_citterio\Python Scripts\Area_ratio_per_image\Area_ratio_3_INT_Normal.csv')

num_nuclei=len(dj.axes[0])
for i in range(1,2):                  #ciclo a colonna che vuol dire a immagine

	for j in range(0,num_nuclei):

		if(not math.isnan(dj[str(i)].iloc[j])):

			data.append(dj[str(i)].iloc[j])
			#area_ratio.append(dj[str(i)].iloc[j])

print(len(data))

data=np.asarray(data)
data=(data)/(np.nanmean(data))

df = pd.DataFrame(data=data)
print(df.head())
hg = pd.qcut(df[0], q=10, precision=0)

print(hg.value_counts())
bin_counts=hg.value_counts()
results, bin_edges = pd.qcut(df[0],q=10,retbins=True)
results_table = pd.DataFrame(zip(bin_edges, bin_counts),columns=['bins', 'counts'])
print(results_table)