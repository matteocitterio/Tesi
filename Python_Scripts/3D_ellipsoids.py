#Used to create an illustrative image with ellipsoids randomly orientated.

import numpy as np
from numpy import linalg
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as colors
from matplotlib.colors import LightSource

fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111, projection='3d')
ls = LightSource(270, 45)
# number of ellipsoids 
ellipNumber = 30

#set colour map so each ellipsoid as a unique colour
norm = colors.Normalize(vmin=0, vmax=ellipNumber)
cmap = cm.jet
m = cm.ScalarMappable(norm=norm, cmap=cmap)

rads=[]
x_maxs=[]
y_maxs=[]
z_maxs=[]

#compute and plot each ellipsoid iteratively
for indx in range(0,ellipNumber):
    # your ellispsoid and center in matrix form
    A = np.array([[np.random.random_sample(),0,0],
                  [0,np.random.random_sample(),0],
                  [0,0,np.random.random_sample()]])
    center = [indx*np.random.random_sample(),indx*np.random.random_sample(),indx*np.random.random_sample()]

    # find the rotation matrix and radii of the axes
    U, s, rotation = linalg.svd(A)
    radii = 1.0/np.sqrt(s) * 1 #reduce radii by factor 0.3 
    x_maxs.append(center[0]+radii[0])
    y_maxs.append(center[1]+radii[1])
    z_maxs.append(center[2]+radii[2])
    rads.append(max(radii))

    # calculate cartesian coordinates for the ellipsoid surface
    u = np.linspace(0.0, 2.0 * np.pi, 60)
    v = np.linspace(0.0, np.pi, 60)
    x = radii[0] * np.outer(np.cos(u), np.sin(v))
    y = radii[1] * np.outer(np.sin(u), np.sin(v))
    z = radii[2] * np.outer(np.ones_like(u), np.cos(v))

    for i in range(len(x)):
        for j in range(len(x)):
            [x[i,j],y[i,j],z[i,j]] = np.dot([x[i,j],y[i,j],z[i,j]], rotation) + center


    ax.plot_surface(x, y, z,lightsource=ls, rstride=2, cstride=2,  color='red', linewidth=1, alpha=1, shade=True)

print(max(rads))
ax.axes.xaxis.set_ticklabels([])
ax.axes.yaxis.set_ticklabels([])
ax.axes.zaxis.set_ticklabels([])
xaxis = np.linspace(int(min(x_maxs))-3, int(max(x_maxs))+3, 35)
yaxis = np.linspace(int(min(y_maxs))-3, int(max(y_maxs))+3, 35)
X, Y = np.meshgrid(xaxis, yaxis)
print(max(x_maxs),max(y_maxs))
Z = X*0+max(z_maxs)/2
print(Z.shape)

xaxis = np.linspace(int(min(x_maxs))-3, int(max(x_maxs))+3, 35)
zaxis = np.linspace(int(min(z_maxs))-3, int(max(z_maxs))+3, 35)
X, z = np.meshgrid(xaxis, zaxis)
y=X*0+2*max(y_maxs)/3

# plot the plane
ax.plot_surface(X, Y, Z, alpha=0.3)
ax.plot_surface(X, y, z,rstride=10,cstride=10,linewidth=5, color='red', alpha=0.3)

plt.show()
