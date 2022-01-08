import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import pandas as pd
import math
eta=0
colors=[]
columns=['9-norm','9-tum','13-tum','3-norm','3-tum','5-tum','7-tum','8-tum','12-norm','12-tum']
img_ns=[7,64,11,13,6,6,18,24,2,9]
real_img_n=[7, 18, 9, 64, 11, 13, 6, 6, 24, 2]
img_ns_string=['9_INT_Normal','9_INT_Tumor','13_INT_Tumor','3_INT_Normal','3_INT_Tumor','5_INT_Tumor','7_INT_Tumor','8_INT_Tumor','12_INT_Normal','12_INT_Tumor']
quant='AR'
divisors=[60,200,170,850,140,100,44,26,135,28]

df=pd.read_csv(r'D:\matteo_citterio\Python Scripts\AR_per_image\AR_'+img_ns_string[eta]+'.csv')

num_nuclei=len(df.axes[0])

data=[]

for i in range(1,real_img_n[eta]):                  #ciclo a colonna che vuol dire a immagine

    for j in range(0,num_nuclei):

        if(not math.isnan(df[str(i)].iloc[j])):

            data.append(df[str(i)].iloc[j]-1)

gamma = stats.gamma
a, loc, scale = 3, 0, 2
size = 20000
#y = gamma.rvs(a, loc, scale, size=size)

y=data

x = np.linspace(0, max(y), 100)
# fit
param = gamma.fit(y, floc=0)
pdf_fitted = gamma.pdf(x, *param)
plt.plot(x, pdf_fitted, color='r')

# plot the histogram
plt.hist(y, bins=30,density=True)

plt.show()