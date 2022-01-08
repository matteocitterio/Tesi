import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import warnings
import os
import sys
import scipy.spatial as sp
# import TiffImagePlugin

from pandas.core.common import SettingWithCopyWarning
from shapely.geometry import MultiPoint, Point, Polygon, LineString,MultiPolygon
from shapely.ops import polygonize,unary_union
from scipy.spatial import Voronoi, voronoi_plot_2d,ConvexHull, convex_hull_plot_2d
from PIL import Image
from matplotlib.path import Path
from functions import *

#omette i warnings per le copie di dataframes

warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)
warnings.simplefilter(action="ignore", category=RuntimeWarning)

######################################################################################################################################

#per avere in maniera veloce una list di tutti i nomi dei csv che vado poi ad aprire con pandas

#path=str(input('Please enter the path of the folder with the 16-bit images and the Results csv:') or os.getcwd())

temp_path='D:\matteo_citterio\\risultati_segment\ROI_13_INT_Tumor\\'

# fondamentale=True

ch_question=yes_or_no('Membrane or convex_hull contour? (YES for convex_hull)')

Results_='Results_'
csv=".csv"
txt=".txt"
standard='16-bit_'
tail2='_processed'
stats=False
# =False
val=1
moltiplicator=False
path=temp_path

filenames=[]
posizionis_filename=[]
titles=[]
indexis=['<Area>','<Perimetro>','<Aspect Ratio>','<Shape Factor>','<AR std>']

for i in range(1,val+1):										#creo le liste di nomi e le directories a seconda del # di images

	filenames.append(path+'/'+Results_+str(i)+csv)				#path dei vari results.csv
	posizionis_filename.append(path+'/'+"contour"+str(i-1)+csv)	#path dei vari cell.csv
	titles.append(standard+str(i)+tail2)						#16-bit_X_processed
	

#############################################################################################################################################

#leggo e riempio una list di df di pandas con dentro i csv

dfs=[]							#creo una lista di dataframes, ciascuno riempito con i risultati di stardist
posizionis=[]

k=0

for i in filenames:					#leggo i .csv e i .txt

	dfs.append(pd.read_csv(i))
	# print(dfs[k])
	posizionis.append(pd.read_csv(posizionis_filename[k],index_col=False))
	# print(posizionis[k])

	k+=1

#############################################################################################################################################

#creo una data frame che darò poi come output in csv dove saranno messi tutti i valori statistici calcolati

output_frame=pd.DataFrame(index=indexis,columns=titles)	#preparo una dataframe in cui metterò i dati totali in uscita

#############################################################################################################################################

#riempio una list con vari valori per fare statistiche e istogrammi

j=0

av_area=[]
av_perimetro=[]
av_shape_factor=[]
av_ar=[]
av_std=[]

for i in filenames:								#ciclo su tutti i nomi delle immagini

	num_nuclei=len(dfs[j].axes[0])						#numero di ROI=numero di righe

#	#creo una nuova colonna nella df dove mettere lo shape factor

	dfs[j]['Shape_Factor']=0
	
	for i in range(num_nuclei):

		dfs[j]['Shape_Factor'].iloc[i]=dfs[j]['Perim.'].iloc[i]/np.sqrt(dfs[j]['Area'].iloc[i])

#	#creo e riempio le list per le statistiche
	
	if stats==True:				#se viene richiesto di calcolare le statistiche
	
		#creo una figura su plt 4x4 per le quattro grandezze che guardo
	
		fig, ax = plt.subplots(2,2,tight_layout=True,figsize=(16,9))
		fig.suptitle(titles[j], fontsize=16)
	
		run_stats(num_nuclei, dfs[j], 'Area',av_area,ax,0,0,'red',title='area',end_bool=True)
		ax[0,0].set_xlabel('area [pixel^2]')
		run_stats(num_nuclei, dfs[j], 'Perim.',av_perimetro,ax,0,1,'green',title='perimetro',end_bool=True)
		ax[0,1].set_xlabel('perimetro[pixel]')
		run_stats(num_nuclei, dfs[j], 'AR',av_ar,ax,1,0,'yellow',title='AR loglog',end_bool=False,std_list=av_std)
		ax[1,0].set_xlabel('AR')
		run_stats(num_nuclei, dfs[j], 'Shape_Factor',av_shape_factor,ax,1,1,'brown',title='SF log-log',end_bool=False)
		ax[1,1].set_xlabel('SF')
		
		fig.savefig(os.getcwd()+'/risultati/infiltrating/'+titles[j]+'/'+titles[j]+'.pdf')
		
		#riempio la df di output e la mando fuori come csv e creo una list di list dei valori calcolati


			
#############################################################################################################################################

#prendo i contorni delle membrane e vedo quali nuclei ci sono dentro

	rois=[]
	rois_val=[]
	nuclei=create_nuclei_array(num_nuclei,dfs[j])			#crea un array con le pos dei cdm dei nuclei dalla df di stardist
	labels=rois_detections(nuclei,0,dfs[j],posizionis[j],rois,rois_val)	#modifica rois e ritorna la struttura labels

#############################################################################################################################################
	
	#faccio il voronoi per ogni roi
	
	fig2, ax = plt.subplots()
	vor=[]
	rois_restricted=[]
	membranes_restricted=[]
	
	voronoi_on_single_roi(labels,vor,rois_restricted,membranes_restricted,rois)	#riempie le prime tre liste in input
	if ch_question:
		voronoi_vertices_inside_convex_hull(vor,membranes_restricted)				#calcola ch e plotta solo le regioni del voronoi
	else:
		volumes=voronoi_vertices_inside_membrane(vor,membranes_restricted)			#plotta solo vertice all'interno della roi della membrana
																					#ad esso interno

	temp_list_of_things=[]

	for l in range(0,len(volumes)):
		
		for eta in range(0,len(volumes[l])):

			temp_list_of_things.append(volumes[l][eta])

	fig,ax=plt.subplots()
	print(len(temp_list_of_things))
	ax.hist(temp_list_of_things,bins=int(len(temp_list_of_things)/13))
	plt.show()	

	#tolgo ticks degli assi e labels

	plt.tick_params(
		axis='both',          # changes apply to the x-axis	
		which='both',      # both major and minor ticks are affected
		bottom=False,      # ticks along the bottom edge are off
		top=False,         # ticks along the top edge are off
		left=False,
		labelleft=False,
		labelbottom=False) # labels along the bottom edge are off
	
	plt.scatter(nuclei[:,0],nuclei[:,1],marker='o',c='tab:red',s=6)
	plt.xlim(left=np.min(nuclei[:,0]), right=np.max(nuclei[:,0])), plt.ylim(bottom=np.min(nuclei[:,1]), top=np.max(nuclei[:,1]))
	
	title='nuclei_map_'+str(j+1)+'.jpg'
	image_index=j+1
	img2=mpimg.imread('/home/matteo/Scrivania/tesi/fluorescenze/infiltrating/'+str(image_index)+'.tif')
	mpimg.imsave('temp.jpg',arr=img2,origin='lower')

	img=Image.open('temp.jpg')
	a,b=img.size
	fig2.set_size_inches(a/77.48, b/77 )				#abbastanza patologiche, magari da settare le dimensioni
	fig2.savefig(os.getcwd()+'/risultati/infiltrating/'+titles[j]+'/'+'Voronoi-'+titles[j]+'.png',transparent=True,bbox_inches='tight', pad_inches=0) #salvo senza bordi
	
#############################################################################################################################################

	#sovrappongo le immagini

	background=Image.open(os.getcwd()+'/'+'temp.jpg')
	overlay = Image.open(os.getcwd()+'/risultati/infiltrating/'+titles[j]+'/'+'Voronoi-'+titles[j]+'.png')
	
	background.paste(overlay, (0, 0), overlay)
	background.save(os.getcwd()+'/risultati/infiltrating/'+titles[j]+'/'+"tessellation.tiff")
	background.save(os.getcwd()+'/risultati/infiltrating/'+'tuttoinsieme/'+str(j)+".tiff")
	#le rigiro come le originali

	# img2=mpimg.imread(os.getcwd()+'/'+titles[j]+'/'+"tessellation.tiff")
	# mpimg.imsave(os.getcwd()+'/'+titles[j]+'/'+"tessellation.tiff",arr=img2,format='tiff',origin='lower')

#############################################################################################################################################
	
	j+=1									#aumento il contatore
	print('Image ', j, 'processed, # nuclei identified: ', num_nuclei)	#per capire a che immagine sono nell'esecuzione

computed=[av_area,av_perimetro,av_ar,av_shape_factor,av_std]				#list di lists
j=0

for k in indexis:

	output_frame.loc[k]=computed[j]
	j+=1

print('\n')
print(output_frame)
print('\n')
output_frame.to_csv("risultati_"+standard+".csv")
#############################################################################################################################################

