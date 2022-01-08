import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import warnings
import os
import sys
import scipy.spatial as sp
import time
import math
from numpy.linalg import eig, inv

from pandas.core.common import SettingWithCopyWarning, is_full_slice
from shapely.geometry import MultiPoint, Point, Polygon, LineString,MultiPolygon
from shapely.ops import polygonize,unary_union
from scipy.spatial import Voronoi, voronoi_plot_2d,ConvexHull, convex_hull_plot_2d
from PIL import Image
from matplotlib.path import Path
from matplotlib.patches import Ellipse
import matplotlib.colors as colors
from scipy.integrate import quad
import scipy.stats as stats
from scipy.integrate import dblquad
import skimage.io
from skimage.measure import EllipseModel
import matlab.engine

def fitEllipse(x,y):
	x = x[:,np.newaxis]
	y = y[:,np.newaxis]
	D =  np.hstack((x*x, x*y, y*y, x, y, np.ones_like(x)))
	S = np.dot(D.T,D)
	C = np.zeros([6,6])
	C[0,2] = C[2,0] = 2; C[1,1] = -1
	E, V =  eig(np.dot(inv(S), C))
	n = np.argmax(np.abs(E))
	a = V[:,n]
	return a

def ellipse_center(a):
	b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
	num = b*b-a*c
	x0=(c*d-b*f)/num
	y0=(a*f-b*d)/num
	return np.array([x0,y0])


def ellipse_angle_of_rotation( a ):
	b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
	return 0.5*np.arctan(2*b/(a-c))


def ellipse_axis_length( a ):
	b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
	up = 2*(a*f*f+c*d*d+g*b*b-2*b*d*f-a*c*g)
	down1=(b*b-a*c)*( (c-a)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
	down2=(b*b-a*c)*( (a-c)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
	res1=np.sqrt(up/down1)
	res2=np.sqrt(up/down2)
	return np.array([res1, res2])

def ellipse_angle_of_rotation2( a ):
	b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
	if b == 0:
		if a > c:
			return 0
		else:
			return np.pi/2
	else: 
		if a > c:
			return np.arctan(2*b/(a-c))/2
		else:
			return np.pi/2 + np.arctan(2*b/(a-c))/2


def yes_or_no(question):					#utile per gestire lo stream, ritorna True se in tastiera do y e False  
	reply = str(input(question+' (y/n): ')).lower().strip()	#viceversa
	if reply[0] == 'y':
		return True
	if reply[0] == 'n':
		return False
	else:
		return yes_or_no("Uhhhh... please enter ")

def voronoi_vertices_in_path(p,vor,ax=0,disegno=False):

	list_of_regions_inside_path=[]
	nuclei_positions_x=[]
	nuclei_positions_y=[]

	print(ax)

	for i in range(0,(vor.points[:,0]).size):

		region=vor.regions[vor.point_region[i]]

	# for regions in vor.regions:											#ciclo sulle regioni del voronoi
		region=np.asarray(region)										#lo trasformo da lista ad array

		if np.all(region>=0) and np.all(region.size!=0):

			grid=p.contains_points(vor.vertices[region])
			grid=np.asarray(grid)

			if np.all(grid==True):

				list_of_regions_inside_path.append(region)
				nuclei_positions_x.append(vor.points[i][0])
				nuclei_positions_y.append(vor.points[i][1])
				if (disegno):
					ax.scatter(vor.points[i][0],vor.points[i][1],c='orange')
	
	if (disegno):

		for i in list_of_regions_inside_path:

			i=np.asarray(i)		
			ax.plot(vor.vertices[i, 0],vor.vertices[i, 1], 'r-',lw=5)
			ax.plot([vor.vertices[i,0][-1],vor.vertices[i,0][0]],[vor.vertices[i,1][-1],vor.vertices[i,1][0]], 'r-',lw=5)
	
	return list_of_regions_inside_path,nuclei_positions_x,nuclei_positions_y
				
def run_stats(num_nuclei, dataframe, dataframe_kw,av_observable,ax,pos_x_ax,pos_y_ax,histo_color,title,end_bool,std_list=[]):

#parametri:
#-)num_nuclei: numero dei nuclei segmentati in totale di cui si vogliono calcolare le proprietà
#-) dataframe: df di pandas dei dati
#-)dataframe_kw: str del titolo della particolare osservabile nella df
#av_observable: list dove metto poi i valori medi
#ax: asse in cui plottare
#pos_x_ax: posizione all'interno del subplot
#pos_y_ax: vedi sopra
#histo_color:str, colore dell'istogramma
#title:titolo del subplot
#end_bool: definisce alcuni parametri di ax
#std_list: per calcolare la deviazione standard degli aspect ratio

	observable_list=[]

	for i in range(num_nuclei):
	
			observable_list.append(dataframe[dataframe_kw].iloc[i])

	#converto le list in np.array perchè è easy fare statistica

	observable_array = np.array(observable_list)

	#calcolo valori statistici e li metto in lists (che avranno quindi solo # di immagini di valori):

	best_observable=np.average(observable_array)
	std_observable=np.std(observable_array)
	if (dataframe_kw=='AR'):

		std_list.append(std_observable)

	av_observable.append(best_observable)

	#faccio gli istogrammi per le varie statistiche

	num_bins=int(num_nuclei/13)

	ax[pos_x_ax,pos_y_ax].hist(observable_array, bins=num_bins,color=histo_color)
	ax[pos_x_ax,pos_y_ax].set_ylabel('Counts')
	ax[pos_x_ax,pos_y_ax].set_title(title)
		
	if end_bool:
		
		ax[pos_x_ax,pos_y_ax].axvline(best_observable,lw=2, color='black',ls='--', label= 'average')
		ax[pos_x_ax,pos_y_ax].legend()
	else:
		
		ax[pos_x_ax,pos_y_ax].set_xscale('log')
		ax[pos_x_ax,pos_y_ax].set_yscale('log')

def create_nuclei_array(num_nuclei,df):				#crea un array con le pos dei cdm dei nuclei dalla df di stardist

	nuclei_val=np.zeros((num_nuclei,5))
	nuclei=np.zeros((num_nuclei,2))				#array vuoto dove poi metto le posizionivo dei centri di massa.

	for i in range(0,num_nuclei):				#riempio l'array con le posizioni dalla df

		nuclei[i][0]=df['XM'].iloc[i]
		nuclei[i][1]=df['YM'].iloc[i]
		nuclei_val[i][0]=df['Area'].iloc[i]
		nuclei_val[i][1]=df['AR'].iloc[i]
		nuclei_val[i][2]=df['Circ.'].iloc[i]
		nuclei_val[i][3]=df['Perim.'].iloc[i]
		nuclei_val[i][4]=df['Round'].iloc[i]

	return nuclei,nuclei_val

def create_nuclei_quantity_array(num_nuclei,df,string_kw):

	quantity=np.zeros((num_nuclei))

	for i in range(0,num_nuclei):

		quantity[i]=df[string_kw].iloc[i]

	return quantity

def get_distance(x1,y1,x2,y2):

	return ((x1-x2)**2+(y1-y2)**2)**0.5
	
def min_distance(x,y,cx,cy):
#x and y are arrays representing a line, cx,cy are the coordinates of the center of a circle and r is the ray of the last one
#x and y have the same size
	distances=[]
	for i in range(0,len(x)):

		distances.append(get_distance(x[i],y[i],cx,cy))

	distances=np.asarray(distances)

	return np.amin(distances)

def intersection(x,y,cx,cy,r):
#x and y are arrays representing a line, cx,cy are the coordinates of the center of a circle and r is the ray of the last one
#x and y have the same size

	return min_distance(x,y,cx,cy)<=r

def box_intersection(xv1,yv1,xv2,yv2,xh1,yh1,xh2,yh2,cx,cy,r):
#x and y are arrays representing a line, cx,cy are the coordinates of the center of a circle and r is the ray of the last one
#x and y have the same size

	return (min_distance(xv1,yv1,cx,cy)<=r) or (min_distance(xv2,yv2,cx,cy)<=r) or  (min_distance(xh1,yh1,cx,cy)<=r) or  (min_distance(xh2,yh2,cx,cy)<=r)

def get_vertex_distances(cx,cy):

	distances=np.zeros((4))
	distances[0]=get_distance(0,0,cx,cy)
	distances[1]=get_distance(0,1460,cx,cy)
	distances[0]=get_distance(1936,1460,cx,cy)
	distances[0]=get_distance(1936,0,cx,cy)	
	return distances

def get_angoli_buoni(normal_dist,r,vertex_minor_distances):

	angles=np.zeros((4))
	for i in range(0,4):

		angles[i]=np.arccos(normal_dist[i]/r)	

def pair_correlation_function(data,data_size,first_image_side,second_image_side,num_bins):

	fig,ax =plt.subplots(3,2, figsize=(16,9))

	a=get_distances_matrix(data,data_size)
	bins=np.linspace(0,((first_image_side**2)+(second_image_side**2))**0.5,num_bins)
	media_ar=np.average(data[:,2])
	quad_ar=data[:,2]**2
	media_quad_ar=np.average(quad_ar)
	media_perim=np.average(data[:,4])
	quad_perim=data[:,4]**2
	media_quad_perim=np.average(quad_perim)
	media_circ=np.average(data[:,3])
	quad_circ=data[:,3]**2
	media_quad_circ=np.average(quad_circ)
	media_area=np.average(data[:,5])
	quad_area=data[:,5]**2
	media_quad_area=np.average(quad_area)
	media_round=np.average(data[:,6])
	quad_round=data[:,6]**2
	media_quad_round=np.average(quad_round)

	digitized= np.digitize(a[0,:], bins)

	matrix_4_histo=np.zeros((num_bins,data_size))

	j=0
	for i in range(1,len(digitized)):

		if(a[0,i]!=0):

			matrix_4_histo[digitized[i]-1,j]+=1
		else:
			print(j)
			j+=1

	medie=[]

	for i in range(0,num_bins):

		medie.append(matrix_4_histo[i,:].mean())

	bin_means_ar= [a[1,:][digitized == i].mean() for i in range(1, len(bins))]
	bin_means_perim=[a[2,:][digitized == i].mean() for i in range(1, len(bins))]
	bin_means_circ=[a[3,:][digitized == i].mean() for i in range(1, len(bins))]
	bin_means_area=[a[4,:][digitized == i].mean() for i in range(1, len(bins))]
	bin_means_round=[a[5,:][digitized == i].mean() for i in range(1, len(bins))]

	bins_revised=np.zeros((num_bins-1))

	for i in range(1,len(bins)):

		bins_revised[i-1]=(bins[i]+bins[i-1])/2
	ax[0][0].plot(bins_revised,bin_means_ar,'-o',c='red')
	ax[0][0].axhline(media_ar**2,label='$(<ar>)^2$')								#quadrato della media, valore atteso in +infty
	ax[0][0].axhline(media_quad_ar,c='green',label='$<ar^2>$')
	ax[0][0].set_title('Aspect Ratio')
	ax[0][0].set_ylabel('C(r)')
	ax[0][0].legend()
	ax[0][1].plot(bins,medie,'o',c='red')
	ax[0][1].axhline(1, label='1 count')								#quadrato della media, valore atteso in +infty
	ax[0][1].axhline(15,c='green',label='15 counts')
	ax[0][1].set_ylabel('#counts')
	ax[0][1].set_title('<conteggi> per r')
	ax[1][0].plot(bins_revised,bin_means_area,'-o',c='red')
	ax[1][0].axhline(media_area**2,label='$(<area>)^2$')								#quadrato della media, valore atteso in +infty
	ax[1][0].axhline(media_quad_area,c='green',label='$<area^2>$')
	ax[1][0].set_title('Area')
	ax[1][0].set_ylabel('C(r)')
	ax[1][0].legend()
	ax[1][1].plot(bins_revised,bin_means_circ,'-o',c='red')
	ax[1][1].axhline(media_circ**2,label='$(<circ>)^2$')								#quadrato della media, valore atteso in +infty
	ax[1][1].axhline(media_quad_circ,c='green',label='$<circ^2>$')
	ax[1][1].set_title('Circularity')
	ax[1][1].legend()
	ax[2][0].plot(bins_revised,bin_means_perim,'-o',c='red')
	ax[2][0].axhline(media_perim**2,label='$(<perim>)^2$')								#quadrato della media, valore atteso in +infty
	ax[2][0].axhline(media_quad_perim,c='green',label='$<perim^2>$')
	ax[2][0].set_title('Perimeter')
	ax[2][0].set_xlabel('r [pixel]')
	ax[2][0].set_ylabel('C(r)')
	ax[2][0].legend()
	ax[2][1].plot(bins_revised,bin_means_round,'-o',c='red')
	ax[2][1].axhline(media_round**2,label='$(<round>)^2$')								#quadrato della media, valore atteso in +infty
	ax[2][1].axhline(media_quad_round,c='green',label='$<round^2>$')
	ax[2][1].set_title('Roundness')
	ax[2][1].set_xlabel('r [pixel]')
	ax[2][1].set_ylabel('C(r)')
	ax[2][1].legend()

	plt.tight_layout()

	plt.show()
				
def create_circularity_array(num_nuclei,df):

	Circ=np.zeros((num_nuclei))

	for i in range(0,num_nuclei):

		Circ[i]=df['Circ.'].iloc[i]

	return Circ

def create_function_circularity(num_nuclei,df):

	temp_nuc=create_nuclei_array(num_nuclei,df)
	temp_ar=create_circularity_array(num_nuclei,df)

	funct=np.zeros((num_nuclei,3))

	for i in range(0,num_nuclei):

		funct[i][0]=temp_nuc[i][0]		#x
		funct[i][1]=temp_nuc[i][1]		#y
		funct[i][2]=temp_ar[i]			#ar

	return funct

def pair_correlation_function_for_iterations(data,data_size,first_image_side,second_image_side,num_bins,kw):

	if (kw=='AR'):							#capisco la particolare kw

		indicator=2
	elif(kw=='Perim.'):

		indicator=4
	elif(kw=='Circ.'):

		indicator=3
	elif(kw=='Area'):

		indicator=5
	elif(kw=='Round'):

		indicator=6
	
	image_diagonal=((first_image_side**2)+(second_image_side**2))**0.5				#diagonale dell'immagine
	a=get_distances_and_quantities_vectors(data,data_size,name='ok')	#matrice delle distanze e vettore con grandezze
	bins=np.linspace(0,image_diagonal,num_bins)
	media=data[:,indicator].mean()
	media_quad=(data[:,indicator]**2).mean()
	digitized= np.digitize(a[0,:], bins)					#ritorna indici che servono per binnare poi il mio array

	bin_rays=[len(a[indicator-1,:][digitized == i]) for i in range(1, num_bins+1)]			#conteggi a raggio
	bin_means= [a[indicator-1,:][digitized == i].mean() for i in range(1, num_bins+1)]		#C(r)

	return bins, bin_means, bin_rays, media, media_quad

def pololifit(data_size,side_1,side_2,num_bins,):

	alpha=np.random.random(data_size)*side_1
	beta=np.random.random(data_size)*side_2
	gamma=np.random.random(data_size)+1

	punti=np.zeros((data_size,3))
	h=0
	print('Creating points')
	for i in range(0,len(alpha)):

		punti[i][0]=alpha[i]
		punti[i][1]=beta[i]
		punti[i][2]=gamma[h]
		h+=1

	distances_matrix=np.zeros((data_size,data_size))
	h=0
	distances_vector=[]
	value_vector_ar=[]

	for i in range(0,data_size):

		# print(i,'/',data_size)

		for j in range(h,data_size):

			distances_matrix[i][j]=get_distance(punti[i][0],punti[i][1],punti[j][0],punti[j][1])
			distances_vector.append(distances_matrix[i][j])
			value_vector_ar.append(punti[i][2]*punti[j][2])

		h+=1

	print('binning')

	distances_vector,value_vector_ar=np.asarray(distances_vector),np.asarray(value_vector_ar)		#array con distanze reciproche e valori
	a=np.zeros((2,len(distances_vector)))
	a[0,:]=distances_vector
	a[1,:]=value_vector_ar
	bins=np.linspace(0,((side_1**2)+(side_2**2))**0.5,num_bins)		
	digitized= np.digitize(a[0,:], bins)	
	inmezzo=[]

	for i in range(1,num_bins+1):

		# print(i,'/',num_bins+1)
		inmezzo.append(len(a[1,:][digitized == i]))


	inmezzo=np.asarray(inmezzo)
	max=inmezzo.max()

	return np.polyfit(bins,inmezzo,deg=10)

def fit_equation(poly_coeffs,x):

	return poly_coeffs[0]*(x**10)+poly_coeffs[1]*(x**9)+poly_coeffs[2]*(x**8)+poly_coeffs[3]*(x**7)+poly_coeffs[4]*(x**6)+poly_coeffs[5]*(x**5)+poly_coeffs[6]*(x**4)+poly_coeffs[7]*(x**3)+poly_coeffs[8]*(x**2)+poly_coeffs[9]*(x)+poly_coeffs[10]

def density_profile_function(data,data_size,side1,side2,num_bins):

	smallest_side=np.minimum(side1,side2)
	image_diagonal=((side1**2)+(side2**2))**0.5				#diagonale dell'immagine
	print('getting distances')
	a=get_distances_and_quantities_vectors(data,data_size,'ad casum',simple=True)	#matrice delle distanze e vettore con grandezze
	bins=np.linspace(0,image_diagonal,num_bins)
	
	digitized= np.digitize(a[0,:], bins)					#ritorna indici che servono per binnare poi il mio array
	
	bin_rays=[len(a[0,:][digitized == i]) for i in range(1, num_bins+1)]			#conteggi a raggio

	dr=(image_diagonal)/num_bins
	print('dr',dr)
	numer_of=int((smallest_side/(dr*2))-1)
	print('nuovo numero calcoalto',numer_of)
	print('numero input',num_bins)
	print(smallest_side)

	bin_rays_not_modi=[]

	for j in range(0,numer_of):


		print(j,'/',numer_of)

		bin_rays_not_modi.append(bin_rays[j])
		bin_rays[j]=bin_rays[j]/analytic_function(side1,side2,(bins[j]+2*dr),1,data_size,dr=dr/2)
	

	return bins, bin_rays,numer_of,bin_rays_not_modi

def density_matlab():

	eng = matlab.engine.start_matlab()
	eng.density	(nargout=0)
	eng.quit()

def density_check():

	eng = matlab.engine.start_matlab()
	eng.density_check	(nargout=0)
	eng.quit()

def density_profile_function_amara(data,data_size,side1,side2,num_bins):

	smallest_side=np.minimum(side1,side2)
	image_diagonal=((side1**2)+(side2**2))**0.5				#diagonale dell'immagine
	simul_data_size=1000
	alpha=np.random.random(simul_data_size)*side1        #uniforme reali
	beta=np.random.random(simul_data_size)*side2
	gamma=np.random.random(simul_data_size)+1

	punti=np.zeros((data_size,3))
	h=0
	print('Creating points')
	for i in range(0,len(alpha)):

		punti[i][0]=alpha[i]
		punti[i][1]=beta[i]
		h+=1

	print('points created')	
	print('getting real distances')
	a=get_distances_and_quantities_vectors(data,data_size,'real',simple=True)	#matrice delle distanze e vettore con grandezze
	print('getting simulated distances')
	b=get_distances_and_quantities_vectors(punti,simul_data_size,'simu',simple=True)
	bins=np.linspace(0,image_diagonal,num_bins)
	
	print('digitizing real data...')
	digitized_real= np.digitize(a[0,:], bins)					#ritorna indici che servono per binnare poi il mio array
	print('digitizing simulation data...')
	digitized_simu= np.digitize(b[0,:], bins)
	step=image_diagonal/num_bins
	num_iterations=int((smallest_side/(2*step)))
	simul_bin_rays=np.zeros((num_iterations))
	bin_rays=np.zeros((num_iterations))                           #il primo elemento sarà per definiziione nullo: numero di elementi ad una distanza zero (escluso se stesso)
	cum_real=0
	cum_simu=0

	xs=np.zeros((num_iterations))

	area=side1*side2
	density=simul_data_size/area

	for i in range(1,num_iterations):                           #perchè digitized ritorna il bin a cui appartiene, bins invece ha gli estremi

		print(i,'/',num_iterations)
		xs[i]=bins[i]
		bin_rays[i]=cum_real+len(a[0,:][digitized_real == i])      #[i-1] semplicemente per come si contano gli elementi
		cum_real+=len(a[0,:][digitized_real == i])
		simul_bin_rays[i]=cum_simu+len(b[0,:][digitized_simu == i])
		cum_simu+=len(b[0,:][digitized_simu == i])

	simul_bin_rays=simul_bin_rays*data_size/(simul_data_size)
	density_profile=bin_rays/simul_bin_rays
	

	return xs,density_profile

def g_r(img_folder,img_N,random_points,check=False):

	df=pd.read_csv('D:\matteo_citterio\\risultati_segment\\'+img_folder+'\Results_'+str(img_N)+'.csv')
	posizioni=pd.read_csv('D:\matteo_citterio\\risultati_segment\\'+img_folder+'\contour'+str(img_N-1)+'.csv')
	num_nuclei=len(df.axes[0])

	img=skimage.io.imread('D:\matteo_citterio\SCAN - Copia\\'+img_folder+'\\'+str(img_N)+'.tif')
	img_star=Image.open(r'D:\matteo_citterio\SCAN - Copia\\'+img_folder+'\\'+str(img_N)+'.tif')
	a,b=img_star.size
	diagonal=((a**2)+(b**2))**0.5
	print(a,b)

	rois=[]
	rois_val=[]
	nuclei=create_nuclei_array(num_nuclei,df)[0]	

	random_size_number=random_points                                #numero dei nuclei fake da generare

	random_nuclei=np.zeros((random_size_number,2))
	random_nuclei2=np.zeros((10000,2))
	random_size_number2=random_points

	for i in range(0,random_size_number):                   #genero punti random

		temp=np.zeros((2))
		temp[0]=np.random.uniform(0,1)*a
		temp[1]=np.random.uniform(0,1)*b
		random_nuclei[i]=temp

	for i in range(0,random_size_number2):                   #genero punti random

		temp=np.zeros((2))
		temp[0]=np.random.uniform(0,1)*a
		temp[1]=np.random.uniform(0,1)*b
		random_nuclei2[i]=temp

	labels=rois_detections(random_nuclei,0,df,posizioni,rois,rois_val)	#modifica rois e ritorna la struttura labels: vede quali nuclei simulati sono nel contorno
	total_points=[]
	for i in range(len(rois)):                               #total_points ha dentro tutti i nuclei fake nel contorno
		for j in range(len(rois[i])):
			total_points.append(rois[i][j])

	print('len:', len(total_points))
	print(((len(total_points)**2)-len(total_points))/2)

	rois2=[]
	labels2=rois_detections(random_nuclei2,0,df,posizioni,rois2,rois_val)	#modifica rois e ritorna la struttura labels: vede quali nuclei simulati sono nel contorno
	total_points2=[]
	for i in range(len(rois2)):                               #total_points ha dentro tutti i nuclei fake nel contorno
		for j in range(len(rois2[i])):
			total_points2.append(rois2[i][j])

	punti_random=np.zeros((len(total_points),3))             #per distances vector pnti fake
	for i in range(0,len(total_points)):

		punti_random[i][0]=total_points[i][0]
		punti_random[i][1]=total_points[i][1]

	punti_random2=np.zeros((len(total_points2),3))             #per distances vector pnti fake
	for i in range(0,len(total_points2)):

		punti_random2[i][0]=total_points2[i][0]
		punti_random2[i][1]=total_points2[i][1]

	punti_real=np.zeros((num_nuclei,3))                      #per distance vector punti real
	for i in range(0,num_nuclei):
		punti_real[i][0]=nuclei[i][0]*3.33
		punti_real[i][1]=nuclei[i][1]*3.33

	num_bins=500

	distance_vector_simu=get_distances_vector(punti_random,len(total_points),'simu')
	distance_vector_simu2=get_distances_vector(punti_random2,len(total_points2),'simu2')
	distance_vector_real=get_distances_vector(punti_real,num_nuclei,'real')

	step=diagonal/num_bins
	#num_iterations=int((max(a,b)/(2*step)))
	num_iterations=int(diagonal/step)
	print('n bins',num_bins)
	bins=np.linspace(0,diagonal,num_bins)

	digitized_simu= np.digitize(distance_vector_simu[0,:], bins)					#ritorna indici che servono per binnare poi il mio array AL DATO BINNAGGIO
	digitized_simu2= np.digitize(distance_vector_simu2[0,:], bins)
	digitized_real= np.digitize(distance_vector_real[0,:], bins)

	simu_bin_rays=np.zeros((int(num_iterations)))                         #il primo elemento sarà per definiziione nullo: numero di elementi ad una distanza zero (escluso se stesso)
	simu_bin_rays2=np.zeros((int(num_iterations))) 
	real_bin_rays=np.zeros((int(num_iterations)))    
	cum=0
	cum2=0
	cum_real=0
	xs=np.zeros((int(num_iterations)))                           #serve per il plot

	for i in range(1,int(num_iterations)):                       #perchè digitized ritorna il bin a cui appartiene, bins invece ha gli estremi

		print(i,'/',int(num_iterations))
		xs[i]=bins[i]
		simu_bin_rays[i]=cum+len(distance_vector_simu[0,:][digitized_simu == i])         #[i-1] semplicemente per come si contano gli elementi
		cum+=len(distance_vector_simu[0,:][digitized_simu == i])
		simu_bin_rays2[i]=cum2+len(distance_vector_simu2[0,:][digitized_simu2 == i])         #[i-1] semplicemente per come si contano gli elementi
		cum2+=len(distance_vector_simu2[0,:][digitized_simu2 == i])
		real_bin_rays[i]=cum_real+len(distance_vector_real[0,:][digitized_real == i])
		cum_real+=len(distance_vector_real[0,:][digitized_real == i])

	norm_simu_1=((len(total_points)**2)-len(total_points))/2
	norm_simu_2=((len(total_points2)**2)-len(total_points2))/2
	norm_real=((num_nuclei**2)-num_nuclei)/2
	for i in range(0,int(num_iterations)):
		simu_bin_rays[i]=simu_bin_rays[i]/norm_simu_1
		simu_bin_rays2[i]=simu_bin_rays2[i]/norm_simu_2
		real_bin_rays[i]=real_bin_rays[i]/norm_real

	quoz=np.zeros((int(num_iterations)))
	quoz2=np.zeros((int(num_iterations)))

	for i in range(0,int(num_iterations)):

		if(simu_bin_rays[i]==0):
			quoz[i]=0
			quoz2[i]=0

		else:

			quoz[i]=real_bin_rays[i]/simu_bin_rays[i]
			quoz2[i]=simu_bin_rays2[i]/simu_bin_rays[i]

	if (check==True):

		fig,ax=plt.subplots(1,4)

		ax[0].plot(xs,simu_bin_rays,'o',c='red',label='simulation')
		ax[0].plot(xs,real_bin_rays,'o',c='green',label='real')
		ax[0].set_xlabel('R')
		ax[0].set_ylabel('av. #counts per point')
		ax[0].legend()

		ax[1].plot(xs,quoz,'-',c='blue')
		ax[1].set_xlabel('R')
		ax[1].set_ylabel('g(R)')
		ax[1].axhline(1,ls='--',color='red')

		ax[2].plot(xs,quoz2,'-',c='blue')
		ax[2].set_title('control')
		ax[2].axhline(1,ls='--',color='red')

		ax[3].plot(xs,simu_bin_rays,'o',c='red',label='simulation1')
		ax[3].plot(xs,simu_bin_rays2,'o',c='blue',label='simulation2')

		ax[3].legend()

		fig,ax=plt.subplots()

		ax.imshow(img)
		ax.scatter(punti_random[:,0],punti_random[:,1],c='k',s=5)
		ax.scatter(punti_real[:,0],punti_real[:,1])

		plt.show()

	return xs,quoz

def g_r_and_nematic(img_folder,img_N,random_points,multi=10/3,check=False):

	df=pd.read_csv('D:\matteo_citterio\\risultati_segment\\'+img_folder+'\Results_'+str(img_N)+'.csv')
	posizioni=pd.read_csv('D:\matteo_citterio\\risultati_segment\\'+img_folder+'\contour'+str(img_N-1)+'.csv')
	num_nuclei=len(df.axes[0])

	img=skimage.io.imread('D:\matteo_citterio\SCAN - Copia\\'+img_folder+'\\'+str(img_N)+'.tif')
	img_star=Image.open(r'D:\matteo_citterio\SCAN - Copia\\'+img_folder+'\\'+str(img_N)+'.tif')
	a,b=img_star.size
	diagonal=((a**2)+(b**2))**0.5
	print(a,b)

	rois=[]
	rois_val=[]
	nuclei=create_nuclei_array(num_nuclei,df)[0]	

	random_size_number=random_points                                #numero dei nuclei fake da generare

	random_nuclei=np.zeros((random_size_number,2))
	random_nuclei2=np.zeros((10000,2))
	random_size_number2=random_points

	for i in range(0,random_size_number):                   #genero punti random

		temp=np.zeros((2))
		temp[0]=np.random.uniform(0,1)*a
		temp[1]=np.random.uniform(0,1)*b
		random_nuclei[i]=temp

	for i in range(0,random_size_number2):                   #genero punti random

		temp=np.zeros((2))
		temp[0]=np.random.uniform(0,1)*a
		temp[1]=np.random.uniform(0,1)*b
		random_nuclei2[i]=temp

	labels=rois_detections(random_nuclei,0,df,posizioni,rois,rois_val)	#modifica rois e ritorna la struttura labels: vede quali nuclei simulati sono nel contorno
	total_points=[]
	for i in range(len(rois)):                               #total_points ha dentro tutti i nuclei fake nel contorno
		for j in range(len(rois[i])):
			total_points.append(rois[i][j])

	#print('len:', len(total_points))
	#x,print(((len(total_points)**2)-len(total_points))/2)

	rois2=[]
	labels2=rois_detections(random_nuclei2,0,df,posizioni,rois2,rois_val)	#modifica rois e ritorna la struttura labels: vede quali nuclei simulati sono nel contorno
	total_points2=[]
	for i in range(len(rois2)):                               #total_points ha dentro tutti i nuclei fake nel contorno
		for j in range(len(rois2[i])):
			total_points2.append(rois2[i][j])

	punti_random=np.zeros((len(total_points),3))             #per distances vector pnti fake
	for i in range(0,len(total_points)):

		punti_random[i][0]=total_points[i][0]
		punti_random[i][1]=total_points[i][1]

	punti_random2=np.zeros((len(total_points2),3))             #per distances vector pnti fake
	for i in range(0,len(total_points2)):

		punti_random2[i][0]=total_points2[i][0]
		punti_random2[i][1]=total_points2[i][1]

	punti_real=np.zeros((num_nuclei,3))                      #per distance vector punti real
	for i in range(0,num_nuclei):
		punti_real[i][0]=nuclei[i][0]*multi
		punti_real[i][1]=nuclei[i][1]*multi
		punti_real[i][2]=df['Angle'].iloc[i]

	num_bins=500

	distance_vector_simu=get_distances_vector(punti_random,len(total_points),'simu')
	distance_vector_simu2=get_distances_vector(punti_random2,len(total_points2),'simu2')
	distance_vector_real=get_distances_and_angles(punti_real,num_nuclei,'real')

	step=diagonal/num_bins
	#num_iterations=int((max(a,b)/(2*step)))
	num_iterations=int(diagonal/step)
	#print('n bins',num_bins)
	bins=np.linspace(0,diagonal,num_bins)

	digitized_simu= np.digitize(distance_vector_simu[0,:], bins)					#ritorna indici che servono per binnare poi il mio array AL DATO BINNAGGIO
	digitized_simu2= np.digitize(distance_vector_simu2[0,:], bins)
	digitized_real= np.digitize(distance_vector_real[0,:], bins)

	simu_bin_rays=np.zeros((int(num_iterations)))                         #il primo elemento sarà per definiziione nullo: numero di elementi ad una distanza zero (escluso se stesso)
	simu_bin_rays2=np.zeros((int(num_iterations))) 
	real_bin_rays=np.zeros((int(num_iterations)))
	nem_energy=np.zeros((num_iterations))     
	cum=0
	cum2=0
	cum_real=0
	xs=np.zeros((int(num_iterations)))                           #serve per il plot

	for i in range(1,int(num_iterations)):                       #perchè digitized ritorna il bin a cui appartiene, bins invece ha gli estremi

		print(i,'/',int(num_iterations))
		xs[i]=bins[i]
		simu_bin_rays[i]=cum+len(distance_vector_simu[0,:][digitized_simu == i])         #[i-1] semplicemente per come si contano gli elementi
		cum+=len(distance_vector_simu[0,:][digitized_simu == i])
		simu_bin_rays2[i]=cum2+len(distance_vector_simu2[0,:][digitized_simu2 == i])         #[i-1] semplicemente per come si contano gli elementi
		cum2+=len(distance_vector_simu2[0,:][digitized_simu2 == i])
		real_bin_rays[i]=cum_real+len(distance_vector_real[0,:][digitized_real == i])
		nem_energy[i]=np.mean(distance_vector_real[1,:][digitized_real == i])
		cum_real+=len(distance_vector_real[0,:][digitized_real == i])

	norm_simu_1=((len(total_points)**2)-len(total_points))/2
	norm_simu_2=((len(total_points2)**2)-len(total_points2))/2
	norm_real=((num_nuclei**2)-num_nuclei)/2
	for i in range(0,int(num_iterations)):
		simu_bin_rays[i]=simu_bin_rays[i]/norm_simu_1
		simu_bin_rays2[i]=simu_bin_rays2[i]/norm_simu_2
		real_bin_rays[i]=real_bin_rays[i]/norm_real

	quoz=np.zeros((int(num_iterations)))
	quoz2=np.zeros((int(num_iterations)))

	for i in range(0,int(num_iterations)):

		if(simu_bin_rays[i]==0):
			quoz[i]=0
			quoz2[i]=0

		else:

			quoz[i]=real_bin_rays[i]/simu_bin_rays[i]
			quoz2[i]=simu_bin_rays2[i]/simu_bin_rays[i]

	if (check==True):

		fig,ax=plt.subplots(1,4)
		fig.suptitle(img_folder+' img: '+str(img_N))

		ax[0].plot(xs,simu_bin_rays,'o',c='red',label='simulation')
		ax[0].plot(xs,real_bin_rays,'o',c='green',label='real')
		ax[0].set_xlabel('R')
		ax[0].set_ylabel('av. #counts per point')
		ax[0].legend()

		ax[1].plot(xs,quoz,'-',c='blue')
		ax[1].set_xlabel('R')
		ax[1].set_ylabel('g(R)')
		ax[1].axhline(1,ls='--',color='red')

		ax[2].plot(xs,quoz2,'-',c='blue')
		ax[2].set_title('control')
		ax[2].axhline(1,ls='--',color='red')

		ax[3].plot(xs,nem_energy,'-',c='purple',label='nem_energy')
		ax[3].axhline(2/np.pi,ls='--',color='red',label='random orient.')
		ax[3].axhline(1,ls='--',color='blue',label='alignment')
		ax[3].axhline(0,ls='--',color='green',label='perpendicularity')

		ax[3].legend()

		fig.savefig('D:\matteo_citterio\plot_al_fly\\nematics\\'+img_folder+'\\'+str(img_N)+'_nem.pdf')

		fig,ax=plt.subplots()
		fig.suptitle(img_folder+' img: '+str(img_N))

		ax.imshow(img)
		ax.scatter(punti_random[:,0],punti_random[:,1],c='k',s=5)
		ax.scatter(punti_real[:,0],punti_real[:,1])

		fig.savefig('D:\matteo_citterio\plot_al_fly\\nematics\\'+img_folder+'\\'+str(img_N)+'_control.pdf')

		#plt.show()

	return xs,quoz,nem_energy

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
	new_cmap = colors.LinearSegmentedColormap.from_list(
		'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
		cmap(np.linspace(minval, maxval, n)))
	return new_cmap

def voronoi_and_nuclei_per_folder(header,img_folder,img_start,img_N,disegno=False):

	ars=[]
	nuke_ars=[]
	areas=[]
	nuke_areas=[]

	for i in range(img_start,img_N+1):

		print(i,'/',img_N)	

		df=pd.read_csv('D:\matteo_citterio\\risultati_segment\\'+img_folder+'\Results_'+str(i)+'.csv')
		posizioni=pd.read_csv('D:\matteo_citterio\\risultati_segment\\'+img_folder+'\contour'+str(i-1)+'.csv')
		multi_df=pd.read_csv('multipliers.csv')
		pixel_sizes_df=pd.read_csv('pixel_sizes.csv')
		multi=multi_df[img_folder].iloc[i-1]
		coeff_sistemazione=(multi)
		coeff_area_calibration=(pixel_sizes_df[img_folder].iloc[i-1])**2
		num_nuclei=len(df.axes[0])
		ax=0

		nuclei=create_nuclei_array(num_nuclei,df)[0]*coeff_sistemazione			#crea un array con le pos dei cdm dei nuclei dalla df di stardist

		rois=[]
		rois_val=[]
		labels=rois_detections(nuclei,0,df,posizioni,rois,rois_val)	#modifica rois e ritorna la struttura labels

		vor=[]
		rois_restricted=[]
		membranes_restricted=[]

		voronoi_on_single_roi(labels,vor,rois_restricted,membranes_restricted,rois)	#riempie le prime tre liste in input
		volumes=voronoi_vertices_inside_membrane(vor,membranes_restricted,ax,disegno)	

		df['XM']=df['XM']*multi
		df['YM']=df['YM']*multi

		temp_list_of_things=[]

		for l in range(0,len(volumes)):
		
			for eta in range(0,len(volumes[l])):

				temp_list_of_things.append(volumes[l][eta])

		for eta in range(0,len(volumes[1])):          #le macchie

			for j in range(0,len(volumes[1][eta])):   

				
				numero=len(volumes[1][eta][j])
				points=np.zeros((numero,2))
				x,y=np.zeros((numero)),np.zeros((numero))
				for l in range(0,numero):

					x[l]=volumes[1][eta][j][l]
					y[l]=volumes[2][eta][j][l]
					points[l,0]=volumes[1][eta][j][l]
					points[l,1]=volumes[2][eta][j][l]

				ell = EllipseModel()
				ell.estimate(points)

				volumes[3][eta][j]=np.around(volumes[3][eta][j],decimals=2)
				volumes[4][eta][j]=np.around(volumes[4][eta][j],decimals=2)

				print(j,'/',len(volumes[1][eta]))

				

				if(ell.params is not None):

					df['XM']=np.around(df['XM'],decimals=2)
					df['YM']=np.around(df['YM'],decimals=2)
					areas.append(volumes[0][eta][j]*coeff_area_calibration)
					nuke_areas.append(df.loc[(df['XM']==volumes[3][eta][j]) & (df['YM']==volumes[4][eta][j])]['Area']*multi*coeff_area_calibration)
					nuke_ars.append(df.loc[(df['XM']==volumes[3][eta][j]) & (df['YM']==volumes[4][eta][j])]['AR'])
					xc, yc, alpha, beta, theta = ell.params
					ars.append(max(alpha,beta)/min(alpha,beta))

	frame_voronoi_ar_medi=pd.read_csv('AR_voronoi_mean.csv')
	frame_voronoi_area_medi=pd.read_csv('Area_ratio_voronoi_mean.csv')

	print(header)

	print(frame_voronoi_ar_medi[header].iloc[i-1]*5)

	limit_ar=frame_voronoi_ar_medi[header].iloc[i-1]*5
	limit_area=frame_voronoi_area_medi[header].iloc[i-1]*5
	fluctuations_ar,fluctuations_area=[],[]
	aspect_ratios=np.zeros((2,len(ars)))
	
	for i in range(0,len(ars)):

		if(ars[i]<=limit_ar):
			aspect_ratios[0,i]=ars[i]
			aspect_ratios[1,i]=nuke_ars[i]

		

	print(aspect_ratios[1,:])

	bins=np.linspace(min(aspect_ratios[0,:]),max(aspect_ratios[0,:]),80)
	digitized_ars= np.digitize(aspect_ratios[0,:], bins)					#ritorna indici che servono per binnare poi il mio array AL DATO BINNAGGIO
	xs=[]
	bins_step=bins[1]-bins[0]
	bins_ars=np.zeros((80))

	real_areas=np.zeros((2,len(areas)))

	for i in range(0,len(ars)):

		if(areas[i]<=limit_area):

			real_areas[0,i]=areas[i]
			real_areas[1,i]=nuke_areas[i]
		
		

	bins_areas=np.linspace(min(real_areas[0,:]),max(real_areas[0,:]),80)
	digitized_areas= np.digitize(real_areas[0,:], bins_areas)					#ritorna indici che servono per binnare poi il mio array AL DATO BINNAGGIO
	xs_areas=[]
	bins_step_areas=bins_areas[1]-bins_areas[0]
	bins_delle_aree=np.zeros((80))

	for i in range(1,81):                       #perchè digitized ritorna il bin a cui appartiene, bins invece ha gli estremi

		# print(img_folder,':',img_N,':',i,'/',1000)
		xs_areas.append(bins_areas[i-1]+bins_step_areas/2)
		bins_delle_aree[i-1]=np.average(real_areas[1,:][digitized_areas == i]) 
		fluctuations_area.append(np.nanstd(real_areas[1,:][digitized_areas == i]))

	for i in range(1,81):                       #perchè digitized ritorna il bin a cui appartiene, bins invece ha gli estremi

		# print(img_folder,':',img_N,':',i,'/',1000)
		xs.append(bins[i-1]+bins_step/2)
		bins_ars[i-1]=np.average(aspect_ratios[1,:][digitized_ars == i])
		fluctuations_ar.append(np.nanstd(aspect_ratios[1,:][digitized_ars == i]))


	# temp1,temp2=[],[]
	# for eta in range(0,len(ars)):
	# 	if(ars[eta]<4):
	# 		temp1.append(ars[eta])
	# 		temp2.append(nuke_ars[eta])

	print(len(nuke_areas))
	print(len(nuke_ars))

	fig,ax=plt.subplots()	
	ax.plot(ars,nuke_ars,'o')
	ax.set_title('Aspect ratios')
	ax.set_xlabel('voronoi AR')
	ax.set_ylabel('Nuclear AR')
	print('ars',stats.spearmanr(ars,nuke_ars))

	fig,ax=plt.subplots()
	ax.plot(xs,bins_ars,'o')
	ax.set_title('binned AR')
	ax.set_xlabel('Voronoi AR')
	ax.set_ylabel('Nuclear AR')

	fig,ax=plt.subplots()
	ax.plot(xs_areas,bins_delle_aree,'o')
	ax.set_title('binned Area')
	ax.set_xlabel('Voronoi Area $[\mu m^{2}]$')
	ax.set_ylabel('Nuclear Area $[\mu m^{2}]$')

	fig,aj=plt.subplots()
	aj.plot(areas,nuke_areas,'o')
	aj.set_title('Areas')
	aj.set_xlabel('voronoi Area $[\mu m^{2}]$')
	aj.set_ylabel('Nuclear Area $[\mu m^{2}]$')
	print('areas',stats.spearmanr(areas,nuke_areas))

	print('binned ar',stats.spearmanr(xs,bins_ars))
	print('binned areas',stats.spearmanr(xs_areas,bins_delle_aree))

	return xs,bins_ars,xs_areas,bins_delle_aree,fluctuations_ar,fluctuations_area

def voronoi_and_quantities_per_folder(img_folder,img_start,img_N,disegno=False):

	ars=[]
	nuke_ars=[]
	areas=[]
	nuke_areas=[]

	for i in range(img_start,img_N+1):

		print(i,'/',img_N)	

		df=pd.read_csv('D:\matteo_citterio\\risultati_segment\\'+img_folder+'\Results_'+str(i)+'.csv')
		posizioni=pd.read_csv('D:\matteo_citterio\\risultati_segment\\'+img_folder+'\contour'+str(i-1)+'.csv')
		multi_df=pd.read_csv('multipliers.csv')
		pixel_sizes_df=pd.read_csv('pixel_sizes.csv')
		multi=multi_df[img_folder].iloc[i-1]
		coeff_sistemazione=(multi)
		coeff_area_calibration=(pixel_sizes_df[img_folder].iloc[i-1])**2
		
		num_nuclei=len(df.axes[0])

		ax=0

		if (disegno):
			fig,ax=plt.subplots()
			img=skimage.io.imread('D:\matteo_citterio\SCAN - Copia\\'+img_folder+'\\'+str(i)+'.tif')
			img_star=Image.open(r'D:\matteo_citterio\SCAN - Copia\\'+img_folder+'\\'+str(i)+'.tif')
			a,b=img_star.size
			diagonal=((a**2)+(b**2))**0.5
			#print(a,b)
			ax.imshow(img)

		rois=[]
		rois_val=[]
		nuclei=create_nuclei_array(num_nuclei,df)[0]*coeff_sistemazione	

		labels=rois_detections(nuclei,0,df,posizioni,rois,rois_val)	#modifica rois e ritorna la struttura labels

		vor=[]
		rois_restricted=[]
		membranes_restricted=[]

		voronoi_on_single_roi(labels,vor,rois_restricted,membranes_restricted,rois)	#riempie le prime tre liste in input
		volumes=voronoi_vertices_inside_membrane(vor,membranes_restricted,ax=ax,disegno=disegno)

		for eta in range(0,len(volumes[1])):          #le macchie

			for j in range(0,len(volumes[1][eta])):   

				areas.append(volumes[0][eta][j]*coeff_area_calibration)
				numero=len(volumes[1][eta][j])
				points=np.zeros((numero,2))
				x,y=np.zeros((numero)),np.zeros((numero))
				for l in range(0,numero):

					x[l]=volumes[1][eta][j][l]
					y[l]=volumes[2][eta][j][l]
					points[l,0]=volumes[1][eta][j][l]
					points[l,1]=volumes[2][eta][j][l]

				ell = EllipseModel()
				ell.estimate(points)

				if(ell.params is not None):

					

					xc, yc, alpha, beta, theta = ell.params
					# volumes[3][eta][j]=np.around(volumes[3][eta][j],decimals=2)
					# volumes[4][eta][j]=np.around(volumes[4][eta][j],decimals=2)
					# df['XM']=np.around(df['XM'],decimals=2)
					# df['YM']=np.around(df['YM'],decimals=2)
					# print(df['XM'])
					# nuke_ars.append(df.loc[(df['XM']==volumes[3][eta][j]) & (df['YM']==volumes[4][eta][j])]['AR'].item())
					# print(df.loc[(df['XM']==volumes[3][eta][j]) & (df['YM']==volumes[4][eta][j])]['AR'].item())
					# nuke_areas.append(df.loc[(df['XM']==volumes[3][eta][j]) & (df['YM']==volumes[4][eta][j])]['Area'].item())
					#print(i,j,max(alpha,beta)/min(alpha,beta))
					ars.append(max(alpha,beta)/min(alpha,beta))

	plt.show()

	return ars,areas,nuke_ars,nuke_areas

def g_r_and_nematic_AR(img_folder,img_N,random_points,portion=1,check=False,num_bins=500,csv=True):

	multi_df=pd.read_csv('multipliers.csv')
	pixel_sizes_df=pd.read_csv('pixel_sizes.csv')

	multi=multi_df[img_folder].iloc[img_N-1]

	df=pd.read_csv('D:\matteo_citterio\\risultati_segment\\'+img_folder+'\Results_'+str(img_N)+'.csv')
	posizioni=pd.read_csv('D:\matteo_citterio\\risultati_segment\\'+img_folder+'\contour'+str(img_N-1)+'.csv')
	num_nuclei=len(df.axes[0])

	img=skimage.io.imread('D:\matteo_citterio\SCAN - Copia\\'+img_folder+'\\'+str(img_N)+'.tif')
	img_star=Image.open(r'D:\matteo_citterio\SCAN - Copia\\'+img_folder+'\\'+str(img_N)+'.tif')
	a,b=img_star.size
	diagonal=((a**2)+(b**2))**0.5
	print(a,b)

	rois=[]
	rois_val=[]
	nuclei=create_nuclei_array(num_nuclei,df)[0]	

	random_size_number=random_points                                #numero dei nuclei fake da generare

	random_nuclei=np.zeros((random_size_number,2))
	random_nuclei2=np.zeros((10000,2))
	random_size_number2=random_points

	for i in range(0,random_size_number):                   #genero punti random

		temp=np.zeros((2))
		temp[0]=np.random.uniform(0,1)*a
		temp[1]=np.random.uniform(0,1)*b
		random_nuclei[i]=temp

	for i in range(0,random_size_number2):                   #genero punti random

		temp=np.zeros((2))
		temp[0]=np.random.uniform(0,1)*a
		temp[1]=np.random.uniform(0,1)*b
		random_nuclei2[i]=temp

	labels=rois_detections(random_nuclei,0,df,posizioni,rois,rois_val)	#modifica rois e ritorna la struttura labels: vede quali nuclei simulati sono nel contorno
	total_points=[]
	for i in range(len(rois)):                               #total_points ha dentro tutti i nuclei fake nel contorno
		for j in range(len(rois[i])):
			total_points.append(rois[i][j])

	rois2=[]
	labels2=rois_detections(random_nuclei2,0,df,posizioni,rois2,rois_val)	#modifica rois e ritorna la struttura labels: vede quali nuclei simulati sono nel contorno
	total_points2=[]
	for i in range(len(rois2)):                               #total_points ha dentro tutti i nuclei fake nel contorno
		for j in range(len(rois2[i])):
			total_points2.append(rois2[i][j])

	punti_random=np.zeros((len(total_points),3))             #per distances vector pnti fake
	for i in range(0,len(total_points)):

		punti_random[i][0]=total_points[i][0]
		punti_random[i][1]=total_points[i][1]

	punti_random2=np.zeros((len(total_points2),3))             #per distances vector pnti fake
	for i in range(0,len(total_points2)):

		punti_random2[i][0]=total_points2[i][0]
		punti_random2[i][1]=total_points2[i][1]

	punti_real=np.zeros((num_nuclei,4))                      #per distance vector punti real
	for i in range(0,num_nuclei):
		punti_real[i][0]=nuclei[i][0]*multi
		punti_real[i][1]=nuclei[i][1]*multi
		punti_real[i][2]=df['Angle'].iloc[i]
		punti_real[i][3]=df['AR'].iloc[i]

	num_bins=num_bins

	distance_vector_simu=get_distances_vector(punti_random,len(total_points),'simu')
	distance_vector_simu2=get_distances_vector(punti_random2,len(total_points2),'simu2')
	distance_vector_real,mean_AR=get_distances_and_angles_with_AR(punti_real,num_nuclei,'real')

	print(num_bins)

	step=diagonal/num_bins
	#num_iterations=int((max(a,b)/(2*step)))
	num_iterations=int((diagonal*portion)/step)
	print('n bins',num_bins)
	bins=np.linspace(0,diagonal,num_bins)

	digitized_simu= np.digitize(distance_vector_simu[0,:], bins)					#ritorna indici che servono per binnare poi il mio array AL DATO BINNAGGIO
	digitized_simu2= np.digitize(distance_vector_simu2[0,:], bins)
	digitized_real= np.digitize(distance_vector_real[0,:], bins)

	simu_bin_rays=np.zeros((int(num_iterations)))                         #il primo elemento sarà per definiziione nullo: numero di elementi ad una distanza zero (escluso se stesso)
	simu_bin_rays2=np.zeros((int(num_iterations))) 
	real_bin_rays=np.zeros((int(num_iterations)))
	nem_energy=np.zeros((num_iterations))     
	ar_correlations=np.zeros((num_iterations))    
	just_nem=np.zeros((num_iterations))    
	cum=0
	cum2=0
	cum_real=0
	xs=np.zeros((int(num_iterations)))                           #serve per il plot

	for i in range(1,int(num_iterations)):                       #perchè digitized ritorna il bin a cui appartiene, bins invece ha gli estremi

		print(img_folder,':',img_N,':',i,'/',int(num_iterations))
		xs[i]=bins[i]
		simu_bin_rays[i]=len(distance_vector_simu[0,:][digitized_simu == i])         #[i-1] semplicemente per come si contano gli elementi
		cum+=len(distance_vector_simu[0,:][digitized_simu == i])
		simu_bin_rays2[i]=len(distance_vector_simu2[0,:][digitized_simu2 == i])         #[i-1] semplicemente per come si contano gli elementi
		cum2+=len(distance_vector_simu2[0,:][digitized_simu2 == i])
		real_bin_rays[i]=len(distance_vector_real[0,:][digitized_real == i])
		nem_energy[i]=np.nanmean(distance_vector_real[1,:][digitized_real == i])
		ar_correlations[i]=np.nanmean(distance_vector_real[2,:][digitized_real == i])
		just_nem[i]=np.nanmean(distance_vector_real[3,:][digitized_real == i])
		cum_real+=len(distance_vector_real[0,:][digitized_real == i])

	ar_correlations[np.isnan(ar_correlations)] = 0
	nem_energy[np.isnan(nem_energy)] = 0
	just_nem[np.isnan(just_nem)] = 0


	norm_simu_1=((len(total_points)**2)-len(total_points))/2
	norm_simu_2=((len(total_points2)**2)-len(total_points2))/2
	norm_real=((num_nuclei**2)-num_nuclei)/2
	for i in range(0,int(num_iterations)):
		simu_bin_rays[i]=simu_bin_rays[i]/norm_simu_1
		simu_bin_rays2[i]=simu_bin_rays2[i]/norm_simu_2
		real_bin_rays[i]=real_bin_rays[i]/norm_real

	quoz=np.zeros((int(num_iterations)))
	quoz2=np.zeros((int(num_iterations)))

	for i in range(0,int(num_iterations)):

		if(simu_bin_rays[i]==0):
			quoz[i]=0
			quoz2[i]=0

	#ROI_3_INT_Tumor

		else:

			quoz[i]=real_bin_rays[i]/simu_bin_rays[i]
			quoz2[i]=simu_bin_rays2[i]/simu_bin_rays[i]

	# print(ar_correlations)
	coeff_sistemazione=(multi)*(pixel_sizes_df[img_folder].iloc[img_N-1])
	xs=xs*coeff_sistemazione

	shift=np.count_nonzero(np.isnan(quoz))
	quoz = quoz[~np.isnan(quoz)]

	maxi=np.argsort(quoz)

	quad=np.zeros((num_nuclei))
	for i in range(num_nuclei):
		quad[i]=(df['AR'].iloc[i]**2)


		
	if (check==True):

		fig,ax=plt.subplots(4,1)
		ax[0].plot(xs,quoz,'-',c='blue')
		ax[0].set_xlabel('R [$\mu m$]')
		ax[0].set_ylabel('g(R)')
		ax[0].axhline(1,ls='--',color='red')
		
		ax[0].axvline(xs[maxi[-1]+shift],c='red')



		ax[1].plot(xs,ar_correlations,'-',c='purple',label='nem_energy')
		ax[1].axhline(np.average(quad),ls='--',color='red',label='$<AR^{2}>$')
		#ax[3].axhline(max_AR**2,ls='--',color='blue',label='alignment')
		ax[1].axhline(np.average(df['AR'])**2,ls='--',color='green',label='(<AR>)^{2}')
		ax[1].axvline(xs[maxi[-1]+shift],c='red')
		ax[1].set_ylabel('$<AR_{i}AR_{j}>$')
		ax[1].legend()

		ax[2].plot(xs,nem_energy,'-',c='purple',label='nem_energy')
		ax[2].axhline(2/np.pi*mean_AR,ls='--',color='red',label='random orient. ?')
		#ax[3].axhline(max_AR**2,ls='--',color='blue',label='alignment')
		ax[2].axhline(0,ls='--',color='green',label='perpendicularity')
		ax[2].axvline(xs[maxi[-1]+shift],c='red')
		ax[2].legend()
		ax[2].set_ylabel('$<AR_{i}AR_{j}|cos{\\theta}|>$')

		ax[3].plot(xs,just_nem,'-',c='purple',label='nem_energy')
		ax[3].axhline(2/np.pi,ls='--',color='red',label='random orient')
		ax[3].axhline(1,ls='--',color='blue',label='alignment')
		ax[3].axhline(0,ls='--',color='green',label='perpendicularity')
		ax[3].axvline(xs[maxi[-1]+shift],c='red')
		ax[3].legend()
		ax[3].set_ylabel('$<|cos{\\theta}|>$')
		ax[3].set_xlabel('R [$\mu m$]')

		fig,ax=plt.subplots(1,3)
		fig.suptitle(img_folder+' img: '+str(img_N))

		ax[0].plot(xs,simu_bin_rays,'o',c='red',label='simulation')
		ax[0].plot(xs,real_bin_rays,'o',c='green',label='real')
		ax[0].set_xlabel('R [$\mu m$]')
		ax[0].set_ylabel('av. #counts per point')
		ax[0].legend()

		ax[1].plot(xs,quoz,'-',c='blue')
		ax[1].set_xlabel('R [$\mu m$]')
		ax[1].set_ylabel('g(R)')
		ax[1].axhline(1,ls='--',color='red')

		ax[2].plot(xs,quoz2,'-',c='blue')
		ax[2].set_title('control')
		ax[2].axhline(1,ls='--',color='red')

		# ax[3].plot(xs,nem_energy,'-',c='purple',label='nem_energy')
		

		# max_AR=max(df['AR'])

		# ax[3].axhline(2/np.pi*mean_AR,ls='--',color='red',label='random orient. ?')
		# #ax[3].axhline(max_AR**2,ls='--',color='blue',label='alignment')
		# ax[3].axhline(0,ls='--',color='green',label='perpendicularity')

		# ax[3].legend()

		#fig.savefig('D:\matteo_citterio\plot_al_fly\\nematics\\'+img_folder+'\\'+str(img_N)+'_nem.pdf')
		fig,ax=plt.subplots()
		fig.suptitle(img_folder+' img: '+str(img_N))

		ax.imshow(img)
		ax.scatter(punti_random[:,0],punti_random[:,1],c='k',s=5)
		ax.scatter(punti_real[:,0],punti_real[:,1])

		#fig.savefig('D:\matteo_citterio\plot_al_fly\\nematics\\'+img_folder+'\\'+str(img_N)+'_control.pdf')

		plt.show()

	if (csv==True):

		columns=['xs','quoz','nem_energy_AR','nem_energy','autocorrelations','max','quad','quad_top']
		indexis=np.arange(xs.size)+1
		indexis.tolist()
		output_frame=pd.DataFrame(index=indexis,columns=columns)
		output_frame['xs']=xs
		output_frame['quad']=np.average(quad)
		output_frame['quad_top']=np.average(df['AR'])**2
		output_frame['quoz']=quoz
		output_frame['nem_energy_AR']=nem_energy
		output_frame['nem_energy']=just_nem
		output_frame['autocorrelations']=ar_correlations
		output_frame['max']=xs[maxi[-1]+shift]
		output_frame.to_csv('D:\matteo_citterio\plot_al_fly\\nematics\\'+img_folder+'\\'+str(img_N)+'_corr.csv')

	return xs,quoz,nem_energy,just_nem,ar_correlations,xs[maxi[-1]+shift]

def get_nematic_angle(theta1,theta2):

	return (min(abs(theta1-theta2),abs(min(theta1,theta2)-max(theta1,theta2)+180))/180)*np.pi

def get_position_in_get_distances(get_distances_array,i,j,num_nuclei):			#i,j sono numeri ordinali, tipo la prima particella con la seconda i=1,j=1

	
	i=i-1															#li trasformo in indici di python per array
	j=j-1

	if (i==j):

		return -1

	if (i>j):

		t=i
		i=j
		j=t

	return (i*num_nuclei-(i*(i+1)/2)+(j-i))-1						#ritorna l'indice di python dell'array

def get_indices_in_get_distance_per_particle(particle_no,get_distances_array,num_nuclei):		#particle no è ordinale va trasfrormato in indice


	i=particle_no

	indices=np.zeros((2,num_nuclei-1),dtype=int)

	for j in range(1,i):

		indices[0,j-1]=int(get_position_in_get_distances(get_distances_array,i,j,num_nuclei))
		indices[1,j-1]=j
	for j in range(i+1,num_nuclei+1):

		indices[0,j-2]=int(get_position_in_get_distances(get_distances_array,i,j,num_nuclei))		#in indices va messo j-2 per ordinali + i==j
		indices[1,j-2]=j

	return indices								#array dove nella prima riga ci sono gli indici da mettere dentro il vector delle distanze, nella seconda ci sono i numeri della particolare particella

def get_nearest_neighbors(get_distances_array, particle_no,num_nuclei):

	indices=get_indices_in_get_distance_per_particle(particle_no,get_distances_array,num_nuclei)
	distances=np.zeros((2,num_nuclei-1))		#array dove nella prima riga ci sono le distanze nella seconda ci sono i numeri della particolare particella

	for i in range(0,num_nuclei-1):

		distances[0,i]=get_distances_array[0,indices[0,i]]
		distances[1,i]=indices[1,i]				#ci sono dentro numeri ordinali, non indici informatici di array
	
	y=np.argsort(distances[0,:])
	j=np.zeros((6))
	for i in range(0,6):
		j[i]=distances[1,y[i]]

	return j									#ritorna gli ordinali delle particelle

def compute_energy(get_distances_array,num_nuclei,df,img_no):

	E=0

	for i in range(0,num_nuclei):

		print('Energy for',i,'/',num_nuclei,'img',img_no)
		NN=get_nearest_neighbors(get_distances_array,i,num_nuclei)
		
		for j in range(0,NN.size):

			theta=get_nematic_angle(df['Angle'].iloc[int(NN[j]-1)],df['Angle'].iloc[i])

			E+=-math.cos(theta)
		
	return E

def compute_energy_Voronoi_NN(df):

	E=0
	NNs=NN_through_Voronoi(df)

	for i in range(0,len(df.axes[0])):

		e=0

		for j in range(0,len(NNs[i])):

			angle=math.radians(abs(df['Angle'].iloc[NNs[i][j]]-df['Angle'].iloc[i]))
			# print(angle)
			# angle=get_nematic_angle(df['Angle'].iloc[NNs[i][j]],df['Angle'].iloc[i])
			# print(angle)
			# e+=-math.cos(theta)
			e+=(abs(np.cos(angle)))
		e=(e/len(NNs[i]))
		
		E+=e
	
	return E/len(df.axes[0])

def compute_energy_Voronoi_NN_with_metrics(df,max_len,coeff):

	E=0
	NNs=NN_through_Voronoi(df)

	nuclei_with_neighbours=0

	for i in range(0,len(df.axes[0])):

		e=0
		counter=0

		for j in range(0,len(NNs[i])):

			if (get_distance(df['XM'].iloc[i],df['YM'].iloc[i],df['XM'].iloc[NNs[i][j]],df['YM'].iloc[NNs[i][j]])*coeff<=max_len):
				counter+=1
				angle=math.radians(abs(df['Angle'].iloc[NNs[i][j]]-df['Angle'].iloc[i]))
				e+=(abs(np.cos(angle)))
		if (counter!=0):
			e=(e/counter)
			nuclei_with_neighbours+=1
			E+=e

	return E/nuclei_with_neighbours

def create_nuclei_function(num_nuclei,df):

	temp_nuc=create_nuclei_array(num_nuclei,df)[0]
	temp_ar=create_nuclei_quantity_array(num_nuclei,df,'AR')  
	temp_area=create_nuclei_quantity_array(num_nuclei,df,'Area')
	temp_circ=create_nuclei_quantity_array(num_nuclei,df,'Circ.')
	temp_perim=create_nuclei_quantity_array(num_nuclei,df,'Perim.')
	temp_round=create_nuclei_quantity_array(num_nuclei,df,'Round')

	funct=np.zeros((num_nuclei,7))

	for i in range(0,num_nuclei):

		funct[i][0]=temp_nuc[i][0]		#x
		funct[i][1]=temp_nuc[i][1]		#y
		funct[i][2]=temp_ar[i]			#ar
		funct[i][3]=temp_circ[i]		#circ
		funct[i][4]=temp_perim[i]		#perim
		funct[i][5]=temp_area[i]		#area
		funct[i][6]=temp_round[i]		#round

	return funct

def quads(R):

	qs=np.zeros((10))

	# qs[0]=-0.5*(R**2)*((((R**2)-1))**0.5)*(np.pi)
	# qs[1]=-(5/2)*R**3
	qs[2]=R*quad(lambda y: -(np.arcsin(((2*R*y-y**2)**0.5)/R))*(R-(((R**4)-(R**2)*(R-y**2))**0.5)),0,R-((((R**4)-(R**2))**0.5)/R))[0]
	# qs[3]=-R*((R**2)-R*(((R**2)-1)**0.5))*((-18*(R**2)-44+18*R*((R**2)-1))**0.5)*(53/12)*(np.arctan(3*(R-(((R**2)-1)**0.5))/((18*R*(((R**2)-1)**0.5)-44)**0.5)))
	# qs[4]=(R**2)*(quad(lambda y: (((R**4)-(R**2)*(R-y**2))**0.5)*np.arcsin((((R**4)-(R**2)*(R-y**2))**0.5)),0,((((R**4)-(R**2))**0.5)/R)))[0]

	return qs[2]#+qs[3]+qs[4]+qs[0]+qs[1]

def analytic_function(a,b,R,q,num_points,dr=1,bool=True):

	area=a*b
	rho_zero=num_points/area

	if (bool):

		q=quads(R)
	# print('shape a',a.shape)
	# print('a',a)

	# print('b',b)

	# print('type a',type(a))
	# print('type b',type(b))
	# print('type R',type(R))
	# print('type q',type(q))
	# print('type num_points',type(num_points))
	# print('type dr',type(dr))

	# print('dr=',dr)
	# print('R=',R)

	return (1/area)*((-4*(a+b)*(R**2))+2*np.pi*area*R+((np.pi*q))*R**3)*dr*rho_zero*num_points

def get_max_ray(function_array,num_nuclei):

	rays=[]

	for i in range(0,num_nuclei):

		rif_x=function_array[i][0]
		rif_y=function_array[i][1]

		for j in range(0,num_nuclei):

			rays.append(((rif_x-function_array[j][0])**2+(rif_y-function_array[j][1])**2)**0.5)

	rays=np.asarray(rays)
	return(np.amax(rays))					#massimo raggio

def get_distances_vector(function_array,num_nuclei,name):

	h=1
	distances_vector=[]

	for i in range(0,num_nuclei):

		print(i,'/',num_nuclei, name)

		for j in range(h,num_nuclei):

			distances_vector.append(get_distance(function_array[i][0],function_array[i][1],function_array[j][0],function_array[j][1]))

		h+=1

	distances_vector=np.asarray(distances_vector)

	a=np.zeros((1,len(distances_vector)))
	a[0,:]=distances_vector

	return a

def get_distances_and_angles(function_array,num_nuclei,name):
#dato un certo num_nuclei di punti, ritorna la matrice (simmetrica) delle distanze relative tra questi

	h=1														#variabile per prendere la diagonale
	distances_vector=[]										#vettore distanze
	cos_vector=[]										#vettore aspect ratios

	for i in range(0,num_nuclei):

		print(i,'/',num_nuclei, name)

		for j in range(h,num_nuclei):

			distances_vector.append(get_distance(function_array[i][0],function_array[i][1],function_array[j][0],function_array[j][1]))
			angle=math.radians(abs(function_array[i][2]-function_array[j][2]))
			cos_vector.append(abs(np.cos(angle)))

		h+=1

	distances_vector,cos_vector=np.asarray(distances_vector),np.asarray(cos_vector)

	a=np.zeros((2,len(distances_vector)))
	a[0,:]=distances_vector
	a[1,:]=cos_vector

	return a

def get_distances_and_angles_with_AR(function_array,num_nuclei,name):
#dato un certo num_nuclei di punti, ritorna la matrice (simmetrica) delle distanze relative tra questi

	h=1														#variabile per prendere la diagonale
	distances_vector=[]										#vettore distanze
	cos_vector=[]										#vettore aspect ratios
	ar_prod_vec=[]
	just_nem_parameter=[]

	for i in range(0,num_nuclei):

		print(i,'/',num_nuclei, name)

		for j in range(h,num_nuclei):

			distances_vector.append(get_distance(function_array[i][0],function_array[i][1],function_array[j][0],function_array[j][1]))
			angle=math.radians(abs(function_array[i][2]-function_array[j][2]))
			magnitude=function_array[i][3]*function_array[j][3]
			ar_prod_vec.append(magnitude)
			cos_vector.append(abs(np.cos(angle))*magnitude)
			just_nem_parameter.append(abs(np.cos(angle)))

		h+=1

	distances_vector,cos_vector,ar_prod_vec,just_nem_parameter=np.asarray(distances_vector),np.asarray(cos_vector),np.asarray(ar_prod_vec),np.asarray(just_nem_parameter)

	a=np.zeros((4,len(distances_vector)))
	a[0,:]=distances_vector
	a[1,:]=cos_vector
	a[2,:]=ar_prod_vec
	a[3,:]=just_nem_parameter


	return a,np.average(ar_prod_vec)

def get_distances_and_quantities_vectors(function_array,num_nuclei,name,simple=False):
#dato un certo num_nuclei di punti, ritorna la matrice (simmetrica) delle distanze relative tra questi

	h=1														#variabile per prendere la diagonale
	distances_vector=[]										#vettore distanze
	value_vector_ar=[]										#vettore aspect ratios
	value_vector_circ=[]									#vettore circularity
	value_vector_perim=[]
	value_vector_area=[]
	value_vector_round=[]

	for i in range(0,num_nuclei):

		print(i,'/',num_nuclei, name)

		for j in range(h,num_nuclei):

			distances_vector.append(get_distance(function_array[i][0],function_array[i][1],function_array[j][0],function_array[j][1]))
			value_vector_ar.append(function_array[i][2]*function_array[j][2])

			if (simple == False):

				value_vector_circ.append(function_array[i][3]*function_array[j][3])
				value_vector_perim.append(function_array[i][4]*function_array[j][4])
				value_vector_area.append(function_array[i][5]*function_array[j][5])
				value_vector_round.append(function_array[i][6]*function_array[j][6])

		h+=1

	distances_vector,value_vector_ar=np.asarray(distances_vector),np.asarray(value_vector_ar)

	if (simple == False):

		value_vector_area,value_vector_circ=np.asarray(value_vector_area),np.asarray(value_vector_circ)
		value_vector_perim,value_vector_round=np.asarray(value_vector_perim),np.asarray(value_vector_round)

	a=np.zeros((6,len(distances_vector)))
	a[0,:]=distances_vector
	a[1,:]=value_vector_ar

	if (simple == False):

		a[2,:]=value_vector_perim
		a[3,:]=value_vector_circ
		a[4,:]=value_vector_area
		a[5,:]=value_vector_round

	return a

def corr_function(function_array_ars,num_nuclei,ray,distances_matrix,dr):

	# start_time = time.time()

	tot_corr=0								#correlazione totale
	counts_in_r=[]

	for i in range(0,num_nuclei):			#somme su n, il punto su cui calcolo il raggio

		rif_ar=function_array_ars[i]		#prendo l'ar del punto di rif
		corr_val_in_r_for_xi=0				#dato un certo riferimento, la somma dei prodotti in r
		nuc_number_in_r_for_xi=0			#il numero totale di vicini al raggio r incontrati

		for j in range(0, num_nuclei):	#chekko se gli altri nuclei sono alla distanza giusta

			if (distances_matrix[i][j]>=ray) and (distances_matrix[i][j]<=ray+dr)	:		#se il raggio sta in r+/-dr

				corr_val_in_r_for_xi+=(rif_ar*function_array_ars[j])						#per il singolo xi, il prodotto tra le due funzioni
				nuc_number_in_r_for_xi+=1													#i nuclei che vanno bene per il xi

		counts_in_r.append(num_nuclei)

		if(nuc_number_in_r_for_xi!=0):														#se ci sono dei nuclei per l'xi

			tot_corr+=corr_val_in_r_for_xi/nuc_number_in_r_for_xi							#contributo alla corr per xi

	# print("--- %s seconds ---" % (time.time() - start_time))

	plt.hist(counts_in_r,bins=10)

	return(tot_corr/num_nuclei)	

def rois_detections(nuclei,nuclei_val, nuclei_df, membranes_df, rois,rois_val):  #cambia la lista delle rois e ritorna la struttura labels

#num_nuclei: numero dei nuclei segmentati
#nuclei_df: dataframe dei nuclei segmentati(stardist)
#membranes_df: dataframe dall'elaborazione sulle membrane
#rois is an empty list that will be filled with a region of interest in each element. Every ROI contains different nuclei found to be contained there

	temp_list=[]
	
	temp_df=membranes_df.drop_duplicates(subset = ["Label"])	#ho fatto collassare una colonna in modo da avere tutti i label
	
	for i in range(0,len(temp_df.axes[0])):
	
		temp_list.append(temp_df['Label'].iloc[i])
	
	number_of_rois=len(temp_list)
	labels=np.empty((number_of_rois,3), dtype=object)
	
	for i in range(0, number_of_rois):		#labels ha l'etichetta della singola roi + il numero di punti+ le coordinate
	
		labels[i][0]=temp_list[i]
		labels[i][1]=(membranes_df['Label'].value_counts())[temp_list[i]]
		temp_df=membranes_df.loc[membranes_df['Label'] == temp_list[i]]
		borders=np.empty((len(temp_df.axes[0]),2))
		
		for k in range (0, len(temp_df.axes[0])):
		
			borders[k][0]=temp_df['X'].iloc[k]
			borders[k][1]=temp_df['Y'].iloc[k]
		
		labels[i][2]=borders

	for i in range(0,number_of_rois):				#creo un poligono per ogni roi e vedo che nuclei ho dentro
	
		p=Path(labels[i][2])					#capisco quali nuclei sono dentro i punti membrana
		grid = p.contains_points(nuclei)
		single_roi=[]
		val_on_single_roi=[]
		
		for k in range (0,len(grid)):
		
			if grid[k]==True:
			
				single_roi.append(nuclei[k])
				#val_on_single_roi.append(nuclei_val[k])
		
		rois.append(single_roi)				#rois list of different roi, label[i][2] array of different membranes
		rois_val.append(val_on_single_roi)

	return labels
	
def voronoi_on_single_roi(labels,vor,rois_restricted,membranes_restricted,rois):

#labels: struttura calcolata con rois_detection function
#liste che vado a riempire con la funzione: vor conterrà tutti gli oggetti voronoi di scipy.spatial, rois_restricted contiene le roi su #cui effettivamente si riesce a fare un voronoi e membranes restricted le membrane che racchiudono i nuclei su cui effettivamente si fa #un voronoi
#rois: list of regions of interest found with the rois_detection function
#temporary_list contiene al primo posto il voronoi, al secondo il numero della regione, al terzo gli indici della regione e al quarto il punto, poi si ripetono al numero di regione

	for i in range(0,len(rois)):
		
		if len(rois[i])>=4:				#ogni rois[i] ha dentro una collezione di nuclei dentro la membrana
		
			temporary_list=[]
			vor.append(Voronoi(rois[i]))		#voronoi per ogni roi
			temporary_list.append(Voronoi(rois[i]))
			temp=np.zeros((len(rois[i]),2))
			for k in range(0,len(rois[i])):
				
				temp[k][0]=rois[i][k][0]
				temp[k][1]=rois[i][k][1]
	
			rois_restricted.append(temp)
			membranes_restricted.append(labels[i][2]) #mi salvo le membrane che contano
			# for k, reg in enumerate(vor[i].regions):

			# 	print ('Region:', k)
			# 	print('Indices of vertices of Voronoi region:',reg)
			# 	temporary_list.append(k)
			# 	temporary_list.append(reg)
			# 	print ('Associated point:', np.where(vor[i].point_region == i)[0][0], '\n')
						
def voronoi_vertices_inside_convex_hull(vor,membranes_restricted):

#liste con voroi e membrane segmentate su cui si è fatto il voronoi

	for i in range(0,len(vor)):								#ciclo sulle roi_restricted
	
		temp=[]
		hull = ConvexHull(membranes_restricted[i])					#faccio il convex hull della roi
		
		for k in range(0,membranes_restricted[i][hull.vertices,0].size):		#prendo i punti del convex hull
		
			single_ch_point=np.zeros(2)
			single_ch_point[0]=membranes_restricted[i][hull.vertices,0][k]
			single_ch_point[1]=membranes_restricted[i][hull.vertices,1][k]
			temp.append(single_ch_point)
		
		for simplex in hull.simplices:							#plotto il convex hull della singola roi

			plt.plot(membranes_restricted[i][simplex, 0], membranes_restricted[i][simplex, 1], 'y-',lw=3)
	
		temp=np.asarray(temp)
#		plt.plot(temp[:,0],temp[:,1],marker='o',markerfacecolor='red')			#plotta i punti del convex hull
		p=Path(temp)									#path del covex hull della single roi
		
		list=voronoi_vertices_in_path(p,vor[i])

def voronoi_volumes(vor, list_of_regions_inside_path):

	vol = []
	xs=[]
	ys=[]
	for i in list_of_regions_inside_path:

		i=np.asarray(i)		
		hull=ConvexHull(vor.vertices[i])
		vol.append(hull.volume)	

		xs.append(vor.vertices[i][:,0])
		ys.append(vor.vertices[i][:,1])

		cx = np.mean(vor.vertices[i][:,0])
		cy = np.mean(vor.vertices[i][:,1])

		# xs.append(cx)
		# ys.append(cy)

		# print(cx,cy)

	vol=np.asarray(vol)
	# xs=np.asarray(xs)
	# ys=np.asarray(ys)
	
	# print(vol.shape)
	# print(xs.shape)

	return vol,xs,ys

def voronoi_vertices_inside_membrane(vor,membranes_restricted,ax=0,disegno=False):

	print(ax)

	vol = []
	xs=[]
	ys=[]								#ogni elemento è la lista dei volumi per una sola roi
	nuke_x,nuke_y=[],[]

	for i in range(0,len(vor)):								#ciclo sulle roi_restricted

		p=Path(membranes_restricted[i])
		list,temp_x,temp_y=voronoi_vertices_in_path(p,vor[i],ax=ax,disegno=disegno)
		
		if (disegno):
			ax.plot(membranes_restricted[i][:,0],membranes_restricted[i][:,1],marker='o',markerfacecolor='red')
		temp=voronoi_volumes(vor=vor[i],list_of_regions_inside_path=list)

		temp_x,temp_y=np.asarray(temp_x),np.asarray(temp_y)
		nuke_x.append(temp_x)
		nuke_y.append(temp_y)
		vol.append(temp[0])
		xs.append(temp[1])
		ys.append(temp[2])

	return(vol,xs,ys,nuke_x,nuke_y)

def voronoi_cell_aspect_ratio(vor,list_of_regions_inside_path):

	for i in list_of_regions_inside_path:

		i=np.asarray(i)		
		hull=ConvexHull(vor.vertices[i])						#convex hull della singola cella di voronoi
		points[hull.vertices[0],0], points[hull.vertices[0],1], 'ro'
		eigvals, eigvecs = np.linalg.eig(np.cov(xy))
	return vol

def NN_through_Voronoi(df):

	num_nuclei=len(df.axes[0])
	points=np.zeros((num_nuclei,2))

	for i in range(0,num_nuclei):					#riempio i punti

		points[i,0]=df['XM'].iloc[i]
		points[i,1]=df['YM'].iloc[i]

	vor=Voronoi(points)								#faccio un voronoi

	NNs=[]

	for i in range(0,num_nuclei):

		current_point=points[i,:]
		region_for_current_point=vor.point_region[i]
		indexis_of_vertices_for_current_region=(vor.regions[vor.point_region[i]])

		#remove duplicate points at infinity:

		if (-1 in indexis_of_vertices_for_current_region):

			try:
				while (True):
		   	 		indexis_of_vertices_for_current_region.remove(-1)
			except ValueError:
	   			pass

			indexis_of_vertices_for_current_region.append(-1)

		indexis_of_vertices_for_current_region=set(indexis_of_vertices_for_current_region)							#a questo punto lui è pronto per essere confrontato

		NN=[]

		for eta in range(0,num_nuclei):																				#ciclo su tutte le altre voronoi cells

			if(i!=eta):
				
				comparer=vor.regions[vor.point_region[eta]]
				if (-1 in indexis_of_vertices_for_current_region):
																		#tolgo ripetizioni infinito
					try:
						while (True):
		   	 				comparer.remove(-1)
					except ValueError:
	   					pass
					comparer.append(-1)
			
				comparer=set(comparer)

				if(len(set.intersection(comparer,indexis_of_vertices_for_current_region))>=2):

					NN.append(eta)																					#lui è NN per i

		NNs.append(NN)										
	
	return NNs																										#ogni elemento è una lista: alla posizione i è presente una lista dei NN all'elemento i della df
