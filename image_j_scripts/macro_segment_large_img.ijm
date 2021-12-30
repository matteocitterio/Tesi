
number_of_images=7;

for (j=3; j<number_of_images; j++){
	print('This iteration is: ', j);
   	//wait(0.5); {just a delay to see the answer in the Info window}

	roiManager("reset");

	hello=j+1;
	title_temp=toString(hello);
	save_temp=toString(j);
	title=title_temp+'.tif';
	new_path="ROI_13_INT_Tumor";						//da cambiare quando guardi altre immagini da una altra cartella
	//old_path_image="infiltrating";
	path_processed="D:\\matteo_citterio\\risultati_segment\\"+ new_path+ "\\";
	path_image="D:\\matteo_citterio\\SCAN - Copia\\"+new_path+ "\\";

	print(title);

	open(path_image+title);
	selectWindow(title);
	close("\\Others");
	img_title = getTitle();
	selectImage(title);
	rename("nuc");
	nuc_ID = getImageID();
	selectWindow("nuc");
	
	run("Scale...", "x=0.3 y=0.3 interpolation=Bilinear average create title=downScale");
	run("Command From Macro", "command=[de.csbdresden.stardist.StarDist2D], args=['input':'downScale', 'modelChoice':'Versatile (H&E nuclei)', 'normalizeInput':'true', 'percentileBottom':'1.0', 'percentileTop':'99.8', 'probThresh':'0.692478', 'nmsThresh':'0.3', 'outputType':'Both', 'nTiles':'6', 'excludeBoundary':'2', 'roiPosition':'Automatic', 'verbose':'false', 'showCsbdeepProgress':'true', 'showProbAndDist':'false'], process=[false]");
	//selectWindow("Label Image");
	//close();
	roiManager("Measure");
	selectWindow("Results");  
	saveAs("Results", path_processed+"Results_"+hello+".csv");
	selectWindow("Results");
	run("Clear Results");
	close();
	roiManager("reset");

// you have made an RGB with gray and magenta,
// which is a bad idea because now hte cell staining is in all the channels
// try to keep independant colors (Red, Green, Blue)
// here is a trick
//	run("Duplicate...", "title=["+img_title+"]");
	//run("Split Channels");

	//if(new_path=="infiltrating"){

	//	imageCalculator("Subtract create 32-bit", title_temp+".tif (blue)",title_temp+".tif (red)");
		//selectImage("Result of "+title_temp+".tif (blue)");
	}
	
	

	//selectImage(title+" (green)");
	//rename("cell");
	//cell_ID = getImageID();
	//selectImage(cell_ID);

print('finito');
