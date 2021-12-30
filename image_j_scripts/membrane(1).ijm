
number_of_images=8;

for (j=0; j<number_of_images; j++){
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

// you have made an RGB with gray and magenta,
// which is a bad idea because now hte cell staining is in all the channels
// try to keep independant colors (Red, Green, Blue)
// here is a trick
	run("Duplicate...", "title=["+img_title+"]");
	run("Split Channels");

	selectImage(title+" (green)");
	rename("green");
	//cell_ID = getImageID();
	selectImage(title+" (red)");
	rename("red");
	//cell_ID = getImageID();
	selectImage(title+" (blue)");
	rename("blue");
	//cell_ID = getImageID();

	selectWindow("green");
	//run("Gaussian Blur...", "sigma=5");
	setThreshold(0, 225);
	run("Convert to Mask");
	selectWindow("red");
	setThreshold(0, 225);
	run("Convert to Mask");
	selectWindow("blue");
	setThreshold(0, 225);
	run("Convert to Mask");
	imageCalculator("Add create","red","green");
	selectWindow("Result of red");
	imageCalculator("Add create", "Result of red","blue");
	selectWindow("Result of Result of red");
	
//setAutoThreshold("Default dark");
	//run("Convert to Mask");



// substract the Voronoi to the Cell mask
//imageCalculator("Subtract create", "cell", voronoi_img);

//select cell_mask and detect cells 
//run("Analyze Particles...", "add");

	setBatchMode(true);
//run("Invert LUT");
//setAutoThreshold("Default dark");
//setOption("BlackBackground", false);
//run("Convert to Mask");
	str = split(getTitle(), ".");
	str = str[0]+"_roi-";
	run("Analyze Particles...", "size=250-Infinity show=Overlay clear add");
	roiManager("Measure");
	run("Clear Results");
	n = roiManager("count");
	for ( i=0; i<n; i++ ) { 
		roiManager("select", i);
		outline2results(str+(i+1));
	}

	setBatchMode(false);

	selectWindow("Results");
	saveAs("Results", path_processed+"contour"+save_temp+".csv");
	
	//exit();
	function outline2results(lbl) {
		nR = nResults;
		Roi.getCoordinates(x, y);
		for (i=0; i<x.length; i++) {
			setResult("Label", i+nR, lbl);
			setResult("X", i+nR, x[i]);
			setResult("Y", i+nR, y[i]);
		}
	}

	//selectWindow(title);
	//roiManager("Combine");

	selectWindow("blue");
	close();
	selectWindow("red");
	close();
	selectWindow("green");
	close();
	selectWindow(title);
	close();
	selectWindow("Result of red");
	close();
	selectWindow("Result of Result of red");
	close();
	//selectImage(img_title);
	//roiManager("Show All");
}

print('finito');
