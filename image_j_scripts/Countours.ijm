//This script segments an image using a simple binary threshold (0-254) and the countors of the detected object are stored as (x,y) coordinates. These last are used in the normalization of the pair correlation function g(r)
//The code takes in input a '.tif' H&E image saved as 'X.tif' where X is an integer number

//PARAMETERS
//number_of_images= number of images to analyze [TYPE INT]
//folder_name= name of the folder containing the considered images [TYPE STRING]
//path_processed= pathname of the folder where the results will be stored [TYPE STRING]
//path_image= pathname of the folder folder_name

number_of_images= ... ;

for (j=0; j<number_of_images; j++){
	print('This iteration is: ', j);

	roiManager("reset");

	title=toString(j+1)+'.tif';
	folder_name= ... ;						
	path_processed= ... + folder_name+ "\\";
	path_image= ... +folder_name+ "\\";

	print(title);

	open(path_image+title);							//this simply opens the image and says to Fiji to use the opened image
	selectWindow(title);
	close("\\Others");
	img_title = getTitle();

	run("Duplicate...", "title=["+img_title+"]");				//this duplicate the image and splits it into the RGB channels
	run("Split Channels");

	selectImage(title+" (green)");						//every channel is renamed as 'color' 
	rename("green");
	selectImage(title+" (red)");
	rename("red");
	selectImage(title+" (blue)");
	rename("blue");

	selectWindow("green");							//every channel is thresholded. Depending on the particular set of images
	setThreshold(0, 225);							//a different threshold should be applied
	run("Convert to Mask");
	
	selectWindow("red");
	setThreshold(0, 225);
	run("Convert to Mask");
	
	selectWindow("blue");
	setThreshold(0, 225);
	run("Convert to Mask");
	
	imageCalculator("Add create","red","green");				//the segmented images are then summed together
	selectWindow("Result of red");
	imageCalculator("Add create", "Result of red","blue");
	selectWindow("Result of Result of red");
	

	setBatchMode(true);

	str = split(getTitle(), ".");
	str = str[0]+"_roi-";
	run("Analyze Particles...", "size=250-Infinity show=Overlay clear add"); //this measures the object
	roiManager("Measure");
	run("Clear Results");
	n = roiManager("count");
	for ( i=0; i<n; i++ ) { 
		roiManager("select", i);
		outline2results(str+(i+1));
	}

	setBatchMode(false);

	selectWindow("Results");
	saveAs("Results", path_processed+"contour"+toString(j)+".csv");
	
	function outline2results(lbl) {						//function used to save the segmentation results as (x,y) coordinates
		nR = nResults;
		Roi.getCoordinates(x, y);
		for (i=0; i<x.length; i++) {
			setResult("Label", i+nR, lbl);
			setResult("X", i+nR, x[i]);
			setResult("Y", i+nR, y[i]);
		}
	}

	selectWindow("blue");							//closing all the windows
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

}

print('Done');
