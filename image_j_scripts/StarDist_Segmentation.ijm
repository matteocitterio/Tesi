//script used to automatically segment the images through Fiji and StardDist. StraDist needs to be properly install in Fiji as its plugin as well as CSBdeep library.
//Images are here considered to be H&E saved as '.tif' extension
// The rescaling (***) is not mandatory: sometimes StarDist has an hardtime in segmenting objects which are too big: by scaling down the image this problem is most of the times solved
//Since normally this kind of images are quite big, StarDist here works by dividing the image into small pieces (called tiles) in order not to overload the GPU. Here the number of tiles is set to six but it could be changed.
//Results of the measurements are then stored as a .csv file. This files contains all the quantities selected in 'Set_measurements' in the Fiji pop-up menu. This could be changed as well.

//PARAMETERS
//number_of_images= Number of images to analyze [TYPE INT]
//folder_path= name of the folder containing the images [TYPE STRING]
//path_processed= pathname in which the results of the segmentation will be saved [TYPE STRING]
//path_image= pathname (which contains folder path) of the images which are going to be analyzed [TYPE STRING]


number_of_images= ...;

for (j=0; j<number_of_images; j++){					//iteration over the N images
	print('This iteration is: ', j);

	roiManager("reset");

	title=toString(j+1)+'.tif';
	folder_path= ...						
	path_processed= ... + folder_path+ "\\";
	path_image= ... +folder_path+ "\\";

	print(title);

	open(path_image+title);						//opening the image
	selectWindow(title);
	close("\\Others");
	img_title = getTitle();
	selectImage(title);
	rename("nuc");
	nuc_ID = getImageID();
	selectWindow("nuc");
	
	//run("Scale...", "x=0.3 y=0.3 interpolation=Bilinear average create title=downScale");   //*** NOT MANDATORY ***
	run("Command From Macro", "command=[de.csbdresden.stardist.StarDist2D], args=['input':'downScale', 'modelChoice':'Versatile (H&E nuclei)', 'normalizeInput':'true', 'percentileBottom':'1.0', 'percentileTop':'99.8', 'probThresh':'0.692478', 'nmsThresh':'0.3', 'outputType':'Both', 'nTiles':'6', 'excludeBoundary':'2', 'roiPosition':'Automatic', 'verbose':'false', 'showCsbdeepProgress':'true', 'showProbAndDist':'false'], process=[false]");		//StarDist command, this might take a while

	roiManager("Measure");						//measurements of the segmented objects
	selectWindow("Results");  
	saveAs("Results", path_processed+"Results_"+toString(j+1)+".csv");
	selectWindow("Results");
	run("Clear Results");
	close();
	roiManager("reset");

	}

print('Done');
