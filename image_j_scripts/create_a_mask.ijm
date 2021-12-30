number_of_images=922

for (j=0; j<number_of_images; j++){

	roiManager("Select", j+1);
	run("Create Mask");

	}