
//script used to quickly extract the pixel size out of a folder of imagers in ImageJ. Images needed to be '.tif' in order to be opened

//PARAMETERS:
//N_images= number of images inside the considered folder [TYPE INT]
//Folder_path= path of the considered folder containing N_images [TYPE STRING]
//Output_pathname= pathname of the output filename (ending in .csv) [TYPE STRING]
//output_folder_name= name of the folder in which the csv file will be saved [TYPE STRING]

number_of_images=N_Images;
print('new cycle     ')
path_image= Folder_path;
resultFilename = Output_pathname+".csv";
f = File.open(resultFilename);
for (j=0; j<number_of_images; j++){
	
	roiManager("reset");
	pixel_size=0.1214;
	title=toString(j+1)+'.tif';			//name of the image which is needed to open

	open(path_image+title);
	selectWindow(title);
	close("\\Others");
	img_title = getTitle();
	getPixelSize(unit, pw, ph, pd);
	print(unit,pw,ph,pd);
	print(f,d2s(pw,6));
	selectWindow(title);
	close();

}
File.close(f);
print('Done');
