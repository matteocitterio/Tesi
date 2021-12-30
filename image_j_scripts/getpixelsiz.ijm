
number_of_images=8;
print('new cycle     ')
new_path="Cropped_13_Tumor";						//da cambiare quando guardi altre immagini da una altra cartella
path_processed="D:\\matteo_citterio\\risultati_segment\\"+ new_path+ "\\";
path_image="C:\\Users\\stefa\\Documents\\MEGAsync Downloads\\"+new_path+ "\\";
resultFilename = "D:\\matteo_citterio\\risultati_segment\\"+new_path+"pixel_sizes.csv";
f = File.open(resultFilename);
for (j=0; j<number_of_images; j++){
   	//wait(0.5); {just a delay to see the answer in the Info window}
	
	roiManager("reset");
	pixel_size=0.1214;
	hello=j+1;
	title_temp=toString(hello);
	save_temp=toString(j);
	title=title_temp+'.tif';

	//print(title);

	open(path_image+title);
	selectWindow(title);
	close("\\Others");
	img_title = getTitle();
	getPixelSize(unit, pw, ph, pd);
	print(unit,pw,ph,pd);
	print(f,d2s(pw,6));


// you have made an RGB with gray and magenta,
// which is a bad idea because now hte cell staining is in all the channels
// try to keep independant colors (Red, Green, Blue)
// here is a trick

	
	selectWindow(title);
	close();

	//selectImage(img_title);
	//roiManager("Show All");
}
File.close(f);
print('finito');
