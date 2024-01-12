//setup









//Channels must be 2 for MYH and 4 for DAPI








dir = getDir("Choose a Directory"); //user chooses directory with images in individual folders 
dirOut = getDirectory("Choose a Directory"); //user chooses directory to save the results in a new folder called "counts"
print(dir);

File.makeDirectory(dirOut +"/counts");
newPathL = dirOut + "/counts/";
run("Set Measurements...", "area mean shape redirect=None decimal=3");


subdir = getFileList(dir);
Array.show(subdir);

//check if ROi manager already has old data. Delete if it does
if (roiManager("Count") > 0){
	roiManager("Deselect");
	roiManager("Delete");
	run("Clear Results");
		close("*");
	}
run("Close All");

for (i = 0; i < subdir.length; i++){ //loop should iterate from 0:subdir.length

	// To add more channels: 
	//file# = dir+subdir[i]+"Image_CH#.tif"; 
	//open(file#);
	localdir = dir+subdir[i];
	print(localdir);
	dir1 = getFileList(localdir);
	for (j = 0; j < dir1.length; j++){
		fileName = substring(dir1[j], 0, dir1[j].length-5);
		print("File name: " + fileName);

		if (endsWith(dir1[j], "CH2.tif")){
			file1 = localdir+fileName+"2.tif";
			print(file1);
			open(file1);
		}

		if (endsWith(dir1[j], "CH4.tif")){
			file1 = localdir+fileName+"4.tif";
			print(file1);
			open(file1);
		}
		
    	
	
	}
	analyzeChannels(newPathL + "results.xlsx", dirOut, subdir[i]); // Edit this file name
		print(i+" is done");
		
			
					
}

function analyzeChannels(output, outLoc, filepath) {
	
	// Make a mask from DAPI channel (Ch4):
	// 		Note: adjust file number if needed 
	names = substring(dir1[0], 0, dir1[0].length-5);
	print(dir);
	print(subdir[i]);
	folderName = replace(subdir[i], "/", "");
	folderNameExcel = substring(subdir[i], 0, 10);
	selectWindow(names+"4.tif");
	run("Duplicate...", "duplicate");
	run("Split Channels");
	selectWindow("C2-"+names+"4-1.tif");
	close();
	selectWindow("C1-"+names+"4-1.tif");
	close();
	selectWindow("C3-"+names+"4-1.tif");
	run("Duplicate...", "duplicate");
	saveAs("PNG", newPathL + folderName + "Nuclei");
	close();
	run("8-bit");
	
	//mask and analyze the nuclei
	setAutoThreshold("Intermodes dark no-reset");
	setOption("BlackBackground", true);
	run("Convert to Mask");
	run("Watershed");
	run("Analyze Particles...", "size=40-Infinity pixel circularity=0.40-1.00 show=Outlines add");
	close();
	if (roiManager("Count") > 0){
		roiManager("Show None");
		roiManager("Show All");
		roiManager("Measure");
	
	
		run("Read and Write Excel","file=[" + newPathL + folderNameExcel+ "results.xlsx] sheet="+subdir[i]); //stack_results
	// Clears ROI Manager before making new mask: 
		roiManager ("Delete");
		run("Clear Results");
	}
		selectWindow("C3-"+names+"4-1.tif");
	run("Duplicate...", "duplicate");
	saveAs("PNG", newPathL + folderName + "Total nuclei");
	close();
	//MAKE myhc MASK from Ch2
	selectWindow(names+"2.tif");
	run("Duplicate...", "duplicate");
	run("Split Channels");
	selectWindow("C3-"+names+"2-1.tif");
	close();
	selectWindow("C1-"+names+"2-1.tif");
	close();
	selectWindow("C2-"+names+"2-1.tif");
	run("Duplicate...", "duplicate");
	saveAs("PNG", newPathL + folderName + "MYHC");
	close();
	run("8-bit");
	setAutoThreshold("Triangle dark no-reset");
	run("Convert to Mask");
	run("Fill Holes");
	run("Duplicate...", "duplicate");
	saveAs("PNG", newPathL + folderName + "MYHC mask");
	close();
	
	//Make double mask		
	print("File name: " + names);
		imageCalculator("AND create", "C2-"+names+"2-1.tif","C3-"+ names + "4-1.tif");
	selectWindow("Result of C2-" +names +"2-1.tif");
	run("Duplicate...", "duplicate");
	saveAs("PNG", newPathL + folderName + "combined mask");
	close();
	run("Analyze Particles...", "size=30-Infinity pixel circularity=0.40-1.00 show=Outlines add");
	if (roiManager("Count") > 0){
		roiManager("Show None");
		roiManager("Show All");
		roiManager("Measure");
		run("Read and Write Excel","file=[" + newPathL + folderNameExcel+ "results.xlsx] sheet="+subdir[i]); //stack_results
	// Clears ROI Manager before making new mask: 
		roiManager ("Delete");
	}
	run("Clear Results");
	
	close("*");
}




