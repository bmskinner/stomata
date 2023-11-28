// Test stomata detection using FFT bandpass filtering

function detectStomata(imageId) {
	run("Clear Results");
	selectImage(imageId);
	run("Duplicate...", "title=original");
	selectImage(imageId);
	roiManager("reset");

	run("Bandpass Filter...", "filter_large=100 filter_small=50 suppress=None tolerance=5 autoscale saturate");
	setAutoThreshold("Default");

	setThreshold(0, 146, "raw");
	setThreshold(60, 255);
	setOption("BlackBackground", true);
	run("Convert to Mask");
	run("Convert to Mask");
	run("Analyze Particles...", "size=500-10000 circularity=0.50-1.00 show=Nothing display exclude add");

	for (i=0; i<nResults; i++) {
		stomata_x = getResult("XM", i); // centre of mass
		stomata_y = getResult("YM", i);
		print(outFileHandle, imageName+"\t"+stomata_x+"\t"+stomata_y);	
	}

	run("Clear Results");
	selectImage("original");
	run("From ROI Manager");
	run("Flatten");
	saveAs("tiff", outDir+"/"+imageId+".tiff");
	close("original");
	close(imageId);
}

function analyseFromDirectory(dir) {

	resetLog();

	var outDir = dir +"//output"; // global to allow for tiff saving within analysis function
	if(!File.exists(outDir)){
		File.makeDirectory(outDir);
	}

	outFile = outDir + "/CoMs.txt";
	if(File.exists(outFile)){
		File.delete(outFile);
	}

	var outFileHandle = File.open(outFile);
	print(outFileHandle, "image\tx\ty");	
	list = getFileList(dir);
	count = 1;

	for (i=0; i<list.length; i++) {

		if (
		  endsWith(list[i], ".jpg") ||
		  endsWith(list[i], ".tif") ||
		  endsWith(list[i], ".tiff")
		) {
			print("image " + (count++) + ": " + dir + list[i]);
			open(dir + list[i]);
			imageId = getImageID();
			imageName = getTitle();
			detectStomata(imageName);

			close(imageName);
			close("Drawing of "+imageName);
		} 
	}
}

function resetLog() {
	if (isOpen("Log")) {
		selectWindow("Log");
		run("Close");
	}
}



macro "analyseFromDirectory"

{
	requires("1.47");
	setBatchMode(true)
	saveSettings();
	run("Colors...", "background=black"); //set background colour to black
	dir = getDirectory("Choose a Directory ");
	analyseFromDirectory(dir);
	restoreSettings();
	setBatchMode(false);
}