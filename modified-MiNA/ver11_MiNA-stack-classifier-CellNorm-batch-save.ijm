////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//      MACRO: Mitochondrial Network Analysis (MiNA) - V6.0           //
//      AUTHOR: Andrew Valente     & Kurt Weiss (BOC)                          //
//      EMAIL: valentaj94@gmail.com                                           //
//      LAST EDITED: October 6th,2017(AV), June 28,2024 (KRW)                         //
////////////////////////////////////////////////////////////////////////////////

//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//		current runs on nd2, cna change to .tif in if statement
//		the 'summary' shows the cell area (approximate), after removing nuclei

//Instructions
//-NOTE RUNS WELL ON MAX PROJ OF 3 CONFOCAL SlICES 200NM APART, 100x, deconvolution optional, mito561, nuclei-405>>This V11 has line 113 modified to run-z-projection on substack
//Consider first testing with origianl macro to confirm improved performance of this more complicated macro. To do this, assign weka=false in line 482 
//1. run the Train-data-for mina macro - this opens .nd2, does max projection, does any preprocessing steps on the mito(561) channel, 
//   saves a tif.
//2. Save all those .tif files into folder and drag folder into fiji which should open as a stack
//3. Run Plugins>segmentation>weka segmenation (2D)
//4. Train model on many images and save model as classifier1.model (or update name within the Mina-v4-classifier macro, line 115). 
//   Also 'save data' as ARFF to allow later modification of labels. 
//5. Place the .nd2 images you wish to image along with the classifier1.model and the Mina-v7 macro in same folder. 
//6. Open and run the macro. 
//7. Check output files for accruacy: for mito-mask, modify training data and thresholding; for cell-nuc-mask modify crudeSegment(); for overlay modify analyzeSkeleton() (line 195)*However
// analyze skeleton is largely built-in Fiji with few parameters, better to modify the training data and thresholding to get accurate mito. 
//-Ensure any preprocessing is done to both the training and testing images.
//-Ensure correct model is referenced in line 115. 
//-Macro optionally pauses with Red/white threshold image - use ctrl+atl+t to open the threshold and you can toggle between 0.1-0.9, 
//    currently set at default 0.5 (see line 490 to change). Alternativly retrain model to avoid ambiguous regions. 
//-Note optional features can be included in the random forest during weka training>advanced when training. 
//TODO
// make separate crudeSegment() functions for nuclei and cytoplasm

//Global Variables--------------------------------------------------------------
//Output Arrays...
var filenameARRAY = newArray();
var individualsARRAY = newArray();
var networksARRAY = newArray();
var meanLengthARRAY = newArray();
var medianLengthARRAY = newArray();
var sdLengthARRAY = newArray();
var meanBranchesARRAY = newArray();
var medianBranchesARRAY = newArray();
var sdBranchesARRAY = newArray();
var mitoAreaARRAY = newArray();

//Preprocessing Modifiers
var CLAHE = false;
var CLAHE_BLOCKSIZE = 127;
var CLAHE_HISTOGRAM = 256;
var CLAHE_MAXSLOPE = 3;
var MED = true;
var MED_RADIUS = 2;
var TOPHAT = false;
var UNSHARP = true;
var UNSHARP_RADIUS = 2;
var UNSHARP_STRENGTH = 0.6;

	
	////////////////preprocess/batch open	
   dir = getDirectory("Choose a Directory ");
   setUp();
   //setBatchMode(true);
   count = 0;
   countFiles(dir);
  // File.makeDirectory(dir+ "_output");
   n = 0;
   processFiles(dir);
   function countFiles(dir) {
      list = getFileList(dir);
      for (i=0; i<list.length; i++) {
          if (endsWith(list[i], "/")){
              countFiles(""+dir+list[i]);
          }
          else{
              count++;
     		 }
      }
  }

   function processFiles(dir) {
      list = getFileList(dir);
      for (i=0; i<list.length; i++) {//was <list.length but this repeats process for each file
         
          if (endsWith(list[i], "/")){
              processFiles(""+dir+list[i]);   
          }
          else {
             showProgress(n++, count);
             path = dir+list[i];
             processFile(path);
          }
      }
  }

  function processFile(path) {
       if (endsWith(path, ".nd2"))  {  //can change to .tif
       		run("Bio-Formats", "open=["+path+"] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
			imageName = getTitle();
			selectImage(imageName); 
			//run("Z Project...", "projection=[Max Intensity]");
			run("Z Project...", "start=2 stop=4 projection=[Max Intensity]");
			run("Split Channels"); 
			selectImage("C1-MAX_"+imageName);
			rename("dapi"); 
			selectImage("C2-MAX_"+imageName);
			rename("mito");
			//Generate path to model
			path_to_model= dir+"classifier1.model";
			mina(); //run Mina 
      }
  }

	
//mina Main
function mina(){

	//Produce GUI to set preprocessing preferences
	//setUp();//run setup here if all images have different preprocessing parameters
    //Duplicate region of interest and collect general information
    selectWindow("mito");
    showStatus("MiNA: Getting image information...");
	run("Duplicate...", "title=Original");
	selectWindow("Original");
	run("Grays");
	getDimensions(width, height, channels, slices, frames);
	getPixelSize(unit, pixelWidth, pixelHeight);
	
	//Generate nucleus mask, cell mask, subtract, find area
	selectWindow("dapi");
	crudeSegment(); //check output image 'cell-nuc-mask' for accuracy and modify crudeSegment as needed, can add waitforuser() to validate
	//waitForUser;
	selectWindow("Original");
	run("Duplicate...", "title=cells");
	crudeSegment();		//check output image 'cell-nuc-mask' for accuracy and modify crudeSegment as needed, can add waitforuser() to validate
	imageCalculator("Subtract create", "cells","dapi");
	run("Analyze Particles...", "size=100-Infinity pixel summarize");
	selectWindow("Summary"); 
	IJ.renameResults("Summary","Results");
	cellArea=getResult("Total Area", 0);
	
	//continue origianl macro
	selectWindow("Original"); 
	run("Duplicate...", "title=raw");
	selectWindow("Original"); 
	//Preprocess image including weka segmentation
	preprocessing();

	//Produce binary duplicate image for skeleton and binary outline...
	selectWindow("segmented");
	run("Duplicate...", "title=TestSkeleton");
	selectWindow("TestSkeleton");
	run("32-bit");
	run("Make Binary");
	run("8-bit");
    run("Duplicate...", "title=Outline");

    //Calculate the mitochondrial footprint
	getStatistics(area, Mean, min, max);
	mitoArea = pow(pixelWidth, 2.0) * parseFloat(width) * parseFloat(height) * (Mean / parseFloat(max)) ;

    //Skeletonize the binary image and overlay it onto the original
    selectWindow("TestSkeleton");
	run("Skeletonize (2D/3D)");
	run("Green");
	selectWindow("raw");
	run("Add Image...", "image=TestSkeleton x=0 y=0 opacity=100 zero");

    //Create an outline and overlay it...
    selectWindow("Outline");
    run("Find Edges");
    run("Magenta");
    selectWindow("raw");
    run("Add Image...", "image=Outline x=0 y=0 opacity=50 zero");
	//Add a scale bar to it.
	size = toString(round(0.25*parseFloat(width)*parseFloat(pixelWidth)));
	run("Scale Bar...", "width=" + size + " height=4 font=14 color=White background=Black location=[Lower Right] bold overlay");

	//Ask the user if this is acceptable.
	//Dialog.create("Skeleton Proofing");
	//Dialog.addMessage("Is the skeleton produced acceptable? Click OK if so.");
	//Dialog.show()

	//Analyze the skeleton using the Analyze Skeleton plugin.
	showStatus("MiNA: Processing skeleton...");
	selectWindow("TestSkeleton");
	run("Analyze Skeleton (2D/3D)", "prune=none show display");
	
	//save and close: 
	selectImage("raw");
	//rename(imageName); //redundant as imageName is in the 'path' variable?
	saveAs("tiff", path+"-overlay");
	selectImage("Result of cells");
	saveAs("tiff", path+"cell-nuc-mask");
	close();
	selectImage("segmented");
	saveAs("tiff", path+"_mito-mask");
	close(); 
    close("Outline");
	close("TestSkeleton");
	close("TestSkeleton-labeled-skeletons");
	close("cells");
	close("mito");
	close("dapi"); 
	
	
	//Collect relevant output from the tables and close these tables.
	selectWindow("Results"); rows = nResults;
	BranchCounts = newArray(rows);
	for (i=0; i<rows; i++) {
		BranchCounts[i] = getResult("# Branches", i);
	}
	
	//count only skeleton >3 px ie slab voxels>1 (2 endpoint plus 1 slab)
	SlabCounts = newArray(rows);
	for (i=0; i<rows; i++) {
		SlabCounts[i] = getResult("# Slab voxels", i);
	}
	
	selectWindow("Results"); run("Close");

	IJ.renameResults("Branch information", "Results");
	selectWindow("Results"); rows = nResults;
	BranchLengths = newArray(rows);
	for (i=0; i<rows; i++) {
		BranchLengths[i] = getResult("Branch length", i);
	}
	run("Close");

	//Feature counts
	individuals = parseFloat(countIndividuals(BranchCounts));
	networks = parseFloat(countNetworks(BranchCounts));
	slabs = parseFloat(countSlabs(SlabCounts));

	//Size descriptors for lengths
	meanLength = mean(BranchLengths);
	medianLength = median(BranchLengths);
	sdLength = sd(BranchLengths);

	//Strip non networked branches
	 //and count lg mito are linear or circular, zero branches, >1 slab voxel
	networkBranchCounts = newArray(0);
	lgMito = 0;
	for (i=0; i<BranchCounts.length; i++) {
		if((BranchCounts[i]==1)&&(SlabCounts[i]>1)){
			lgMito = lgMito+1; 
		}
		
		if (BranchCounts[i]>1) {
			networkBranchCounts = Array.concat(networkBranchCounts, BranchCounts[i]);
		}
	}

	meanBranches = mean(networkBranchCounts);
	medianBranches = median(networkBranchCounts);
	sdBranches = sd(networkBranchCounts);

	Measurement = newArray("circle or line mito (>3 vox & non-networked)", "branches (>3 voxel, all networks+circle/line)","Individuals (branches>1 voxel)",
			       "Networks (>1 branch, exclusive of circle or line mito)",
			       "Mean Branch Length",
			       "Median Branch Length",
			       "Length Standard Deviation",
			       "Mean Network Size (Branches)",
			       "Median Network Size (Branches)",
			       "Network Size Standard Deviation",
			       "Mitochondrial Footprint", "cell area");

	Value = newArray(lgMito, slabs,individuals,
			 networks,
			 meanLength,
			 medianLength,
			 sdLength,
			 meanBranches,
			 medianBranches,
			 sdBranches,
			 mitoArea, cellArea);

	Units = newArray("Counts","Counts","Counts",
			 "Counts",
			 unit,
			 unit,
			 unit,
			 "Counts",
			 "Counts",
			 "Counts",
			 unit+" squared", unit+" squared");

	Array.show("MiNA Output", Measurement, Value, Units);
	//save
	Table.rename("MiNA Output", "MiNA"+imageName);
	saveAs("Results", path+".csv");
	
	//print("next okay close all, get measuremnts and save");
	//waitForUser;
	close("*");
}

//crude segment for cell and nuclei


//Count Individuals...
function countIndividuals(data) {
	entries = data.length;
	total = 0.0;
	for (i=0; i<entries; i++) {
		if (data[i] <= 1) {  
			total = total + 1;
		}
		else {
		}
	}
	return(total);
}

//Count Networks...
function countNetworks(data) {
	entries = data.length;
	total = 0.0;
	for (i=0; i<entries; i++) {
		if (data[i] > 1) {
			total = total + 1;
		}
		else {
		}
	}
	return(total);
}

//Count Slabs
function countSlabs(data) {
	entries = data.length;
	total = 0.0;
	for (i=0; i<entries; i++) {
		if (data[i] > 1) {
			total = total + 1;
		}
		else {
		}
	}
	return(total);
}

//Count Branches in Networks...
function countNetworkBranches(data) {
	entries = data.length;
	total = 0.0;
	for (i=0; i<entries; i++) {
		if (data[i] > 1) {   
			total = total + data[i];
		}
		else {
		}
	}
	return(total);
}

//Mean...
function mean(data) {
	entries = data.length;
	total = 0.0;
	for (i=0; i<entries; i++) {
		total = total + data[i];
	}
	ave = total/parseFloat(entries);
	return(ave);
}

//Median...
function median(data) {
	entries = data.length;
	sorted = Array.sort(data);
	if (entries%2 == 0) {
		iF = (entries/2);
		med = parseFloat((sorted[(iF-1)] + sorted[iF])/2.0);
	}
	else {
		iF = ((entries-1)/2);
		med = sorted[iF];
	}

	return(med);
}

//Standard deviation...
function sd(data) {
	entries = data.length;
	N = parseFloat(entries);
	u = mean(data);
	sse = 0.0;
	for (i=0; i<entries; i++) {
		sse = sse + (pow((data[i] - u),2)/N);
	}
	std = sqrt(sse);
	return(std);
}

//Preprocessing set up GUI...
function setUp() {
	//Produce GUI to set preprocessing preferences
    Dialog.create("Preprocessing: ensure all steps have been applied to training images");
    Dialog.addMessage("Processing to Apply");
    Dialog.addCheckbox("CLAHE", false);
	Dialog.addSlider("Blocksize: ", 1, 256, 127);
	Dialog.addSlider("Histogram Bins: ", 1, 256, 256);
	Dialog.addSlider("Maximum Slope: ", 0, 10, 3);
	Dialog.addCheckbox("Median Filter", false);
	Dialog.addSlider("Radius: ", 0, 20, 2);
	Dialog.addCheckbox("Unsharp Mask", false);
	Dialog.addSlider("Radius: ", 0, 20, 2);
	Dialog.addSlider("Mask Strength: ", 0, 0.9, 0.6);
	Dialog.addCheckbox("Tophat (Iannetti et al., 2016)" , false);

    Dialog.show();

	CLAHE = Dialog.getCheckbox();
	CLAHE_BLOCKSIZE = Dialog.getNumber();
	CLAHE_HISTOGRAM = Dialog.getNumber();
	CLAHE_MAXSLOPE = Dialog.getNumber();

	MED = Dialog.getCheckbox();
	MED_RADIUS = Dialog.getNumber();

	UNSHARP = Dialog.getCheckbox();
	UNSHARP_RADIUS = Dialog.getNumber();
	UNSHARP_STRENGTH = Dialog.getNumber();
	TOPHAT = Dialog.getCheckbox();

}

//crude segmentation for cell area and nuclear area thanks to mountain_man on image.sc forum
function crudeSegment(){
		run("Variance...", "radius=0.5");
		run("8-bit");	//convert to 8-bit grayscale
		run("Gray Morphology", "radius=2 type=circle operator=erode");	//gray morphology to smooth dark edges
		run("Find Edges");		//essentially the derivative of pixel intensty vs position - highlights contrast
		run("Gray Morphology", "radius=4 type=circle operator=dilate");  // smooth edges to fill in minor gaps
		run("Auto Threshold", "method=Li white");
		run("Make Binary");		//convert to binary image
		run("Options...", "iterations=5 count=1 black pad do=Dilate");		//Dilate to fill gaps in edge
		run("Fill Holes");	//Fill center of defined perimeter
		run("Options...", "iterations=5 count=1 black pad do=Erode");			//Erode to undo previous dilation
		run("Erode");	//Erode to undo enlarging effect of find edges
		run("Open");	//Erode/Dilate to remove remaining small defects (dots) and smooth edges
		run("Convert to Mask");
}

//Preprocessing...
function preprocessing() {
	//CAUTION: ensure any preprocessing was also done on the trained images
	//Apply unsharp mask to image if selected
	selectWindow("Original");
	if (UNSHARP == true) {
		run("Unsharp Mask...", "radius="+toString(UNSHARP_RADIUS)+" mask="+toString(UNSHARP_STRENGTH));
	}

    //Apply contrast limited adaptive histogram equalization if selected
	if (CLAHE == true) {
		run("Enhance Local Contrast (CLAHE)",
                    "blocksize="+toString(CLAHE_BLOCKSIZE)+" histogram="+toString(CLAHE_HISTOGRAM)+" maximum="+toString(CLAHE_MAXSLOPE)+" mask=*None* fast_(less_accurate)");
	}

    //Apply median filtering if selected
	if (MED == true) {
		run("Median...", "radius="+toString(MED_RADIUS));
	}

	//Apply tophat filter if selected
	if (TOPHAT == true) {
		run("Convolve...", "text1=[0 0 -1 -1 -1 0 0 \n0 -1 -1 -1 -1 -1 0\n-1 -1 3 3 3 -1 -1\n-1 -1 3 4 3 -1 -1\n-1 -1 3 3 3 -1 -1\n0 -1 -1 -1 -1 -1 0\n0 0 -1 -1 -1 0 0 \n] normalize");
	}
	
	//call trainable weka 	
	weka=true;
	if (weka == true){
		
	//run("Trainable Weka Segmentation 3D");
	run("Trainable Weka Segmentation"); 
	wait(2000);
	call("trainableSegmentation.Weka_Segmentation.loadClassifier", path_to_model);
	wait(3000); 
	//USE PROBABILITY can set upper bound from 0.5 = 50% probability
	call("trainableSegmentation.Weka_Segmentation.getProbability");
	call("ij.plugin.frame.ThresholdAdjuster.setMode", "Red");
	setAutoThreshold("Default no-reset");
	//run("Threshold...");
	setThreshold(-1000000000000000000000000000000.0000, 0.5000);//confidence interval parameter
	//can comment out next line after first run to optimize the thresholding, use those values in previous line to apply to all cases
	waitForUser("type ctrl-shift-T to bring up threshold window. \n Ensure only 'dont reset range is checked' this choose apply \n Then convet to mask (default,default,black-backgrd) \n After acceptable confidence threshold is set, hard-code it line 491 comment out line 493");
	run("Convert to Mask", " ");
	run("Duplicate...", "title=segmented duplicate channels=2");//weirdly seems like ch2 gives better skeletons

	
	//USE AUTO-PROBABILITY RESULT
	// 50% probability seems default for get result command, alternatively try probability command set to 50-90% as above
	//call("trainableSegmentation.Weka_Segmentation.getResult");
	//setOption("ScaleConversions", false);
	//selectImage("Classified image");
	//setAutoThreshold("Default dark no-reset");
	//run("Threshold...");
	//setOption("BlackBackground", false);
	//run("Convert to Mask", "background=Dark calculate black");
	//run("Invert", "stack");
	//setOption("ScaleConversions", false);
	//run("8-bit");
	//rename("segmented");
	//selectWindow("segmented"); 
		
	}
	//set weka false and do basic thresholding for fast run as done in original macro, suggest to use median and tophat filters
	else{
	run("32-bit");
	run("Make Binary", "calculate black");//could modify method from defual to minError to get 'whole cell'
	run("8-bit");
	rename("segmented");
	}

}
