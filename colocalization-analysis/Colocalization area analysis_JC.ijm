// reference: https://www.youtube.com/watch?v=4umlxVsjY04
// Channel 2 (green) and channel 3 (red) colocalization area measurement
// Generating colocalized area in yellow and save pictures in jpg

// Set the folder path
folder = getDirectory("Choose a Directory");

// Get the list of files in the folder
list = getFileList(folder);

// Define Bio-formats options to avoid the pop-up
options = "open=[open] autoscale color_mode=Default rois_import=[ROI manager]";

for (i = 0; i < list.length; i++) {
    // Check if the file has the ".nd2" extension
    if (endsWith(list[i], ".nd2")) {
        // Open each nd2 file with Bio-formats options
        run("Bio-Formats Importer", "open=[" + folder + list[i] + "]" + " " + options);
        
        // Split the channels
        run("Split Channels");

        // Duplicate channel 1 and rename it as "dup_ch1"
        selectWindow("C2-" + list[i]);
        run("Duplicate...", "title=Dup_ch2");
        
        // Duplicate channel 2 and rename it as "dup_ch2"
        selectWindow("C3-" + list[i]);
        run("Duplicate...", "title=Dup_ch3");
        
        // Use Median Filter with a radius of 1 pixel for better segmentation.
        selectImage("Dup_ch2");
        run("Median...", "radius=1");
        selectImage("Dup_ch3");
        run("Median...", "radius=1");

        // Gray LUT
        selectImage("Dup_ch2");
        run("Grays");
        selectImage("Dup_ch3");
        run("Grays");
        
        // Run threshold and wait for the user to adjust the threshold manually
        selectImage("Dup_ch2");
        setAutoThreshold("Otsu dark no-reset");
        call("ij.plugin.frame.ThresholdAdjuster.setMode", "Red");
        waitForUser("Adjust the threshold", "Please adjust the threshold, then click OK to continue.");
        run("Create Selection");
        roiManager("Add");
        
        selectImage("Dup_ch3");
        setAutoThreshold("Otsu dark no-reset");
        call("ij.plugin.frame.ThresholdAdjuster.setMode", "Red");
        waitForUser("Adjust the threshold", "Please adjust the threshold, then click OK to continue.");
        run("Create Selection");
        roiManager("Add");
        
        // Set ROI 
        roiManager("Select", 0);
        roiManager("Rename", "Ch2_Green");
        roiManager("Select", 1);
        roiManager("Rename", "Ch3_Red");
        roiManager("Select", newArray(0, 1));
        roiManager("AND");
        roiManager("Add");
        roiManager("Select", 2);
        roiManager("Set Fill Color", "yellow");
        roiManager("Rename", "Merged");
        
        // Save co-localization images
        selectWindow("C3-" + list[i]);
		roiManager("Select", 2); 
		run("Add Selection...");
		saveAs("Jpeg", folder + list[i] + "_colocal.jpg"); 
        
        // Measurement
        roiManager("Select", newArray(0,1,2));
        roiManager("Measure");
        
        // Save the results
        saveAs("Results", folder + list[i] + "_results.csv");
        
        
        // Close all images and clear the ROI Manager for the next file
        roiManager("Deselect");
        roiManager("Delete");
        run("Close All");
    }
}