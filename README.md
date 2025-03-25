To use my pipeline to the Ole Rømer Observatory, the folder structure is important. You should have a main folder, which you can name what you want e.g. Observations or M101, in that folder you should have one folder, which is named after the date of the observations, where all the raw observations are placed a folder called "Raw", this is all the observations made so both light images and calibration images. The pipeline will automatically make the folders "Master" and "Calibrated", which is where the master frames of the calibration images and the calibrated light images will be placed respectively. You can see the folder structure below for a clearer overview ;)

If no calibration images were taken, or if you only got bias and darks, do not worry, I have made sure you can still get calibrated light images. You just need to place either your calibration images from a different night or mine calibration images into folders called "Bias", "Flats", and "Darks" where bias, flat, and dark calibration images go respectively. This should be placed in the main folder, so NOT the folder for the observation night, and it should work just fine. The code will automatically find the folders and use those images ;)

The log folder is a folder where two txt files for a given observation night will be placed. The folders will be named "log_cal_"date".txt" and "log_raw_"date".txt", which will give some basic information for all the images, both before and after calibration. Some of the information is the weighted mean and the standard deviation, but also the filter, image type, and object can be viewed. Hopeful you shouldn't need to look at the logs, but they are there just in case. If you want to find the bias frame or flat frame that was used in the calibration of an image, you need to look at the header :)

-> Observations 
      -> Observation night (e.g. 2025-02-16) 
              -> Raw
              -> Master 
              -> Calibrated
       -> Bias
       -> Flats
       -> Darks
       -> Log

Another thing you need to remember is that the code uses the image type, which is written in the header of the fits file directly at the Ole Rømer Observatory, so when making observations you need to change the image type before taking images, you do that in the SIPS program just above where you change filters. The code will NOT work if you have forgotten to change from flats to lights when taking observations of the object you are interested in.
