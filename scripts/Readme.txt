############################################################

In this Readme, the procedure to run the codes is explained.

This folder contains matlab codes to extract detect laughter segments in the given audio file.

##########################################################################################

1) Call_Detect_Laughter_video.sh: Execute this shell script, if the input is in video format and we have to detect the laughter segments.
   To execute this code: sh Call_Detect_Laughter_video.sh /path/to/videofiles /path/where/audiofiles/to/be/saved
   Note: The output of the laughter detections will also be saved in the same folder, where audio files are saved.

##########################################################################################

2) Call_Detect_Laughter_audio.sh: Execute this shell script, if audio files are available and laughter detection need to be done.
   To execute this code: Call_Detect_Laughter_audio.sh /path/to/audiofiles/ /path/where/laughterDetectionOutput/to/be/saved

##########################################################################################

Following are matlab codes which are called from the script files

1) Audio_from_video.py: This is a python code to extract audio from the given video file.

2) Detect_laughter.m: This code calls all other codes and also writes the detected laughter regions to a text (.txt) file.

3) gci_location.m: This code is used to extract the Glottal closure instant locations

4) mean_smooth.m: This code is used to smooth the feature contours obtained

5) voice_decision.m: This code is used to obtain voiced/unvoiced segments in the given audio signal

6) zff_method.m and zffsig.m which are used for the feature extraction.
