############################################################

In this Readme, the procedure to run the codes is explained.

This folder contains matlab codes to extract detect laughter segments in the given audio file.

1) call_Detect_laugh.sh: This is a shell script which is used to call the main matlab code i.e., Detect_laughter.m.

2)Audio_from_video.py: This is a python code to extract audio from the given video file.

3) Detect_laughter.m: This code calls all other codes and also writes the detected laughter regions to a text (.txt) file.

4) gci_location.m: This code is used to extract the Glottal closure instant locations

5) mean_smooth.m: This code is used to smooth the feature contours obtained

6) voice_decision.m: This code is used to obtain voiced/unvoiced segments in the given audio signal

7) zff_method.m and zffsig.m which are used for the feature extraction.
