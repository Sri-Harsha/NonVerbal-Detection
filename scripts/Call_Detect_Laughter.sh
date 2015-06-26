#!/bin/bash

# This code is used to call the python code to extract audio from video files and then detect laughtert regions in the given audio file

data_dir="$1"; # enter full path of the folder containing audio files
out_dir="$2";  # path where the output has to be saved

# call the python code to extract audio from the given video file

python Audio_from_video.py $data_dir $out_dir

# call matlab code to detect laughter segments in the given audio file

matlab -nodesktop -nosplash -nodisplay -nojvm -r "Detect_laughter $out_dir $out_dir"
