#!/bin/bash

data_dir='/home/sriharsha/Desktop/scripts/audio_files/'; # enter full path of the folder containing audio files
out_dir='/home/sriharsha/Desktop/scripts/audio_files/';  # path where the output has to be saved

matlab -nodesktop -nosplash -nodisplay -nojvm -r "Detect_laughter $data_dir $out_dir"
