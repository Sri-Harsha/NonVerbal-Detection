#!/bin/bash

data_dir="$1"; # enter full path of the folder containing audio files
out_dir="$2";  # path where the output has to be saved

python Audio_from_video.py $data_dir $out_dir

matlab -nodesktop -nosplash -nodisplay -nojvm -r "Detect_laughter $out_dir $out_dir"
