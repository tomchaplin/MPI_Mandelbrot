#!/bin/bash
cd $*
ffmpeg -r 10 -i file%04d.bmp -c:v libx264 -preset slow -crf 21 output.mp4

