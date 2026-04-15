#!/bin/bash
cd "C:/Users/mic23/prototype-embed/chameleon-research-loop"
C:/Users/mic23/miniconda3/envs/chameleon/python.exe experiments/iter_18_premmff_cv.py > experiments/iter_18_output.txt 2>&1
echo "Exit code: $?"
