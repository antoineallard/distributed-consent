#!/bin/bash

# Antoine Allard
# antoineallard.info
# Octobre 2021


module load python
module load scipy-stack

python compile_single_layer_results.py

git add ../results/single_layer/*.pkl
git commit -m "updated single_layer pkl files"
git push origin main