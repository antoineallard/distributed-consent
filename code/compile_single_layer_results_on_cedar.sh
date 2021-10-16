#!/bin/bash

# Antoine Allard
# antoineallard.info
# Octobre 2021


module load python
module load scipy-stack

python compile_single_layer_results.py

git add ../results/single_layer/*.json
git commit -m "updated single_layer json files"
git push origin main