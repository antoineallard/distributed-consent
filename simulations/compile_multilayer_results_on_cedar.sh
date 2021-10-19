#!/bin/bash

# Antoine Allard
# antoineallard.info
# Octobre 2021


module load python
module load scipy-stack

python compile_multilayer_results.py

git add ../results/multilayer/*.pkl
git commit -m "updated multilayer pkl files"
git push origin main