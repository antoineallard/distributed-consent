#!/bin/bash

# Antoine Allard
# antoineallard.info
# Octobre 2021


module load python
module load scipy-stack

python compile_multilayer_results.py

# git add ../results/multilayer/*.json
# git commit -m "updated multilayer json files"
# git push origin main