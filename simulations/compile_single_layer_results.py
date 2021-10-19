# -*- coding: utf-8 -*-
# @author: Antoine Allard <antoineallard.info>


# Packages
import glob
import os
import pandas as pd
import numpy as np


for results_filename in glob.glob("../results/single_layer/*.dat"):

    print("Compiling " + results_filename.rsplit(".", 1)[-2] + '.json')

    header = open(results_filename, 'r').readline().replace('#', ' ').split()
    df = pd.read_table(results_filename, names=header, comment="#", delimiter=r"\s+")

    df["counter_culture"] = (df["NObsGCNbType0"] + df["NObsGCNbType1"] + df["NObsGCNbType2"]) / df["NbVertices"]
    df["herd_immunity"] = (df["ObsNbType1"] / (df["ObsNbType1"] + df["NObsNbType1"]))
    df["obs_comp"] = (df["ObsNbType0"] + df["ObsNbType1"] + df["ObsNbType2"]) / df["NbVertices"]

    pt = df.pivot_table(columns=["ObsDepth", "PrivProfFrac", "AppCoverage", "AdoptionRate"], values = ["counter_culture", "herd_immunity", "obs_comp"], aggfunc = [np.mean])
    # pt.to_json(results_filename.rsplit(".", 1)[-2] + '.json')
    pt.to_pickle(results_filename.rsplit(".", 1)[-2] + '.pkl', compression="xz")
