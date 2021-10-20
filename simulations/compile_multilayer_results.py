# -*- coding: utf-8 -*-
# @author: Antoine Allard <antoineallard.info>


# Packages
import glob
import os
import pandas as pd
import numpy as np


for results_filename in glob.glob("../results/multilayer/*.dat"):

    print("Updating " + results_filename.rsplit(".", 1)[-2] + '.pkl')

    header = open(results_filename, 'r').readline().replace('#', ' ').split()
    df = pd.read_table(results_filename, names=header, comment="#", delimiter=r"\s+")

    df["counter_culture"] = (df["NObsGCNbType0"] + df["NObsGCNbType1"] + df["NObsGCNbType2"]) / df["NbVertices"]
    df["herd_immunity"] = (df["ObsNbType1"] / (df["ObsNbType1"] + df["NObsNbType1"]))
    df["herd_immunity2"] = (df["ObsNbType2"] / (df["ObsNbType2"] + df["NObsNbType2"]))
    df["obs_comp"] = (df["ObsNbType0"] + df["ObsNbType1"] + df["ObsNbType2"]) / df["NbVertices"]

    pt = df.pivot_table(columns=["ObsDepth", "PrivProfFrac", "AppCoverage", "PassportAdop", "AdoptionRate"], values = ["counter_culture", "herd_immunity", "herd_immunity2", "obs_comp"], aggfunc = [np.mean])
    # pt.to_json(results_filename.rsplit(".", 1)[-2] + '.json')
    pt.to_pickle(results_filename.rsplit(".", 1)[-2] + '.pkl', compression="xz")
