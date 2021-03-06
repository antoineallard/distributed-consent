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

    if len(df['Name'].unique()) > 1:
        print(df['Name'].unique())

    df["CounterCulture"] = (df["NObsGCNbType0"] + df["NObsGCNbType1"] + df["NObsGCNbType2"]) / df["NbVertices"]
    df["HerdImmunityType1"] = (df["ObsNbType1"] / (df["ObsNbType1"] + df["NObsNbType1"]))
    df["HerdImmunityType2"] = (df["ObsNbType2"] / (df["ObsNbType2"] + df["NObsNbType2"]))
    df["ObsCompRelSize"] = (df["ObsNbType0"] + df["ObsNbType1"] + df["ObsNbType2"]) / df["NbVertices"]

    pt = df.pivot_table(columns=["ObsDepth", "PrivProfFrac", "AppCoverage", "PassportAdop", "AdoptionRate"],
                        values = ["CounterCulture", "HerdImmunityType1", "HerdImmunityType2", "ObsCompRelSize"],
                        aggfunc = [np.mean])
    pt.to_pickle(results_filename.rsplit(".", 1)[-2] + '.pkl', compression="xz")
