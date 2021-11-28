# -*- coding: utf-8 -*-
# @author: Antoine Allard <antoineallard.info>


# Packages
import glob
import os
import pandas as pd
import numpy as np


for results_filename in glob.glob("../results/single_layer/*.dat"):

    print("Updating " + results_filename.rsplit(".", 1)[-2] + '.pkl')

    header = open(results_filename, 'r').readline().replace('#', ' ').split()
    df = pd.read_table(results_filename, names=header, comment="#", delimiter=r"\s+", on_bad_lines='warn')

    if len(df['NbVertices'].unique()) > 1:
        d = df['NbVertices'].value_counts().to_dict()
        del d[max(d, key=lambda k: d[k])]
    df.drop(df[df['NbVertices'].isin(d.keys())].index, inplace=True)

    if len(df['Name'].unique()) > 1:
        d = df['Name'].value_counts().to_dict()
        del d[max(d, key=lambda k: d[k])]
    df.drop(df[df['Name'].isin(d.keys())].index, inplace=True)

    df["CounterCulture"] = (df["NObsGCNbType0"] + df["NObsGCNbType1"] + df["NObsGCNbType2"]) / df["NbVertices"]
    df["HerdImmunityType1"] = (df["ObsNbType1"] / (df["ObsNbType1"] + df["NObsNbType1"]))
    df["HerdImmunityType2"] = (df["ObsNbType2"] / (df["ObsNbType2"] + df["NObsNbType2"]))
    df["ObsCompRelSize"] = (df["ObsNbType0"] + df["ObsNbType1"] + df["ObsNbType2"]) / df["NbVertices"]

    pt = df.pivot_table(columns=["ObsDepth", "PrivProfFrac", "AppCoverage", "AdoptionRate"],
                        values = ["CounterCulture", "HerdImmunityType1", "HerdImmunityType2", "ObsCompRelSize"],
                        aggfunc = [np.mean])
    pt.to_pickle(results_filename.rsplit(".", 1)[-2] + '.pkl', compression="xz")
