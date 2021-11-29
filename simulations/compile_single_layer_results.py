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
    header.append("NaN?")
    df = pd.read_table(results_filename, names=header, comment="#", delimiter=r"\s+", on_bad_lines='warn')

    if len(df['NbVertices'].unique()) > 1:
        d = df['NbVertices'].value_counts().to_dict()
        del d[max(d, key=lambda k: d[k])]
        row_idx_to_remove = df[df['NbVertices'].isin(d.keys())].index
        df.drop(row_idx_to_remove, inplace=True)

    if len(df['Name'].unique()) > 1:
        d = df['Name'].value_counts().to_dict()
        del d[max(d, key=lambda k: d[k])]
        row_idx_to_remove = df[df['Name'].isin(d.keys())].index
        df.drop(row_idx_to_remove, inplace=True)

    if len(df['NaN?'].unique()) > 1:
        row_idx_to_remove = df[~df['NaN?'].isna()].index
        df.drop(row_idx_to_remove, inplace=True)
    df.drop(columns=['NaN?'], inplace=True)

    df["CounterCulture"] = (df["NObsGCNbType0"].astype('float64') + df["NObsGCNbType1"].astype('float64') + df["NObsGCNbType2"].astype('float64')) / df["NbVertices"].astype('float64')
    df["HerdImmunityType1"] = (df["ObsNbType1"].astype('float64') / (df["ObsNbType1"].astype('float64') + df["NObsNbType1"].astype('float64')))
    df["HerdImmunityType2"] = (df["ObsNbType2"].astype('float64') / (df["ObsNbType2"].astype('float64') + df["NObsNbType2"].astype('float64')))
    df["ObsCompRelSize"] = (df["ObsNbType0"].astype('float64') + df["ObsNbType1"].astype('float64') + df["ObsNbType2"].astype('float64')) / df["NbVertices"].astype('float64')

    pt = df.pivot_table(columns=["ObsDepth", "PrivProfFrac", "AppCoverage", "AdoptionRate"],
                        values = ["CounterCulture", "HerdImmunityType1", "HerdImmunityType2", "ObsCompRelSize"],
                        aggfunc = [np.mean])
    pt.to_pickle(results_filename.rsplit(".", 1)[-2] + '.pkl', compression="xz")
