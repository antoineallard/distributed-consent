# -​*- coding: utf-8 -*​-
# @author: Antoine Allard <antoine.allard.1@gmail.com>

# Packages
from pylab import *
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec
import cycler
import glob
import os
import seaborn as sns
import pandas as pd
import numpy as np
matplotlib.use('agg')


def lowerbound(x):
    return np.percentile(x, 5)


def upperbound(x):
    return np.percentile(x, 95)

def frequency(x):
    return len(x)



# Global parameters for the figures.
plt.style.use('seaborn-deep')

# new_axes_prop_cycle = cycler.cycler('linestyle', ['-', '--', '-.', ':']) * \
#                       plt.rcParams["axes.prop_cycle"]
# plt.rcParams["axes.prop_cycle"] = new_axes_prop_cycle

plt.rcParams["text.usetex"] = True
plt.rcParams["font.size"] = 24
plt.rcParams["legend.frameon"] = False
plt.rcParams["legend.fancybox"] = False

fig = plt.figure(constrained_layout=True, figsize=(18, 5.5))
gs = matplotlib.gridspec.GridSpec(nrows=1,
                                  ncols=3,
                                  figure=fig,
                                  width_ratios=[1, 1, 1],
                                  height_ratios=[1],
                                  wspace=0.2, hspace=0.05
                                  )
ax0 = fig.add_subplot(gs[0, 1])
ax1 = fig.add_subplot(gs[0, 2])
ax2 = fig.add_subplot(gs[0, 0])

linestyles = {0.333333: "-",
              0.50000: "--",
              0.800000: "-."}
colors = {0.333333: plt.rcParams["axes.prop_cycle"].by_key()['color'][0],
          0.500000: plt.rcParams["axes.prop_cycle"].by_key()['color'][1],
          0.800000: plt.rcParams["axes.prop_cycle"].by_key()['color'][2]}


# for name, sub_df in df0.loc[quantity_to_plot].groupby(level=0):


for names in glob.glob("../../../../data_ethics_data/Facebook100/*.txt.tar.xz"):

    name = os.path.basename(names).split(".")[0]

    graph_properties_filename = "../../../../data_ethics_data/results/"\
                                "importance_of_passports/" + name + ".dat"

    if os.path.isfile(graph_properties_filename): # and name not in ["Caltech36"]:

        header = open(graph_properties_filename, 'r')\
            .readline().replace('#', ' ').split()
        df = pd.read_table(graph_properties_filename, names=header,
                           comment="#", delimiter=r"\s+")

        if df["NbVertices"][0].item() < 2000:
            continue


        df["counter_culture"] = (df["NObsGCNbType0"] + df["NObsGCNbType1"]
                                 + df["NObsGCNbType2"]) / df["NbVertices"]
        df["herd_immunity"] = (df["ObsNbType1"]
                               / (df["ObsNbType1"] + df["NObsNbType1"]))
        df["obs_comp"] = (df["ObsNbType0"] + df["ObsNbType1"]
                          + df["ObsNbType2"]) / df["NbVertices"]

        df0 = df.pivot_table(columns=["ObsDepth", "PrivProfFrac", "AppCoverage", "PassportAdop", "AdoptionRate"], values = ["counter_culture", "herd_immunity", "obs_comp"], aggfunc = [np.mean, frequency])

        # nb_vertices = df.iloc[0]["NbVertices"]
        # if nb_vertices < 1e4:
        #     print(name, end=" ")

        print("{:15} {:7d} {:10.2f}".format(name, df.iloc[0]["NbVertices"], df0["frequency"].mean()))

        quantity_to_plot = "counter_culture"
        for obs_depth, sub_df in df0.loc[quantity_to_plot].groupby(level=0):
            if obs_depth in [2]:
                
                for priv_prof_frac, sub_df in df0.loc[quantity_to_plot, obs_depth].groupby(level=0):
                    if priv_prof_frac in [0.333333]:

                        for app_coverage, sub_df in df0.loc[quantity_to_plot, obs_depth, priv_prof_frac].groupby(level=0):
                            if app_coverage in [0.002500]:
                
                                for passport_adopt, sub_df in df0.loc[quantity_to_plot, obs_depth, priv_prof_frac, app_coverage].groupby(level=0):
                                    if passport_adopt in [0.500000]:

                                        ax2.plot(df0.loc[quantity_to_plot, obs_depth, priv_prof_frac, app_coverage, passport_adopt]['mean'].index.values * priv_prof_frac,
                                            df0.loc[quantity_to_plot, obs_depth, priv_prof_frac, app_coverage, passport_adopt]['mean'].values,
                                            linewidth=1.5, color='k', alpha=0.35)
                
                                    if passport_adopt in [0.900000]:

                                        ax0.plot(df0.loc[quantity_to_plot, obs_depth, priv_prof_frac, app_coverage, passport_adopt]['mean'].index.values * priv_prof_frac,
                                            df0.loc[quantity_to_plot, obs_depth, priv_prof_frac, app_coverage, passport_adopt]['mean'].values,
                                            linewidth=1.5, color='k', alpha=0.35)
                
                                    if passport_adopt in [1.000000]:

                                        ax1.plot(df0.loc[quantity_to_plot, obs_depth, priv_prof_frac, app_coverage, passport_adopt]['mean'].index.values * priv_prof_frac,
                                            df0.loc[quantity_to_plot, obs_depth, priv_prof_frac, app_coverage, passport_adopt]['mean'].values,
                                            linewidth=1.5, color='k', alpha=0.35)



ax0.set_xlabel("Adoption of distributed consent")
ax1.plot(np.NaN, np.NaN, color=None, linewidth = 0, label = r'Consent passport adoption = 100\%')
ax2.plot(np.NaN, np.NaN, color=None, linewidth = 0, label = r'Consent passport adoption = 50\%')
ax0.plot(np.NaN, np.NaN, color=None, linewidth = 0, label = r'Consent passport adoption = 90\%')
ax0.plot(np.arange(0,0.34,0.01), (2/3)*np.arange(0,0.34,0.01), linestyle = '--', color='black', linewidth = 2, label = 'Fraction of prevented data flow')
ax1.plot(np.arange(0,0.34,0.01), (2/3)*np.arange(0,0.34,0.01), linestyle = '--', color='black', linewidth = 2, label = 'Fraction of prevented data flow')
ax2.plot(np.arange(0,0.34,0.01), (2/3)*np.arange(0,0.34,0.01), linestyle = '--', color='black', linewidth = 2, label = 'Fraction of prevented data flow')

ax2.set_ylabel("Largest unobserved component")
# ax0.set_aspect('equal', adjustable='box')
# ax1.set_aspect('equal', adjustable='box')
# ax2.set_aspect('equal', adjustable='box')
ax2.legend(loc="upper left", ncol=1, prop={"size": "x-small"})
ax1.legend(loc="upper left", ncol=1, prop={"size": "x-small"})
ax0.legend(loc="upper left", ncol=1, prop={"size": "x-small"})

ax0.set_ybound(0,0.4)
ax1.set_ybound(0,0.4)
ax2.set_ybound(0,0.4)

ax1.set_xlabel("Adoption of distributed consent")
# # ax1.legend(loc="upper left", ncol=1, prop={"size": "x-small"})

ax2.set_xlabel("Adoption of distributed consent")
# # ax2.legend(loc="lower left", ncol=1, prop={"size": "x-small"})

fig.savefig("../pdf/importance_of_passports.pdf")
# fig.savefig("../png/emergence_of_components.png")
# fig.savefig("../svg/emergence_of_components.svg")
