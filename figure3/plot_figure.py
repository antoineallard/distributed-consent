# -​*- coding: utf-8 -*​-
# @author: Antoine Allard <antoine.allard.1@gmail.com>
# @author: Laurent-Hébert Dufresne

# Packages
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


for names in glob.glob("../Facebook100/*.txt.tar.xz"):

    name = os.path.basename(names).split(".")[0]

    graph_properties_filename = "results/" + name + ".dat"

    if os.path.isfile(graph_properties_filename):

        header = open(graph_properties_filename, 'r').readline().replace('#', ' ').split()
        df = pd.read_table(graph_properties_filename, names=header, comment="#", delimiter=r"\s+")

        if df["NbVertices"][0].item() < 2000:
            continue

        df["counter_culture"] = (df["NObsGCNbType0"] + df["NObsGCNbType1"] + df["NObsGCNbType2"]) / df["NbVertices"]
        df["herd_immunity"] = (df["ObsNbType1"] / (df["ObsNbType1"] + df["NObsNbType1"]))
        df["obs_comp"] = (df["ObsNbType0"] + df["ObsNbType1"] + df["ObsNbType2"]) / df["NbVertices"]

        df0 = df.pivot_table(columns=["ObsDepth", "PrivProfFrac", "AppCoverage", "PassportAdop", "AdoptionRate"], values = ["counter_culture", "herd_immunity", "obs_comp"], aggfunc = [np.mean, frequency])

        print("{:15} {:7d} {:10.2f}".format(name, df.iloc[0]["NbVertices"], df0["frequency"].loc['obs_comp'].mean()))

        quantity_to_plot = "counter_culture"
        obs_depth = 2
        priv_prof_frac = 0.333333
        app_coverage = 0.002500

        passport_adopt = 0.500000
        ax2.plot(df0['mean'].loc[quantity_to_plot].loc[obs_depth, priv_prof_frac, app_coverage, passport_adopt].index.values * priv_prof_frac,
                 df0['mean'].loc[quantity_to_plot].loc[obs_depth, priv_prof_frac, app_coverage, passport_adopt].values,
                 linewidth=1.5, color='k', alpha=0.35)

        passport_adopt = 0.900000
        ax0.plot(df0['mean'].loc[quantity_to_plot].loc[obs_depth, priv_prof_frac, app_coverage, passport_adopt].index.values * priv_prof_frac,
                 df0['mean'].loc[quantity_to_plot].loc[obs_depth, priv_prof_frac, app_coverage, passport_adopt].values,
                 linewidth=1.5, color='k', alpha=0.35)

        passport_adopt = 1.000000
        ax1.plot(df0['mean'].loc[quantity_to_plot].loc[obs_depth, priv_prof_frac, app_coverage, passport_adopt].index.values * priv_prof_frac,
                 df0['mean'].loc[quantity_to_plot].loc[obs_depth, priv_prof_frac, app_coverage, passport_adopt].values,
                 linewidth=1.5, color='k', alpha=0.35)


ax0.plot(np.NaN, np.NaN, color=None, linewidth = 0, label = r'Consent passport adoption = 90\%')
ax0.plot(np.arange(0,0.34,0.01), (2/3)*np.arange(0,0.34,0.01), linestyle = '--', color='black', linewidth = 2, label = 'Fraction of prevented data flow')
ax0.set_xlabel("Adoption of distributed consent")
ax0.legend(loc="upper left", ncol=1, prop={"size": "x-small"})
ax0.set_ybound(0,0.4)

ax1.plot(np.NaN, np.NaN, color=None, linewidth = 0, label = r'Consent passport adoption = 100\%')
ax1.plot(np.arange(0,0.34,0.01), (2/3)*np.arange(0,0.34,0.01), linestyle = '--', color='black', linewidth = 2, label = 'Fraction of prevented data flow')
ax1.legend(loc="upper left", ncol=1, prop={"size": "x-small"})
ax1.set_xlabel("Adoption of distributed consent")
ax1.set_ybound(0,0.4)

ax2.plot(np.NaN, np.NaN, color=None, linewidth = 0, label = r'Consent passport adoption = 50\%')
ax2.plot(np.arange(0,0.34,0.01), (2/3)*np.arange(0,0.34,0.01), linestyle = '--', color='black', linewidth = 2, label = 'Fraction of prevented data flow')
ax2.legend(loc="upper left", ncol=1, prop={"size": "x-small"})
ax2.set_xlabel("Adoption of distributed consent")
ax2.set_ybound(0,0.4)
ax2.set_ylabel("Largest unobserved component")

fig.savefig("figure3.pdf", bbox_inches='tight')
