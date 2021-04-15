# -​*- coding: utf-8 -*​-
# @author: Antoine Allard <antoineallard.info>

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

        df["counter_culture"] = (df["NObsGCNbType0"] + df["NObsGCNbType1"] + df["NObsGCNbType2"]) / df["NbVertices"]
        df["herd_immunity"] = (df["ObsNbType1"] / (df["ObsNbType1"] + df["NObsNbType1"]))
        df["obs_comp"] = (df["ObsNbType0"] + df["ObsNbType1"] + df["ObsNbType2"]) / df["NbVertices"]

        df0 = df.pivot_table(columns=["ObsDepth", "PrivProfFrac", "AppCoverage", "AdoptionRate"], values = ["counter_culture", "herd_immunity", "obs_comp"], aggfunc = [np.mean, frequency])

        print("{:15} {:7d} {:10.2f}".format(name, df.iloc[0]["NbVertices"], df0["frequency"].loc['obs_comp'].mean()))

        obs_depth = 2
        priv_prof_frac = 0.333333
        app_coverage = 0.01

        quantity_to_plot = "counter_culture"
        ax0.plot(df0['mean'].loc[quantity_to_plot].loc[obs_depth, priv_prof_frac, app_coverage].index.values * priv_prof_frac,
                 df0['mean'].loc[quantity_to_plot].loc[obs_depth, priv_prof_frac, app_coverage].values,
                 linewidth=1.5, color="#f29e78", alpha=0.35)

        quantity_to_plot = "herd_immunity"
        ax1.plot(df0['mean'].loc[quantity_to_plot].loc[obs_depth, priv_prof_frac, app_coverage].index.values * priv_prof_frac,
                 df0['mean'].loc[quantity_to_plot].loc[obs_depth, priv_prof_frac, app_coverage].values / df0['mean'].loc[quantity_to_plot].loc[obs_depth, priv_prof_frac, app_coverage].values[0],
                 linewidth=1.5, color="#f278cc", alpha=0.35)

        quantity_to_plot = "obs_comp"
        ax2.plot(df0['mean'].loc[quantity_to_plot].loc[obs_depth, priv_prof_frac, app_coverage].index.values * priv_prof_frac,
                 df0['mean'].loc[quantity_to_plot].loc[obs_depth, priv_prof_frac, app_coverage].values,
                 linewidth=1.5, color="#78ccf2", alpha=0.35)


ax0.set_xlabel("Adoption of distributed consent")
ax0.plot(np.arange(0,0.34,0.01), (2/3)*np.arange(0,0.34,0.01), linestyle = '--', color='black', linewidth = 2, label = 'Fraction of prevented data flow')
ax0.set_ylabel("Largest unobserved component")
ax0.legend(loc="upper left", ncol=1, prop={"size": "x-small"})

ax1.set_xlabel("Adoption of distributed consent")
ax1.set_ylabel("Fraction of observed type 1")

ax2.set_xlabel("Adoption of distributed consent")
ax2.set_ylabel("Fraction of observed individuals")

fig.savefig("figure2.pdf")
