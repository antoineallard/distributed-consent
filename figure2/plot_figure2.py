# -​*- coding: utf-8 -*​-
# @author: Antoine Allard <antoine.allard.1@gmail.com>

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


for names in glob.glob("../Facebook100/*.txt.tar.xz"):

    name = os.path.basename(names).split(".")[0]

    graph_properties_filename = "results/emergence_of_components/" + name + ".dat"

    if os.path.isfile(graph_properties_filename): # and name not in ["Caltech36"]:

        header = open(graph_properties_filename, 'r')\
            .readline().replace('#', ' ').split()
        df = pd.read_table(graph_properties_filename, names=header,
                           comment="#", delimiter=r"\s+")

        df["counter_culture"] = (df["NObsGCNbType0"] + df["NObsGCNbType1"]
                                 + df["NObsGCNbType2"]) / df["NbVertices"]
        df["herd_immunity"] = (df["ObsNbType1"]
                               / (df["ObsNbType1"] + df["NObsNbType1"]))
        df["obs_comp"] = (df["ObsNbType0"] + df["ObsNbType1"]
                          + df["ObsNbType2"]) / df["NbVertices"]

        df0 = df.pivot_table(columns=["ObsDepth", "PrivProfFrac", "AppCoverage", "AdoptionRate"], values = ["counter_culture", "herd_immunity", "obs_comp"], aggfunc = [np.mean, frequency])

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
                            if app_coverage in [0.01]:

                                ax0.plot(df0.loc[quantity_to_plot, obs_depth, priv_prof_frac, app_coverage]['mean'].index.values * priv_prof_frac,
                                         df0.loc[quantity_to_plot, obs_depth, priv_prof_frac, app_coverage]['mean'].values,
                                         # color=sns.color_palette("Blues")[4], alpha=0.3)
                                         linewidth=1.5, color="#f29e78", alpha=0.35)

        quantity_to_plot = "herd_immunity"
        for obs_depth, sub_df in df0.loc[quantity_to_plot].groupby(level=0):
            if obs_depth in [2]:

                for priv_prof_frac, sub_df in df0.loc[quantity_to_plot, obs_depth].groupby(level=0):
                    if priv_prof_frac in [0.333333]:

                        for app_coverage, sub_df in df0.loc[quantity_to_plot, obs_depth, priv_prof_frac].groupby(level=0):
                            if app_coverage in [0.01]:

                                ax1.plot(df0.loc[quantity_to_plot, obs_depth, priv_prof_frac, app_coverage]['mean'].index.values * priv_prof_frac,
                                         df0.loc[quantity_to_plot, obs_depth, priv_prof_frac, app_coverage]['mean'].values / df0.loc[quantity_to_plot, obs_depth, priv_prof_frac, app_coverage]['mean'].values[0],
                                         # color=sns.color_palette("Greens")[4], alpha=0.3)
                                         linewidth=1.5, color="#f278cc", alpha=0.35)


        quantity_to_plot = "obs_comp"
        for obs_depth, sub_df in df0.loc[quantity_to_plot].groupby(level=0):
            if obs_depth in [2]:

                for priv_prof_frac, sub_df in df0.loc[quantity_to_plot, obs_depth].groupby(level=0):
                    if priv_prof_frac in [0.333333]:

                        for app_coverage, sub_df in df0.loc[quantity_to_plot, obs_depth, priv_prof_frac].groupby(level=0):
                            if app_coverage in [0.01]:

                                ax2.plot(df0.loc[quantity_to_plot, obs_depth, priv_prof_frac, app_coverage]['mean'].index.values * priv_prof_frac,
                                         df0.loc[quantity_to_plot, obs_depth, priv_prof_frac, app_coverage]['mean'].values,
                                         # color=sns.color_palette("Reds")[4], alpha=0.3)
                                         linewidth=1.5, color="#78ccf2", alpha=0.35)


ax0.set_xlabel("Adoption of distributed consent")
ax0.plot(np.arange(0,0.34,0.01), (2/3)*np.arange(0,0.34,0.01), linestyle = '--', color='black', linewidth = 2, label = 'Fraction of prevented data flow')
# ax0.set_ylabel("Fraction of unobserved individuals")
ax0.set_ylabel("Largest unobserved component")
# ax0.set_aspect('equal', adjustable='box')
# ax1.set_aspect('equal', adjustable='box')
# ax2.set_aspect('equal', adjustable='box')
ax0.legend(loc="upper left", ncol=1, prop={"size": "x-small"})

ax1.set_xlabel("Adoption of distributed consent")
ax1.set_ylabel("Fraction of observed type 1")
# # ax1.legend(loc="upper left", ncol=1, prop={"size": "x-small"})

ax2.set_xlabel("Adoption of distributed consent")
ax2.set_ylabel("Fraction of observed individuals")
# # ax2.legend(loc="lower left", ncol=1, prop={"size": "x-small"})


fig.savefig("figure2.pdf")
# fig.savefig("../png/emergence_of_components.png")
# fig.savefig("../svg/emergence_of_components.svg")
