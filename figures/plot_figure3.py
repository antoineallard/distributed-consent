# -*- coding: utf-8 -*-
# @author: Antoine Allard <antoineallard.info>
# @author: Laurent-Hébert Dufresne

# Packages
import glob
import numpy as np
import matplotlib
import matplotlib.gridspec
import matplotlib.pyplot as plt
import os
import pandas as pd


# Global parameters for the figures.
matplotlib.use('agg')
plt.style.use('seaborn-deep')
plt.rcParams["text.usetex"] = True
plt.rcParams["font.size"] = 24
plt.rcParams["legend.frameon"] = False
plt.rcParams["legend.fancybox"] = False


fig = plt.figure(constrained_layout=True, figsize=(18, 11))

gs = matplotlib.gridspec.GridSpec(nrows=2,
                                  ncols=3,
                                  figure=fig,
                                  width_ratios=[1, 1, 1],
                                  height_ratios=[1, 1],
                                  wspace=0.2, hspace=0.05
                                  )

ax0 = fig.add_subplot(gs[0, 0])
ax1 = fig.add_subplot(gs[0, 1])
ax2 = fig.add_subplot(gs[0, 2])

ax3 = fig.add_subplot(gs[1, 0])
ax4 = fig.add_subplot(gs[1, 1])
ax5 = fig.add_subplot(gs[1, 2])


for filename in glob.glob("../results/multilayer/*.pkl"):

    name = os.path.basename(filename).rsplit(".", 1)[-2]

    # Ignores networks with fewer than 2000 nodes.
    if name not in ['Swarthmore42', 'Haverford76', 'Simmons81', 'Reed98', 'Caltech36']:

        pt = pd.read_pickle(filename, compression="xz")

        # =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
        # First row.

        # valeurs utilisées dans la version actuelle du papier
        quantity_to_plot = "CounterCulture"
        obs_depth = 2
        priv_prof_frac = 0.333333
        app_coverage = 0.0025

        passport_adopt = 0.500000
        ax0.plot(pt['mean', obs_depth, priv_prof_frac, app_coverage, passport_adopt].loc[quantity_to_plot].index.values * priv_prof_frac,
                 pt['mean', obs_depth, priv_prof_frac, app_coverage, passport_adopt].loc[quantity_to_plot].values,
                 linewidth=1.5, color='k', alpha=0.35)

        passport_adopt = 0.900000
        ax1.plot(pt['mean', obs_depth, priv_prof_frac, app_coverage, passport_adopt].loc[quantity_to_plot].index.values * priv_prof_frac,
                 pt['mean', obs_depth, priv_prof_frac, app_coverage, passport_adopt].loc[quantity_to_plot].values,
                 linewidth=1.5, color='k', alpha=0.35)

        passport_adopt = 1.000000
        ax2.plot(pt['mean', obs_depth, priv_prof_frac, app_coverage, passport_adopt].loc[quantity_to_plot].index.values * priv_prof_frac,
                 pt['mean', obs_depth, priv_prof_frac, app_coverage, passport_adopt].loc[quantity_to_plot].values,
                 linewidth=1.5, color='k', alpha=0.35)


        # # =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
        # # Second row.

        # quantity_to_plot = "CounterCulture"
        # obs_depth = 2
        # priv_prof_frac = 0.333333  # options: 0.333333 0.500000 0.66666
        # app_coverage = 0.0005      # options: 0.01 0.025 0.005 0.0025 0.001 0.0005 0.0001

        # passport_adopt = 0.500000
        # ax0.plot(pt['mean', obs_depth, priv_prof_frac, app_coverage, passport_adopt].loc[quantity_to_plot].index.values * priv_prof_frac,
        #          pt['mean', obs_depth, priv_prof_frac, app_coverage, passport_adopt].loc[quantity_to_plot].values,
        #          linewidth=1.5, color='k', alpha=0.35)

        # passport_adopt = 0.900000
        # ax1.plot(pt['mean', obs_depth, priv_prof_frac, app_coverage, passport_adopt].loc[quantity_to_plot].index.values * priv_prof_frac,
        #          pt['mean', obs_depth, priv_prof_frac, app_coverage, passport_adopt].loc[quantity_to_plot].values,
        #          linewidth=1.5, color='k', alpha=0.35)

        # passport_adopt = 1.000000
        # ax2.plot(pt['mean', obs_depth, priv_prof_frac, app_coverage, passport_adopt].loc[quantity_to_plot].index.values * priv_prof_frac,
        #          pt['mean', obs_depth, priv_prof_frac, app_coverage, passport_adopt].loc[quantity_to_plot].values,
        #          linewidth=1.5, color='k', alpha=0.35)


ax0.plot(np.NaN, np.NaN, color=None, linewidth = 0, label = r'Consent passport adoption = 50\%')
ax0.plot(np.arange(0,0.34,0.01), (2/3)*np.arange(0,0.34,0.01), linestyle = '--', color='black', linewidth = 2, label = 'Fraction of prevented data flow')
ax0.legend(loc="upper left", ncol=1, prop={"size": "x-small"})

ax1.plot(np.NaN, np.NaN, color=None, linewidth = 0, label = r'Consent passport adoption = 90\%')
ax1.plot(np.arange(0,0.34,0.01), (2/3)*np.arange(0,0.34,0.01), linestyle = '--', color='black', linewidth = 2, label = 'Fraction of prevented data flow')
ax1.legend(loc="upper left", ncol=1, prop={"size": "x-small"})

ax2.plot(np.NaN, np.NaN, color=None, linewidth = 0, label = r'Consent passport adoption = 100\%')
ax2.plot(np.arange(0,0.34,0.01), (2/3)*np.arange(0,0.34,0.01), linestyle = '--', color='black', linewidth = 2, label = 'Fraction of prevented data flow')
ax2.legend(loc="upper left", ncol=1, prop={"size": "x-small"})

ax0.set_xlabel("Adoption of distributed consent")
ax1.set_xlabel("Adoption of distributed consent")
ax2.set_xlabel("Adoption of distributed consent")

ax0.set_ylabel("Largest unobserved component")

ax0.set_ybound(0,0.4)
ax1.set_ybound(0,0.4)
ax2.set_ybound(0,0.4)

fig.savefig("figure3.pdf", bbox_inches='tight')
fig.savefig("figure3.svg", bbox_inches='tight')