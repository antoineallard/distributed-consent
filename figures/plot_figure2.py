# -*- coding: utf-8 -*-
# @author: Antoine Allard <antoineallard.info>

# Packages
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec
import glob
import os
import numpy as np
matplotlib.use('agg')


# Global parameters for the figures.
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


for filename in glob.glob("../results/single_layer/*.pkl"):

    name = os.path.basename(filename).rsplit(".", 1)[-2]

    # Ignores networks with fewer than 2000 nodes.
    if name not in ['Swarthmore42', 'Haverford76', 'Simmons81', 'Reed98', 'Caltech36']:

        pt = pd.read_pickle(filename, compression="xz")

        # valeurs utilisées dans la version actuelle du papier
        obs_depth = 2
        priv_prof_frac = 0.333333
        app_coverage = 0.01

        quantity_to_plot = "ObsCompRelSize"
        ax0.plot(pt['mean', obs_depth, priv_prof_frac, app_coverage].loc[quantity_to_plot].index.values * priv_prof_frac,
                 pt['mean', obs_depth, priv_prof_frac, app_coverage].loc[quantity_to_plot].values,
                 linewidth=1.5, color="#78ccf2", alpha=0.35)

        quantity_to_plot = "CounterCulture"
        ax1.plot(pt['mean', obs_depth, priv_prof_frac, app_coverage].loc[quantity_to_plot].index.values * priv_prof_frac,
                 pt['mean', obs_depth, priv_prof_frac, app_coverage].loc[quantity_to_plot].values,
                 linewidth=1.5, color="#f29e78", alpha=0.35)

        # quantity_to_plot = "HerdImmunityType1"
        # ax2.plot(pt['mean', obs_depth, priv_prof_frac, app_coverage].loc[quantity_to_plot].index.values * priv_prof_frac,
        #          pt['mean', obs_depth, priv_prof_frac, app_coverage].loc[quantity_to_plot].values / pt['mean', obs_depth, priv_prof_frac, app_coverage].loc[quantity_to_plot].values[1],
        #          linewidth=1.5, color="#f278cc", alpha=0.35)

        quantity_to_plot = "HerdImmunityType2"
        ax2.plot(pt['mean', obs_depth, priv_prof_frac, app_coverage].loc[quantity_to_plot].index.values * priv_prof_frac,
                 pt['mean', obs_depth, priv_prof_frac, app_coverage].loc[quantity_to_plot].values / pt['mean', obs_depth, priv_prof_frac, app_coverage].loc[quantity_to_plot].values[1],
                 linewidth=1.5, color="#f278cc", alpha=0.35)


        obs_depth = 2
        priv_prof_frac = 0.333333  # options: 0.333333 0.500000 0.66666
        app_coverage = 0.0005      # options: 0.01 0.025 0.005 0.001 0.0005 0.0001

        quantity_to_plot = "ObsCompRelSize"
        ax3.plot(pt['mean', obs_depth, priv_prof_frac, app_coverage].loc[quantity_to_plot].index.values * priv_prof_frac,
                 pt['mean', obs_depth, priv_prof_frac, app_coverage].loc[quantity_to_plot].values,
                 linewidth=1.5, color="#78ccf2", alpha=0.35)

        quantity_to_plot = "CounterCulture"
        ax4.plot(pt['mean', obs_depth, priv_prof_frac, app_coverage].loc[quantity_to_plot].index.values * priv_prof_frac,
                 pt['mean', obs_depth, priv_prof_frac, app_coverage].loc[quantity_to_plot].values,
                 linewidth=1.5, color="#f29e78", alpha=0.35)

        # quantity_to_plot = "HerdImmunityType1"
        # ax5.plot(pt['mean', obs_depth, priv_prof_frac, app_coverage].loc[quantity_to_plot].index.values * priv_prof_frac,
        #          pt['mean', obs_depth, priv_prof_frac, app_coverage].loc[quantity_to_plot].values / pt['mean', obs_depth, priv_prof_frac, app_coverage].loc[quantity_to_plot].values[1],
        #          linewidth=1.5, color="#f278cc", alpha=0.35)

        quantity_to_plot = "HerdImmunityType2"
        ax5.plot(pt['mean', obs_depth, priv_prof_frac, app_coverage].loc[quantity_to_plot].index.values * priv_prof_frac,
                 pt['mean', obs_depth, priv_prof_frac, app_coverage].loc[quantity_to_plot].values / pt['mean', obs_depth, priv_prof_frac, app_coverage].loc[quantity_to_plot].values[1],
                 linewidth=1.5, color="#f278cc", alpha=0.35)


# Il faudra peut-être ajuster cette courbe.
ax1.plot(np.arange(0,0.34,0.01), (2/3)*np.arange(0,0.34,0.01), linestyle = '--', color='black', linewidth = 2, label = 'Fraction of prevented data flow')
ax1.legend(loc="upper left", ncol=1, prop={"size": "x-small"})


ax0.set_xlabel("Adoption of distributed consent")
ax1.set_xlabel("Adoption of distributed consent")
ax2.set_xlabel("Adoption of distributed consent")
ax3.set_xlabel("Adoption of distributed consent")
ax4.set_xlabel("Adoption of distributed consent")
ax5.set_xlabel("Adoption of distributed consent")


ax0.set_ylabel("Fraction of observed individuals")
ax3.set_ylabel("Fraction of observed individuals")
ax1.set_ylabel("Largest unobserved component")
ax4.set_ylabel("Largest unobserved component")
# ax2.set_ylabel("Fraction of observed type 1")
# ax5.set_ylabel("Fraction of observed type 1")
ax2.set_ylabel("Fraction of observed type 2")
ax5.set_ylabel("Fraction of observed type 2")

ax0.set_ylim(bottom=-0.025, top=1.025)
ax1.set_ylim(bottom=-0.025, top=1.025)
ax2.set_ylim(bottom=-0.025, top=1.025)
ax3.set_ylim(bottom=-0.025, top=1.025)
ax4.set_ylim(bottom=-0.025, top=1.025)
ax5.set_ylim(bottom=-0.025, top=1.025)


fig.savefig("figure2.pdf", bbox_inches='tight')
