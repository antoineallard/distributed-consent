# -*- coding: utf-8 -*-
# @author: Antoine Allard <antoineallard.info>

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


fig = plt.figure(constrained_layout=True, figsize=(24, 11))

gs = matplotlib.gridspec.GridSpec(nrows=2,
                                  ncols=4,
                                  figure=fig,
                                  width_ratios=[1, 1, 1, 1],
                                  height_ratios=[1, 1],
                                  wspace=0.05, hspace=0.05
                                  )

ax0 = fig.add_subplot(gs[0, 0])
ax1 = fig.add_subplot(gs[0, 1])
ax2 = fig.add_subplot(gs[0, 2])
ax3 = fig.add_subplot(gs[0, 3])

ax4 = fig.add_subplot(gs[1, 0])
ax5 = fig.add_subplot(gs[1, 1])
ax6 = fig.add_subplot(gs[1, 2])
ax7 = fig.add_subplot(gs[1, 3])


for filename in glob.glob("../results/single_layer/*.pkl"):

    name = os.path.basename(filename).rsplit(".", 1)[-2]

    # Ignores networks with fewer than 2000 nodes.
    if name not in ['Swarthmore42', 'Haverford76', 'Simmons81', 'Reed98', 'Caltech36']:

        pt = pd.read_pickle(filename, compression="xz")

        # =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
        # First row.

        obs_depth = 2
        priv_prof_frac = 0.333333
        app_coverage = 0.01

        quantity_to_plot = "ObsCompRelSize"
        ax0.plot(pt['mean', obs_depth, priv_prof_frac, app_coverage].loc[quantity_to_plot].index.values,# * priv_prof_frac,
                 pt['mean', obs_depth, priv_prof_frac, app_coverage].loc[quantity_to_plot].values,
                 linewidth=1.5, color="#78ccf2", alpha=0.35)

        quantity_to_plot = "CounterCulture"
        ax1.plot(pt['mean', obs_depth, priv_prof_frac, app_coverage].loc[quantity_to_plot].index.values,# * priv_prof_frac,
                 pt['mean', obs_depth, priv_prof_frac, app_coverage].loc[quantity_to_plot].values,
                 linewidth=1.5, color="k", alpha=0.35)

        quantity_to_plot = "HerdImmunityType1"
        ax2.plot(pt['mean', obs_depth, priv_prof_frac, app_coverage].loc[quantity_to_plot].index.values,# * priv_prof_frac,
                  pt['mean', obs_depth, priv_prof_frac, app_coverage].loc[quantity_to_plot].values / pt['mean', obs_depth, priv_prof_frac, app_coverage].loc[quantity_to_plot].values[1],
                  linewidth=1.5, color="#f29e78", alpha=0.35)

        quantity_to_plot = "HerdImmunityType2"
        ax3.plot(pt['mean', obs_depth, priv_prof_frac, app_coverage].loc[quantity_to_plot].index.values,# * priv_prof_frac,
                 pt['mean', obs_depth, priv_prof_frac, app_coverage].loc[quantity_to_plot].values / pt['mean', obs_depth, priv_prof_frac, app_coverage].loc[quantity_to_plot].values[1],
                 linewidth=1.5, color="#f278cc", alpha=0.35)


        # =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
        # Second row.

        obs_depth = 2
        priv_prof_frac = 0.66666  # options: 0.333333 0.500000 0.66666
        app_coverage = 0.01      # options: 0.01 0.025 0.005 0.001 0.0005 0.0001

        quantity_to_plot = "ObsCompRelSize"
        ax4.plot(pt['mean', obs_depth, priv_prof_frac, app_coverage].loc[quantity_to_plot].index.values,# * priv_prof_frac,
                 pt['mean', obs_depth, priv_prof_frac, app_coverage].loc[quantity_to_plot].values,
                 linewidth=1.5, color="#78ccf2", alpha=0.35)

        quantity_to_plot = "CounterCulture"
        ax5.plot(pt['mean', obs_depth, priv_prof_frac, app_coverage].loc[quantity_to_plot].index.values,# * priv_prof_frac,
                 pt['mean', obs_depth, priv_prof_frac, app_coverage].loc[quantity_to_plot].values,
                 linewidth=1.5, color="k", alpha=0.35)

        quantity_to_plot = "HerdImmunityType1"
        ax6.plot(pt['mean', obs_depth, priv_prof_frac, app_coverage].loc[quantity_to_plot].index.values,# * priv_prof_frac,
                  pt['mean', obs_depth, priv_prof_frac, app_coverage].loc[quantity_to_plot].values / pt['mean', obs_depth, priv_prof_frac, app_coverage].loc[quantity_to_plot].values[1],
                  linewidth=1.5, color="#f29e78", alpha=0.35)

        quantity_to_plot = "HerdImmunityType2"
        ax7.plot(pt['mean', obs_depth, priv_prof_frac, app_coverage].loc[quantity_to_plot].index.values,# * priv_prof_frac,
                 pt['mean', obs_depth, priv_prof_frac, app_coverage].loc[quantity_to_plot].values / pt['mean', obs_depth, priv_prof_frac, app_coverage].loc[quantity_to_plot].values[1],
                 linewidth=1.5, color="#f278cc", alpha=0.35)


ax0.plot(np.NaN, np.NaN, color=None, linewidth = 0, label = r'Taste for privacy = 33\%')
ax1.plot(np.NaN, np.NaN, color=None, linewidth = 0, label = r'Taste for privacy = 33\%')
ax2.plot(np.NaN, np.NaN, color=None, linewidth = 0, label = r'Taste for privacy = 33\%')
ax3.plot(np.NaN, np.NaN, color=None, linewidth = 0, label = r'Taste for privacy = 33\%')
ax4.plot(np.NaN, np.NaN, color=None, linewidth = 0, label = r'Taste for privacy = 66\%')
ax5.plot(np.NaN, np.NaN, color=None, linewidth = 0, label = r'Taste for privacy = 66\%')
ax6.plot(np.NaN, np.NaN, color=None, linewidth = 0, label = r'Taste for privacy = 66\%')
ax7.plot(np.NaN, np.NaN, color=None, linewidth = 0, label = r'Taste for privacy = 66\%')
ax1.plot(np.arange(0,1.0,0.01), (2/9)*np.arange(0,1.0,0.01), linestyle = '--', color='black', linewidth = 2, label = 'Fraction of prevented data flow')
ax5.plot(np.arange(0,1.0,0.01), (2/9)*np.arange(0,1.0,0.01), linestyle = '--', color='black', linewidth = 2, label = 'Fraction of prevented data flow')
ax0.legend(loc="lower left", ncol=1, handlelength=0, prop={"size": "x-small"})
ax1.legend(loc="upper left", ncol=1, handlelength=1,  prop={"size": "x-small"})
ax2.legend(loc="lower left", ncol=1, handlelength=0,  prop={"size": "x-small"})
ax3.legend(loc="lower left", ncol=1, handlelength=0,  prop={"size": "x-small"})
ax4.legend(loc="lower left", ncol=1, handlelength=0,  prop={"size": "x-small"})
ax5.legend(loc="upper left", ncol=1, handlelength=1,  prop={"size": "x-small"})
ax6.legend(loc="lower left", ncol=1, handlelength=0,  prop={"size": "x-small"})
ax7.legend(loc="lower left", ncol=1, handlelength=0,  prop={"size": "x-small"})


ax0.annotate(r'a.', xy= (0.99,0.99), xycoords='axes fraction', horizontalalignment='right', verticalalignment='top', fontsize='small')
ax1.annotate(r'b.', xy= (0.99,0.99), xycoords='axes fraction', horizontalalignment='right', verticalalignment='top', fontsize='small')
ax2.annotate(r'c.', xy= (0.99,0.99), xycoords='axes fraction', horizontalalignment='right', verticalalignment='top', fontsize='small')
ax3.annotate(r'd.', xy= (0.99,0.99), xycoords='axes fraction', horizontalalignment='right', verticalalignment='top', fontsize='small')
ax4.annotate(r'e.', xy= (0.99,0.99), xycoords='axes fraction', horizontalalignment='right', verticalalignment='top', fontsize='small')
ax5.annotate(r'f.', xy= (0.99,0.99), xycoords='axes fraction', horizontalalignment='right', verticalalignment='top', fontsize='small')
ax6.annotate(r'g.', xy= (0.99,0.99), xycoords='axes fraction', horizontalalignment='right', verticalalignment='top', fontsize='small')
ax7.annotate(r'h.', xy= (0.99,0.99), xycoords='axes fraction', horizontalalignment='right', verticalalignment='top', fontsize='small')


ax0.set_xlabel("Distributed consent / taste for privacy")
ax1.set_xlabel("Distributed consent / taste for privacy")
ax2.set_xlabel("Distributed consent / taste for privacy")
ax3.set_xlabel("Distributed consent / taste for privacy")
ax4.set_xlabel("Distributed consent / taste for privacy")
ax5.set_xlabel("Distributed consent / taste for privacy")
ax6.set_xlabel("Distributed consent / taste for privacy")
ax7.set_xlabel("Distributed consent / taste for privacy")


ax0.set_ylabel("Fraction of observed individuals")
ax4.set_ylabel("Fraction of observed individuals")
ax1.set_ylabel("Largest unobserved component")
ax5.set_ylabel("Largest unobserved component")
ax2.set_ylabel("Relative fraction of observed type 1")
ax6.set_ylabel("Relative fraction of observed type 1")
ax3.set_ylabel("Relative fraction of observed type 2")
ax7.set_ylabel("Relative fraction of observed type 2")


ax0.set_ylim(bottom=-0.025, top=1.025)
ax1.set_ylim(bottom=-0.025, top=1.025)
ax2.set_ylim(bottom=0.8, top=1.005)
ax3.set_ylim(bottom=0.25, top=1.02)
ax4.set_ylim(bottom=-0.025, top=1.025)
ax5.set_ylim(bottom=-0.025, top=1.025)
ax6.set_ylim(bottom=0.8, top=1.005)
ax7.set_ylim(bottom=0.25, top=1.02)


fig.savefig("figure3.pdf", bbox_inches='tight')
fig.savefig("figure3.svg", bbox_inches='tight')