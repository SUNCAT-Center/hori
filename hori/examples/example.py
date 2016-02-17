#!/usr/bin/env python
import os
import matplotlib.pyplot as plt
from matplotlib import rc
from hori.thermo import (AdsorbateThermodynamics,
                         ProtonElectronThermodynamics,
                         GasThermodynamics)
from hori.common import hbond_dict
from hori.pathway import Pathway
from hori.plots import PlotStates
from hori.io import *
import scipy


### OPTIONS ####
save_fig = True  # set to true to save pdf files

# font sizes
font_size = 11
legend_font_size = 11


# CHANGE TO YOUR HORI PATH
set_directory(os.environ['SCRATCH'] + '/usr/pylib/hori/data_QE_BEEF')


####################
## Set up pathway ##
####################


# Define the reaction pathway.  Note that steps are added to the pathway via the method path.addstep
# Defines '*' as the starting point, and 0. as the reference free energy associated with state '*'
path = Pathway('*', 0.)
# - for reactants, + for products. a_X defines X* while, while pe is a proton-electron pair
path.addstep('*', 'CHO*', ['-s', '-CO', '-pe', '+a_CHO'])
path.addstep('CHO*', 'CH2O*', ['-a_CHO', '-pe', '+a_CH2O'])
path.addstep('CH2O*', 'CH3O*', ['-a_CH2O', '-pe', '+a_CH3O'])
path.addstep('CH3O*', 'CH3OH', ['-a_CH3O', '-pe', '+s', '+CH3OH'])

# Thermodynamics
temperature = 298.15  # K

at = AdsorbateThermodynamics(['CHO', 'CH2O', 'CH3O', ''],           # list of adsorbates
                             hbond_dict=hbond_dict(),    # hydrogen bond dictionary (hbond_dict() is the default)
                             temperature=temperature,    # temperature
                             vibrationsmode='generic',  # vibrations - generic uses the same for all surfaces
                             surface='MoS2-vac-no-H2O')              # surface name

voltage = 0  # V vs RHE
# define proton-electron thermodynamics
pet = ProtonElectronThermodynamics(voltage=voltage, temperature=temperature, functional='BEEF')

#
gt = GasThermodynamics([('H2', temperature, 101345),
                        ('CO', temperature, 101345),
                        ('CH3OH', temperature, 101345),
                        ],
                       functional='BEEF')


##########################
## Plotting the pathway ##
##########################

figsize = (6, 4)

fig = plt.figure(4, figsize=figsize)
ax = fig.add_subplot(111)

labels = {
    '*': '$\mathrm{*+CO}$',
    'CHO*': '$\mathrm{CHO*}$',
    'CH2O*': '$\mathrm{CH_2O*}$',
    'CH3O*': '$\mathrm{CH_3O*}$',
    'CH3OH': '$\mathrm{CH_3OH}$',
}

rxn_path = ['*', 'CHO*', 'CH2O*', 'CH3O*', 'CH3OH']

# plot colors
plot1color = 'black'
plot2color = '#e14f4f'

# calculate free energies
# you have to specify the thermodynamics modules to get the pathway to
# calculate free energies as shown
path.calculate_Gs(BG=at.G, mu=gt.G, mu_pe=pet.G)
# plot

# create a PlotStates object, that contains all the plotting information

ps = PlotStates(ax,                     # matplotlib ax object
                path.G,                 # pathway to plot
                halfwidth=0.25,         # customizing the appearance
                textwidth=0.,
                fontsize=font_size,
                color=plot1color,
                dashcolor=plot1color,
                textcolor=plot1color,
                textposition='above',
                text_vspace=0.02)


# starting from the first state, plot each successive step
loc = 0.
prevstate = '*'
for state in rxn_path:
    ps.plotstate(state, loc, labels[state], special_text_vspace=0.05)
    if loc != 0:
        # PlotStates.connect draws a line between two states (prevstate and state)
        # and needs two locations (loc-1 and loc).  Note it can also take a
        # keyword "barrier" to draw a parabolic barrier in between states, but
        # default just draws a line
        ps.connect(prevstate, loc - 1, state, loc, linestyle='-')
    prevstate = state
    loc += 1


# find limiting potential
UL = path.find_limiting_potential(path=rxn_path)
pet.set_voltage(UL[0])  # set voltage to limiting potential

# you have to specify the thermodynamics modules to get the pathway to
# calculate free energies as shown
path.calculate_Gs(BG=at.G, mu=gt.G, mu_pe=pet.G)


# same as before, but now the path object has the updated free energies at the limiting potential.
# then just replot as before:
ps = PlotStates(ax, path.G, halfwidth=0.25, textwidth=0.,
                fontsize=font_size, color=plot2color, dashcolor=plot2color, textcolor=plot2color, textposition='above',
                text_vspace=-0.5)
loc = 1.
prevstate = '*'
for state in rxn_path[1:]:
    ps.plotstate(state, loc, labels[state], special_text_vspace=0.05)
    # PlotStates.connect draws a line between two states (prevstate and state)
    # and needs two locations (loc-1 and loc).  Note it can also take a
    # keyword "barrier" to draw a parabolic barrier in between states, but
    # default just draws a line
    ps.connect(prevstate, loc - 1, state, loc, linestyle='-')
    prevstate = state
    loc += 1


# add labels
plt.subplots_adjust(wspace=1, hspace=0, top=0.92, bottom=0.13, left=0.12, right=0.96)
ax.set_xticks(scipy.arange(-1, len(rxn_path[1:]) + 2, 1))
xtick_labels = ['', '0']
for i in np.arange(1, len(rxn_path[1:]) + 1):
    xtick_labels.append(str(i))

ax.set_xticklabels(xtick_labels)
xmax = len(xtick_labels) - 1.5
plt.xlim([-0.5, xmax])
plt.xlabel('H$^+$ + e$^-$ transfers')
plt.ylabel('Free energy (eV)')

# generate legend
ax1 = fig.add_subplot(111)
ax1.plot([1], color=plot1color, linewidth=2, label="$0 \,\mathrm{V}$ vs. RHE")
ax1.plot([1], color=plot2color, linewidth=2, label="$%s \,\mathrm{V}$ vs. RHE" % round(UL[0], 2))
ax1.legend(loc=3,  ncol=1, shadow=False, prop={'size': legend_font_size})

if save_fig == True:
    fig.savefig('pathway.pdf', bbox_inches='tight')

plt.show()
