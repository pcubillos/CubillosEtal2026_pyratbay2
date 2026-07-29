import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
from scipy.ndimage import gaussian_filter1d as gaussf

import pyratbay.io as io
import pyratbay.constants as pc
import pyratbay.spectrum as ps
import chemcat as cat

matplotlib.rcParams['axes.labelpad'] = 1.0

def ax_move(ax, dx=0.0, dy=0.0):
    x0, y0, width, height = ax.get_position().bounds
    ax.set_position([x0+dx, y0+dy, width, height])


def ax_connect(ax1, ax2, fig):
   x1, y1, w1, h1 = ax1.get_position().bounds   # bounds of ax1
   x2, y2, w2, h2 = ax2.get_position().bounds   # bounds of ax2

   # Coordinates of the gap polygon (in figure coords)
   # ax1 right side
   ax1_right_bottom = (x1 + w1, y1)
   ax1_right_top    = (x1 + w1, y1 + h1)

   # ax2 left side
   ax2_left_top    = (x2, y2 + h2)
   ax2_left_bottom = (x2, y2)

   # Define polygon vertices in order
   verts = [
       ax1_right_bottom,
       ax1_right_top,
       ax2_left_top,
       ax2_left_bottom
   ]

   # Create and add polygon patch
   poly = plt.Polygon(
       verts,
       closed=True,
       transform=fig.transFigure,
       color='0.9',
       ec='none',
       zorder=-100,
   )
   fig.patches.append(poly)


def main():
    j = 0
    with np.load('WASP69b_transmission_spectra.npz') as d:
        spectrum = d['spectra'][j]
        contribution_spectrum = d['contribution_spectra'][j]
        wl = d['wl']
        molecs = d['molecs']
        models = d['models']
        quenched_vmr = d['quenched_vmr'][j]
    nmol = len(molecs)
    molecs[list(molecs).index('potassium_vdw')] = 'K'

    atm = io.read_atm(models[j])
    species, pressure, temp, vmr = atm[1:5]

    # Nudge SO2 to show well in plot
    quenched_vmr[:,list(species).index('SO2')] = 2.75e-6

    ## Calculate equal CH4 -- CO boundary
    net = cat.Network(pressure, temp, species)
    ntemps = 101
    nlayers = len(pressure)
    temp_array = np.linspace(400, 1800, ntemps)
    c_to_o = 0.59
    metallicity = np.log10(3.0)

    i_CH4 = list(species).index('CH4')
    i_CO = list(species).index('CO')
    i_balance = np.zeros(ntemps, int)
    for i in range(ntemps):
        t_iso = np.tile(temp_array[i], nlayers)
        c_vmr = net.thermochemical_equilibrium(
            temperature=t_iso,
            e_ratio={'C_O': c_to_o},
            metallicity=metallicity,
        )
        methane_to_co = c_vmr[:,i_CH4] / c_vmr[:,i_CO]
        if np.any(methane_to_co>1):
            i_balance[i] = np.where(methane_to_co>1)[0][0]
        else:
            i_balance[i] = nlayers-1
    press_rec_059 = gaussf(pressure[i_balance], sigma=3.0)

    # Bin spectra
    bin_wl = ps.constant_resolution_spectrum(0.6, 12.5, 100.0)
    nbins = len(bin_wl)
    bin_depth = ps.bin_spectrum(bin_wl, wl, spectrum) / pc.percent

    bin_contrib_depths = np.zeros((nmol, nbins))
    for j in range(nmol):
        depth = ps.bin_spectrum(bin_wl, wl, contribution_spectrum[j])
        bin_contrib_depths[j] = depth / pc.percent


    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # The plot

    # Species to show
    colors = {
        'H2O': 'xkcd:blue',
        'CH4': 'xkcd:goldenrod',
        'CO': 'xkcd:green',
        'CO2': 'red',
        'SO2': 'deepskyblue',
        'NH3': 'magenta',
        'H2S': 'xkcd:indigo',
        'K': '0.6',
        'C2H2': 'darkorange',
        'HCN': 'green',
    }
    nmol = len(species)
    p_ran = 1e2, 3e-8

    savefile = '../plots/fig_WASP69b_model_atmosphere_3x_solar_032.png'

    fs = 10.0
    args = dict(fontsize=fs-1, weight='bold', ha='center')
    y1 = 0.18
    y2 = 0.985 - y1
    y3 = 0.85 - y1

    x0  = 0.054
    dx0 = 0.122
    x1  = 0.18
    dx1 = 0.23
    x2  = 0.47
    dx2 = 0.529


    fig = plt.figure(1)
    plt.clf()
    fig.set_size_inches(9.0, 2.25)
    # TP
    ax = plt.axes([x0, y1, dx0, y2])
    ax.plot(temp, pressure, c='k', lw=1.5)
    ax.plot(temp_array, press_rec_059, c='red', lw=1.25, dashes=(3,1))
    ax.text(490, 3e-5, 'CO', color='r', fontsize=fs-1.5, rotation=-70)
    ax.text(395, 5e-4, r'CH$_4$', color='r', fontsize=fs-1.5, rotation=-70)
    ax.set_yscale('log')
    ax.set_yticks(np.logspace(2, -7, 4))
    ax.set_ylim(p_ran)
    ax.set_xlim(425, 1675)
    ax.tick_params(which='both', direction='in', labelsize=fs-1)
    ax.set_ylabel('pressure (bar)', fontsize=fs, labelpad=0.0)
    ax.set_xlabel('temperature (K)', fontsize=fs)
    # VMRs
    ax = plt.axes([x1, y1, dx1, y2])
    for j,mol in enumerate(species):
        if mol not in colors:
            continue
        col = colors[mol]
        alpha = 1.0
        if mol in ['NH3', 'SO2']:
            alpha = 0.45
            q_vmr = quenched_vmr[:,j]
            ax.plot(q_vmr, pressure, lw=1.75, c=col, dashes=(6,1))
        ax.plot(vmr[:,j], pressure, lw=1.75, c=col, alpha=alpha)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xticks(10.0**np.arange(-3, -25, -3))
    ax.set_yticks(np.logspace(2, -7, 4))
    ax.set_yticklabels([])
    ax.set_xlim(2e-16, 0.5e-2)
    ax.set_ylim(p_ran)
    h1, = ax.plot(1,1, lw=1.5, c='k', label='equilibrium')
    h3, = ax.plot(1,1, lw=1.5, c='k', dashes=(4,1), label='disequilibrium')
    ax.tick_params(which='both', top=True, direction='in', labelsize=fs-1)
    ax.legend(fontsize=fs-1, loc='upper left', framealpha=0.7)
    ax.set_xlabel('volume mixing ratio', fontsize=fs)
    plt.savefig(savefile, dpi=300)
    # Spectra
    ax = plt.axes([x2, y1, dx2, y3])
    ax.clear()
    ax.plot(bin_wl, bin_depth, c='black', lw=1.35, zorder=300)
    for mol,col in colors.items():
        if mol not in molecs:
            ax.fill_between(
                bin_wl, 0.0, 0.1, color=col, label=mol,
                fc=to_rgba(col, alpha=0.20),
            )
            continue
        j = list(molecs).index(mol)
        cont = bin_contrib_depths[j]
        ax.fill_between(
            bin_wl, 1.5, cont, color=col, label=mol,
            fc=to_rgba(col, alpha=0.20),
        )
        ax.plot(bin_wl, cont, color=col, lw=1.0, zorder=150)
    ax.set_xscale('log')
    ax.set_ylabel('transit depth (%)', fontsize=fs)
    ax.legend(loc=(0.0,1.01), fontsize=fs-2.5, ncols=5, labelspacing=0.15, borderpad=0.3)
    ax.text(0.81, 1.680, 'K', color=colors['K'], **args)
    ax.text(1.17, 1.665, r'H$\mathbf{_2}$O', color=colors['H2O'], **args)
    ax.text(3.33, 1.715, r'CH$\mathbf{_4}$', color=colors['CH4'], **args)
    ax.text(3.87, 1.685, r'SO$\mathbf{_2}$', color=colors['SO2'], **args)
    ax.text(4.66, 1.742, r'CO$\mathbf{_2}$', color=colors['CO2'], **args)
    ax.text(4.87, 1.695, 'CO', color=colors['CO'], **args)
    ax.text(10.7, 1.68, r'NH$\mathbf{_3}$', color=colors['NH3'], **args)
    ax.set_xlabel(r'wavelength ($\mathrm{\mu}$m)', fontsize=fs)
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.set_xticks([0.7, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0])
    ax.set_xlim(0.74, 12.0)
    ax.set_ylim(1.54, 1.76)
    ax.tick_params(which='both', right=True, direction='in', labelsize=fs-1)
    plt.savefig(savefile, dpi=300)


if __name__ == '__main__':
    main()
