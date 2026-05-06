import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
from matplotlib.colors import to_rgba
from scipy.ndimage import gaussian_filter1d as gaussf

import pyratbay.io as io
import pyratbay.constants as pc
import pyratbay.spectrum as ps
import chemcat as cat

matplotlib.rcParams['axes.labelpad'] = 1.0
font_dir = ['/home/pcubillos/tmp/fonts']
for font in font_manager.findSystemFonts(font_dir):
    font_manager.fontManager.addfont(font)

# Set font family globally
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial'] + matplotlib.rcParams['font.sans-serif']


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


def plot():
    with np.load('WASP69b_transmission_spectra.npz') as d:
        spectra = d['spectra']
        contribution_spectra = d['contribution_spectra']
        wl = d['wl']
        molecs = d['molecs']
        models = d['models']
        quenched_vmr = d['quenched_vmr']
    nmodels = len(models)
    nmol = len(molecs)
    molecs[list(molecs).index('potassium_vdw')] = 'K'

    atms = [io.read_atm(model) for model in models]
    species = atms[0][1]
    pressure = atms[0][2]
    temps = np.array([atm[3] for atm in atms])
    vmrs = np.array([atm[4] for atm in atms])

    ## Calculate equal CH4 -- CO boundary
    net = cat.Network(pressure, temps[0], species)
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
        vmr = net.thermochemical_equilibrium(
            temperature=t_iso,
            e_ratio={'C_O': c_to_o},
            metallicity=metallicity,
        )
        methane_to_co = vmr[:,i_CH4] / vmr[:,i_CO]
        if np.any(methane_to_co>1):
            i_balance[i] = np.where(methane_to_co>1)[0][0]
        else:
            i_balance[i] = nlayers-1
    press_rec_059 = gaussf(pressure[i_balance], sigma=3.0)

    # Bin spectra
    bin_wl = ps.constant_resolution_spectrum(0.6, 12.5, 100.0)
    nbins = len(bin_wl)
    bin_depths = np.zeros((nmodels, nbins))
    for k in range(nmodels):
        bin_depths[k] = ps.bin_spectrum(bin_wl, wl, spectra[k]) / pc.percent

    bin_contrib_depths = np.zeros((nmodels, nmol, nbins))
    for k in range(nmodels):
        for j in range(nmol):
            bin_depth = ps.bin_spectrum(bin_wl, wl, contribution_spectra[k,j])
            bin_contrib_depths[k,j] = bin_depth / pc.percent



    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # The plot:

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
    d_ran = [
        (1.54, 1.77),
        (1.53, 1.82),
    ]
    fs = 9.5
    dashes = [
        (),
        (6,1),
    ]
    q_dash = (6,1,2,1)

    c_o_labs = ['C/O = 0.59', 'C/O = 1.05']
    savefile = '../plots/WASP69b_forward_model_atmosphere_3x_036.png'


    fs = 9.75
    fig = plt.figure(1)
    plt.clf()
    fig.set_size_inches(9.0, 4.5)
    plt.subplots_adjust(0.05, 0.08, 0.995, 0.985, hspace=0.3, wspace=0.26)
    # TP
    ax = plt.subplot(331)
    for i in range(nmodels):
        lab = c_o_labs[i]
        ax.plot(temps[i], pressure, c='k', lw=1.5, dashes=dashes[i], label=lab)
    ax.plot(temp_array, press_rec_059, c='red', lw=1.25, dashes=(3,1))
    ax.text(490, 1e-5, 'CO', color='r', fontsize=fs, rotation=-45)
    ax.text(420, 5e-4, 'CH4', color='r', fontsize=fs, rotation=-45)
    ax.set_yscale('log')
    ax.set_yticks(np.logspace(2, -7, 4))
    ax.set_ylim(p_ran)
    ax.set_xlim(425, 1650)
    ax.tick_params(which='both', direction='in', labelsize=fs-1)
    ax.set_ylabel('pressure (bar)', fontsize=fs)
    ax.set_xlabel('temperature (K)', fontsize=fs)
    ax.legend(fontsize=fs-1, loc='upper right')
    # VMRs
    vaxes = []
    for i in range(nmodels):
        ax = plt.subplot(3,3,4+3*i)
        for j,mol in enumerate(species):
            if mol not in colors:
                continue
            col = colors[mol]
            alpha = 1.0
            if mol in ['NH3', 'SO2']:
                alpha = 0.45
                q_vmr = quenched_vmr[k,:,j]
                ax.plot(q_vmr, pressure, lw=1.75, c=col, dashes=q_dash)
            ax.plot(vmrs[i,:,j], pressure, lw=1.75, c=col, alpha=alpha)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_yticks(np.logspace(2, -7, 4))
        ax.set_xticks(10.0**np.arange(-3, -25, -3))
        ax.set_xlim(1e-16, 1e-2)
        ax.set_ylim(p_ran)
        h1, = ax.plot(1,1, lw=1.5, c='k', label='equilibrium')
        h3, = ax.plot(1,1, lw=1.5, c='k', dashes=(5,1,2,1),label='disequilibrium')
        ax.tick_params(which='both', top=True, direction='in', labelsize=fs-1)
        ax.set_ylabel('pressure (bar)', fontsize=fs)
        if i==0:
            ax.legend(fontsize=fs-1, loc='upper left', labelspacing=0.15, framealpha=0.7)
            ax_move(ax, dy=-0.02)
        else:
            ax.set_xlabel('volume mixing ratio', fontsize=fs)
        vaxes.append(ax)
    # Spectra
    for i in range(nmodels):
        ax = plt.subplot(2,3,(2+3*i,3+3*i))
        ax.plot(bin_wl, bin_depths[i], c='black', zorder=300)
        for mol,col in colors.items():
            if mol not in molecs:
                ax.fill_between(
                    bin_wl, 0.0, 0.1, color=col, label=mol,
                    fc=to_rgba(col, alpha=0.25),
                )
                continue
            j = list(molecs).index(mol)
            cont = bin_contrib_depths[i,j]
            ax.fill_between(
                bin_wl, 1.5, cont, color=col, label=mol,
                fc=to_rgba(col, alpha=0.25),
            )
            ax.plot(bin_wl, cont, color=col, lw=1.15, zorder=150)
        ax.text(
            0.01, 0.9, c_o_labs[i], fontsize=fs, weight='bold',
            transform=ax.transAxes,
        )
        ax.set_xscale('log')
        ax.set_ylabel('transit depth (%)', fontsize=fs)
        if i == 0:
            args = dict(fontsize=fs, weight='bold', ha='center')
            ax.legend(
                loc=(0.34,1.02), fontsize=fs-2, ncols=5,
                labelspacing=0.2, framealpha=0.4,
            )
            ax.text(0.81, 1.680, 'K', color=colors['K'], **args)
            ax.text(1.15, 1.665, 'H2O', color=colors['H2O'], **args)
            ax.text(3.40, 1.715, 'CH4', color=colors['CH4'], **args)
            ax.text(3.90, 1.685, 'SO2', color=colors['SO2'], **args)
            ax.text(4.64, 1.745, 'CO2', color=colors['CO2'], **args)
            ax.text(4.82, 1.695, 'CO', color=colors['CO'], **args)
            ax.text(10.8, 1.688, 'NH3', color=colors['NH3'], **args)
            ax_move(ax, dy=-0.065)
        else:
            ax.set_xlabel('wavelength (um)', fontsize=fs)
        ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
        ax.set_xticks([0.7, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0])
        ax.set_xlim(0.725, 12.5)
        ax.set_ylim(d_ran[i])
        ax.tick_params(which='both', right=True, direction='in', labelsize=fs-1)
        ax_connect(vaxes[i], ax, fig)
    plt.savefig(savefile, dpi=300)

