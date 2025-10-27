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
    molecs[list(molecs).index('sodium_vdw')] = 'Na'

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

    # Species to showcase:
    colors = {
        'H2O': 'xkcd:blue',
        'CH4': 'xkcd:goldenrod',
        'CO': 'xkcd:green',
        'CO2': 'red',
        'SO2': 'deepskyblue',
        'NH3': 'magenta',
        'H2S': 'xkcd:indigo',
        'Na': '0.6',
        'K': '0.0',
        'C2H2': 'darkorange',
    }
    nmol = len(species)

    p_ran = 1e2, 3e-8
    d_ran = [
        (1.53, 1.76),
        (1.52, 1.82),
    ]
    fs = 9.5
    dashes = [
        (),
        (6,1),
    ]
    q_dash = (6,1,2,1)
    quench_text = 'NH$_3$ quenching pressure'

    c_o_labs = ['C/O = 0.59', 'C/O = 1.10']
    savefile = '../plots/WASP69b_forward_model_atmosphere_3x_036.png'


    fs = 9.5
    fig = plt.figure(1)
    plt.clf()
    fig.set_size_inches(9.0, 4.5)
    plt.subplots_adjust(0.055, 0.075, 0.995, 0.985, hspace=0.3, wspace=0.26)
    # TP
    ax = plt.subplot(331)
    for i in range(nmodels):
        ax.plot(temps[i], pressure, c='k', lw=1.5,dashes=dashes[i])
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
    # VMRs
    for i in range(nmodels):
        ax = plt.subplot(3,3,4+3*i)
        for j,mol in enumerate(species):
            if mol not in colors:
                continue
            col = colors[mol]
            ax.plot(vmrs[i,:,j], pressure, lw=1.75, c=col)
            if mol in ['NH3', 'SO2']:
                q_vmr = quenched_vmr[k,:,j]
                ax.plot(q_vmr, pressure, lw=1.75, c=col, dashes=q_dash)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_yticks(np.logspace(2, -7, 4))
        ax.set_xticks(10.0**np.arange(-3, -25, -3))
        ax.set_xlim(1e-16, 1e-2)
        ax.set_ylim(p_ran)
        h1, = ax.plot(1,1, lw=1.5, c='k', label='equilibrium')
        h3, = ax.plot(1,1, lw=1.5, c='k', dashes=(5,1,2,1),label='quenched')
        ax.axhspan(6.0, 0.5, color='pink', alpha=0.5, zorder=50)
        ax.tick_params(which='both', top=True, direction='in', labelsize=fs-1)
        ax.set_ylabel('pressure (bar)', fontsize=fs)
        if i==0:
            ax.legend(fontsize=fs-1, loc='upper left', labelspacing=0.15, framealpha=0.7)
            ax.text(2e-16, 3.0, quench_text, fontsize=fs-1, zorder=51)
            ax_move(ax, dy=-0.03)
        else:
            ax.set_xlabel('volume mixing ratio', fontsize=fs)
    # Spectra
    for i in range(nmodels):
        ax = plt.subplot(2,3,(2+3*i,3+3*i))
        ax.plot(bin_wl, bin_depths[i], c='black', zorder=300)
        for mol,col in colors.items():
            if mol not in molecs:
                continue
            j = list(molecs).index(mol)
            cont = bin_contrib_depths[i,j]
            mol = mol.replace('Na', 'Na+K')
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
            ax.legend(
                loc=(0.4,1.02), fontsize=fs-2, ncols=4,
                labelspacing=0.2, framealpha=0.4,
            )
            #ax.text(0.17, 0.71, 'CH4',  fontsize=fs, transform=ax.transAxes)
            #ax.text(0.275,0.485, 'SO2', fontsize=fs, transform=ax.transAxes)
            #ax.text(0.325, 0.54, 'CO2', fontsize=fs, transform=ax.transAxes)
            #ax.text(0.39, 0.45, 'CO',   fontsize=fs, transform=ax.transAxes)
            #ax.text(0.45, 0.47, 'H2O',  fontsize=fs, transform=ax.transAxes)
            #ax.text(0.88, 0.44, 'NH3',  fontsize=fs, transform=ax.transAxes)
            ax_move(ax, dy=-0.065)
        else:
            ax.set_xlabel('wavelength (um)', fontsize=fs)
        ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
        ax.set_xticks([0.7, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0])
        ax.set_xlim(0.6, 12.5)
        ax.set_ylim(d_ran[i])
        ax.tick_params(which='both', right=True, top=True, direction='in')
    plt.savefig(savefile, dpi=300)



