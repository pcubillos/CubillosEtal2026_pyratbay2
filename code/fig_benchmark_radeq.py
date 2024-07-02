import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import scipy.interpolate as si
import pyratbay.io as io
import pyratbay.constants as pc
import pyratbay.spectrum as ps
import os

import matplotlib.font_manager as font_manager

# Add every font at the specified location
font_dir = ['/home/pcubillos/tmp/fonts']
for font in font_manager.findSystemFonts(font_dir):
    font_manager.fontManager.addfont(font)


# Set font family globally
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']



def read_helios(folder):
    path, root = os.path.split(folder)

    tp_file = f'{folder}/{root}_tp.dat'
    tp = np.loadtxt(tp_file, unpack=True, skiprows=3, usecols=(0,1,2,3,4))
    h_temp = tp[1]
    h_press = tp[2] / pc.bar  # bar
    h_alt = tp[3] / pc.rjup

    # Spectral fluxes given in [erg s^-1 cm^-3].
    flux_file = f'{folder}/{root}_TOA_flux_eclipse.dat'
    data = np.loadtxt(flux_file, unpack=True, skiprows=4)
    h_wl = data[1]
    h_dwl = data[3]
    h_fstar = data[4]    # erg s^-1 cm^-3
    h_fplanet = data[5]  # erg s^-1 cm^-3
    h_fpfs = data[6]

    return (
        h_press, h_temp, h_alt,
        h_wl, h_dwl, h_fstar, h_fplanet, h_fpfs,
    )


def read_pyratbay(folder):
    files = [
        'WASP107b_radeq_benchmark_100K_Tirr.atm',
        'WASP107b_radeq_benchmark_350K_Tirr.atm',
        'WASP39b_radeq_benchmark_01x_solar_no_heat.atm',
        'WASP39b_radeq_benchmark_50x_solar_no_heat.atm',
        'WASP121b_radeq_benchmark_with_TiO_VO_no_heat.atm',
        'WASP121b_radeq_benchmark_no_TiO_VO_no_heat.atm',
    ]
    nfiles = len(files)

    pressures = []
    temps = []
    vmrs = []
    spectra = []
    for i in range(nfiles):
        afile = files[i]
        if not os.path.exists(f'{folder}/{afile}'):
            continue
        atm = io.read_atm(f'{folder}/{afile}')
        species = atm[1]
        press = atm[2]
        temp = atm[3]
        vmr = atm[4]
        temps.append(temp)
        pressures.append(press)
        vmrs.append(vmr)

        wl, spec = io.read_spectrum(afile.replace('atm', 'dat'), wn=False)
        spectra.append(spec)

    return press, temps, vmrs, species, wl, spectra


def read_stars():
    wn, spec = io.read_spectrum('WASP107b_radeq_benchmark_100K_Tirr.dat')

    s_files = [
        '../inputs/phoenix_WASP107b_4400K_m00_logg4.5.dat',
        '../inputs/phoenix_WASP39b_5500K_m00_logg4.5.dat',
        '../inputs/phoenix_WASP121b_6400K_m00_logg4.5.dat',
    ]

    fluxes = []
    for i in range(len(s_files)):
        starflux, starwn, star_temps = io.read_spectra(s_files[i])
        sinterp = si.interp1d(starwn, starflux)
        fluxes.append(sinterp(wn))

    return 1e4/wn, fluxes


def main():

    h_files = [
        '../inputs/helios_v01/wasp-107b-f0.667-solar-tint_100',
        '../inputs/helios_v01/wasp-107b-f0.667-solar-tint_350',
        '../inputs/helios_v01/wasp-39b-f0.667-solar-no_tio_vo-tint_0',
        '../inputs/helios_v02/wasp-39b-f0.667-mh_+1.7dex-no_tio_vo-tint_0',
        '../inputs/helios_v02/wasp-121b-f0.667-solar-tint_0',
        '../inputs/helios_v02/wasp-121b-f0.667-solar-no_tio_vo-tint_0',
    ]
    nhelios = len(h_files)

    h_temps = []
    h_fstar = []
    h_fplanet = []
    h_fpfs = []
    for i in range(nhelios):
        data = read_helios(h_files[i])
        h_press = data[0]
        h_temps.append(data[1])
        h_fstar.append(data[5])
        h_fplanet.append(data[6])
        h_fpfs.append(data[7])

    press, temps, vmrs, species, wl, spectra = read_pyratbay('.')


    labels = [
        r'$T_{\rm int} = 100$ K',
        r'$T_{\rm int} = 350$ K',
        '1x solar',
        '50x solar',
        'with TiO/VO',
        'no TiO/VO',
    ]
    planets_short = [
        'WASP-107 b',
        'WASP-39 b',
        'WASP-121 b',
    ]
    planets = [
        r'WASP-107 b ($T_{\rm eq} = 750 {\rm K}$)',
        r'WASP-39 b ($T_{\rm eq} = 1200 {\rm K}$)',
        r'WASP-121 b ($T_{\rm eq}=2300 {\rm K}$)',
    ]
    nplanets = len(planets)

    rstars = np.array([0.67, 0.94, 1.458]) * pc.rsun
    rplanets = np.array([0.943, 1.280, 1.753]) * pc.rjup

    # The stars:
    wl, starfluxes = read_stars()
    bin_wl = ps.constant_resolution_spectrum(0.15, 33.0, 150)
    nbins = len(bin_wl)
    rs_smaxis = [
        0.670 * pc.rsun / (0.055 * pc.au),
        0.940 * pc.rsun / (0.048 * pc.au),
        1.458 * pc.rsun / (0.026 * pc.au),
    ]
    beta_irr = 0.667
    fstar = [
        ps.bin_spectrum(bin_wl, wl, starflux * beta_irr * rs_a**2)
        for starflux, rs_a in zip(starfluxes, rs_smaxis)
    ]


    bin_spectra = np.zeros((2*nplanets, nbins))
    for i in range(2*nplanets):
        starflux = starfluxes[i//2]
        rprs = (rplanets/rstars)[i//2]
        bin_spectra[i] = ps.bin_spectrum(bin_wl, wl, spectra[i]/starflux*rprs**2)

    fs = 11.0
    lw = 1.5
    hlw = 1.85
    dashes = (7,1)


    sed_cols = [
        'cornflowerblue',
        'xkcd:green',
        'salmon',
    ]

    # SED at top of atmosphere
    fig = plt.figure(0)
    fig.set_size_inches(5.0, 3.75)
    plt.subplots_adjust(0.15, 0.11, 0.995, 0.99)
    plt.clf()
    ax = plt.subplot(111)
    for i in reversed(range(nplanets)):
        ax.plot(bin_wl, fstar[i], c=sed_cols[i], label=planets_short[i])
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.set_xticks([0.3, 1.0, 3.0, 10, 30])
    ax.set_xlim(0.15, 30.0)
    ax.set_ylim(1e-2, 4e5)
    ax.set_ylabel(
        'SED at top of atmosphere\n(erg s$^{-1}$ cm$^{-2}$ cm)',
        fontsize=fs, labelpad=0,
    )
    ax.set_xlabel('wavelength (um)', fontsize=fs)
    ax.tick_params(which='both', direction='in', labelsize=fs-1)
    ax.legend(loc='lower right', fontsize=fs)
    plt.savefig('../plots/benchmark_radeq_sed.png', dpi=300)

    lw = 1.25
    hlw = 1.95
    colors = [
        'mediumblue',
        'cornflowerblue',
        'limegreen',
        'darkgreen',
        'darkred',
        'tomato',
    ]

    temp_ranges = [
        (450, 3200),
        (250, 2600),
        (1500, 4300),
    ]
    leg_loc = [
        'upper right',
        'lower left',
        'upper left',
    ]

    depth_max = 4250, 5000, 8750
    nrows = 4
    molecs = {
        'H2O': 'xkcd:blue',
        'CH4': 'gold',
        'CO': 'xkcd:green',
        'CO2': 'red',
        'NH3': 'chocolate',
        'N2': 'darkorange',
        'Na': '0.25',
        'K': '0.65',
        'TiO': 'cornflowerblue',
        'VO': 'green',
        'H2': 'xkcd:violet',
        'He': 'yellowgreen',
        'H': 'xkcd:pink',
        'e': 'deepskyblue',
    }
    vmr_range = (1e-10, 2)

    dx = 0.25
    x0 = 0.062
    deltax = 0.08
    x1 = x0+dx+deltax
    x2 = x1+dx+0.075
    dy = 0.19
    y0 = 0.73
    deltay = 0.32
    dx2 = 0.27
    dy2 = 0.12

    fig = plt.figure(1)
    fig.set_size_inches(9.0, 10.0)
    plt.subplots_adjust(0.07, 0.05, 0.99, 0.97, wspace=0.2, hspace=0.44)
    plt.clf()
    # Temperatures
    for i in range(nplanets):
        j = 2*i
        ax = plt.axes([x0, y0-i*deltay, dx, dy])
        ax.plot(temps[j], press, lw=lw, c=colors[j], label=labels[j])
        ax.plot(temps[j+1], press, lw=lw, c=colors[j+1], label=labels[j+1])
        ax.plot(h_temps[j], h_press, lw=hlw, c=colors[j], dashes=dashes)
        ax.plot(h_temps[j+1], h_press, lw=hlw, c=colors[j+1], dashes=dashes)
        ax.set_yscale('log')
        ax.tick_params(which='both', direction='in', top=True, labelsize=fs-1)
        ax.set_ylim(np.amax(press), np.amin(press))
        ax.set_xlim(temp_ranges[i])
        ax.set_title(planets[i], fontsize=fs+1, weight='bold')
        ax.legend(loc=leg_loc[i], fontsize=fs-1, handlelength=1.5)
        ax.set_xlabel('temperature (K)', fontsize=fs)
        ax.set_ylabel('pressure (bar)', fontsize=fs)
        if i == 0:
            legend = ax.get_legend()
            code_handles = [
                Line2D([], [], color='k', lw=1.5, dashes=(), label='Pyrat Bay'),
                Line2D([], [], color='k', lw=1.5, dashes=(6,2), label='HELIOS'),
            ]
            ax.legend(handles=code_handles, loc=(1.35,1.05))
            ax.add_artist(legend)
    # Spectra
    for i in range(nplanets):
        ax = plt.axes([x1, y0-i*deltay, dx, dy])
        ax.plot(bin_wl, bin_spectra[2*i]/pc.ppm, c=colors[2*i])
        ax.plot(bin_wl, bin_spectra[2*i+1]/pc.ppm, c=colors[2*i+1])
        ax.set_xscale('log')
        ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
        ax.set_xticks([0.3, 1.0, 3.0, 10])
        ax.set_ylim(0.0, depth_max[i])
        ax.set_xlim(0.5, 20.0)
        ax.tick_params(which='both', top=True, direction='in', labelsize=fs-1)
        ax.set_xlabel('wavelength (um)', fontsize=fs)
        ax.set_ylabel(r'$F_{\rm p} / F_{\rm s}$ (ppm)', fontsize=fs)
    # VMRs
    for i in range(nplanets):
        for j in range(2):
            ax = plt.axes([x2, 0.812-i*deltay-j*(dy2+0.005), dx2, dy2])
            for mol,col in molecs.items():
                for ion in ['', '+', '-']:
                    name = mol+ion
                    label = mol if (name==mol or name=='e-') else None
                    dash = (6,1.4) if ion=='+' else (2,1) if ion=='-' else ()
                    if name not in species:
                        continue
                    k = list(species).index(name)
                    ax.loglog(
                        vmrs[2*i+j][:,k], press, label=label, c=col,
                        lw=2.0, dashes=dash,
                    )
            ax.set_yscale('log')
            ax.set_xlim(vmr_range)
            ax.set_ylim(np.amax(press), np.amin(press))
            ax.set_yticks(np.logspace(-8,0, 3))
            ax.tick_params(
                axis='both', which='both', direction='in', labelsize=fs-1,
                length=6.0, width=1.25,
            )
            ax.text(
                0.03, 0.93, labels[2*i+j],
                fontsize=fs-1, transform=ax.transAxes, va='top',
                bbox=dict(boxstyle='Round', facecolor='w', alpha=0.78),
            )
            ax.set_ylabel('pressure (bar)', fontsize=fs)
            if i == 0 and j == 0:
                ax.legend(
                    loc=(-0.4, 1.08), ncols=5, fontsize=fs-1,
                    labelspacing=0.1,
                    columnspacing=1.0,
                    handlelength=1.5,
                )
            if j == 0:
                ax.set_xticklabels([])
        ax.set_xlabel('volume mixing ratio', fontsize=fs)
    plt.savefig('../plots/benchmark_radeq_atmospheres.png', dpi=300)



