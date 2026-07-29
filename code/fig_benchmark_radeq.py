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
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial']


def main():
    # Read helios data
    h_files = [
        '../inputs/helios/wasp107b_tint_000K',
        '../inputs/helios/wasp107b_tint_350K',
        '../inputs/helios/wasp39b_01x_solar',
        '../inputs/helios/wasp39b_50x_solar',
        '../inputs/helios/wasp121b_with_tio',
        '../inputs/helios/wasp121b_no_tio',
    ]

    h_temps = []
    h_fpfs = []
    for i,file in enumerate(h_files):
        h_press, h_temp = np.loadtxt(f'{file}_atm.dat', unpack=True)
        h_wl, fpfs = np.loadtxt(f'{file}_spectrum.dat', unpack=True)
        h_temps.append(h_temp)
        h_fpfs.append(fpfs)


    # Read pyratbay data
    pb_files = [
        'WASP107b_radeq_benchmark_control',
        'WASP107b_radeq_benchmark_tint_350K',
        'WASP39b_radeq_benchmark_control',
        'WASP39b_radeq_benchmark_50x',
        'WASP121b_radeq_benchmark_TiO_VO',
        'WASP121b_radeq_benchmark_control',
    ]

    temps = []
    vmrs = []
    fplanet = []
    for i,file in enumerate(pb_files):
        u, species, press, temp, vmr, radius = io.read_atm(f'{file}.atm')
        temps.append(temp)
        vmrs.append(vmr)
        wl, spec = io.read_spectrum(f'{file}.dat', wn=False)
        fplanet.append(spec)


    # Read stellar SEDs
    s_files = [
        '../inputs/helios_phoenix_WASP107b_4400K_m00_logg4.5.dat',
        '../inputs/helios_phoenix_WASP39b_5500K_m00_logg4.5.dat',
        '../inputs/helios_phoenix_WASP121b_6400K_m00_logg4.5.dat',
    ]
    starfluxes = []
    for file in s_files:
        star_wl, flux = io.read_spectrum(file, wn=False)
        sinterp = si.interp1d(star_wl, flux)
        starfluxes.append(sinterp(wl))

    # Bin all specta to a common wavelenth
    bin_wl = ps.constant_resolution_spectrum(0.15, 30.0, 150.0)
    nbins = len(bin_wl)

    nplanets = len(s_files)
    bin_spectra = np.zeros((2*nplanets, 2, nbins))
    rstars = np.array([0.67, 0.94, 1.458]) * pc.rsun
    rplanets = np.array([0.943, 1.280, 1.753]) * pc.rjup
    for i in range(2*nplanets):
        starflux = starfluxes[i//2]
        rprs = (rplanets/rstars)[i//2]
        fpfs = fplanet[i]/starflux * rprs**2
        bin_spectra[i,0] = ps.bin_spectrum(bin_wl, wl, fpfs)
        bin_spectra[i,1] = ps.bin_spectrum(bin_wl, h_wl, h_fpfs[i])

    # Contribution functions
    #with np.load('contribution_functions_benchmark.npz') as d:
    #    cf_data = d['cf']
    #    cf_wl = d['wl']
    #    cf_press = d['press']
    #median_cf = np.mean(cf_data, axis=2)


    # The big one
    labels = [
        r'$T_{\rm int} =   0$ K',
        r'$T_{\rm int} = 350$ K',
        '1x solar',
        '50x solar',
        'with TiO/VO',
        'no TiO/VO',
    ]
    planets = [
        r'WASP-107 b ($T_{\rm eq} = 750 {\rm K}$)',
        r'WASP-39 b ($T_{\rm eq} = 1200 {\rm K}$)',
        r'WASP-121 b ($T_{\rm eq}=2300 {\rm K}$)',
    ]

    fs = 11.0
    dashes = (7,1)
    lw = 1.25
    hlw = 1.95

    colors = [
        'blue',
        'cornflowerblue',
        'limegreen',
        'darkgreen',
        'xkcd:red',
        'salmon',
    ]

    temp_ranges = [
        (400, 2400),
        (800, 2250),
        (2400, 4300),
    ]
    leg_loc = [
        'upper right',
        'lower left',
        'upper left',
    ]

    depth_max = 3500, 4300, 8100
    depth_max = 3300, 4200, 7950
    molecs = {
        'H2O': 'xkcd:blue',
        'CO': 'xkcd:green',
        'CO2': 'red',
        'CH4': 'gold',
        'NH3': 'chocolate',
        'Fe': 'darkorange',
        'Na': '0.25',
        'K': '0.65',
        'TiO': 'cornflowerblue',
        'VO': 'green',
        'H2': 'xkcd:violet',
        'He': 'yellowgreen',
        'H': 'xkcd:pink',
        'e': 'deepskyblue',
    }
    vmr_range = (1e-12, 2)
    p_range = 1e2, 1e-8

    dx0 = 0.17
    dx1 = 0.3
    dx2 = 0.308
    x0 = 0.062
    deltax = 0.08
    x1 = x0+dx0+deltax
    x2 = x1+dx1+0.075

    dy = 0.22
    y0 = 0.71
    y1 = 0.815
    deltay = 0.32
    dy2 = 0.125

    fig = plt.figure(2)
    plt.clf()
    fig.set_size_inches(9.0, 10.0)
    plt.subplots_adjust(0.07, 0.05, 0.99, 0.97, wspace=0.2, hspace=0.44)
    # Temperatures
    for i in range(nplanets):
        j = 2*i
        ax = plt.axes([x0, y0-i*deltay, dx0, dy])
        ax.plot(temps[j], press, lw=lw, c=colors[j], label=labels[j])
        ax.plot(temps[j+1], press, lw=lw, c=colors[j+1], label=labels[j+1])
        ax.plot(h_temps[j], h_press, lw=hlw, c=colors[j], dashes=dashes)
        ax.plot(h_temps[j+1], h_press, lw=hlw, c=colors[j+1], dashes=dashes)
        ax.set_yscale('log')
        ax.tick_params(which='both', direction='in', top=True, labelsize=fs-1)
        ax.set_ylim(p_range)
        ax.set_xticks(np.arange(500, 4501, 500))
        ax.set_xlim(temp_ranges[i])
        ax.set_title(planets[i], fontsize=fs+1, weight='bold', pad=8)
        ax.legend(loc=leg_loc[i], fontsize=fs-1, handlelength=1.5)
        ax.set_xlabel('temperature (K)', fontsize=fs)
        ax.set_ylabel('pressure (bar)', fontsize=fs)
        if i == 0:
            legend = ax.get_legend()
            code_handles = [
                Line2D([],[], lw=1.5, color='royalblue',label='Pyrat Bay'),
                Line2D([],[], lw=1.5, color='k', dashes=(6,2), label='HELIOS'),
            ]
            ax.legend(handles=code_handles, loc=(1.5, 1.07))
            ax.add_artist(legend)
        # CF
        #xran = ax.get_xlim()
        #height = 0.075 * (xran[1]-xran[0])
        #ax.axvline(xran[0]+height, lw=0.75, color='0.5', zorder=-10, dashes=(16,2))
        #for k in range(2):
        #    cf = median_cf[j+k] / np.amax(median_cf[j+k]) * height
        #    cf = xran[0] + cf
        #    ax.plot(cf, press, color=colors[j+k], lw=1.75, zorder=10)
        #ax.set_xlim(temp_ranges[i])
    # Spectra
    for i in range(nplanets):
        for j in range(2):
            ax = plt.axes([x1, y1-i*deltay-j*(dy2+0.005), dx1, dy2])
            ax.plot(bin_wl, bin_spectra[2*i+j,1]/pc.ppm, c='0.05', lw=1.0)
            ax.plot(bin_wl, bin_spectra[2*i+j,0]/pc.ppm, c=colors[2*i+j])
            ax.set_xscale('log')
            ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
            ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
            ax.set_xticks([0.3, 1.0, 3.0, 10])
            ax.set_ylim(0.0, depth_max[i])
            ax.set_xlim(0.5, 15.0)
            ax.tick_params(which='both', direction='in', right=True, labelsize=fs-1)
            ax.set_xlabel(r'wavelength ($\mathrm{\mu}$m)', fontsize=fs)
            ax.set_ylabel(r'$F_{\rm p} / F_{\rm s}$ (ppm)', fontsize=fs)
            ax.text(
                0.03, 0.93, labels[2*i+j],
                fontsize=fs-1, transform=ax.transAxes, va='top',
                bbox=dict(boxstyle='Round', facecolor='w', alpha=0.78),
            )
    # VMRs
    for i in range(nplanets):
        for j in range(2):
            ax = plt.axes([x2, y1-i*deltay-j*(dy2+0.005), dx2, dy2])
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
            ax.set_ylim(p_range)
            ax.set_yticks(np.logspace(-8,0, 3))
            ax.tick_params(
                axis='both', which='both', direction='in', labelsize=fs-1,
                length=5.0, width=1.0,
            )
            ax.text(
                0.03, 0.93, labels[2*i+j],
                fontsize=fs-1, transform=ax.transAxes, va='top',
                bbox=dict(boxstyle='Round', facecolor='w', alpha=0.78),
            )
            ax.set_ylabel('pressure (bar)', fontsize=fs)
            if i == 0 and j == 0:
                ax.legend(
                    loc=(-0.65, 1.10), ncols=7, fontsize=fs-1,
                    labelspacing=0.1,
                    columnspacing=1.0,
                    handlelength=1.3,
                )
            if j == 0:
                ax.set_xticklabels([])
        ax.set_xlabel('volume mixing ratio', fontsize=fs)
    plt.savefig('../plots/benchmark_radeq_atmospheres.png', dpi=300)


if __name__ == '__main__':
    main()

