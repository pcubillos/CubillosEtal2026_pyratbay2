# pip install stsynphot
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pyratbay.spectrum as ps
import pyratbay.constants as pc
import pyratbay.io as io
import os
from astropy import units as u
from stsynphot.catalog import grid_to_spec


def remake_from_helios():
    """
    synphot is equivalent to pysynphot (which I originally used 
    to create my PHOENIX SEDs)
    """
    resolution = 750.0
    bin_wl = ps.constant_resolution_spectrum(0.1, 50.0, resolution)

    # Synphot models
    teff = 4700.0
    log_gstar = 4.5
    Fe_H_star = 0.15
    filename = f'inputs/synphot_phoenix_WASP69b_{teff:.0f}K_m00_logg4.5.dat'
    sp = grid_to_spec('phoenix', teff, Fe_H_star, log_gstar)
    sp_wl = sp.waveset.to('micron').value
    sp_flux = sp(sp.waveset, flux_unit=u.mJy).value * 1e-26 * pc.c
    bin_flux = ps.bin_spectrum(bin_wl, sp_wl, sp_flux)
    io.write_spectrum(bin_wl, bin_flux, filename, 'emission')



def plot_SEDs():
    # System parameters
    teff = 4700.0
    log_gstar = 4.5
    resolution = 150.0
    bin_wl = ps.constant_resolution_spectrum(0.1, 50.0, resolution)
    nbin = len(bin_wl)

    filename = f'inputs/synphot_phoenix_WASP69b_{teff:.0f}K_m00_logg4.5.dat'
    wl, flux = io.read_spectrum(filename, wn=False)
    s_sed = ps.bin_spectrum(bin_wl, wl, flux)

    flux, wn, t, g = ps.read_kurucz('inputs/fp02k2odfnew.pck', teff, log_gstar)
    wl = 1.0 / wn / pc.um
    k_sed = ps.bin_spectrum(bin_wl, wl, flux, gaps='ignore')

    # The plot
    fs = 11
    colors = ['salmon', 'xkcd:green', 'royalblue']
    fig = plt.figure(4)
    plt.clf()
    fig.set_size_inches(5.0, 3.15)
    plt.subplots_adjust(0.11, 0.13, 0.995, 0.995)
    ax = plt.subplot(111)
    ax.plot(bin_wl, k_sed, c='royalblue', label='Kurucz')
    ax.plot(bin_wl, s_sed, c='salmon', label='synphot', alpha=0.95)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.set_xticks([0.3, 1.0, 3.0, 10, 30])
    ax.set_xlim(0.15, 30.0)
    ax.set_ylim(1e0, 3e6)
    ax.set_ylabel(
        'WASP-69 SED (erg s$^{-1}$ cm$^{-2}$ cm)',
        fontsize=fs, labelpad=0,
    )
    ax.set_xlabel(r'wavelength ($\mathrm{\mu}$m)', fontsize=fs)
    ax.legend(loc='lower right', fontsize=fs)
    ax.tick_params(which='both', direction='in', labelsize=fs-1)
    plt.savefig('plots/WASP69_sed.png', dpi=300)
