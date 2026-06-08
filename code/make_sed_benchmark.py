import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pyratbay.spectrum as ps
import pyratbay.constants as pc
import pyratbay.io as io
import os
from astropy import units as u
from stsynphot.catalog import grid_to_spec
import matplotlib.font_manager as font_manager


font_dir = ['/home/pcubillos/tmp/fonts']
for font in font_manager.findSystemFonts(font_dir):
    font_manager.fontManager.addfont(font)
# Set font family globally
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial'] + matplotlib.rcParams['font.sans-serif']


def remake_from_helios():
    """
    synphot is equivalent to pysynphot (which I originally used 
    to create my PHOENIX SEDs)
    """
    resolution = 750.0
    bin_wl = ps.constant_resolution_spectrum(0.1, 50.0, resolution)

    # System parameters
    planets = [
        'WASP107b',
        'WASP39b',
        'WASP121b',
    ]
    Fe_H_star = 0.0
    log_gstar = 4.5
    beta_irr = 0.667
    teff = 4400.0, 5500.0, 6400.0
    rstar = np.array([0.670, 0.940, 1.458]) * pc.rsun
    smaxis = np.array([0.055, 0.048, 0.026]) * pc.au
    rs_smaxis = rstar / smaxis

    # HELIOS models
    h_files = [
        'inputs/helios_v03/wasp-107b-f0.667-solar-no_tio_vo-tint_0',
        'inputs/helios_v03/wasp-39b-f0.667-solar-no_tio_vo-tint_0',
        'inputs/helios_v03/wasp-121b-f0.667-solar-no_tio_vo-tint_0',
    ]

    nplanets = len(planets)
    for i in range(nplanets):
        h_folder = h_files[i]
        path, root = os.path.split(h_folder)
        flux_file = f'{h_folder}/{root}_TOA_flux_eclipse.dat'
        data = np.loadtxt(flux_file, unpack=True, skiprows=4)
        h_wl = data[1]
        # Spectral flux
        # Convert F_lamda to F_nu (erg s-1 cm-2 cm-1  --> erg s-1 cm-2 cm)
        # Convert TOA to surface flux
        h_star = data[4] * (h_wl * pc.um)**2
        h_star = h_star / beta_irr / rs_smaxis[i]**2
        h_flux = ps.bin_spectrum(bin_wl, h_wl, h_star)

        planet = planets[i]
        filename = f'inputs/helios_phoenix_{planet}_{teff[i]:.0f}K_m00_logg4.5.dat'
        io.write_spectrum(bin_wl, h_flux, filename, 'emission')

    # Synphot models
    nplanets = len(planets)
    for i in range(nplanets):
        planet = planets[i]
        filename = f'inputs/synphot_phoenix_{planet}_{teff[i]:.0f}K_m00_logg4.5.dat'
        sp = grid_to_spec('phoenix', teff[i], Fe_H_star, log_gstar)
        sp_wl = sp.waveset.to('micron').value
        sp_flux = sp(sp.waveset, flux_unit=u.mJy).value * 1e-26 * pc.c
        bin_flux = ps.bin_spectrum(bin_wl, sp_wl, sp_flux)
        io.write_spectrum(bin_wl, bin_flux, filename, 'emission')



def plot_sed():
    # System parameters
    planets = [
        'WASP-121 b',
        'WASP-39 b',
        'WASP-107 b',
    ]
    nmodels = len(planets)

    resolution = 150.0
    bin_wl = ps.constant_resolution_spectrum(0.1, 50.0, resolution)
    nbin = len(bin_wl)

    beta_irr = 0.667
    teff = 6400.0, 5500.0, 4400.0
    rstar = np.array([1.458, 0.940, 0.670]) * pc.rsun
    smaxis = np.array([0.026, 0.048, 0.055]) * pc.au
    rs_smaxis = rstar / smaxis

    h_sed = np.zeros((nmodels,nbin))
    s_sed = np.zeros((nmodels,nbin))
    for i in range(nmodels):
        planet = planets[i].replace('-', '').replace(' ', '')
        filename = f'phoenix_{planet}_{teff[i]:.0f}K_m00_logg4.5.dat'
        wl, flux = io.read_spectrum(f'inputs/helios_{filename}', wn=False)
        h_sed[i] = ps.bin_spectrum(bin_wl, wl, flux) * beta_irr*rs_smaxis[i]**2


    # The plot
    fs = 12
    colors = ['salmon', 'xkcd:green', 'royalblue']
    fig = plt.figure(4)
    plt.clf()
    fig.set_size_inches(5.0, 3.15)
    plt.subplots_adjust(0.15, 0.135, 0.995, 0.995)
    ax = plt.subplot(111)
    for i in range(nmodels):
        planet = planets[i].replace('-', '')
        ax.plot(bin_wl, h_sed[i], c=colors[i], label=planets[i], alpha=0.95)
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
    ax.set_xlabel(r'wavelength ($\mathrm{\mu}$m)', fontsize=fs)
    ax.legend(loc='lower right', fontsize=fs)
    ax.tick_params(which='both', direction='in', labelsize=fs-1)
    plt.savefig('plots/benchmark_radeq_sed.png', dpi=300)



