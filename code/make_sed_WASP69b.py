import matplotlib.pyplot as plt
import pyratbay.spectrum as ps
import pyratbay.constants as pc
import pyratbay.io as io
from astropy import units as u
from stsynphot.catalog import grid_to_spec
import gen_tso.pandeia_io as pandeia


def make_sed():
    """
    Here are three other ways to generate SEDs
    """
    # OPTION 1: Synphot models
    teff = 4700.0
    log_gstar = 4.5
    Fe_H_star = 0.15
    filename = f'inputs/synphot_phoenix_WASP69b_{teff:.0f}K_m00_logg4.5.dat'
    sp = grid_to_spec('phoenix', teff, Fe_H_star, log_gstar)
    sp_wl = sp.waveset.to('micron').value
    sp_flux = sp(sp.waveset, flux_unit=u.mJy).value * 1e-26 * pc.c
    resolution = 750.0
    bin_wl = ps.constant_resolution_spectrum(0.1, 50.0, resolution)
    bin_flux = ps.bin_spectrum(bin_wl, sp_wl, sp_flux)
    io.write_spectrum(bin_wl, bin_flux, filename, 'emission')

    # OPTION 2: Use Pandeia package to get a PHOENIX SED model
    # Closest SED to WASP-69 is an K2V model (Teff=4750K, logg=4.5)
    scene = pandeia.make_scene(
        sed_type='phoenix',
        sed_model='k2v',
    )
    sed_wl, flux = pandeia.extract_sed(scene, wl_range=(0.35,12.0))
    # Convert flux from mJy to erg s-1 cm-2 cm-1
    sed_flux = flux * pc.c / 1e26
    bin_sed_flux = ps.bin_spectrum(bin_wl, sed_wl, sed_flux, gaps='interpolate')
    starspec_file = 'inputs/phoenix_WASP69_K2V_4750K.dat'
    io.write_spectrum(
        bin_wl,
        bin_sed_flux,
        starspec_file,
        type='emission',
    )

    # OPTION 3: Fetch a PHOENIX New ERA model
    teff = 4700
    logg = 4.5
    metal = 0.15
    folder = 'inputs/'
    ps.fetch_phoenix(teff, logg, metal, folder)


def compare_seds():
    """
    These are all PHOENIX models
    """
    # Synphot model
    filename = 'inputs/synphot_phoenix_WASP69b_4700K_m00_logg4.5.dat'
    wl, flux = io.read_spectrum(filename, False)

    # New Era model
    p_file = 'inputs/lte04700-4.50-0.0.PHOENIX-NewEra-ACES-COND-2023.HSR.h5'
    p_wl, p_sed = ps.read_phoenix(p_file)
    bin_new_era = ps.bin_spectrum(wl, p_wl, p_sed)

    # Through Pandeia
    starspec_file = 'inputs/phoenix_WASP69_K2V_4750K.dat'
    gt_wl, gt_sed = io.read_spectrum(starspec_file, False)

    # Take a look
    plt.figure(0, (7., 4.))
    plt.clf()
    plt.subplots_adjust(0.07, 0.11, 0.98, 0.95)
    ax = plt.subplot(111)
    ax.plot(gt_wl, gt_sed, color='tomato', alpha=0.85, label='pandeia')
    ax.plot(wl, flux, color='royalblue', alpha=0.85, label='synphot')
    ax.plot(wl, bin_new_era, color='xkcd:green', alpha=0.85, label='new era')
    ax.set_title('WASP-69 SED Spectrum')
    ax.set_xscale('log')
    ax.set_xlim(0.35, 10.5)
    ax.set_xlabel(r'Wavelength ($\mathrm{\mu}$m)', fontsize=12)
    ax.set_ylabel(r'$F_{\rm p}$ (erg s$^{-1}$ cm$^{-2}$ cm)', fontsize=12)
    ax.tick_params(direction='in', which='both', labelsize=11)
    ax.legend(loc='best')


