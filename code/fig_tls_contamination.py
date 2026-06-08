import pyratbay.spectrum as ps
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager


font_dir = ['/home/pcubillos/tmp/fonts']
for font in font_manager.findSystemFonts(font_dir):
    font_manager.fontManager.addfont(font)
# Set font family globally
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial'] + matplotlib.rcParams['font.sans-serif']


def setup():
    """
    Fetch all PHOENIX New-Era models for logg=4.5 and [M/H] = 0.5
    """
    # Leave teff as None to select all Teff models
    # Find models with closest log_g and metallicity
    teff = None
    logg = 4.57
    metal = 0.5
    # Make sure that the output folder already exists
    folder = 'phoenix/'
    ps.fetch_phoenix(teff, logg, metal, folder)


def main():
    # Initialize TLS model
    sed_folder = 'phoenix/'
    teff = 4800.0
    wl = ps.constant_resolution_spectrum(0.3, 12.0, resolution=400.0)
    tls = ps.TransitLightSource(sed_folder, teff, wl, sampling='bin')

    # Evaluate TLS effect for a range of star spot/faculae temperatures
    nspots = 9
    f_spot = 0.01
    t_spots = teff + np.linspace(-1600, 1600, nspots)
    epsilon = [tls.epsilon(t_spot, f_spot) for t_spot in t_spots]


    fs = 12
    fig = plt.figure(1)
    fig.set_size_inches(5.5, 3.25)
    plt.clf()
    plt.subplots_adjust(0.12, 0.14, 0.995, 0.995)
    ax = plt.subplot(111)
    for i,t in enumerate(t_spots):
        col = 'red' if t==teff else plt.cm.viridis(i/(nspots-1))
        label = f'$T$spot = {t_spots[i]:.0f} K'
        plt.plot(tls.wl, epsilon[i], color=col, label=label, alpha=0.75)
    plt.legend(loc='lower right', fontsize=fs-2, ncols=2, handlelength=1.5)
    ax.set_xlim(0.45, 12)
    ax.set_ylim(0.9825, 1.0105)
    ax.set_xlabel(r"wavelength ($\mathrm{\mu}$m)", fontsize=fs)
    ax.set_ylabel(r"TLS contamination $\epsilon$", fontsize=fs)
    ax.set_xscale('log')
    ax.tick_params(which='both', direction='in', labelsize=fs-1)
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.set_xticks([0.5, 0.7, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0])
    plt.savefig('../plots/fig_tls_4800K_star_spots.png', dpi=300)


if __name__ == '__main__':
    main()

