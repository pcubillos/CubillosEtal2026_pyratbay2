import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pyratbay.spectrum as ps
import pyratbay.constants as pc
import pyratbay.io as io


def main():
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
    for i in range(nmodels):
        planet = planets[i].replace('-', '').replace(' ', '')
        filename = f'phoenix_{planet}_{teff[i]:.0f}K_m00_logg4.5.dat'
        wl, flux = io.read_spectrum(f'inputs/helios_{filename}', wn=False)
        h_sed[i] = ps.bin_spectrum(bin_wl, wl, flux) * beta_irr*rs_smaxis[i]**2


    # The plot
    fs = 11.5
    colors = ['salmon', 'xkcd:green', 'royalblue']
    fig = plt.figure(4)
    plt.clf()
    fig.set_size_inches(5.0, 3.15)
    plt.subplots_adjust(0.16, 0.135, 0.995, 0.995)
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


if __name__ == '__main__':
    main()

