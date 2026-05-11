import pickle

import numpy as np
import pyratbay.constants as pc
import pyratbay.io as io
import pyratbay.spectrum as ps
import pyratbay.tools as pt
import matplotlib
import matplotlib.pyplot as plt
import mc3

matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial'] + matplotlib.rcParams['font.sans-serif']


def plot_box(ax, rect, text, lw, fs, bcol, tcol=None):
    if tcol is None:
        tcol = bcol
    box = (
        [rect[0], rect[0], rect[1], rect[1], rect[0]],
        [rect[2], rect[3], rect[3], rect[2], rect[2]],
    )
    midx = 0.48*(rect[0]+rect[1])
    midy = 0.5*(rect[2]+rect[3])
    ax.plot(box[0], box[1], lw=lw, c=bcol)
    ax.text(
        midx, midy, text,
        fontsize=fs, color=tcol, ha='center', va='center',
    )


def main():
    # Input models
    with np.load('WASP69b_transmission_spectra.npz') as d:
        models = d['models']
        model_wl = d['wl']
        true_spectra = d['spectra']
        press = d['press']
        true_temp = d['temp']
        species = d['species']
        quenched_vmr = d['quenched_vmr']

    bin_wl = ps.constant_resolution_spectrum(0.8, 12.0, 150.0)

    # Retrieved models
    root_files = [
        'slant_059_snm/WASP69b_transit_jwst_0.59_slant_snm',
        'iso_110_snm/WASP69b_transit_jwst_1.10_iso',
    ]
    nfiles = len(root_files)

    themes = [
        mc3.plots.Theme('salmon'),
        mc3.plots.Theme('royalblue'),
    ]

    posteriors = []
    for i in range(nfiles):
        with open(f'{root_files[i]}_posteriors_info.pickle', 'rb') as f:
            post_info = pickle.load(f)
        posteriors.append(post_info)
        pressure = post_info['pressure']
        wl = post_info['wl']
        tex_names = post_info['params_texnames']
        npars = len(tex_names)

        post = pt.weighted_to_equal(f'{root_files[i]}.txt')
        indices = np.arange(6,npars)
        if i == 0:
            # Load the whole thing, but then remove the CH4 pars, same shape
            full_posterior = mc3.plots.Posterior(
                post[:,indices], pnames=tex_names[indices], theme=themes[i],
            )
            indices = np.arange(6,npars-2)
            indices[4:] += 2
        posterior = mc3.plots.Posterior(
            post[:,indices], pnames=tex_names[indices], theme=themes[i],
        )
        post_info['posterior'] = posterior
        post_info['pnames'] = post_info['params_names'][indices]

        post_info['bin_depths'] = [
            ps.bin_spectrum(bin_wl, wl, depth)/pc.percent
            for depth in post_info['depth_posterior']
        ]

        cf = post_info['cf_posterior_median']
        median_cf = np.median(cf, axis=1)
        idx_bounds = np.where((np.abs(median_cf-0.5)<=0.475))[0]
        p_bounds = pressure[idx_bounds[0]], pressure[idx_bounds[-1]]
        post_info['p_bounds'] = p_bounds


    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # The plots  C/O = 0.59  vs C/O = 1.05
    lw = 1.25
    tick_params = {
        'which': 'both',
        'right': True,
        'top': True,
        'direction': 'in',
        'labelsize': fs-1,
    }
    bbox = dict(
        boxstyle="round",
        ec='none',
        fc=(1.0, 1.0, 1.0),
        alpha=0.6,
    )
    true_theme = mc3.plots.Theme('xkcd:green')

    savefile = '../plots/WASP69b_simulated_retrieval_v04.png'
    niriss_box  = [0.83, 2.80, 1.56, 1.58]
    #nirspec_box = [2.86, 5.13, 1.58, 1.60]
    nirspec_box = [2.86, 5.13, 1.56, 1.58]
    miri_box    = [5.00, 11.9, 1.56, 1.58]

    niriss_col = 'xkcd:green'
    nirspec_col = 'xkcd:violet'
    miri_col = 'goldenrod'

    labels = [
        'C/O = 0.59',
        'C/O = 1.05',
    ]
    true_cols = [
        'red',
        'xkcd:blue',
    ]

    ranges = [
        (-5.5, -4.7),   # p_ref
        (-5.0, -2.3),   # H2O
        (-3.75, -1.75),  # CO
        (-7.15, -4.35),   # CO2
        (-5.75, -2.85),   # CH4
        (-7.5, -5.25),    # K
        (-6.25, -5.1),   # NH3
        (-6.25, -5.1),  # SO2
        (-7.5, -3.6),     # H2S
        (-10, -6.),   # HCN
        (-10, -6.),    # C2H2
        (-3.5, 2.0),    # p_cloud
    ]

    fs = 12.0
    tick_params = {
        'which': 'both',
        'right': True,
        'top': True,
        'direction': 'in',
        'labelsize': fs-1,
    }

    ret_species = 'H2O CO CO2 CH4 K NH3 SO2 H2S HCN C2H2'.split()
    ret_labels = [f'$\\log\\ X_{{\\rm {mol}}}$' for mol in ret_species]

    niriss_col = 'green'
    nirspec_col = 'limegreen'
    miri_col = 'forestgreen'

    x0 = 0.067
    y0 = 0.4
    height = 0.98 - y0

    fig = plt.figure(10)
    fig.set_size_inches(9.5, 4.25)
    plt.clf()
    # Spectrum
    ax = plt.axes([x0, y0, 0.674, height])
    for i in [0,1]:
        post_info = posteriors[i]
        bin_depth = post_info['bin_depths']  # best-fit
        ax.plot(bin_wl, bin_depth[0], c=themes[i].light_color, alpha=0.9)
        obs_wl = post_info['band_wl']
        obs_depth = post_info['data'] / pc.percent
        obs_err = post_info['uncert'] / pc.percent
        ax.errorbar(
            obs_wl, obs_depth, obs_err, fmt='o',
            c=true_cols[i], ms=4.0, lw=lw, label=labels[i], zorder=10-i,
            mew=1.0, mfc='w', elinewidth=1.0, ecolor=themes[i].color,
        )
    ax.set_xscale('log')
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.set_xticks([1.0, 1.4, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0])
    ax.tick_params(**tick_params)
    ax.set_xlim(0.82, 12.0)
    ax.set_ylim(1.55, 1.815)
    ax.set_ylabel('transit depth (%)', fontsize=fs)
    ax.set_xlabel(r'wavelength ($\mathrm{\mu}$m)', fontsize=fs)
    ax.legend(loc='upper left', fontsize=fs-0.5)
    ymin, ymax = ax.get_ylim()
    ax.set_ylim(ymin, ymax)
    plot_box(ax, niriss_box, 'NIRISS / SOSS', 1.25, fs-1.25, niriss_col)
    plot_box(ax, nirspec_box, 'NIRSpec / G395H', 1.25, fs-1.25, nirspec_col)
    plot_box(ax, miri_box, 'MIRI / LRS', 1.25, fs-1.25, miri_col)

    # Temperature
    ax = plt.axes([0.805, y0, 0.19, height])
    ax.clear()
    for i in [0,1]:
        tp = posteriors[i]['temperature_posterior']
        p_bounds = posteriors[i]['p_bounds']
        ax.axhspan(
            p_bounds[0], p_bounds[-1], edgecolor=themes[i].color,
            alpha=0.5, zorder=-10, hatch='//', facecolor='none',
        )
        ax.fill_betweenx(
            pressure, tp[1], tp[2], facecolor=themes[i].color,
            edgecolor='none', alpha=0.45,
        )
        ax.plot(
            true_temp[i], pressure, c=themes[i].color,
            lw=1.5, dashes=(6,1.0), zorder=10,
        )
    ax.text(
        0.53, 0.525, 'probed\nregion',
        transform=ax.transAxes, fontsize=fs, bbox=bbox,
    )
    ax.set_yscale('log')
    ax.set_ylabel('pressure (bar)', fontsize=fs)
    ax.set_xlabel('temperature (K)', fontsize=fs)
    ax.set_yticks(np.logspace(-7, 1, 5))
    ax.tick_params(**tick_params)
    ax.set_ylim(1e2, 1e-7)
    ax.set_ylim(1e2, 5e-8)

    # Histograms
    #plt.clf()
    h_rect = [x0, 0.11, 0.99, 0.28]
    xmargin = 0.012
    positions = np.arange(12)
    positions[0] = 11
    positions[11] = 12
    h_axes = [
        mc3.plots.subplot(h_rect, xmargin, pos, nx=12, ny=1)
        for pos in positions
    ]
    post = posteriors[1]['posterior']
    post.theme = themes[1]
    post.ranges = ranges
    hfig = post.plot_histogram(
        axes=h_axes, show_estimates=False, show_texts=False,
    )
    hfig.fontsize = fs - 1
    post = posteriors[0]['posterior']
    post.theme = themes[0]
    post.ranges = ranges
    hfig.overplot(post)

    for j,ax in enumerate(h_axes):
        lab = ax.get_xlabel()
        if lab not in ret_labels:
            continue
        bx = h_axes[j].twinx()
        idx = ret_labels.index(lab)
        index = list(species).index(ret_species[idx])
        xran = ax.get_xlim()
        if ret_species[idx] == 'CH4':
            ax.set_ylim(-2, -1)
            for i in [0,1]:
                imol = list(posteriors[i]['species']).index(ret_species[idx])
                post_vmr = np.log10(posteriors[i]['vmr_posterior'][:,:,imol])
                bx.fill_betweenx(
                    pressure, post_vmr[3], post_vmr[4],
                    color=themes[i].light_color, alpha=0.12,
                )
                bx.fill_betweenx(
                    pressure, post_vmr[1], post_vmr[2],
                    color=themes[i].light_color, alpha=0.5,
                )
        for i in [0,1]:
            bx.yaxis.set_label_position("left")
            bx.yaxis.tick_left()
            true_vmr = np.log10(quenched_vmr[i,:,index])
            col = true_cols[i]
            bx.plot(true_vmr, pressure, color=col, lw=1.25, dashes=(9,1))
            ax.set_xlim(xran)
            bx.set_yscale('log')
            bx.set_ylim(10, 1e-7)
            bx.tick_params(direction='in', which='both', labelsize=fs-1)
            ax.set_yticks([])
            bx.set_yticks(np.logspace(-6, 0, 3))
            if j != 1:
                bx.set_yticklabels([])
            else:
                bx.set_ylabel('pressure (bar)', fontsize=fs-1, loc='top')
            bx.axhspan(
                p_bounds[0], p_bounds[-1], edgecolor='0.8',
                alpha=0.5, zorder=-10, hatch='//', facecolor='none',
            )
    plt.savefig(savefile, dpi=300)

