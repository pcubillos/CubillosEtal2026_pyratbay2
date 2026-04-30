import pickle

import numpy as np
import pyratbay.constants as pc
import pyratbay.io as io
import pyratbay.spectrum as ps
import pyratbay.tools as pt
import matplotlib
import matplotlib.pyplot as plt
import mc3
import gen_tso.pandeia_io as jwst


#font_dir = ['/home/pcubillos/tmp/fonts']
#for font in font_manager.findSystemFonts(font_dir):
#    font_manager.fontManager.addfont(font)
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial'] + matplotlib.rcParams['font.sans-serif']


def find_pos(ret_names, pnames, nout):
    indices = np.tile(nout, len(ret_names))
    for i,pname in enumerate(ret_names):
        if pname in pnames:
            indices[i] = list(pnames).index(pname)
    return indices


def main():
    # Input models
    with np.load('WASP69b_transmission_spectra.npz') as d:
        species = d['species']
        quenched_vmr = d['quenched_vmr']

    bin_wl = ps.constant_resolution_spectrum(0.8, 12.0, 150.0)

    # Retrieved models
    root_files = [
        'iso_059_nm/WASP69b_transit_jwst_0.59_iso_nm',
        'iso_059_snm/WASP69b_transit_jwst_0.59_iso_snm',
        'slant_059_snm/WASP69b_transit_jwst_0.59_slant_snm',
    ]
    nfiles = len(root_files)

    themes = [
        mc3.plots.Theme('xkcd:blue', alpha_dark=0.75),
        mc3.plots.Theme('tomato', alpha_dark=0.75),
        mc3.plots.Theme('xkcd:green', alpha_dark=0.75),
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
        posterior = mc3.plots.Posterior(
            post[:,indices], pnames=tex_names[indices], theme=themes[i],
        )
        post_info['posterior'] = posterior
        print(post[:,indices].shape)

        post_info['bin_depths'] = [
            ps.bin_spectrum(bin_wl, wl, depth, gaps='ignore')/pc.percent
            for depth in post_info['depth_posterior']
        ]

        cf = post_info['cf_posterior_median']
        median_cf = np.median(cf, axis=1)
        idx_bounds = np.where((np.abs(median_cf-0.5)<=0.475))[0]
        p_bounds = pressure[idx_bounds[0]], pressure[idx_bounds[-1]]
        post_info['p_bounds'] = p_bounds

        model_idx = 0 if '0.59' in root_files[i] else 1
        log_vmr = np.log10(quenched_vmr[model_idx])
        ret_species = 'H2O CO CO2 CH4 K NH3 SO2 HCN H2S C2H2'.split()
        median_vmr = [
            np.median(log_vmr[idx_bounds, list(species).index(spec)])
            for spec in ret_species
        ]
        median_vmr = np.clip(median_vmr, -10.0, 1)
        post_info['median_vmr'] = median_vmr

    post_info = posteriors[0]
    obs_wl = post_info['band_wl']
    obs_depths = post_info['data'] / pc.percent
    obs_errors = post_info['uncert'] / pc.percent

    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # The plots

    fs = 9.5
    lw = 1.25
    ms = 2.5
    tick_params = {
        'which': 'both',
        'right': True,
        'direction': 'in',
        'labelsize': fs-1,
    }

    savefile = '../plots/WASP69b_simulated_retrieval_methane_v07.png'
    # Happy with this
    savefile = '../plots/WASP69b_simulated_retrieval_methane.png'
    labels = [
       'isobaric CH$_4$: G395H + LRS',
       'isobaric CH$_4$: SOSS + G395H + LRS',
       'slanted CH$_4$: SOSS + G395H + LRS',
    ]

    ranges = {
        'H2O':  (-3.6, -1.75),
        'CO':   (-4.0, -1.25),
        'CO2':  (-5.7, -3.4),
        'CH4':  (-5.3, -3.7),
        'K':    (-7.8, -5.8),
        'NH3':  (-6.3, -4.0),
        'SO2':  (-6.1, -4.5),
        'H2S':  (-5.4, -2.4),
        'HCN':  (-9.0, -5.5),
        'C2H2': (-10,  -4.25),
    }

    all_pars = list(set(np.concatenate(
        [post['params_names'][6:] for post in posteriors]
    )))
    npars = len(all_pars)
    pnames = [
        'log_H2O', 'log_CO', 'log_CO2',
        'log_CH4', 'log_SO2', 'log_C2H2',
    ]

    fig = plt.figure(1)
    fig.set_size_inches(4.5, 5.75)
    plt.clf()
    plt.subplots_adjust(0.11, 0.0, 0.99, 0.9, hspace=0.0)
    # main
    ax = plt.subplot(211)
    for i in [0, 1, 2]:
        post_info = posteriors[i]
        bin_depths = post_info['bin_depths']
        ax.plot(
            bin_wl, bin_depths[0],
            c=themes[i].color, alpha=0.85, label=labels[i],
        )
    ax.errorbar(
        obs_wl, obs_depths, obs_errors, fmt='o',
        color='0.15', ms=ms, lw=lw, zorder=-10, mew=0.0, elinewidth=1.0,
    )
    ax.set_xscale('log')
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.set_xticks([0.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0])
    ax.set_yticks(np.arange(1.65, 1.86, 0.03))
    ax.tick_params(**tick_params)
    ax.set_xlim(3.0, 12.0)
    ax.set_ylim(1.64, 1.85)
    ax.set_ylabel('transit depth (%)', fontsize=fs)
    ax.set_xlabel('wavelength (um)', fontsize=fs)
    ax.legend(
        loc=(0.0, 1.01), fontsize=fs-0.5, labelspacing=0.25,
        borderpad=0.4, handletextpad=0.5)
    # insets
    bx = plt.axes([0.18, 0.72, 0.35, 0.17])
    cx = plt.axes([0.63, 0.72, 0.35, 0.17])
    bx.clear()
    cx.clear()
    for i in [0, 1, 2]:
        post_info = posteriors[i]
        bin_depths = post_info['bin_depths']
        bx.plot(bin_wl, bin_depths[0], c=themes[i].color, alpha=0.85)
        cx.plot(bin_wl, bin_depths[0], c=themes[i].color, alpha=0.85)
    bx.errorbar(
        obs_wl, obs_depths, obs_errors, fmt='o',
        color='0.15', ms=ms, lw=lw, zorder=10, mew=0.0, elinewidth=1.0,
    )
    bx.set_xticks(np.arange(2.8, 3.7, 0.2))
    bx.set_yticks(np.arange(1.68, 1.80, 0.03))
    bx.set_xlim(3.15, 3.64)
    bx.set_ylim(1.67, 1.741)
    bx.tick_params(**tick_params)
    bx.plot([3.34, 3.39, 3.45], [1.693, 1.712, 1.712], lw=0.75, c='0.5')
    bx.text(3.39, 1.715, r'$p \approx 1.5$ mbar', fontsize=fs-1)
    cx.errorbar(
        obs_wl, obs_depths, obs_errors, fmt='o',
        color='0.15', ms=ms, lw=lw, zorder=10, mew=0.0, elinewidth=1.0,
    )
    cx.set_xticks(np.arange(6.0, 9.0, 0.5))
    cx.set_yticks(np.arange(1.69, 1.80, 0.03))
    cx.set_xlim(6.80, 8.6)
    cx.set_ylim(1.685, 1.752)
    cx.tick_params(**tick_params)
    cx.plot([7.6, 7.745, 8.0], [1.73, 1.74, 1.74], lw=0.75, c='0.5')
    cx.text(7.745, 1.743, r'$p \approx 20$ $\mathrm{\mu}$bar', fontsize=fs-1)

    # Histograms
    h_rect = [0.11, 0.055, 0.99, 0.37]
    xmargin = 0.02
    positions = np.arange(npars)
    h_axes = [
        mc3.plots.subplot(h_rect, xmargin, pos+1, nx=3, ny=2, ymargin=0.065)
        for pos in positions
    ]
    for i in range(len(pnames), npars):
        h_axes[i].set_position([1.5, 0.5, 0.1, 0.1])
    # Isobaric
    post = posteriors[0]['posterior']
    locs = find_pos(posteriors[0]['params_names'][6:], pnames, npars-1)
    hfig = post.plot_histogram(
        axes=np.array(h_axes)[locs], show_estimates=False, show_texts=False,
    )
    hfig.fontsize = fs
    for ax in h_axes:
        ax.set_xlabel(ax.get_xlabel(), labelpad=0)
    # Non-isobaric (TBD skip CH4, overplot as VMR(p) boundaries)
    post = posteriors[1]['posterior']
    locs = find_pos(posteriors[1]['params_names'][6:], pnames, npars-1)
    axes = [h_axes[j].twinx() for j in locs]
    hfig = post.plot_histogram(
        axes=axes, show_estimates=False, show_texts=False,
    )
    for ax in axes:
        ax.set_yticks([])
    post = posteriors[2]['posterior']
    locs = find_pos(posteriors[2]['params_names'][6:], pnames, npars-1)
    axes = [h_axes[j].twinx() for j in locs]
    hfig = post.plot_histogram(
        axes=axes, show_estimates=False, show_texts=False,
    )
    for ax in axes:
        ax.set_yticks([])

    # Overplot true VMR profile
    p_bounds = posteriors[2]['p_bounds']
    for j in range(len(pnames)):
        mol = pnames[j][4:]
        bx = h_axes[j].twinx()
        h_axes[j].set_yticks([])
        bx.yaxis.set_label_position("left")
        bx.yaxis.tick_left()
        imol = list(species).index(mol)
        true_vmr = np.log10(quenched_vmr[0][:,imol])
        bx.plot(true_vmr, pressure, color='k', lw=1.25, dashes=(9,1))
        bx.set_yscale('log')
        bx.tick_params(direction='in', which='both', labelsize=fs-1)
        if j%3 == 0:
            bx.set_ylabel('pressure (bar)', fontsize=fs)
        else:
            bx.set_yticklabels([])
        bx.axhspan(
            p_bounds[0], p_bounds[-1], edgecolor='0.5',
            alpha=0.3, zorder=-20, hatch='//', facecolor='none',
        )
        if mol == 'CH4':
            i = 2
            vmr_post = np.log10(posteriors[i]['vmr_posterior'])
            post_species = posteriors[i]['species']
            imol = list(post_species).index(mol)
            bx.fill_betweenx(
                pressure, vmr_post[3,:,imol], vmr_post[4,:,imol],
                color=themes[i].light_color, alpha=0.2,
            )
            bx.fill_betweenx(
                pressure, vmr_post[1,:,imol], vmr_post[2,:,imol],
                color=themes[i].light_color, alpha=0.6,
            )
        bx.set_yticks(np.logspace(0, -6, 3))
        bx.set_ylim(3, 1e-7)
        bx.set_xlim(ranges[mol])
    plt.savefig(savefile, dpi=300)
