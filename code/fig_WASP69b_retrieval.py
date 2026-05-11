import pickle
import numpy as np
import pyratbay.constants as pc
import pyratbay.spectrum as ps
import pyratbay.tools as pt
import matplotlib
import matplotlib.pyplot as plt
import mc3

matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = (
    ['Arial'] + matplotlib.rcParams['font.sans-serif']
)


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


def get_axes(posterior, ranges):
    pnames = posterior['params_names']
    indices = np.zeros(len(pnames), int)
    pranges = []
    for i,pname in enumerate(pnames):
        pname = pname[4:]
        if pname in ranges:
            indices[i] = list(ranges).index(pname)
            pranges.append(ranges[pname])
        else:
            indices[i] = -1
            pranges.append(None)
    return indices, pranges


def main():
    # Input models
    with np.load('WASP69b_transmission_spectra.npz') as d:
        models = d['models']
        true_temp = d['temp']
        species = d['species']
        quenched_vmr = d['quenched_vmr'][0]

    bin_wl = ps.constant_resolution_spectrum(0.8, 12.0, 150.0)

    # Retrieved models
    root_files = [
        'iso_059_nm/WASP69b_transit_jwst_0.59_iso_nm',
        'iso_059_snm/WASP69b_transit_jwst_0.59_iso_snm',
        'slant_059_nm/WASP69b_transit_jwst_0.59_slant_nm',
        'slant_059_snm/WASP69b_transit_jwst_0.59_slant_snm',
    ]
    nruns = len(root_files)

    blue = mc3.plots.alphatize('royalblue', 0.5, 'xkcd:blue')
    themes = [
        mc3.plots.Theme('xkcd:goldenrod', alpha_dark=0.75),
        mc3.plots.Theme('tomato'),
        mc3.plots.Theme('xkcd:green', alpha_dark=0.75),
        mc3.plots.Theme(blue),
    ]

    posteriors = []
    for i in range(nruns):
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
        post_info['params_names'] = post_info['params_names'][indices]

        post_info['bin_depths'] = [
            ps.bin_spectrum(bin_wl, wl, depth, gaps='ignore')/pc.percent
            for depth in post_info['depth_posterior']
        ]

        cf = post_info['cf_posterior_median']
        median_cf = np.median(cf, axis=1)
        idx_bounds = np.where((np.abs(median_cf-0.5)<=0.475))[0]
        p_bounds = pressure[idx_bounds[0]], pressure[idx_bounds[-1]]
        post_info['p_bounds'] = p_bounds

    post_info = posteriors[1]
    obs_wl = post_info['band_wl']
    obs_depths = post_info['data'] / pc.percent
    obs_errors = post_info['uncert'] / pc.percent

    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # The plot
    savefile = '../plots/fig_WASP69b_simulated_retrieval.png'

    titles = [
       'constant CH$_4$',
       'non-isobaric CH$_4$',
    ]
    labels = [
       'fit: G395H + LRS',
       'fit: SOSS + G395H + LRS',
       'fit: G395H + LRS',
       'fit: SOSS + G395H + LRS',
    ]

    bbox = dict(
        boxstyle="round",
        ec='none',
        fc=(1.0, 1.0, 1.0),
        alpha=0.6,
    )

    lw = 1.25
    x0 = 0.062
    x1 = 0.72
    y0 = 0.38
    dy = 0.105
    delta_y = 0.14
    y2 = 0.2
    xmargin = 0.0075

    fs = 12.0
    tick_params = {
        'which': 'both',
        'right': True,
        'direction': 'in',
        'labelsize': fs-1,
    }

    ranges = {
        'H2O':  (-3.6, -1.75),
        'CO':   (-3.9, -1.3),
        'CO2':  (-5.7, -3.4),
        'CH4':  (-5.75, -3.8),
        'K':    (-8.0, -2.),
        'NH3':  (-6.3, -4.20),
        'SO2':  (-6.1, -4.7),
        'H2S':  (-6.5, -2.51),
        'HCN':  (-10.0, -4.1),
        'C2H2': (-10,  -4.25),
    }

    leg_args = dict(
        fontsize=fs-1, framealpha=0.75, labelspacing=0.25,
        borderpad=0.4, handletextpad=0.5, handlelength=1.25,
    )


    legs = []
    fig = plt.figure(10)
    fig.set_size_inches(9.5, 6.0)
    plt.clf()
    # spectra
    ax = plt.axes([x0, y0, x1-x0, 0.99-y0])
    bx = plt.axes([0.12, 0.77, 0.22, 0.21])
    cx = plt.axes([0.44, 0.77, 0.22, 0.21])
    for i in range(nruns):
        post_info = posteriors[i]
        bin_depths = post_info['bin_depths']
        legs += ax.plot(
            bin_wl, bin_depths[0], lw=1.0,
            c=themes[i].color, alpha=1.0, label=labels[i],
        )
        bx.plot(bin_wl, bin_depths[0], c=themes[i].color)
        cx.plot(bin_wl, bin_depths[0], c=themes[i].color)
    ax.errorbar(
        obs_wl, obs_depths, obs_errors, fmt='o',
        color='k', mfc='w', ms=3.0, lw=lw, zorder=-10, mew=1.0, elinewidth=1.0,
    )
    ax.set_xscale('log')
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.set_xticks([1.0, 1.4, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0])
    ax.tick_params(**tick_params)
    ax.set_xlim(0.82, 12.0)
    ax.set_ylim(1.56, 1.89)
    ax.set_ylabel('transit depth (%)', fontsize=fs)
    ax.set_xlabel(r'wavelength ($\mathrm{\mu}$m)', fontsize=fs)
    leg1 = ax.legend(handles=legs[0:2], title=titles[0], loc=(0.28, 0.015), **leg_args)
    ax.legend(handles=legs[2:4], title=titles[1], loc=(0.64, 0.015), **leg_args)
    ax.add_artist(leg1)
    # insets
    ax.plot(
        [3.15, 3.64, 3.64, 3.15, 3.15],
        [1.67, 1.67, 1.741, 1.741, 1.67],
        c='0.5', lw=1, dashes=(5,1), zorder=-10,
    )
    ax.plot(
        [6.8, 8.6, 8.6, 6.8, 6.8],
        [1.685, 1.685, 1.754, 1.754, 1.685],
        c='0.5', lw=1, dashes=(5,1), zorder=-10,
    )
    bx.errorbar(
        obs_wl, obs_depths, obs_errors, fmt='o',
        color='k', mfc='w', ms=3.5, lw=lw, zorder=10, mew=1.0, elinewidth=1.0,
    )
    bx.set_xticks(np.arange(2.8, 3.7, 0.2))
    bx.set_yticks(np.arange(1.68, 1.80, 0.03))
    bx.set_xlim(3.15, 3.64)
    bx.set_ylim(1.67, 1.744)
    bx.tick_params(which='both', direction='in', labelsize=fs-2)
    bx.plot([3.34, 3.39, 3.45], [1.698, 1.714, 1.714], lw=0.75, c='0.5')
    bx.text(3.39, 1.717, r'$p \approx 2$ mbar', fontsize=fs-1)

    cx.errorbar(
        obs_wl, obs_depths, obs_errors, fmt='o',
        color='k', mfc='w', ms=3.5, lw=lw, zorder=10, mew=1.0, elinewidth=1.0,
    )
    cx.set_xticks(np.arange(6.0, 9.0, 0.5))
    cx.set_yticks(np.arange(1.69, 1.80, 0.03))
    cx.set_xlim(6.80, 8.6)
    cx.set_ylim(1.685, 1.754)
    cx.tick_params(which='both', direction='in', labelsize=fs-2)
    cx.plot([7.6, 7.8, 8.02], [1.72, 1.737, 1.737], lw=0.75, c='0.5')
    cx.text(7.8, 1.741, r'$p \approx 50$ $\mathrm{\mu}$bar', fontsize=fs-1)

    # Temperature
    ax = plt.axes([0.78, y0+0.05, 0.215, 0.5])
    ax.set_yscale('log')
    ax.set_ylim(1e1, 1e-7)
    ax.set_xlim(450, 1600)
    for i in range(nruns):
        tp = posteriors[i]['temperature_posterior']
        p_bounds = posteriors[i]['p_bounds']
        ax.plot(
            tp[0], pressure, color=themes[i].color, lw=1.25,
        )
        ax.fill_betweenx(
            pressure, tp[1], tp[2], facecolor=themes[i].color,
            edgecolor='none', alpha=0.4,
        )
    ax.plot(
        true_temp[0], pressure, c='k',
        lw=1.5, dashes=(6,1.0), zorder=10,
    )
    ax.axhspan(
        p_bounds[0], p_bounds[-1], edgecolor='0.5',
        alpha=0.5, zorder=-10, hatch='//', facecolor='none',
    )
    ax.text(
        0.5, 0.5, 'probed\nregion', color='0.5',
        transform=ax.transAxes, fontsize=fs, bbox=bbox,
    )
    ax.set_ylabel('pressure (bar)', fontsize=fs, labelpad=0)
    ax.set_xlabel('temperature (K)', fontsize=fs)
    ax.set_yticks(np.logspace(-7, 1, 5))
    ax.tick_params(**tick_params)

    # Histograms
    for k in [0, 1]:
        y = y2 - delta_y*k
        h_rect = [x0, y, 0.99, y+dy]
        positions = 1 + np.arange(11)
        h_axes = np.array([
            mc3.plots.subplot(h_rect, xmargin, pos, nx=10, ny=1)
            for pos in positions
        ])
        h_axes[-1].set_position([1.5, 0.5, 0.1, 0.1])

        posterior = posteriors[2*k]
        indices, pranges = get_axes(posterior, ranges)
        axes = h_axes[indices]
        posterior['posterior'].ranges = pranges
        hfig = posterior['posterior'].plot_histogram(
            axes=axes, show_estimates=False, show_texts=False,
        )
        hfig.fontsize = fs - 1

        posterior = posteriors[2*k+1]
        indices, pranges = get_axes(posterior, ranges)
        axes = [h_axes[j].twinx() for j in indices]
        posterior['posterior'].ranges = pranges
        hfig = posterior['posterior'].plot_histogram(
            axes=axes, show_estimates=False, show_texts=False,
        )
        hfig.fontsize = fs - 1
        for ax in axes:
            ax.set_yticks([])

        p_bounds = posteriors[2]['p_bounds']
        for j,ax in enumerate(h_axes[:-1]):
            mol = list(ranges)[j]
            ax.set_xlim(ranges[mol])
            ax.set_yticks([])
            ax.tick_params(direction='in', which='both', labelsize=fs-1)
            xlabel = '' if k==0 else f'$\\log\\ X_{{\\rm {mol}}}$'
            ax.set_xlabel(xlabel, fontsize=fs-1, labelpad=0)
            bx = ax.twinx()
            index = list(species).index(mol)
            if k==1 and mol=='CH4':
                for i in [2, 3]:
                    imol = list(posteriors[i]['species']).index(mol)
                    post_vmr = np.log10(posteriors[i]['vmr_posterior'][:,:,imol])
                    bx.fill_betweenx(
                        pressure, post_vmr[3], post_vmr[4],
                        color=themes[i].light_color, alpha=0.12,
                    )
                    bx.fill_betweenx(
                        pressure, post_vmr[1], post_vmr[2],
                        color=themes[i].light_color, alpha=0.5,
                    )
            bx.yaxis.set_label_position("left")
            bx.yaxis.tick_left()
            true_vmr = np.log10(quenched_vmr[:,index])
            bx.plot(true_vmr, pressure, color='k', lw=1.25, dashes=(9,1))
            bx.set_yscale('log')
            bx.set_ylim(10, 1e-7)
            bx.tick_params(direction='in', which='both', labelsize=fs-1)
            bx.set_yticks(np.logspace(-6, 0, 3))
            if j != 0:
                bx.set_yticklabels([])
            elif k==0:
                bx.set_ylabel('pressure (bar)      ', fontsize=fs-1, loc='top')
            bx.axhspan(
                p_bounds[0], p_bounds[-1], edgecolor='0.5',
                alpha=0.35, zorder=-10, hatch='//', facecolor='none',
            )
    plt.savefig(savefile, dpi=300)

