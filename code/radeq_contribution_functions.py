import numpy as np
import pyratbay as pb
import pyratbay.io as io
import pyratbay.spectrum as ps


def make_obs_file():
    """
    Generate low-res NIR obserbing bands for radeq contribution functions.
    """
    obs_file = 'obs_generic_nir.dat'
    bin_wl = ps.constant_resolution_spectrum(1.0, 10.0, 50.0)
    half_widths = 0.5*np.ediff1d(bin_wl)
    wl = 0.5*(bin_wl[1:] + bin_wl[:-1])
    names = ['nir' for _ in wl]
    io.write_observations(
        obs_file, names, wl, half_widths,
    )


def contributions():
    """
    Generate radeq contribution functions.
    """
    cfg_files = [
        'cf_benchmark_WASP107b_control.cfg',
        'cf_benchmark_WASP107b_tint_350K.cfg',
        'cf_benchmark_WASP39b_control.cfg',
        'cf_benchmark_WASP39b_50x.cfg',
        'cf_benchmark_WASP121b_control.cfg',
        'cf_benchmark_WASP121b_TiO_VO.cfg',
    ]

    cf_data = []
    temps = []
    for cfg_file in cfg_files:
        pyrat = pb.Pyrat(cfg_file)
        pyrat.eval([])
        band_cf = pyrat.band_contribution()
        cf_data.append(band_cf)
        temps.append(pyrat.atm.temp)

    np.savez(
        'contribution_functions_benchmark.npz',
        cf=np.array(cf_data),
        press=pyrat.atm.press,
        temps=np.array(temps),
        wl=pyrat.obs.band_wl,
        models=cfg_files,
    )


if __name__ == '__main__':
    contributions()
