import pickle
import numpy as np
import pyratbay as pb
import pyratbay.io as io
import pyratbay.constants as pc
import pyratbay.spectrum as ps

import gen_tso.pandeia_io as jwst
import gen_tso.catalogs as cat


def main():
    """
    Take radeq profiles and generate transmission spectra.
    """
    pyrat = pb.Pyrat('transit_WASP69b.cfg')
    pyrat.log.verb = 0
    pyrat.spec.specfile = None
    wl = pyrat.spec.wl
    press = pyrat.atm.press
    nwave = len(wl)
    nlayers, nspecies = np.shape(pyrat.atm.vmr)

    radeq_models = [
        'WASP69b_radeq_3x_0.59_cto_032_beta.atm',
        'WASP69b_radeq_3x_1.05_cto_032_beta.atm',
    ]
    nmodels = len(radeq_models)

    temp = []
    vmr = []
    for i in range(nmodels):
       atm = io.read_atm(radeq_models[i])
       species = list(atm[1])
       pressure = atm[2]
       temp.append(atm[3])
       vmr.append(atm[4])
    vmr = np.array(vmr)

    # Quenched abundances:
    p_quench = 2.0
    idx_quench = np.where(press >= p_quench)[0][0]

    quenched_vmr = np.copy(vmr)
    for k in range(nmodels):
        imol = species.index('NH3')
        vmr_quench = vmr[k,idx_quench,imol]
        quenched_vmr[k,press<p_quench,imol] = vmr_quench
        # SO2/H2S = 0.2
        imol = species.index('SO2')
        quenched_vmr[k,:,imol] = 2.0e-6

    spectra = np.zeros((nmodels, nwave))
    for k in range(nmodels):
        pyrat.run(temp[k], quenched_vmr[k])
        spectra[k] = np.copy(pyrat.spec.spectrum)

    # Species to showcase:
    molecs = [
        'H2O',
        'CO',
        'CO2',
        'CH4',
        'H2S',
        'SO2',
        'NH3',
        'potassium_vdw',
        'sodium_vdw',
    ]
    nmolecs = len(molecs)

    # Single-species contribution:
    contribution_spectra = np.zeros((nmodels, nmolecs, nwave))
    for i in range(nmolecs-1):
        skips = [mol for mol in molecs if mol!=molecs[i]]
        if 'sodium_vdw' not in skips:
            skips.remove('potassium_vdw')
        for k in range(nmodels):
            pyrat.run(temp[k], quenched_vmr[k], skip=skips)
            contribution_spectra[k,i] = np.copy(pyrat.spec.spectrum)

    np.savez(
        'WASP69b_transmission_spectra.npz',
        spectra=spectra,
        contribution_spectra=contribution_spectra,
        wl=wl,
        molecs=molecs,
        models=radeq_models,
        press=press,
        temp=temp,
        species=species,
        quenched_vmr=quenched_vmr,
    )


def make_jwst_simulation():
    """
    Simulate JWST TSOs for the transit models.
    """
    catalog = cat.Catalog()
    target = catalog.get_target('WASP-69 b')

    sed_type = 'phoenix'
    t_eff = target.teff
    logg_star = target.logg_star
    sed_model = jwst.find_closest_sed(t_eff, logg_star, sed_type)
    norm_band = '2mass,ks'
    norm_mag = target.ks_mag

    transit_dur = target.transit_dur
    t_start = 1.0
    t_settling = 0.75
    t_base = np.max([0.5*transit_dur, 1])
    obs_dur = t_start + t_settling + transit_dur + 2*t_base
    obs_dur = 7.5

    insts = [
        'niriss',
        'nirspec',
        'miri',
    ]
    ninst = len(insts)
    pandos = []
    for inst in insts:
        pando = jwst.PandeiaCalculation(inst)
        pando.set_scene(sed_type, sed_model, norm_band, norm_mag)
        if inst == 'niriss':
            pando.set_config(subarray='substrip96')
        pandos.append(pando)


    with np.load('WASP69b_transmission_spectra.npz') as d:
        spectra = d['spectra']
        wl = d['wl']
        models = d['models']
    nmodels = len(models)
    labels = [model[14:25] for model in models]

    # Run TSO simulations:
    obs_type = 'transit'
    ngroups = [3, 5, 15]
    obs_dur = [7.03, 8.27, 7.93]
    for j in range(nmodels):
        depth_model = wl, spectra[j]
        for i in range(ninst):
            pando = pandos[i]
            inst = insts[i]
            tso_file = f'tso_WASP69b_{labels[j]}_{inst}.pickle'
            print(tso_file)
            tso = pando.tso_calculation(
                obs_type, transit_dur,
                obs_dur[i], depth_model, ngroup=ngroups[i],
            )
            pando.save_tso(tso_file)


def make_obs_files():
    """
    Generate JWST transit spectra from the TSO models.
    """
    # Load models
    insts = [
        'niriss',
        'nirspec',
        'miri',
    ]
    ninst = len(insts)

    labels = [
        '3x_0.59_cto',
        '3x_1.05_cto',
    ]
    nmodels = len(labels)

    tsos = []
    for j in range(nmodels):
        tsos.append([])
        for i in range(ninst):
            inst = insts[i]
            label = labels[j]
            tso_file = f'tso_WASP69b_{label}_{inst}.pickle'
            with open(tso_file, 'rb') as f:
                tso = pickle.load(f)
            tso['label'] = inst
            tsos[j].append(tso)

    wl, spectrum = tso['input_depth']
    bin_wl = ps.constant_resolution_spectrum(0.6, 12.0, resolution=150)
    nwave = len(bin_wl)

    bin_spectra = np.zeros((nmodels, nwave))
    for i in range(nmodels):
        if tsos[i][0] is None:
            continue
        wl, spectrum = tsos[i][0]['input_depth']
        bin_spectra[i] = ps.bin_spectrum(bin_wl, wl, spectrum) / pc.percent

    resolutions = {
        'niriss': 180,
        'nirspec': 180,
        'miri': 80,
    }

    # The simulations with inflated noise
    observations = []
    for j in range(nmodels):
        obs_wl, obs_depths, obs_errors, widths, obs_inst = [], [], [], [], []
        for i in range(ninst):
            inst = insts[i]
            tso = tsos[j][i]
            if tso is None:
                continue
            res = resolutions[inst]
            obs_wave, obs_depth, obs_error, half_width = jwst.simulate_tso(
                tso, resolution=res, noiseless=False, inflation=1.5,
            )
            mask = obs_error < 5*np.median(obs_error)
            if inst == 'miri':
                mask &= obs_wave < 12.0
            obs_wl = np.append(obs_wl, obs_wave[mask])
            widths = np.append(widths, half_width[mask])
            obs_depths = np.append(obs_depths, obs_depth[mask]/pc.percent)
            obs_errors = np.append(obs_errors, obs_error[mask]/pc.percent)
            obs_inst = np.append(obs_inst, [inst for _ in obs_wave[mask]])
        observations.append((obs_wl, widths, obs_depths, obs_errors, obs_inst))

    # Save to file
    for j in range(nmodels):
        obs_wl, widths, obs_depths, obs_errors, inst_labs = observations[j]
        obs_file = f'../inputs/obs_WASP69b_transit_jwst_{labels[j]}_soss_g395h_lrs.dat'
        io.write_observations(
            obs_file, inst_labs,
            obs_wl, widths, obs_depths, obs_errors,
            depth_units='percent',
        )

