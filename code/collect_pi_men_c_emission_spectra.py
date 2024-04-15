import numpy as np
import pyratbay as pb
import pyratbay.atmosphere as pa
import pyratbay.io as io
import pyratbay.constants as pc
import pyratbay.spectrum as ps
import chemcat as cat


def main():
    # The grid
    labs = [
        '0.1x solar',
        '1x solar',
        '10x solar',
        '100x solar',
        '1000x solar',
        '3000x solar',
        'water, C-depleted',
        'water world',
        'comet world',
        'no atmosphere',
    ]
    atms = [
        '0000.1x',
        '0001.0x',
        '0010.0x',
        '0100.0x',
        '1000.0x',
        '3000.0x',
        'water_dep1',
        'water',
        'comet',
    ]
    natms = len(atms)

    atm = io.read_atm('radeq_emission_pi_men_c_0001.0x.atm')
    press = atm[2]
    nlayers = len(press)
    species = atm[1]
    nmol = len(species)

    tp = np.zeros((natms, nlayers))
    vmr = np.zeros((natms, nlayers, nmol))
    for j in range(natms):
        atmfile = f'radeq_emission_pi_men_c_{atms[j]}.atm'
        atm = io.read_atm(atmfile)
        tp[j] = atm[3]
        vmr[j] = atm[4]

    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # Elemental composition
    metal = [
        -1.0, 0.0, 1.0, 2.0, 3.0, 3.5,
         0.0, 0.0, 0.0,
    ]
    e_ratio = [
        None, None, None, None, None, None,
        {'O_H': 0.5},
        {'O_H': 0.5},
        {'O_H': 0.75, 'C_H': 0.125},
    ]
    e_scale = {'C': -1.0}

    net = cat.Network(press, tp[0], species)
    elements = net.elements
    e_abundances = np.zeros((natms+1, len(elements)))
    for j in range(natms):
        scale = None if atms[j] != 'water_dep1'  else e_scale
        _ = net.thermochemical_equilibrium(
            temperature=tp[j],
            metallicity=metal[j],
            e_ratio=e_ratio[j],
            e_scale=scale,
        )
        e_abundances[j] = net.element_rel_abundance

    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # Spectra
    pyrat = pb.Pyrat('pi_men_c_radeq_0100.0x.cfg', log=None, mute=True)
    starflux = pyrat.spec.starflux

    spec_file = f'radeq_emission_pi_men_c_{atms[0]}.dat'
    d = np.loadtxt(spec_file, unpack=True)
    nwave = len(d[0])

    rprs = pyrat.atm.rplanet / pyrat.phy.rstar
    depth = rprs**2.0
    # No-atmosphere model, i.e.: black body with  no heat redistribution
    tstar = pyrat.phy.tstar
    rstar = pyrat.phy.rstar
    smaxis = pyrat.atm.smaxis
    teq, teq_err = pa.equilibrium_temp(tstar, rstar, smaxis, A=0.0, f=0.5)
    # 1363.9


    spectra = np.zeros((natms+1, nwave))
    flux_ratio = np.zeros((natms+1, nwave))
    for j in range(natms):
        spec_file = f'radeq_emission_pi_men_c_{atms[j]}.dat'
        wl, spectra[j] = np.loadtxt(spec_file, unpack=True)
        flux_ratio[j] = spectra[j] / starflux * depth / pc.ppm
    spectra[natms] = ps.bbflux(1e4/wl, teq)
    flux_ratio[natms] = spectra[natms] / starflux * depth / pc.ppm


    np.savez(
        'pi_men_c_eclipse_flux_ratios.npz',
        wl=wl,
        flux_ratio=flux_ratio*pc.ppm,
        rprs=rprs,
        starflux=starflux,
        tp=tp,
        vmr=vmr,
        species=species,
        elements=elements,
        e_abundances=e_abundances,
        atm_labels=labs,
    )


if __name__ == '__main__':
    main()
