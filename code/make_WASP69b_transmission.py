import numpy as np
import pyratbay as pb
import pyratbay.io as io


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
    ]
    nmodels = len(radeq_models)

    temp = []
    vmr = []
    for i in range(nmodels):
       atm = io.read_atm(radeq_models[i])
       species = list(atm[1])
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
    ]
    nmolecs = len(molecs)

    # Single-species contribution:
    contribution_spectra = np.zeros((nmodels, nmolecs, nwave))
    for i in range(nmolecs):
        skips = [mol for mol in molecs if mol!=molecs[i]]
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


if __name__ == '__main__':
    main()
