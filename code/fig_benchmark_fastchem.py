import subprocess

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import chemcat as cat


fc_dict = {
    'H3N': 'NH3',
    'CHN': 'HCN',
    'HO': 'OH',
    'OTi': 'TiO',
    'O2Ti': 'TiO2',
    'OV': 'VO',
    'O2V': 'VO2',
    'OS': 'SO',
    'O2S': 'SO2',
}


def run_fastchem(pressure, temperature, molecules):
    # Generate TP profile:
    with open('fc_atmosphere.dat', 'w') as f:
        f.write('# temperature in K, pressure in bar\n')
        for p,t in zip(pressure, temperature):
            f.write(f'{t:.6e}    {p:.6e}\n')
    # Run fastchem:
    proc = subprocess.Popen('../fastchem/fastchem fc_config.input'.split())
    proc.communicate()


def read_fastchem():
    data = np.loadtxt('fc_chem_output.dat', unpack=True)
    with open('fc_chem_output.dat', 'r') as f:
        header = f.readline()
    header = [h.strip().replace('1','') for h in header.split('\t')]
    header = [fc_dict[h] if h in fc_dict else h for h in header]
    return header, data


def main():
    # Setup atmosphere:
    nlayers = 101
    pressure = np.logspace(-10, 3, nlayers)

    neutrals = 'H2O CH4 CO CO2 NH3 N2 H2 HCN C2H2 C2H4 OH H He C N O'.split()
    ions = 'e- H- H+ H2+ He+'.split()
    alkali = 'Na Na- Na+ K K- K+'.split()
    metals = 'Mg Mg+ Fe Fe+'.split()
    sulfurs = 'S S2 HS H2S SO SO2'.split()
    metal_oxides = 'Ti TiO TiO2 Ti+ TiO+ V VO VO2 V+'.split()
    metal_oxides = 'Ti TiO TiO2 Ti+ V VO VO2 V+'.split()
    molecules = neutrals + ions + alkali + metals + metal_oxides + sulfurs

    temperatures = [300.0, 1400.0, 3000.0]
    temperature = np.tile(temperatures[0], nlayers)

    # Setup network:
    net = cat.Network(
        pressure, temperature, molecules,
        sources='janaf',
        e_source='asplund_2009',
    )

    # Run:
    species = chemcat_species = net.species
    ncases = len(temperatures)
    nspecies = len(chemcat_species)
    vmr_chemcat = np.zeros((ncases, nlayers, nspecies))
    vmr_fastchem = np.zeros((ncases, nlayers, nspecies))

    for i in range(ncases):
        temperature = np.tile(temperatures[i], nlayers)
        run_fastchem(pressure, temperature, chemcat_species)
        fastchem_species, data = read_fastchem()
        fastchem_idx = np.array([
            fastchem_species.index(mol) for mol in chemcat_species])
        vmr_fastchem[i] = data[fastchem_idx].T
        vmr_chemcat[i] = net.thermochemical_equilibrium(temperature)


    # Plot:
    fs = 10
    xlim = 1e-20, 3.0
    lw = 1.5

    dashes = {
        'ion': (1.75, 1.0),
        'cation': (5.0, 1.25),
        'neutral': (),
        'chemcat': (),
    }

    labels = ['ion', 'cation', 'neutral']
    legs = [
        Line2D([], [], color='black', lw=lw, dashes=dashes[label], label=label)
        for label in labels
    ]
    colors = cat.utils.resolve_colors(species)

    plt.figure(1, (5.0,6.5))
    plt.clf()
    plt.subplots_adjust(0.12, 0.06, 0.85, 0.99, hspace=0.06)
    for i in range(ncases):
        ax = plt.subplot(3,1,i+1)
        ax.tick_params(
            which='both', right=True, top=True, direction='in', labelsize=8)
        for j in range(nspecies):
            spec = species[j].replace('+','').replace('-','')
            if species[j].endswith('+'):
                dash = dashes['cation']
                label = None
            elif species[j].endswith('-'):
                dash = dashes['ion']
                label = species[j] if species[j]=='e-' else None
            else:
                dash = dashes['neutral']
                label = species[j]
            plt.loglog(
                vmr_fastchem[i,:,j], pressure, lw=2.5,
                color=colors[spec],
                dashes=dash, label=label,
            )
            plt.loglog(
                vmr_chemcat[i,:,j], pressure, lw=0.5, color='black',
                dashes=dashes['chemcat'],
            )
        if i == 0:
            ax.legend(
                handles=legs,
                loc=(1.01, 0.74),
                fontsize=fs-2,
                handlelength=1.0,
                handletextpad=0.4,
                borderpad=0.3,
            )
        if i == 2:
            ax.legend(
                loc=(1.01, 0.0),
                fontsize=fs-2,
                handlelength=1.0,
                labelspacing=0.2,
            )
            ax.set_xlabel("Volume mixing ratio", fontsize=fs)
        else:
            ax.set_xticklabels([])
        ax.text(
            0.03, 0.05, fr'$T={temperatures[i]:.0f}$ K',
            fontsize=fs-1, transform=ax.transAxes,
            bbox=dict(boxstyle='Round', facecolor='w', alpha=0.75),
        )
        ax.set_ylabel("Pressure (bar)", fontsize=fs)
        ax.set_ylim(np.amax(pressure), np.amin(pressure))
        ax.set_xlim(xlim)
    plt.savefig('../plots/benchmark_chemcat_fastchem.pdf')


if __name__ == '__main__':
    main()
