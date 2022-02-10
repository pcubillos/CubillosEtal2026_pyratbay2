import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

import tea.utils as u
sys.path.append('../code')
import benchmark_woitke


def main():
    # The composition from Woitke et al. (2018):
    species_lists = benchmark_woitke.load_species()

    # Read GGchem data:
    pressure, temperature, tea_species, vmr_tea = u.read_tea(
        'vmr_tea_benchmark_woitke.dat')

    # Read GGchem data:
    with open('Static_Conc.dat', 'r') as f:
        gg_lines = f.readlines()

    n_atoms, n_molecs, n_cond, nlayers = np.array(gg_lines[1].split(), int)
    nspecies = n_molecs + n_atoms + 1
    gg_columns = gg_lines[2].split()
    gg_species = gg_columns[3:nspecies+3]

    temperature = np.zeros(nlayers)
    pressure = np.zeros(nlayers)
    vmr = np.zeros((nlayers,nspecies))
    for i in range(nlayers):
        line = gg_lines[i+3].split()
        temperature[i] = line[0]
        pressure[i] = line[2]
        vmr[i] = line[3:nspecies+3]
    # Need to normalize by total VMR in layer:
    vmr_gg = 10**vmr / np.atleast_2d(np.sum(10**vmr, axis=1)).T

    # Translate molecule names:
    gg_dict = {
        'PCH': 'CHP',
        'HAlO': 'OAlH',
        'AlOCl': 'OAlCl',
        'AlOF': 'OAlF',
        'Al2O2': '(AlO)2',
        'Na2F2': '(NaF)2',
        'AlO2H': 'OAlOH',
        'AlF2O': 'OAlF2',
        'SN': 'NS',
        'P4O10': '(P2O5)2',
        'P4O6': '(P2O3)2',
        'K2F2': '(KF)2',
        'TiOF': 'OTiF',
        'TiOCl': 'OTiCl',
        'H2WO4': 'O2W(OH)2',
    }

    atom_names = 'MG LI AL CL NA CA SI TI FE ZR CR NI MN'.split()

    for i in range(nspecies):
        s = gg_species[i]
        for atom in atom_names:
            s = s.replace(atom, atom.capitalize())
        if s in gg_dict:
            s = gg_dict[s]
        gg_species[i] = s


    # Plot all panels as in Woitke et al. (2018):
    npanels = len(species_lists)
    pnames = [
        'hydrogen', 'lithium', 'carbon',
        'nitrogen', 'oxygen', 'fluorine',
        'sodium', 'magnesium', 'aluminum',
        'silicon', 'phosphorus', 'sulphur',
        'chlorine', 'potassium', 'calcium',
        'titanium', 'vanadium', 'chromium',
        'manganese', 'iron', 'nickel',
        'zirconium', 'tungsten',
    ]

    ylim = [
        (1e-14, 1e1),
        (1e-25, 1e-10),
        (1e-17, 1e-2),

        (1e-18, 1e-3),
        (1e-18, 1e-2),
        (1e-22, 1e-6),

        (1e-20, 1e-4),
        (1e-18, 1e-3),
        (1e-22, 1e-5),

        (1e-18, 1e-4),
        (1e-20, 1e-6),
        (1e-18, 1e-4),

        (1e-20, 1e-6),
        (1e-20, 1e-6),
        (1e-20, 1e-5),

        (1e-20, 1e-6),
        (1e-22, 1e-7),
        (1e-20, 1e-4),

        (1e-7, 1e-6),
        (1e-18, 1e-4),
        (1e-20, 1e-5),

        (1e-23, 1e-8),
        (1e-25, 1e-10),
    ]

    cols = [
        'royalblue', 'deepskyblue', 'mediumseagreen', 'red',
        'cyan', 'magenta', 'gold', 'black',
        'orange', 'cornflowerblue', 'salmon', 'greenyellow',
        'crimson', '0.75', 'darkviolet', 'sienna',
        'skyblue', 'aquamarine', 'dodgerblue', '0.4',
    ]


    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # As in Woitke et al. (2018):
    ny = 3
    for k in range(npanels):
        j = k%9
        if k%9 == 0:
            plt.figure(100+k//9, (8.5,9.0))
            plt.clf()
            plt.subplots_adjust(0.07, 0.05, 0.99, 0.97, wspace=0.18, hspace=0.2)
        ax = plt.subplot(3, ny, j+1)
        for i,name in enumerate(species_lists[k]):
            if name in gg_species:
                plt.loglog(
                    temperature, vmr_gg[:,gg_species.index(name)],
                    label=name, dashes=(), lw=2.0, color=cols[i])
            else:
                print(f'{name} not found in GGchem.')
            if name in tea_species:
                plt.loglog(
                    temperature, vmr_tea[:,tea_species.index(name)],
                    dashes=(10,1), lw=0.75, color='k', alpha=0.7)
            else:
                print(f'{name} not found in TEA.')
        plt.xlim(100, 6000)
        ax.set_ylim(ylim[k])
        plt.legend(loc='upper left', fontsize=6.0, framealpha=0.5)
        if j >= 6 or k>= npanels-3:
            plt.xlabel('Temperature (K)', fontsize=10)
        ax.set_title(pnames[k], fontsize=10)
        ax.tick_params(
            which='both', right=True, top=True, direction='in', labelsize=8)
        if j%ny == 0:
            plt.ylabel('Volume mixing ratio', fontsize=10)
        if (k+1)%9 == 0 or k == npanels-1:
            plt.savefig(
                f'../plots/benchmark_tea_ggchem_woitke2018_{k//9:02d}.png',
                dpi=300,
                facecolor='w')


    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # The real deal:
    nx = 3
    ny = 2
    panels = (
        species_lists[2],
        species_lists[3],
        species_lists[4],
        species_lists[6],
        species_lists[15],
        species_lists[19]+species_lists[0],
    )
    pnames = [
        'carbon', 'nitrogen', 'oxygen',
        'sodium', 'titanium', 'iron+hydrogen+helium',
    ]
    ylim = [
        (1e-17, 1e-2),
        (1e-18, 1e-3),
        (1e-18, 1e-2),
        (1e-20, 1e-4),
        (1e-20, 1e-6),
        (1e-18, 1e+1),
    ]
    cols = [
        'royalblue', 'deepskyblue', 'mediumseagreen', 'red', 'cyan',
        'magenta', 'gold', '0.4', 'darkorange', 'slateblue', '0.75',
        'greenyellow', 'crimson', 'paleturquoise',
    ]

    plt.figure(10, (8.5,5.5))
    plt.clf()
    plt.subplots_adjust(0.07, 0.08, 0.995, 0.96, wspace=0.18, hspace=0.2)
    for k in range(len(panels)):
        ax = plt.subplot(ny, nx, k+1)
        for i,name in enumerate(panels[k]):
            if name in gg_species:
                plt.loglog(
                    temperature, vmr_gg[:,gg_species.index(name)],
                    label=name, dashes=(), lw=2.5, color=cols[i])
            if name in tea_species:
                plt.loglog(
                    temperature, vmr_tea[:,tea_species.index(name)],
                    dashes=(10,2), lw=0.75, color='k', alpha=0.7)
        plt.gca().xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
        ax.set_xticks([100, 300, 1000, 3000])
        ax.set_xlim(100, 6000)
        ax.set_ylim(ylim[k])
        plt.legend(loc='lower left', fontsize=6.0, framealpha=0.5)
        if k >= nx:
            plt.xlabel('Temperature (K)', fontsize=10)
        ax.set_title(pnames[k], fontsize=10)
        ax.tick_params(
            which='both', right=True, top=True, direction='in', labelsize=8)
        if k%nx == 0:
            plt.ylabel('Volume mixing ratio', fontsize=10)
    plt.savefig('../plots/benchmark_tea_ggchem.png', dpi=300, facecolor='w')
    plt.savefig('../plots/benchmark_tea_ggchem.pdf', facecolor='w')


if __name__ == '__main__':
    main()

