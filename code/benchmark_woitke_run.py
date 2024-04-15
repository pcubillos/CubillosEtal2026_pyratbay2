import subprocess
import itertools
import pickle

import numpy as np
import chemcat as cat


def main():
    # The composition (Woitke et al. 2018):
    species_lists = dict(
        hydrogen = 'H2 H He'.split(),
        lithium = '(LiOH)2 LiOH LiCl LiAlF4 LiF LiH Li LiO LiN'.split(),
        carbon = 'CH4 CO CO2 Si(CH3)4 HCN CHP CH2 CN CH C'.split(),
        nitrogen = 'NH3 N2 HCN NH2 PN NO NH N CN'.split(),
        oxygen = (
            'H2O Mg(OH)2 Fe(OH)2 Ca(OH)2 OAlOH AlOH Al2O SiO (NaOH)2 '
            'MgOH CO CaOH OH O').split(),
        fluorine = (
            'OAlF2 TiF3 SiH3F TiF2 OAlF NaAlF4 AlF ZrF4 HF NaF '
            'MgClF TiF CaF2 MgF F').split(),
        sodium = (
            '(NaCl)2 (NaOH)2 NaCl NaOH Na NaAlF4 NaH NaF (NaF)2 Na2 NaCN '
            'NaO').split(),
        magnesium = (
            'Mg(OH)2 MgOH Mg MgCl2 MgH MgClF MgS MgCl MgF MgO MgN').split(),
        aluminum = (
            'OAlOH AlOH Al2O OAlF2 (AlO)2 OAlF OAlCl NaAlF4 AlF AlCl AlH '
            'AlS Al OAlH AlO').split(),
        silicon = (
            'SiH4 SiO SiS SiH3F SiH3 SiO2 Si(CH3)4 SiH2 SiH SiN Si '
            'Si2 SiC').split(),
        phosphorus = '(P2O3)2 PH3 PH2 PN PS PO2 CHP PO PH P2 P CP'.split(),
        sulfur = (
            'H2S NiS CrS SiS HS S2 COS FeS MgS AlS CaS SO S PS NS CS '
            'SO2 S2O').split(),
        chlorine = (
            'NiCl (NaCl)2 (KCl)2 TiOCl2 KCl NaCl TiCl3 ZrCl4 HCl CaCl2 '
            'AlCl MgCl2 FeCl2 CaCl MgCl Cl').split(),
        potassium = '(KCl)2 (KOH)2 KOH KCl K KH KF KCN (KF)2 K2 KO'.split(),
        calcium = 'Ca(OH)2 CaCl2 CaOH CaCl Ca CaF2 CaH CaS CaF CaO Ca2'.split(),
        titanium = (
            'TiO2 TiF3  TiOCl2 TiF2 TiCl3 OTiCl TiO TiCl2 OTiF TiS TiF '
            'TiCl Ti TiH').split(),
        vanadium = 'VO2 VO V VN'.split(),
        chromium = 'Cr CrS CrH CrO CrO2 CrN'.split(),
        manganese = 'MnH Mn MnS MnCl MnF MnO'.split(),
        iron = 'Fe(OH)2 FeCl2 Fe FeS FeH FeO FeCl FeF FeF2'.split(),
        nickel = 'NiS NiH Ni NiCl NiCl2 NiF NiO'.split(),
        zirconium = (
            'ZrCl4 ZrO2 ZrF4 ZrCl3 ZrO ZrF3 ZrCl2 ZrF2 ZrCl ZrH ZrF '
            'Zr ZrN').split(),
        tungsten = (
            'O2W(OH)2 (WO3)3 (WO3)4 W3O8 (WO3)2 WO2Cl2 WO3 WCl2 WO2 WO WCl '
            'W WF').split(),
    )

    # Set up parameter space:
    nlayers = 100
    pressure = np.tile(1.0, nlayers)
    temperature = np.logspace(np.log10(5999.9999), np.log10(100), nlayers)
    all_molecules = species_lists.values()
    molecules = list(set(itertools.chain.from_iterable(all_molecules)))

    # Run chemcat:
    net = cat.Network(
        pressure,
        temperature,
        molecules,
        e_source='asplund_2009',
    )
    vmr_cat = net.thermochemical_equilibrium()
    chemcat_species = list(net.species)

    # Run GGchem:
    proc = subprocess.Popen(
        '../ggchem/ggchem benchmark_ggchem_woitke_neutrals.in'.split(),
        stdout=subprocess.PIPE,
        universal_newlines=True,
    )
    proc.communicate()

    # Read GGchem output:
    with open('Static_Conc.dat', 'r') as f:
        gg_lines = f.readlines()
    n_atoms, n_molecs, n_cond, nlayers = np.array(gg_lines[1].split(), int)
    nspecies = n_molecs + n_atoms + 1
    gg_columns = gg_lines[2].split()
    ggchem_species = gg_columns[3:nspecies+3]

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
        s = ggchem_species[i]
        for atom in atom_names:
            s = s.replace(atom, atom.capitalize())
        if s in gg_dict:
            s = gg_dict[s]
        ggchem_species[i] = s

    results = {
        'temperature': temperature,
        'vmr_gg': vmr_gg,
        'species_gg': ggchem_species,
        'vmr_cat': vmr_cat,
        'species_cat': chemcat_species,
        'species_lists': species_lists,
    }

    with open('benchmark_chemcat_ggchem.pickle', 'wb') as handle:
        pickle.dump(results, handle, protocol=4)


if __name__ == '__main__':
    main()

