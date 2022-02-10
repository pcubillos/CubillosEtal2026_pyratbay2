import subprocess
import itertools

import numpy as np
import tea


def load_species():
    # The composition (Woitke et al. 2018):
    hydrogen = 'H2 H He'.split()
    lithium = '(LiOH)2 LiOH LiCl LiAlF4 LiF LiH Li LiO LiN'.split()
    carbon = 'CH4 CO CO2 Si(CH3)4 HCN CHP CH2 CN CH C'.split()
    nitrogen = 'NH3 N2 HCN NH2 PN NO NH N CN'.split()
    oxygen = (
        'H2O Mg(OH)2 Fe(OH)2 Ca(OH)2 OAlOH AlOH Al2O SiO (NaOH)2 '
        'MgOH CO CaOH OH O').split()
    fluorine = (
        'OAlF2 TiF3 SiH3F TiF2 OAlF NaAlF4 AlF ZrF4 HF NaF '
        'MgClF TiF CaF2 MgF F').split()
    sodium = (
        '(NaCl)2 (NaOH)2 NaCl NaOH Na NaAlF4 NaH NaF (NaF)2 Na2 NaCN '
        'NaO').split()
    magnesium = 'Mg(OH)2 MgOH Mg MgCl2 MgH MgClF MgS MgCl MgF MgO MgN'.split()
    aluminum = (
        'OAlOH AlOH Al2O OAlF2 (AlO)2 OAlF OAlCl NaAlF4 AlF AlCl AlH '
        'AlS Al OAlH AlO').split()
    silicon = (
        'SiH4 SiO SiS SiH3F SiH3 SiO2 Si(CH3)4 SiH2 SiH SiN Si Si2 SiC').split()
    phosphorus = '(P2O3)2 PH3 PH2 PN PS PO2 CHP PO PH P2 P CP'.split()
    sulphur = 'H2S NiS CrS SiS HS S2 COS FeS MgS AlS CaS SO S PS NS CS'.split()
    chlorine = (
        'NiCl (NaCl)2 (KCl)2 TiOCl2 KCl NaCl TiCl3 ZrCl4 HCl CaCl2 '
        'AlCl MgCl2 FeCl2 CaCl MgCl Cl').split()
    potassium = '(KCl)2 (KOH)2 KOH KCl K KH KF KCN (KF)2 K2 KO'.split()
    calcium = 'Ca(OH)2 CaCl2 CaOH CaCl Ca CaF2 CaH CaS CaF CaO Ca2'.split()
    titanium = (
        'TiO2 TiF3  TiOCl2 TiF2 TiCl3 OTiCl TiO TiCl2 OTiF TiS TiF '
        'TiCl Ti TiH').split()
    vanadium = 'VO2 VO V VN'.split()
    chromium = 'Cr CrS CrH CrO CrO2 CrN'.split()
    manganese = 'MnH Mn MnS MnCl MnF MnO'.split()
    iron = 'Fe(OH)2 FeCl2 Fe FeS FeH FeO FeCl FeF FeF2'.split()
    nickel = 'NiS NiH Ni NiCl NiCl2 NiF NiO'.split()
    zirconium = (
        'ZrCl4 ZrO2 ZrF4 ZrCl3 ZrO ZrF3 ZrCl2 ZrF2 ZrCl ZrH ZrF Zr ZrN').split()
    tungsten = (
        'O2W(OH)2 (WO3)3 (WO3)4 W3O8 (WO3)2 WO2Cl2 WO3 WCl2 WO2 WO WCl '
        'W WF').split()

    # Put all species in a single list (remove any duplicates):
    molecules_lists = [
        hydrogen, lithium, carbon, nitrogen, oxygen, fluorine,
        sodium, magnesium, aluminum, silicon, phosphorus, sulphur,
        chlorine, potassium, calcium, titanium, vanadium, chromium,
        manganese, iron, nickel, zirconium, tungsten,
    ]
    return molecules_lists


def run(molecules_lists):
    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # tea run:
    nlayers = 100
    pressure = np.tile(1.0, nlayers)
    temperature = np.logspace(np.log10(5999.9999), np.log10(100), nlayers)
    molecules = list(set(itertools.chain.from_iterable(molecules_lists)))

    tea_net = tea.Tea_Network(pressure, temperature, molecules)
    vmr = tea_net.thermochemical_equilibrium()

    tea.utils.write_tea(
        'vmr_tea_benchmark_woitke.dat',
        tea_net.species, tea_net.pressure, tea_net.temperature, vmr)


    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # ggchem run:
    proc = subprocess.Popen(
        '../ggchem/ggchem benchmark_ggchem_woitke_neutrals.in'.split(),
        stdout=subprocess.PIPE,
        universal_newlines=True)
    proc.communicate()


def make_custom_ggchem_dispol():
    """A very convoluted code to setup ggchem for neutral species only"""
    # Run ggchem with all species and ions to capture the screen output:
    proc = subprocess.Popen(
        '../ggchem/ggchem benchmark_ggchem_woitke.in'.split(),
        stdout=subprocess.PIPE,
        universal_newlines=True)
    stdout = proc.communicate()

    # Extract only the lines with the species-source pairs:
    lines = [
        line for line in stdout[0].splitlines()
        if ' -> ' in line]

    ggchem_log = 'ggchem_screen_output.log'
    with open(ggchem_log, 'w') as f:
        f.write('\n'.join(lines))


    # OK, now create a custom dispol file with only the species I want
    # (neutrals only, in this case):
    molecs, sources = np.loadtxt(
        ggchem_log, dtype=str, usecols=(1,3), unpack=True)

    # Sort out ggchem sources into a single file:
    with open('../ggchem/data/dispol_BarklemCollet.dat', 'r') as f:
        nspecies = int(f.readline().split()[0])
        b_lines = f.readlines()[0:2*nspecies]
        b_molecs = [l.split()[0].upper() for l in b_lines[::2]]

    with open('../ggchem/data/dispol_StockKitzmann_withoutTsuji.dat', 'r') as f:
        nspecies = int(f.readline().split()[0])
        s_lines = f.readlines()[0:2*nspecies]
        s_molecs = [l.split()[0].upper() for l in s_lines[::2]]

    with open('../ggchem/data/dispol_WoitkeRefit.dat', 'r') as f:
        nspecies = int(f.readline().split()[0])
        w_lines = f.readlines()[0:2*nspecies]
        w_molecs = [l.split()[0].upper() for l in w_lines[::2]]

    bench_dispol = []
    for mol, source in zip(molecs, sources):
        # Remove ions:
        if '-' in mol or '+' in mol:
            continue

        if mol == 'CNN':
            mol = 'CN2_CNN'
        elif mol == 'NCN':
            mol = 'CN2_NCN'
        elif mol == 'CLOCL':
            mol = 'CL2O_CLOCL'
        elif mol == 'CLCLO':
            mol = 'CL2O_CLCLO'
        elif mol == 'CLO2CL':
            mol = 'CL2O2_CLO2CL'
        elif mol == 'CLOCLO':
            mol = 'CL2O2_CLOCLO'
        elif mol == 'OFO':
            mol = 'FO2_OFO'
        elif mol == 'FOO':
            mol = 'FO2_FOO'

        if source[7] == 'B':
            idx = b_molecs.index(mol.upper())
            bench_dispol += b_lines[idx*2:(idx+1)*2]
        if source[7] == 'S':
            idx = s_molecs.index(mol.upper())
            bench_dispol += s_lines[idx*2:(idx+1)*2]
        if source[7] == 'W':
            idx = w_molecs.index(mol.upper())
            bench_dispol += w_lines[idx*2:(idx+1)*2]

    with open('data/dispol_benchmark_neutrals.dat',  'w') as f:
        f.write(f'{len(bench_dispol)//2}  Woitke Stock Barklem data\n')
        f.writelines(bench_dispol)


if __name__ == '__main__':
    molecules_lists = load_species()
    run(molecules_lists)
