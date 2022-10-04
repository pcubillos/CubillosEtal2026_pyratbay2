import subprocess
import numpy as np


def main():
    """A very convoluted code to setup ggchem for neutral species only"""
    # Run ggchem with all species and ions to capture the screen output:
    proc = subprocess.Popen(
        '../ggchem/ggchem benchmark_ggchem_woitke.in'.split(),
        stdout=subprocess.PIPE,
        universal_newlines=True,
    )
    stdout = proc.communicate()

    # Extract only the lines with the species-source pairs:
    lines = [
        line for line in stdout[0].splitlines()
        if ' -> ' in line
    ]

    ggchem_log = 'ggchem_screen_output.log'
    with open(ggchem_log, 'w') as f:
        f.write('\n'.join(lines))


    # OK, now create a custom dispol file with only the species I want
    # (neutrals only, in this case):
    molecs, sources = np.loadtxt(
        ggchem_log, dtype=str, usecols=(1,3), unpack=True,
    )

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

    with open('data/dispol_benchmark_neutrals_test.dat', 'w') as f:
        f.write(f'{len(bench_dispol)//2}  Woitke Stock Barklem data\n')
        f.writelines(bench_dispol)


if __name__ == '__main__':
    main()
