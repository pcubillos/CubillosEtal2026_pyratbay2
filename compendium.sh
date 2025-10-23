# Define topdir (in your top working directory) to make your life easier:
topdir=`pwd`

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Setup:

# pip install pyratbay==2.0.0
git clone -b radeq https://github.com/pcubillos/pyratbay
cd pyratbay
pip install -e .

# Install ggchem:
cd $topdir
git clone https://github.com/pw31/GGchem ggchem
cd ggchem
git checkout 017a60e
cd $topdir/ggchem/src16
# This makefile worked for me, use what works for you:
cp ../../inputs/ggchem_makefile makefile
make
cd $topdir
cp ggchem/data/dispol_* ggchem/data/DustChem.dat benchmark_tea_ggchem/data/

# Install fastchem:
cd $topdir
git clone https://github.com/exoclime/FastChem fastchem
cd fastchem
git checkout aa68c20
make demo2


# Download line-list data:
# Mol   Source  Line-list   Partition
# -----------------------------------
# H2O   exomol  pokazatel   exomol
# TiO   exomol  toto        exomol
# VO    exomol  hyvo        exomol
# C2H2  exomol  acety       exomol
# NH3   exomol  coyute      tips
# H2S   exomol  ayt2        tips
# SO2   exomol  exoames     tips
# CO    HITEMP  Li          tips
# CO2   ames    ai3000k     states
# CH4   exomol  mm          states
# HCN   exomol  harris      states

cd $topdir/inputs/opacity
wget -i wget_linelists.txt
bzip2 -d 05_HITEMP2019.par.bz2

# Generate partition-function files
cd $topdir/run_opacities
pbay -pf tips NH3
pbay -pf tips H2S
pbay -pf tips SO2
pbay -pf tips CO

pbay -pf exomol ../inputs/opacity/1H2-16O__POKAZATEL.pf
pbay -pf exomol ../inputs/opacity/12C2-1H2__aCeTY.pf
pbay -pf exomol ../inputs/opacity/51V-16O__HyVO.pf
pbay -pf exomol \
    ../inputs/opacity/46Ti-16O__Toto.pf \
    ../inputs/opacity/47Ti-16O__Toto.pf \
    ../inputs/opacity/48Ti-16O__Toto.pf \
    ../inputs/opacity/49Ti-16O__Toto.pf \
    ../inputs/opacity/50Ti-16O__Toto.pf

pbay -pf states 5.0 6000.0 5.0 ../inputs/opacity/12C-1H4__MM.states.bz2
pbay -pf states 5.0 6000.0 5.0 \
    ../inputs/opacity/1H-12C-14N__Harris.states.bz2 \
    ../inputs/opacity/1H-13C-14N__Larner.states.bz2
pbay -pf states 5.0 6000.0 5.0 \
    exomol/12C-16O2__AI3000K.states.bz2 \
    exomol/13C-16O2__AI3000K.states.bz2 \
    exomol/16O-12C-17O__AI3000K.states.bz2 \
    exomol/16O-12C-18O__AI3000K.states.bz2

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Benchmarking:

# Comparison with ggchem for Woitke et al. (2018) network:
cd $topdir/benchmark_chemcat_ggchem
python ../code/fig_benchmark_woitke.py

cd $topdir/benchmark_chemcat_fastchem
python ../code/fig_benchmark_fastchem.py

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Make TLI files:
cd $topdir/run_setup
pbay -c tli_H2O_exomol_pokazatel.cfg
pbay -c tli_HCN_exomol_harris-larner.cfg
pbay -c tli_NH3_exomol_coyute-byte.cfg
pbay -c tli_CO2_hitemp_2010.cfg
pbay -c tli_CO_hitemp_2019.cfg
pbay -c tli_CH4_hitemp_2020.cfg
pbay -c tli_C2H2_exomol_acety.cfg
pbay -c tli_SiO_exomol_siouvenir.cfg

pbay -c tli_CH4_hitemp_2020_extrap.cfg

# Make opacity files:
cd $topdir/run_setup
sh ../inputs/launch_opacity.sh

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# WASP-69b

# Radiative-equilibrium profiles:
cd $topdir/run_radeq_WASP69b
sh ../inputs/launch_wasp69b_radeq.sh

# Transmission spectra:
python make_WASP69b_transmission.py


# Then, repeat on run_radeq_WASP107b folder
# cd $topdir/run_radeq_WASP69b
# cd $topdir/run_radeq_WASP80b
pbay -c radeq_000.1x-solar_CO-0.10.cfg

