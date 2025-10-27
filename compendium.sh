# Define topdir (in your top working directory) to make your life easier:
topdir=`pwd`

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Setup:
pip install pyratbay==2.0.0
pip install gen_tso


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


# Download cross-section data:
# Mol   Source  Line-list
# -----------------------
# H2O   exomol  pokazatel
# H2S   exomol  ayt2
# HCN   exomol  harris
# C2H2  exomol  acety
# CO    HITEMP  Li
# CO2   ames    ai3000k
# CH4   exomol  mm
# NH3   exomol  coyute
# SO2   exomol  exoames
# TiO   exomol  toto
# VO    exomol  hyvo

cd $topdir/inputs
wget -i wget_cross_sections.txt

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Benchmarking:

# Comparison with ggchem for Woitke et al. (2018) network:
cd $topdir/benchmark_chemcat_ggchem
python ../code/fig_benchmark_ggchem.py

cd $topdir/benchmark_chemcat_fastchem
python ../code/fig_benchmark_fastchem.py

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# WASP-69b

# Radiative-equilibrium profiles:
cd $topdir/run_WASP69b
sh ../inputs/launch_wasp69b_radeq.sh

# Forward models
python make_WASP69b_transmission.py
python fig_WASP69b_transmission_spectra.py

# Retrievals
mpirun -n 128 pbay -c ret_WASP69b_transit_jwst_0.59_iso.cfg
mpirun -n 128 pbay -c ret_WASP69b_transit_jwst_0.59_slant.cfg
mpirun -n 128 pbay -c ret_WASP69b_transit_jwst_1.10_iso.cfg
