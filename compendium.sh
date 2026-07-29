# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 0. Install requirements (ideally in a new Python working environment)

pip install 'pyratbay==2.1.0'
pip install gen_tso
pip install mpi4py
pip install pymultinest
# Also, make sure to have multinest and MPI running (for retrieval runs)


# Define topdir (in your top working directory) to make your life easier
topdir=`pwd`


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 1. Benchmark chemcat

# 1.1 Install ggchem
cd $topdir
git clone https://github.com/pw31/GGchem ggchem
cd ggchem
git checkout 017a60e
cd $topdir/ggchem/src16
# This makefile worked for me, use what works for you
cp ../../inputs/ggchem_makefile makefile
make
cd $topdir
cp ggchem/data/dispol_* ggchem/data/DustChem.dat benchmark_chemcat_ggchem/data/

# 1.2 Install fastchem
cd $topdir
git clone https://github.com/exoclime/FastChem fastchem
cd fastchem
git checkout aa68c20
make demo2

# 1.3 Comparison with ggchem for Woitke et al. (2018) network
cd $topdir/benchmark_chemcat_ggchem
python ../code/fig_benchmark_ggchem.py

cd $topdir/benchmark_chemcat_fastchem
python ../code/fig_benchmark_fastchem.py


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 2. Benchmark radiative equilibrium

# 2.1 Download cross-section data
cd $topdir/inputs
wget -i wget_cross_sections.txt

# 2.2 Run radiative-equilibrium simulations
cd $topdir/benchmark_radeq
sh launch_radeq.sh

# 2.3 Make figures
cd $topdir
python code/fig_benchmark_sed.py

cd $topdir/benchmark_radeq
python ../code/radeq_contribution_functions.py
python ../code/fig_benchmark_radeq.py


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 3. VMR and TLS figures

# 3.1 Figure (5) Non-isobaric VMR profiles vs Moses VMR profiles
python code/fig_non_isobaric_vmr.py

# 3.2 Fetch PHOENIX New era SEDs
cd $topdir
python code/tls_contamination.py fetch

# Figure (6) TLS contamination
python code/tls_contamination.py figure


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 4. WASP-69b retrieval simulation

# 4.0.1 Make sure you installed multinest, pymultinest, MPI, and mpi4py
# 4.0.2 Download cross-section data (2.1) if you haven't already
# 4.0.3 See code/make_sed_WASP69b.py for alternaives to generate SEDs

cd $topdir/run_wasp69b
# 4.1 Compute radiative-equilibrium profile
pbay -c radeq_WASP69b_3x_0.59_cto_032_beta.cfg

# 4.2 Computer transmission model
python ../code/make_WASP69b_transmission.py
python ../code/fig_WASP69b_transmission_spectra.py

# 4.3 Simulate JWST observations
python ../code/simulate_WASP69b_jwst.py

# 4.4 Run retrievals (adjust number of CPUs at convenience)
mpirun -n 128 pbay -c ret_WASP69b_transit_jwst_0.59_iso_nm.cfg
mpirun -n 128 pbay -c ret_WASP69b_transit_jwst_0.59_iso_snm.cfg
mpirun -n 128 pbay -c ret_WASP69b_transit_jwst_0.59_slant_nm.cfg
mpirun -n 128 pbay -c ret_WASP69b_transit_jwst_0.59_slant_snm.cfg

# 4.5 Retrieval output figure
python ../code/fig_WASP69b_retrieval.py

