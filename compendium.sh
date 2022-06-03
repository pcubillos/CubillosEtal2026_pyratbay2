# Define topdir (in your top working directory) to make your life easier:
topdir=`pwd`

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Setup:

# Install tea:
pip install tea_chemistry==0.2.0

# Install pyratbay:
# TBD: Updated to 1.X.0 after adding radeq option
pip install pyratbay==1.1.0

# Install ggchem:
cd $topdir
git clone https://github.com/pw31/GGchem ggchem
cd ggchem
git checkout 017a60e
cd $topdir/ggchem/src16
# This ggchem makefile worked for me, use what works for you:
cp $topdir/inputs/ggchem_makefile makefile
make
cp $topdir/ggchem/data/dispol_* $topdir/benchmark_tea_ggchem/data/
cp $topdir/ggchem/data/DustChem.dat $topdir/benchmark_tea_ggchem/data/


# Install fastchem:
cd $topdir
git clone https://github.com/exoclime/FastChem fastchem
cd fastchem
git checkout aa68c20
make demo2


# Download ExoMol/repack data (H2O, HCN, NH3, C2H2):
cd $topdir/inputs/opacity
wget -i wget_repack.txt

# Dowload HITEMP data (CO2, CO, CH4), and HITRAN data (H2O):
cd $topdir/inputs/opacity
wget -i wget_hitemp_CO2.txt
wget -i wget_hitran_H2O.txt
unzip '*.zip'
rm -f *.zip

wget https://hitran.org/hitemp/data/bzip2format/05_HITEMP2019.par.bz2
wget https://hitran.org/hitemp/data/bzip2format/06_HITEMP2020.par.bz2
bzip2 -d 05_HITEMP2019.par.bz2
bzip2 -d 06_HITEMP2020.par.bz2


# Download Kurucz stellar models:
#cd $topdir/inputs
#wget http://kurucz.harvard.edu/grids/gridp00odfnew/fp00k2odfnew.pck


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Benchmarking:

# Comparison with ggchem for Woitke et al. (2018) network:
cd $topdir/benchmark_tea_ggchem
python $topdir/code/benchmark_woitke.py
python $topdir/code/fig_benchmark_woitke.py

cd $topdir/benchmark_fastchem
python $topdir/code/fig_benchmark_fastchem.py

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Generate partition-function files for Exomol species:
cd $topdir/run_setup
pbay -pf exomol ../inputs/opacity/1H2-16O__POKAZATEL.pf
pbay -pf exomol \
    ../inputs/opacity/1H-12C-14N__Harris.pf \
    ../inputs/opacity/1H-13C-14N__Larner.pf
pbay -pf exomol ../inputs/opacity/12C2-1H2__aCeTY.pf
pbay -pf tips NH3 as_exomol

# Make TLI files:
cd $topdir/run_setup
pbay -c tli_H2O_exomol_pokazatel.cfg
pbay -c tli_HCN_exomol_harris-larner.cfg
pbay -c tli_NH3_exomol_coyute-byte.cfg
pbay -c tli_CO2_hitemp_2010.cfg
pbay -c tli_CO_hitemp_2019.cfg
pbay -c tli_CH4_hitemp_2020.cfg
pbay -c tli_C2H2_exomol_acety.cfg

# Make atmospheric files:
cd $topdir/run_setup
pbay -c atmosphere_solar_isothermal.cfg

# Make opacity files:
cd $topdir/run_setup
pbay -c opacity_H2O_0.3-33.0um.cfg
pbay -c opacity_HCN_0.3-33.0um.cfg
pbay -c opacity_NH3_0.3-33.0um.cfg
pbay -c opacity_CO2_0.3-33.0um.cfg
pbay -c opacity_CO_0.3-33.0um.cfg
pbay -c opacity_CH4_0.3-33.0um.cfg
pbay -c opacity_C2H2_0.3-33.0um.cfg



# Radiative-equilibrium runs:
cd $topdir/run_radeq_WASP69b
# Then, repeat on run_radeq_WASP107b folder
# cd $topdir/run_radeq_WASP69b
# cd $topdir/run_radeq_WASP80b
pbay -c radeq_000.1x-solar_CO-0.10.cfg
pbay -c radeq_000.1x-solar_CO-0.55.cfg
pbay -c radeq_000.1x-solar_CO-0.90.cfg
pbay -c radeq_000.1x-solar_CO-0.95.cfg
pbay -c radeq_000.1x-solar_CO-1.05.cfg
pbay -c radeq_000.1x-solar_CO-1.50.cfg
pbay -c radeq_000.1x-solar_CO-2.00.cfg
pbay -c radeq_000.1x-solar_CO-5.00.cfg

pbay -c radeq_001.0x-solar_CO-0.10.cfg
pbay -c radeq_001.0x-solar_CO-0.55.cfg
pbay -c radeq_001.0x-solar_CO-0.90.cfg
pbay -c radeq_001.0x-solar_CO-0.95.cfg
pbay -c radeq_001.0x-solar_CO-1.05.cfg
pbay -c radeq_001.0x-solar_CO-1.50.cfg
pbay -c radeq_001.0x-solar_CO-2.00.cfg
pbay -c radeq_001.0x-solar_CO-5.00.cfg

pbay -c radeq_010.0x-solar_CO-0.10.cfg
pbay -c radeq_010.0x-solar_CO-0.55.cfg
pbay -c radeq_010.0x-solar_CO-0.90.cfg
pbay -c radeq_010.0x-solar_CO-0.95.cfg
pbay -c radeq_010.0x-solar_CO-1.05.cfg
pbay -c radeq_010.0x-solar_CO-1.50.cfg
pbay -c radeq_010.0x-solar_CO-2.00.cfg
pbay -c radeq_010.0x-solar_CO-5.00.cfg

pbay -c radeq_100.0x-solar_CO-0.10.cfg
pbay -c radeq_100.0x-solar_CO-0.55.cfg
pbay -c radeq_100.0x-solar_CO-0.90.cfg
pbay -c radeq_100.0x-solar_CO-0.95.cfg
pbay -c radeq_100.0x-solar_CO-1.05.cfg
pbay -c radeq_100.0x-solar_CO-1.50.cfg
pbay -c radeq_100.0x-solar_CO-2.00.cfg
pbay -c radeq_100.0x-solar_CO-5.00.cfg


