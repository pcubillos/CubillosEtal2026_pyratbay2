import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d as gaussf
import chemcat as cat
import scipy.interpolate as si
import pyratbay.atmosphere as pa


def main():
    # Moses+2011 Fig2
    press = np.array([
       7.97165329e-07, 8.26207410e-07, 8.32082617e-07, 8.32082617e-07,
       8.32082617e-07, 8.32082617e-07, 8.32082617e-07, 8.37999602e-07,
       8.62091306e-07, 8.93182245e-07, 8.99533713e-07, 8.99533713e-07,
       8.99533713e-07, 9.09145701e-07, 9.51998689e-07, 9.72452596e-07,
       9.72452596e-07, 1.04016776e-06, 1.05128250e-06, 1.05128250e-06,
       1.11654704e-06, 1.13650259e-06, 1.14053630e-06, 1.21564112e-06,
       1.25502829e-06, 1.33399141e-06, 1.43589743e-06, 1.53724482e-06,
       1.62099069e-06, 1.78301074e-06, 2.01939614e-06, 2.31774283e-06,
       2.83211356e-06, 3.65877213e-06, 5.38143118e-06, 6.69392092e-06,
       1.19327936e-05, 3.69172164e-05, 1.58588900e-04, 3.69244721e-04,
       8.33790989e-04, 1.96142575e-03, 4.99785691e-03, 1.34577840e-02,
       4.00762880e-02, 8.75836885e-02, 1.88937055e-01, 3.09167536e-01,
       4.46906604e-01, 6.19120923e-01, 8.66862536e-01, 1.18820916e+00,
       1.67550224e+00, 2.50931231e+00, 3.92128674e+00, 6.55446780e+00,
       1.11122242e+01, 2.64194132e+01, 7.78520235e+01, 2.09341416e+02,
       2.76792653e+02, 3.13731368e+02, 3.64395787e+02, 4.16746810e+02,
       4.68434396e+02, 5.16371806e+02, 5.64306489e+02, 6.09390791e+02,
       6.61546427e+02, 7.06816711e+02, 7.64113380e+02, 8.28293363e+02,
       8.89158500e+02, 9.55308392e+02, 1.02131270e+03, 1.09933262e+03,
       1.17310876e+03, 1.26546159e+03, 1.32678746e+03,
    ])
    temp = np.array([
       2478.89302815, 2435.53982963, 2253.86928343, 2390.12219308,
       2344.70455653, 2208.45164688, 2299.28691998, 2163.03401032,
       2117.61637377, 2072.19873722, 1935.94582757, 1981.36346412,
       2026.78110067, 1890.52819102, 1845.11055447, 1799.69291792,
       1754.27528137, 1708.85764482, 1663.44000827, 1618.02237172,
       1572.60473517, 1527.18709862, 1481.76946207, 1436.35182551,
       1390.93418896, 1348.72790045, 1300.09891586, 1250.03629376,
       1208.80487876, 1162.9284782 , 1110.68672707, 1066.52249932,
       1019.33534446,  971.62388788,  927.67429615,  890.30796789,
        881.01799678,  862.43805456,  878.95355876,  914.04900518,
        955.33776568,  992.49765013, 1017.27090643, 1031.72197261,
       1050.30191483, 1085.39736126, 1130.81499781, 1176.23263436,
       1221.65027091, 1267.06790746, 1312.48554401, 1357.90318056,
       1403.32081711, 1448.73845366, 1494.15609022, 1539.57372677,
       1566.41142109, 1560.21810702, 1531.31597467, 1565.72327508,
       1618.84814693, 1662.2603294 , 1707.4813528 , 1755.65157339,
       1810.70325405, 1856.98107311, 1909.79627925, 1958.13853634,
       2000.8281655 , 2046.04918891, 2094.21940949, 2143.42184909,
       2192.39490668, 2236.17410035, 2285.52399981, 2333.69422039,
       2384.61702501, 2432.9592821 , 2475.64891126,
    ])
    temp[0:37] = temp[37]

    # Moses+2011 Fig3
    NH3 = np.array([
       [7.22514548e-14, 1.49768027e-13, 3.16860050e-13, 6.70371996e-13,
        1.41828739e-12, 3.00063119e-12, 6.34835196e-12, 1.34310317e-11,
        2.84156602e-11, 6.01182221e-11, 1.27190451e-10, 2.69093300e-10,
        5.69313210e-10, 1.20448013e-09, 2.54828511e-09, 5.39133594e-09,
        1.14062995e-08, 2.41319907e-08, 5.10553816e-08, 1.08016451e-07,
        2.28527402e-07, 4.83489068e-07, 1.02290438e-06, 2.16413034e-06,
        4.57859038e-06, 9.68679634e-06, 1.50830987e-05, 1.67059008e-05,
        1.67059008e-05, 1.67059008e-05, 1.70077633e-05, 1.81108590e-05,
        2.15033279e-05, 5.14093916e-05, 6.52520591e-05],
       [1.31328610e-04, 1.67413223e-04, 2.06580974e-04, 2.46135425e-04,
        2.90342278e-04, 3.34023730e-04, 3.82358316e-04, 4.37687113e-04,
        4.83771431e-04, 5.40087744e-04, 5.96953882e-04, 6.53235222e-04,
        7.04169057e-04, 7.66711411e-04, 8.30640483e-04, 9.27336167e-04,
        9.99642108e-04, 1.11043930e-03, 1.25217945e-03, 1.45506120e-03,
        1.76872517e-03, 2.32927420e-03, 3.32324258e-03, 4.83725255e-03,
        7.07634999e-03, 1.32293093e-02, 4.49818046e-02, 1.46942017e-01,
        4.93413631e-01, 2.36976145e+00, 9.78227736e+00, 2.94805117e+01,
        9.05493659e+01, 3.60312549e+02, 1.16094193e+03],
    ])

    CH4_diseq = np.array([
       [1.96375957e-10, 4.29742852e-10, 9.40435488e-10, 1.96667336e-09,
        4.50370259e-09, 1.01971285e-08, 2.23150693e-08, 4.88335826e-08,
        1.06865847e-07, 2.33861796e-07, 4.84639807e-07, 1.08246069e-06,
        2.28952105e-06, 5.27286619e-06, 8.93836918e-06, 8.34992560e-06,
        8.63913871e-06, 9.24796398e-06, 9.24796398e-06, 9.24796398e-06,
        9.24796398e-06, 9.24796398e-06, 9.24796398e-06, 9.14359222e-06,
        1.14216999e-05, 2.31909127e-05, 4.83877258e-05, 8.26793866e-05,
        1.81344108e-04, 3.65706751e-04, 4.33582938e-04, 4.53589171e-04],
       [1.06366684e-09, 1.23190339e-09, 1.40366028e-09, 1.66142445e-09,
        1.97724769e-09, 2.32765053e-09, 2.87760702e-09, 3.38756926e-09,
        4.32686387e-09, 5.61751124e-09, 7.29314197e-09, 1.07887405e-08,
        1.76011418e-08, 4.91969348e-08, 9.93719240e-01, 2.96584895e-07,
        6.21756542e-03, 2.32550902e-06, 1.08966946e-05, 5.21818037e-05,
        2.34099624e-04, 1.05022498e-03, 2.40841301e-02, 1.46513766e-01,
        2.41392704e+00, 4.26825353e+00, 6.85815053e+00, 9.97022691e+00,
        1.71003314e+01, 4.37662274e+01, 2.23720665e+02, 1.10688308e+03],
    ])
    isort = np.argsort(CH4_diseq[1])
    CH4_diseq = CH4_diseq[:,isort]


    # Interpolate to common grid
    pressure = np.logspace(3, -8, 101)
    log_press = np.log10(pressure)
    kwargs = dict(bounds_error=False, fill_value='extrapolate')

    log_p = np.log10(press)
    t_interp = si.interp1d(log_p, temp, **kwargs)
    temperature = gaussf(t_interp(np.log10(pressure)), 2.0)

    log_p = np.log10(NH3[1])
    nh3_interp = si.interp1d(log_p, NH3[0], **kwargs)
    vmr_NH3_diseq = gaussf(nh3_interp(log_press), 1.0)

    log_p = np.log10(CH4_diseq[1])
    ch4_interp = si.interp1d(log_p, CH4_diseq[0], **kwargs)
    vmr_CH4_diseq = gaussf(ch4_interp(log_press), 1.0)

    # Compute atmosphere in equilibrium for Moses T(p)
    molecules = '''
        H  He  C  O  N  Na  K  S  Si  Mg Fe  Ti  V Al Cl
        H2  H2O  CH4  CO  CO2  HCN  NH3  N2  OH  C2H2  C2H4
        S2  SH  H2S  SO2  SO  TiO  VO  TiO2  VO2
        e-  H-  H+  H2+  He+  Na-  Na+  Mg+ K-  K+  Fe+  Ti+  V+
        MgOH MgH Mg(OH)2  FeO Fe(OH)2 FeS  SiO SiO2
        MgS Fe- NaOH KOH SiS SH- Si+
        HCl (NaCl)2  NaCl KCl Cl+ Cl- (KCl)2
        AlO- Al+ Al2O AlO2 AlH AlOH OAlOH'''.split()

    net = cat.Network(pressure, temperature, molecules)

    vmr = net.thermochemical_equilibrium()
    iH2O = list(net.species).index('H2O')
    iCO2 = list(net.species).index('CO2')
    vmr_H2O = vmr[:,iH2O]
    vmr_CO2 = vmr[:,iCO2]


    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # The plot
    labs = [
        'H2O (equilibrium)',
        'CO2 (equilibrium)',
        'NH3 (quenching + photochemistry)',
        'CH4 (vertical quenching)',
    ]
    nmodels = len(labs)
    vmr_mos = [
        vmr_H2O,
        vmr_CO2,
        vmr_NH3_diseq*2.5,
        vmr_CH4_diseq,
    ]
    cols = [
        'xkcd:blue',
        'red',
        'orange',
        'xkcd:green',
    ]

    texnames = [
        'm',
        '\\log\\ {\\rm VMR}_{0}',
        '\\log\\ p_{0}',
        '\\log\\ {\\rm VMR}_{\\rm min}',
        '\\log\\ {\\rm VMR}_{\\rm max}',
    ]
    text_cols = (
        'gray  black gray gray gray'.split(),
        'black black gray gray gray'.split(),
        'black black gray gray black'.split(),
        'black black gray black black'.split(),
    )

    vmr_model = pa.vmr_models.SlantVMR('CH4', pressure)
    pars = (
        # m     VMR0   p0    min    max
        [ 0.0,  -3.5,  1.0, -np.inf,  0.0],  # H2O -- iso
        [-0.15, -7.2,  1.0, -np.inf,  0.0],  # CO2 -- slope
        [ 4.0,  -4.8, -2.0, -np.inf, -4.4],  # NH3 -- quenc + photo
        [ 1.2,  -4.0,  1.0, -5.0,    -3.3],  # CH4 -- quench
    )
    vmr = [
        vmr_model(pars[0]),
        vmr_model(pars[1]),
        vmr_model(pars[2]),
        vmr_model(pars[3]),
    ]

    fs = 11.5
    fig = plt.figure(11)
    plt.clf()
    fig.set_size_inches(5.5, 4.6)
    plt.subplots_adjust(0.12, 0.105, 0.99, 0.78)
    ax = plt.subplot(111)
    for i in [0,1,3,2]:
        plt.plot(
            vmr_mos[i], pressure, lw=2.5, dashes=(5,1), c=cols[i], alpha=0.4,
            label=labs[i],
        )
        ax.plot(vmr[i], pressure, label='', lw=2.0, c=cols[i])
        ax.plot(
            10**pars[i][1], 10**pars[i][2],
            'o', ms=7, mec=cols[i], mfc='w', mew=1.5,
        )
    ax.text(
        10**pars[1][1], 10**pars[1][2]/1.5, '(VMR$_0, p_0$)',
        fontsize=fs-2, ha='right', va='bottom',
    )
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_yticks(np.logspace(3, -5, 5))
    ax.set_xlim(1e-11, 3e-3)
    ax.set_ylim(1e3, 0.5e-6)
    ax.set_xlabel('volume mixing ratio', fontsize=fs)
    ax.set_ylabel('pressure (bar)', fontsize=fs)
    ax.legend(loc='upper left', fontsize=fs-1.5, labelspacing=0.2)
    ax.tick_params(which='both', direction='in', labelsize=fs-0.5)
    ax.tick_params(which='major', length=5.0)
    # Over text
    for i in range(nmodels):
        x0 = 0.005 + i/nmodels
        ax.text(
            x0, 0.995, f'{i+1} free parameters', weight='bold',
            fontsize=fs-2.5, transform=fig.transFigure, va='top', color=cols[i],
        )
        for j in range(5):
            ax.text(
                x0, 0.99-0.032*(j+1), f'${texnames[j]} = {pars[i][j]}$',
                fontsize=fs-3, transform=fig.transFigure, color=text_cols[i][j],
                va='top',
            )
    plt.savefig('plots/non_isobaric_vmr_model.png', dpi=300)


if __name__ == '__main__':
    main()

