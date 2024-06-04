'''
Fitting the C10-11 cavity impedance with closed gap relays with resonators

@author: alasheen
'''

# General import
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# Impedance tools imports
script_path = os.path.dirname(os.path.realpath(__file__))
toolbox_path = script_path+'/../../../../../impedance_toolbox'
sys.path.insert(0, toolbox_path)
from impedance_toolbox.impedance_params import ImpedanceParameters

# Loading input impedance and defining the limit in frequency for the
# multi resonators fit
harmonic_list = [8, 16, 21]
cavity_name = 'C10-11_ClosedGapRelays'
data_path = script_path+'/../Measurements/'

for harmonic in harmonic_list:

    loaded_data = np.loadtxt(
        data_path+'/h_'+str(harmonic)+'/f_h'+str(harmonic)+'.csv',
        delimiter=',')
    freq_data = loaded_data[loaded_data < 30e6]

    loaded_data = np.loadtxt(
        data_path+'/h_'+str(harmonic)+'/real_c11_h'+str(harmonic)+'.csv',
        delimiter=',')
    real_data = loaded_data[:len(freq_data)]

    loaded_data = np.loadtxt(
        data_path+'/h_'+str(harmonic)+'/image_c11_h'+str(harmonic)+'.csv',
        delimiter=',')
    imag_data = loaded_data[:len(freq_data)]

    real_data[freq_data < 2e6] = np.interp(freq_data[freq_data < 2e6],
                                           [0, freq_data[freq_data < 2e6][-1]],
                                           [0, real_data[freq_data < 2e6][-1]])
    imag_data[freq_data < 2e6] = np.interp(freq_data[freq_data < 2e6],
                                           [0, freq_data[freq_data < 2e6][-1]],
                                           [0, imag_data[freq_data < 2e6][-1]])

    # Importing an impedance source or a full impedance model
    impedParams = ImpedanceParameters('.')

    impedParams.addImpedanceTable(freq_data,
                                  real_data * 2,  # Two gaps in one cavity
                                  imag_data * 2,  # Two gaps in one cavity
                                  freqArrayInterp=freq_data)

    # Plotting the impedance
    plt.figure('Impedance')
    plt.clf()
    plt.plot(impedParams.freqArray/1e6, np.abs(impedParams.impedance),
             label='Input')
    plt.xlabel('Frequency [MHz]')
    plt.ylabel('Impedance [$\\Omega$]')
    plt.legend(loc='best')
    plt.tight_layout()

    plt.figure('Impedance real')
    plt.clf()
    plt.plot(impedParams.freqArray/1e6, impedParams.impedance.real,
             label='Input')
    plt.xlabel('Frequency [MHz]')
    plt.ylabel('Impedance Re [$\\Omega$]')
    plt.legend(loc='best')
    plt.tight_layout()

    plt.figure('Impedance imag')
    plt.clf()
    plt.plot(impedParams.freqArray/1e6, impedParams.impedance.imag,
             label='Input')
    plt.xlabel('Frequency [MHz]')
    plt.ylabel('Impedance Im [$\\Omega$]')
    plt.legend(loc='best')
    plt.tight_layout()

    # First approximation
    f_start = [23e6]
    Q_start = [10]
    r_start = [60]

    # Setting the fitting parameters, first run is coarse with only one resonator
    fittedParameters = impedParams.fitResonators(r_start,
                                                 f_start,
                                                 Q_start,
                                                 RShuntScale=1e2,
                                                 freqScale=1e6,
                                                 QScale=1,
                                                 freqBound=None,
                                                 ImZoverF=None,
                                                 ImZoverFBound=None,
                                                 ImZoverFScale=1e-6,
                                                 fitResidue='lin_real_and_imag',
                                                 frequencyWindow=[0, 30e6])

    # Plotting the result with one resonator
    plt.figure('Impedance')
    plt.plot(impedParams.resonatorForFit.frequency_array/1e6,
             np.abs(impedParams.resonatorForFit.impedance),
             label='1 Resonator')
    plt.xlabel('Frequency [MHz]')
    plt.ylabel('Impedance [$\\Omega$]')
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig(script_path+'/fitted_h'+str(harmonic)+'.png')

    # Plotting the result with one resonator
    plt.figure('Impedance real')
    plt.plot(impedParams.resonatorForFit.frequency_array/1e6,
             impedParams.resonatorForFit.impedance.real,
             label='1 Resonator')
    plt.xlabel('Frequency [MHz]')
    plt.ylabel('Impedance Re [$\\Omega$]')
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig(script_path+'/fitted_h'+str(harmonic)+'_real.png')

    # Plotting the result with one resonator
    plt.figure('Impedance imag')
    plt.plot(impedParams.resonatorForFit.frequency_array/1e6,
             impedParams.resonatorForFit.impedance.imag,
             label='1 Resonator')
    plt.xlabel('Frequency [MHz]')
    plt.ylabel('Impedance Im [$\\Omega$]')
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig(script_path+'/fitted_h'+str(harmonic)+'_imag.png')

    # Setting the fitting parameters, first run is coarse with only one resonator
    fittedParameters = impedParams.fitResonators(fittedParameters[0],
                                                 fittedParameters[1],
                                                 fittedParameters[2],
                                                 RShuntScale=1e2,
                                                 freqScale=1e6,
                                                 QScale=1,
                                                 freqBound=None,
                                                 ImZoverF=1*2.1e-6 * 2,  # Two gaps in one cavity
                                                 ImZoverFBound=None,
                                                 ImZoverFScale=1e-6,
                                                 fitResidue='lin_real_and_imag',
                                                 frequencyWindow=[0, 30e6])

    # Plotting the result with one resonator
    plt.figure('Impedance')
    plt.plot(impedParams.resonatorForFit.frequency_array/1e6,
             np.abs(impedParams.resonatorForFit.impedance),
             label='1 Resonator and ImZ/f')
    plt.xlabel('Frequency [MHz]')
    plt.ylabel('Impedance [$\\Omega$]')
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig(script_path+'/fitted_h'+str(harmonic)+'.png')

    plt.figure('Impedance real')
    plt.plot(impedParams.resonatorForFit.frequency_array/1e6,
             impedParams.resonatorForFit.impedance.real,
             label='1 Resonator and ImZ/f')
    plt.xlabel('Frequency [MHz]')
    plt.ylabel('Impedance Re [$\\Omega$]')
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig(script_path+'/fitted_h'+str(harmonic)+'_real.png')

    plt.figure('Impedance imag')
    plt.plot(impedParams.resonatorForFit.frequency_array/1e6,
             impedParams.resonatorForFit.impedance.imag,
             label='1 Resonator and ImZ/f')
    plt.xlabel('Frequency [MHz]')
    plt.ylabel('Impedance Im [$\\Omega$]')
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig(script_path+'/fitted_h'+str(harmonic)+'_imag.png')

    # Saving result into text file
    R_S_save = np.array(fittedParameters[0], ndmin=2)
    fr_save = np.array(fittedParameters[1], ndmin=2)
    Q_save = np.array(fittedParameters[2], ndmin=2)
    ImZ_over_F_save = np.array(fittedParameters[3], ndmin=2)

    saved_matrix = np.hstack((fr_save.T, R_S_save.T, Q_save.T, ImZ_over_F_save.T))

    np.savetxt(script_path+'/single_resonator_and_ImZ_over_f_h'+str(harmonic)+'.txt',
               saved_matrix, header='Impedance of %s \nSingle resonator and ImZ/f fit for h=%d\nAuthor: A. Lasheen\nf_r [Hz]\t\t\tR_s [Ohm]\t\t\tQ\t\t\tImZ/f [Hz-1]' %(cavity_name, harmonic))

    # Fitting with multiple resonators
    fittedParameters = impedParams.fitResonators(
        [fittedParameters[0], -0.98*fittedParameters[0]],
        [fittedParameters[1], 0.98*fittedParameters[1]],
        [fittedParameters[2], fittedParameters[2]],
        RShuntScale=1e2, freqScale=1e6, QScale=1,
        freqBound=None, ImZoverF=fittedParameters[3],
        ImZoverFBound=None, ImZoverFScale=1e-6,
        fitResidue='lin_real_and_imag', frequencyWindow=[0, 30e6])

#     n_res_max = 2
#     fittedParameters = impedParams.fitMultiResonators(fittedParameters[0],
#                                                       fittedParameters[1],
#                                                       fittedParameters[2],
#                                                       RShuntScale=1e2,
#                                                       freqScale=1e6,
#                                                       QScale=1,
#                                                       freqBound=None,
#                                                       ImZoverF=fittedParameters[3],
#                                                       ImZoverFBound=None,
#                                                       ImZoverFScale=1e-6,
#                                                       fitResidue='lin_real_and_imag',
#                                                       frequencyWindow=[0, 30e6],
#                                                       n_res_max=n_res_max)

    # Plotting the result with one resonator
    plt.figure('Impedance')
    plt.plot(impedParams.resonatorForFit.frequency_array/1e6,
             np.abs(impedParams.resonatorForFit.impedance),
             label='2 Resonators and ImZ/f')
    plt.xlabel('Frequency [MHz]')
    plt.ylabel('Impedance [$\\Omega$]')
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig(script_path+'/fitted_h'+str(harmonic)+'.png')

    plt.figure('Impedance real')
    plt.plot(impedParams.resonatorForFit.frequency_array/1e6,
             impedParams.resonatorForFit.impedance.real,
             label='2 Resonators and ImZ/f')
    plt.xlabel('Frequency [MHz]')
    plt.ylabel('Impedance Re [$\\Omega$]')
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig(script_path+'/fitted_h'+str(harmonic)+'_real.png')

    plt.figure('Impedance imag')
    plt.plot(impedParams.resonatorForFit.frequency_array/1e6,
             impedParams.resonatorForFit.impedance.imag,
             label='2 Resonators and ImZ/f')
    plt.xlabel('Frequency [MHz]')
    plt.ylabel('Impedance Im [$\\Omega$]')
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig(script_path+'/fitted_h'+str(harmonic)+'_imag.png')

    # Saving result into text file
    R_S_save = np.array(fittedParameters[0], ndmin=2)
    fr_save = np.array(fittedParameters[1], ndmin=2)
    Q_save = np.array(fittedParameters[2], ndmin=2)
    ImZ_over_F_save = fittedParameters[3] * np.ones(R_S_save.shape) / np.prod(R_S_save.shape)

    saved_matrix = np.hstack((fr_save.T, R_S_save.T, Q_save.T, ImZ_over_F_save.T))

    np.savetxt(script_path+'/multi_resonator_and_ImZ_over_f_h'+str(harmonic)+'.txt',
               saved_matrix, header='Impedance of %s \nMulti resonator and ImZ/f fit for h=%d\nThe ImZ/n was divided over the number of resonators and should be summed\nAuthor: A. Lasheen\nf_r [Hz]\t\t\tR_s [Ohm]\t\t\tQ\t\t\tImZ/f [Hz-1]' %(cavity_name, harmonic))

plt.show()