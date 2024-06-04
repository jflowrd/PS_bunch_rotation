'''
Fitting the finement cavity impedance with resonators

@author: alasheen
'''

# General import
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# Impedance tools imports
script_path = os.path.dirname(os.path.realpath(__file__))
toolbox_path = script_path+'/../../../../impedance_toolbox'
sys.path.insert(0, toolbox_path)
from impedance_toolbox.impedance_params import ImpedanceParameters

# Loading input impedance and defining the limit in frequency for the 
# multi resonators fit
data_path = script_path+'/../Measurements/'

loaded_data = np.loadtxt(data_path+'/Zmeas-4-6-2013.csv', skiprows=1,
                         delimiter=';')
freq_data = loaded_data[:, 0]
abs_data = loaded_data[:, 1]*6  # Only one cell, six cells in total !!!
phase_data = loaded_data[:, 2]  # Only one cell, six cells in total !!!

# Importing an impedance source or a full impedance model
impedParams = ImpedanceParameters('.')

impedParams.addImpedanceTable(freq_data,
                              abs_data*np.cos(phase_data*np.pi/180),
                              abs_data*np.sin(phase_data*np.pi/180),
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

plt.figure('Impedance Re/Im')
plt.clf()
plt.plot(impedParams.freqArray/1e6, impedParams.impedance.real,
         label='Input real')
plt.plot(impedParams.freqArray/1e6, impedParams.impedance.imag,
         label='Input imag')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()


# Fitting with a single resonator

Rstart = [80*6]
fstart = [2.e6]
Qstart = [0.5]

RShuntBound = None
QBound = [0.5, 5]
freqBound = None

fittedParameters = impedParams.fitResonators(
    Rstart,
    fstart,
    Qstart,
    RShuntScale=1e2,
    freqScale=1e6,
    QScale=1,
    QBound=QBound,
    freqBound=freqBound,
    fitResidue='lin_real_and_imag',
    disp=True)

# Plotting the result with multiple resonators
plt.figure('Impedance single')
plt.plot(impedParams.freqArray/1e6, np.abs(impedParams.impedance),
         label='Input')
plt.plot(impedParams.resonatorForFit.frequency_array/1e6,
         np.abs(impedParams.resonatorForFit.impedance),
         label='%d Resonators' % (len(impedParams.resonatorForFit.R_S)))
plt.xlabel('Frequency [MHz]')
plt.ylabel('Abs. impedance [$\\Omega$]')
plt.legend(loc='best')
plt.xlim((0.1, 100))
plt.xscale('log')
plt.tight_layout()
plt.savefig(script_path+'/fitted_single_resonators_abs.png')

plt.figure('Impedance Re/Im single')
plt.plot(impedParams.freqArray/1e6, impedParams.impedance.real,
         label='Input real')
plt.plot(impedParams.freqArray/1e6, impedParams.impedance.imag,
         label='Input imag')
plt.plot(impedParams.resonatorForFit.frequency_array/1e6,
         impedParams.resonatorForFit.impedance.real,
         label='Real Fit')
plt.plot(impedParams.resonatorForFit.frequency_array/1e6,
         impedParams.resonatorForFit.impedance.imag,
         label='Imag Fit')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance [$\\Omega$]')
plt.legend(loc='best')
plt.xlim((0.1, 100))
plt.xscale('log')
plt.tight_layout()
plt.savefig(script_path+'/fitted_single_resonators_real_imag.png')

# Saving result into text file

sorted_freq = np.argsort(impedParams.resonatorForFit.frequency_R)

R_S_save = np.array(impedParams.resonatorForFit.R_S[sorted_freq], ndmin=2)
fr_save = np.array(impedParams.resonatorForFit.frequency_R[sorted_freq], ndmin=2)
Q_save = np.array(impedParams.resonatorForFit.Q[sorted_freq], ndmin=2)

saved_matrix = np.hstack((fr_save.T, R_S_save.T, Q_save.T))

np.savetxt(script_path+'/single_resonator_impedance.txt',
           saved_matrix,
           header='Impedance of Finemet cavity (total, 6 cells) \n' +
           'Single resonators fit\n' +
           'Impedance data: M. Paoluzzi\n' +
           'Fit: A. Lasheen\n' +
           'f_r [Hz]\t\t\t\tR_s [Ohm]\t\t\t\tQ')

# Fitting with multiple resonators
n_res_max = 7
RShuntBound = None
QBound = [0.5, 5]
freqBound = None
frequencyWindow = [0.1e6, 65e6]

fittedParameters = impedParams.fitMultiResonators(
    Rstart,
    fstart,
    Qstart,
    RShuntScale=1e2,
    freqScale=1e6,
    QScale=1,
    QBound=QBound,
    freqBound=freqBound,
    fitResidue='lin_real_and_imag',
    frequencyWindow=frequencyWindow,
    n_res_max=n_res_max,
    disp=True)

print(fittedParameters)

impedParams.resonatorForFit.R_S = impedParams.resonatorForFit.R_S[impedParams.resonatorForFit.frequency_R<25e6]
impedParams.resonatorForFit.Q = impedParams.resonatorForFit.Q[impedParams.resonatorForFit.frequency_R<25e6]
impedParams.resonatorForFit.frequency_R = impedParams.resonatorForFit.frequency_R[impedParams.resonatorForFit.frequency_R<25e6]
impedParams.resonatorForFit.n_resonators = len(impedParams.resonatorForFit.R_S)

# Recalculate final result
impedParams.resonatorForFit.imped_calc(impedParams.freqArray)

# Plotting the result with multiple resonators
plt.figure('Impedance')
plt.plot(impedParams.resonatorForFit.frequency_array/1e6,
         np.abs(impedParams.resonatorForFit.impedance),
         label='%d Resonators' % (len(impedParams.resonatorForFit.R_S)))
plt.xlabel('Frequency [MHz]')
plt.ylabel('Abs. impedance [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig(script_path+'/fitted_multi_resonators_abs.png')
plt.xlim((0.1, 100))
plt.xscale('log')
plt.savefig(script_path+'/fitted_multi_resonators_abs_log.png')

plt.figure('Impedance Re/Im')
plt.plot(impedParams.resonatorForFit.frequency_array/1e6,
         impedParams.resonatorForFit.impedance.real,
         label='Real Fit')
plt.plot(impedParams.resonatorForFit.frequency_array/1e6,
         impedParams.resonatorForFit.impedance.imag,
         label='Imag Fit')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig(script_path+'/fitted_multi_resonators_real_imag.png')
plt.xlim((0.1, 100))
plt.xscale('log')
plt.savefig(script_path+'/fitted_multi_resonators_real_imag_log.png')

# Plotting the discrepancy
discrepancy = np.abs(impedParams.resonatorForFit.impedance) - np.abs(impedParams.impedance)

# Plotting the discrepancy
plt.figure('Discrepancy')
plt.clf()
plt.plot(impedParams.freqArray/1e6, discrepancy,
         label='Abs')
plt.hlines(0, impedParams.freqArray[0]/1e6, impedParams.freqArray[-1]/1e6,
           linewidth=0.3)
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()
plt.xlim((0.1, 40))
plt.ylim((-60, 60))
plt.savefig(script_path+'/deviation_multi_resonators.png')

plt.show()

# Saving result into text file

sorted_freq = np.argsort(impedParams.resonatorForFit.frequency_R)

R_S_save = np.array(impedParams.resonatorForFit.R_S[sorted_freq], ndmin=2)
fr_save = np.array(impedParams.resonatorForFit.frequency_R[sorted_freq], ndmin=2)
Q_save = np.array(impedParams.resonatorForFit.Q[sorted_freq], ndmin=2)

saved_matrix = np.hstack((fr_save.T, R_S_save.T, Q_save.T))

np.savetxt(script_path+'/multi_resonators_impedance.txt',
           saved_matrix,
           header='Impedance of Finemet cavity (total, 6 cells) \n' +
           'Four resonators fit\n' +
           'Impedance data: M. Paoluzzi\n' +
           'Fit: A. Lasheen\n' +
           'f_r [Hz]\t\t\t\tR_s [Ohm]\t\t\t\tQ')
