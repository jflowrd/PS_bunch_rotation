'''
Fitting the measured C80-89 cavity impedance with resonators

@author: alasheen
'''

# General import
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import hilbert
from scipy.constants import c, e, m_p

# Impedance tools imports
script_path = os.path.dirname(os.path.realpath(__file__))
toolbox_path = script_path+'/../../../../../../impedance_toolbox'
sys.path.insert(0,toolbox_path)
from impedance_toolbox.impedance_params import ImpedanceParameters

# Loading input impedance and defining the limit in frequency for the 
# multi resonators fit
n_bunches = np.arange(12, 80, 12)
# n_bunches = np.array([12])

cavity_name = 'C80-89'
data_path = script_path+'/../Measurements/'

charge = 39  # e
mass = 120.054e9  # eV
BField = 1.13615  # T
bending_radius = 70.079  # m
momentum = BField*bending_radius*charge*c

t_rev = (2*np.pi*100)/(c*momentum/np.sqrt(momentum**2.+mass**2.))
f_rev = 1/t_rev

for index_bunch in range(len(n_bunches)):

    loaded_data = np.loadtxt(data_path+'/'+cavity_name+'_b'+str(n_bunches[index_bunch])+'.txt')
    freq_data = loaded_data[:,0]
    abs_Z_data = loaded_data[:,1]
    
    new_freq = np.linspace(freq_data[0], freq_data[-1], 1000)
    abs_Z_data = np.interp(new_freq, freq_data, abs_Z_data)
    freq_data = new_freq
    
    # Importing an impedance source or a full impedance model
    impedParams = ImpedanceParameters('.')
     
    impedParams.addImpedanceTable(freq_data,
                                  abs_Z_data,
                                  np.zeros(len(abs_Z_data)),
                                  freqArrayInterp=freq_data)
     
    # Plotting the impedance
    plt.figure('Impedance')
    plt.clf()
    plt.plot(impedParams.freqArray/1e6, np.abs(impedParams.impedance), label='Input')
    plt.xlabel('Frequency [MHz]')
    plt.ylabel('Impedance [$\\Omega$]')
    plt.legend(loc='best')
    plt.tight_layout()
    
    # Getting a first approximation of the resonances of the impedance table 
    f_start = 169*f_rev
    Q_start = 100
    r_start = np.max(abs_Z_data)
    
    # Setting the fitting parameters, first run is coarse with only one resonator
    # and on a narrow bandwidth
    fittedParameters = impedParams.fitResonators(r_start,
                                                 f_start,
                                                 Q_start,
                                                 RShuntScale=1e3,
                                                 freqScale=1e6,
                                                 QScale=1e2,
                                                 freqBound=0.,
                                                 frequencyWindow=[f_start-10e6,
                                                                  f_start+10e6],
                                                 fitResidue='lin_total',
                                                 maxiter=10000)
       
    # Plotting the result with one resonator
    plt.figure('Impedance')
    plt.plot(impedParams.resonatorForFit.frequency_array/1e6,
             np.abs(impedParams.resonatorForFit.impedance),
             label='1 Resonator')
    plt.xlabel('Frequency [MHz]')
    plt.ylabel('Impedance [$\\Omega$]')
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig(script_path+'/fitted_b'+str(n_bunches[index_bunch])+'.png')
    
    # Saving result into text file
    
    sorted_freq = np.argsort(fittedParameters[1])
    
    R_S_save = np.abs(np.array(fittedParameters[0][sorted_freq], ndmin=2))
    fr_save = np.array(fittedParameters[1][sorted_freq], ndmin=2)
    Q_save = np.array(fittedParameters[2][sorted_freq], ndmin=2)
    
    saved_matrix = np.hstack((fr_save.T, R_S_save.T, Q_save.T))
    
    np.savetxt(script_path+'/single_resonator_b'+str(n_bunches[index_bunch])+'.txt',
               saved_matrix, header='Impedance of %s \nSingle resonator fit for %d bunches\nAuthor: A. Lasheen\nf_r [Hz]\tR_s [Ohm]\tQ' %(cavity_name, n_bunches[index_bunch]))
    
    plt.show()
     
    
#     # Fitting with multiple resonators
#     fitFreqMin = 35e6
#     fitFreqMax = 45e6
#     n_res_max = 2
#     fittedParameters = impedParams.fitMultiResonators(np.max(abs_Z_data),
#                                                       freq_data[abs_Z_data==np.max(abs_Z_data)],
#                                                       Q_start,
#                                                       RShuntScale=1e3,
#                                                       freqScale=1e6,
#                                                       QScale=1e2,
#                                                       RShuntBound=0.,
#                                                       freqBound=0.,
#                                                       QBound=None,
#                                                       fitResidue='lin_total',
#                                                       n_res_max=n_res_max,
#                                                       frequencyWindow=[fitFreqMin,
#                                                                        fitFreqMax])
#  
#     # Plotting the result with multiple resonators
#     impedParams.resonatorForFit.imped_calc(impedParams.freqArray)
#     plt.figure('Impedance')
#     plt.plot(impedParams.resonatorForFit.frequency_array/1e6,
#              np.abs(impedParams.resonatorForFit.impedance),
#              label='%d Resonators' %(n_res_max))
#     plt.xlabel('Frequency [MHz]')
#     plt.ylabel('Impedance [$\\Omega$]')
#     plt.legend(loc='best')
#     plt.tight_layout()
#     plt.savefig(script_path+'/fitted_b'+str(n_bunches[index_bunch])+'.png')
#       
#     # Plotting the discrepancy
#     discrepancy = np.abs(impedParams.impedance)-np.abs(impedParams.resonatorForFit.impedance)
#                
#     # Plotting the discrepancy
#     plt.figure('Discrepancy')
#     plt.clf()
#     plt.plot(impedParams.freqArray/1e6, np.abs(discrepancy),
#              label='Abs')
#     plt.hlines(0, impedParams.freqArray[0]/1e6, impedParams.freqArray[-1]/1e6,
#                linewidth=0.3)
#     plt.xlabel('Frequency [MHz]')
#     plt.ylabel('Impedance [$\\Omega$]')
#     plt.legend(loc='best')
#     plt.tight_layout()
#     plt.savefig(script_path+'/deviation_b'+str(n_bunches[index_bunch])+'.png')
#       
#       
#     # Saving result into text file
#       
#     sorted_freq = np.argsort(fittedParameters[1])
#       
#     R_S_save = np.array(fittedParameters[0][sorted_freq], ndmin=2)
#     fr_save = np.array(fittedParameters[1][sorted_freq], ndmin=2)
#     Q_save = np.array(fittedParameters[2][sorted_freq], ndmin=2)
#     print(Q_save)
#     saved_matrix = np.hstack((fr_save.T, R_S_save.T, Q_save.T))
#       
#     np.savetxt(script_path+'/multi_resonators_b'+str(n_bunches[index_bunch])+'.txt',
#                saved_matrix, header='Impedance of %s \nMultiple resonators fit for %d bunches\nAuthor: A. Lasheen\nf_r [Hz]\tR_s [Ohm]\tQ' %(cavity_name, n_bunches[index_bunch]))
#       
     
     
