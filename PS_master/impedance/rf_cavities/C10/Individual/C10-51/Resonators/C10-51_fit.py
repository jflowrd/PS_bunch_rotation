'''
Fitting the C10-51 cavity impedance with resonators

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
sys.path.insert(0,toolbox_path)
from impedance_toolbox.impedance_params import ImpedanceParameters

# Loading input impedance and defining the limit in frequency for the 
# multi resonators fit
harmonic_list = [8, 16, 21]
fitFreqMin_list = [1e6, 3e6, 5e6]
fitFreqMax_list = [25e6, 25e6, 25e6]
cavity_name = 'C10-51'
data_path = script_path+'/../Measurements/'

indexloop = 0
for harmonic in harmonic_list:

    loaded_data = np.loadtxt(data_path+'/measured_h'+str(harmonic)+'.txt')
    freq_data = loaded_data[:,0]
    real_data = loaded_data[:,1]
    imag_data = loaded_data[:,2]
    
    real_data[freq_data<1e6] = 0
    imag_data[freq_data<1e6] = 0
    
    # Importing an impedance source or a full impedance model
    impedParams = ImpedanceParameters('.')
     
    impedParams.addImpedanceTable(freq_data,
                                  real_data,
                                  imag_data,
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
    
    # Getting a first approximation of the resonances of the impedance table 
    r_start, f_start, Q_start = impedParams.fitInitialGuess()
    
    # Taking the max peak
    f_start = f_start[r_start==np.max(r_start)]
    Q_start = Q_start[r_start==np.max(r_start)]
    r_start = r_start[r_start==np.max(r_start)]
    
    # # Setting the fitting parameters, first run is coarse with only one resonator
    # # and on a narrow bandwidth
    # fittedParameters = impedParams.fitResonators(r_start,
    #                                              f_start,
    #                                              Q_start,
    #                                              RShuntScale=1e3,
    #                                              freqScale=1e6,
    #                                              QScale=1,
    #                                              frequencyWindow=[f_start-2e6,
    #                                                               f_start+2e6],
    #                                              fitResidue='lin_real_and_imag')
    #   
    # # Plotting the result with one resonator
    # plt.figure('Impedance')
    # plt.plot(impedParams.resonatorForFit.frequency_array/1e6,
    #          np.abs(impedParams.resonatorForFit.impedance),
    #          label='One-res narrow')
    # plt.xlabel('Frequency [MHz]')
    # plt.ylabel('Impedance [$\\Omega$]')
    # plt.legend(loc='best')
    # plt.tight_layout()
     
    # Setting the fitting parameters, first run is coarse with only one resonator
    # but on all the data
    fittedParameters = impedParams.fitResonators(r_start,
                                                 f_start,
                                                 Q_start,
                                                 RShuntScale=1e3,
                                                 freqScale=1e6,
                                                 QScale=1,
                                                 fitResidue='lin_real_and_imag')
       
    # Plotting the result with one resonator
    plt.figure('Impedance')
    plt.plot(impedParams.resonatorForFit.frequency_array/1e6,
             np.abs(impedParams.resonatorForFit.impedance),
             label='1 Resonator')
    plt.xlabel('Frequency [MHz]')
    plt.ylabel('Impedance [$\\Omega$]')
    plt.legend(loc='best')
    plt.tight_layout()
    
    
    # Saving result into text file
    
    sorted_freq = np.argsort(fittedParameters[1])
    
    R_S_save = np.array(fittedParameters[0][sorted_freq], ndmin=2)
    fr_save = np.array(fittedParameters[1][sorted_freq], ndmin=2)
    Q_save = np.array(fittedParameters[2][sorted_freq], ndmin=2)
    
    saved_matrix = np.hstack((fr_save.T, R_S_save.T, Q_save.T))
    
    np.savetxt(script_path+'/single_resonator_h'+str(harmonic)+'.txt',
               saved_matrix, header='Impedance of %s \nSingle resonator fit for h=%d\nAuthor: A. Lasheen\nf_r [Hz]\t\t\t\tR_s [Ohm]\t\t\t\tQ' %(cavity_name, harmonic))
    
     
    
    # Fitting with multiple resonators
    fitFreqMin = fitFreqMin_list[indexloop]
    fitFreqMax = fitFreqMax_list[indexloop]
    n_res_max = 9
    freqBound=[fitFreqMin_list[indexloop], 50e6]
    QBound=[0.5,10]
    fittedParameters = impedParams.fitMultiResonators(fittedParameters[0],
                                                      fittedParameters[1],
                                                      fittedParameters[2],
                                                      RShuntScale=1e3,
                                                      freqScale=1e6,
                                                      QScale=1,
                                                      freqBound=freqBound,
                                                      QBound=QBound,
                                                      fitResidue='lin_real_and_imag',
                                                      n_res_max=n_res_max,
                                                      frequencyWindow=[fitFreqMin,
                                                                       fitFreqMax])
    
    # Plotting the result with multiple resonators
    impedParams.resonatorForFit.imped_calc(impedParams.freqArray)
    plt.figure('Impedance')
    plt.plot(impedParams.resonatorForFit.frequency_array/1e6,
             np.abs(impedParams.resonatorForFit.impedance),
             label='%d Resonators' %(n_res_max))
    plt.xlabel('Frequency [MHz]')
    plt.ylabel('Impedance [$\\Omega$]')
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig(script_path+'/fitted_h'+str(harmonic)+'.png')
    
    # Plotting the discrepancy
    discrepancy = impedParams.impedance-impedParams.resonatorForFit.impedance
             
    # Plotting the discrepancy
    plt.figure('Discrepancy')
    plt.clf()
    plt.plot(impedParams.freqArray/1e6, discrepancy.real,
             label='Real')
    plt.plot(impedParams.freqArray/1e6, discrepancy.imag,
             label='Imag')
    plt.plot(impedParams.freqArray/1e6, np.abs(discrepancy),
             label='Abs')
    plt.hlines(0, impedParams.freqArray[0]/1e6, impedParams.freqArray[-1]/1e6,
               linewidth=0.3)
    plt.xlabel('Frequency [MHz]')
    plt.ylabel('Impedance [$\\Omega$]')
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig(script_path+'/deviation_h'+str(harmonic)+'.png')
    
    
    # Saving result into text file
    
    sorted_freq = np.argsort(fittedParameters[1])
    
    R_S_save = np.array(fittedParameters[0][sorted_freq], ndmin=2)
    fr_save = np.array(fittedParameters[1][sorted_freq], ndmin=2)
    Q_save = np.array(fittedParameters[2][sorted_freq], ndmin=2)
    
    saved_matrix = np.hstack((fr_save.T, R_S_save.T, Q_save.T))
    
    np.savetxt(script_path+'/multi_resonators_h'+str(harmonic)+'.txt',
               saved_matrix, header='Impedance of %s \nMultiple resonators fit for h=%d\nAuthor: A. Lasheen\nf_r [Hz]\t\t\t\tR_s [Ohm]\t\t\t\tQ' %(cavity_name, harmonic))
    
    
    indexloop += 1
    
    
