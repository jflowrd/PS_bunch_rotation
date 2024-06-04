'''
Created on 9 May 2017

@author: alasheen
'''

# General import
import numpy as np
import matplotlib.pyplot as plt
plt.ioff()

# Impedance tools imports
from impedance_toolbox.impedance_params import ImpedanceParameters


# Importing an impedance source or a full impedance model
impedanceDirectory = '../SPS_impedance'

impedParams = ImpedanceParameters(impedanceDirectory)

impedParams.addResonators([360.5, 29.7119e3], [149e6, 781e6], [1119, 2198], 
                          freqArray=np.linspace(0, 1e9, 50000))

# Plotting the impedance

plt.figure('Impedance')
plt.clf()
plt.plot(impedParams.freqArray/1e9, np.abs(impedParams.impedance))
plt.xlabel('Frequency [GHz]')
plt.ylabel('Impedance [$\\Omega$]')
plt.tight_layout()


plt.show()
    
    
# Getting a first approximation of the resonances of the impedance table
# You can do this manually, especially if the impedance is broadband

r_start, f_start, Q_start = impedParams.fitInitialGuess()


# Setting the fitting parameters, first run is coarse, by imposing strict bouncs
 
fittedParameters = impedParams.fitResonators(r_start, f_start, Q_start, 
                                             RShuntBound = 0.1, freqBound = 0.1, QBound = None,
                                             RShuntScale=1e3, freqScale=1e9, QScale=1e3)
 
 
# Setting the fitting parameters, second run is fine, all parameters are free
  
fittedParameters = impedParams.fitResonators(fittedParameters[0], fittedParameters[1], fittedParameters[2],
                                             RShuntScale=1e3, freqScale=1e9, QScale=1e3)
  
  
# Plotting the fit on top of the input data for the three steps above
   
plt.figure('Impedance')
plt.clf()
plt.plot(impedParams.freqArray/1e9, np.abs(impedParams.impedance))
plt.plot(impedParams.resonatorForFit.frequency_array/1e9, np.abs(impedParams.resonatorForFit.impedance))
plt.xlabel('Frequency [GHz]')
plt.ylabel('Impedance [$\\Omega$]')
plt.tight_layout()

plt.show()


# Fitting only in a given window in frequency (the second peak)
  
fittedParameters = impedParams.fitResonators(r_start[1]*0.8, f_start[1], Q_start[1], frequencyWindow = [0.6e9,1.0e9],
                                             RShuntScale=1e3, freqScale=1e9, QScale=1e3)


# Plotting the fit on top of the input data for the windowed fit
   
plt.figure('Impedance')
plt.clf()
plt.plot(impedParams.freqArray/1e9, np.abs(impedParams.impedance))
plt.plot(impedParams.resonatorForFit.frequency_array/1e9, np.abs(impedParams.resonatorForFit.impedance))
plt.xlabel('Frequency [GHz]')
plt.ylabel('Impedance [$\\Omega$]')
plt.tight_layout()
  
plt.show()
