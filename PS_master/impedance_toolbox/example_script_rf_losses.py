'''
Created on 9 May 2017

@author: alasheen
'''


# General import
import numpy as np
import matplotlib.pyplot as plt
plt.ioff()
from scipy.constants import e

# Impedance tools imports
from impedance_toolbox.machine_params import MachineParameters
from impedance_toolbox.impedance_params import ImpedanceParameters


# Input folders for impedance and machine parameters
machineParamsInput = './beams/LHC/LHC25ns_flat_top.yml'


# Generating the machine parameters and beam, in one object
machineParams = MachineParameters(machineParamsInput)

# machineParams.generateBeamCurrent(1, 2**19)
# machineParams.generateBeamCurrent(1, resolutionTime=1e-10)
machineParams.generateBeamCurrent(1, maxFreq=5e9)

machineParams.generateBeamSpectrum()



# Plotting the beam current

plt.figure('Beam current')
plt.clf()
plt.plot(machineParams.timeArray*1e6, machineParams.beamCurrent)
plt.xlabel('Time [$\\mu s$]')
plt.ylabel('Current [A]')
plt.tight_layout()


# Plotting the beam spectrum

plt.figure('Beam spectrum')
plt.clf()
plt.plot(machineParams.freqArray/1e9, np.abs(machineParams.beamSpectrum))
plt.plot(machineParams.freqArray/1e9, np.abs(machineParams.analyticalSpectrum))
plt.xlabel('Frequency [GHz]')
plt.ylabel('Absolute spectrum [A]')
plt.tight_layout()

plt.figure('Normalized power spectrum')
plt.clf()
plt.plot(machineParams.freqArray/1e9, 10*np.log10(np.abs(machineParams.beamSpectrum)**2./np.max(np.abs(machineParams.beamSpectrum)**2.)))
plt.plot(machineParams.freqArray/1e9, 10*np.log10(np.abs(machineParams.analyticalSpectrum)**2./np.max(np.abs(machineParams.beamSpectrum)**2.)))
plt.xlabel('Frequency [GHz]')
plt.ylabel('Normalized power spectrum [dB]')
plt.tight_layout()


# Importing an impedance source or a full impedance model
impedanceDirectory = '../SPS_impedance'

impedParams = ImpedanceParameters(impedanceDirectory, machineParams)

impedParams.addResonators([360.5, 29.7119e3], [149e6, 781e6], [1119, 2198])

impedParams.inducedVoltageGeneration()


# Plotting the impedance on top of spectrum

plt.figure('Impedance')
plt.clf()
plt.plot(machineParams.freqArray/1e9, np.abs(machineParams.beamSpectrum)/np.max(np.abs(machineParams.beamSpectrum))*np.max(np.abs(impedParams.impedance)))
plt.plot(impedParams.freqArray/1e9, np.abs(impedParams.impedance))
plt.xlabel('Frequency [GHz]')
plt.ylabel('Impedance [$\\Omega$]')
plt.tight_layout()


# Plotting the induced voltage on top of current

plt.figure('Induced Voltage')
plt.clf()
plt.plot(machineParams.timeArray*1e6, machineParams.beamCurrent/np.max(machineParams.beamCurrent)*np.max(np.abs(impedParams.inducedVoltage)))
plt.plot(impedParams.timeArray*1e6, impedParams.inducedVoltage)
plt.xlabel('Time [$\\mu s$]')
plt.ylabel('Induced voltage [V]')
plt.tight_layout()


# Using the induced voltage to compute various things
impedParams.getRFLosses()

print('-- Energy loss per turn: \n- Time: %.2e [eV]\n- Freq: %.2e [eV]\n' %(impedParams.energyLossTime, impedParams.energyLossFreq))

print('-- Energy loss per turn: \n- Time: %.2e [J]\n- Freq: %.2e [J]\n' %(impedParams.energyLossTime*e, impedParams.energyLossFreq*e))

print('-- Power loss per turn: \n- Time: %.2e [W]\n- Freq: %.2e [W]\n' %(impedParams.powerLossTime, impedParams.powerLossTime))

plt.show()




