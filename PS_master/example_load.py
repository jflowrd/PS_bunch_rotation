# coding: utf8
'''
Loading impedance source

@author: alasheen
'''

# General imports
import os
import numpy as np
import matplotlib.pyplot as plt

# Toolbox import
from impedance_toolbox.impedance_toolbox.machine_params import MachineParameters
from impedance_toolbox.impedance_toolbox.impedance_params import ImpedanceParameters

# Scenario import
from impedance_loader import ImpedanceLoader

script_path = os.path.dirname(os.path.realpath(__file__))
toolbox_path = script_path + '/impedance_toolbox'

# Case selector
model = 'trisplit'  # top, trisplit
beam = 'LHC25ns_C1940_end_flat_plateau'
# LHC25ns_C2595_flat_top, LHC25ns_flat_top_after_four_split, LHC25ns_C1940_end_flat_plateau
n_turns = 50
maxFreq = 6e9

fig_path = script_path + '/figures/Beam_' + beam + '_Impedance_' + model + '/'
# try:
#     os.mkdir(fig_path)
# except:
#     pass

# Making machine parameters
machineParamsInput = toolbox_path + '/beams/PS/LHC25ns/' + beam + '.yml'

# Generating the machine parameters and beam, in one object
machineParams = MachineParameters(machineParamsInput)
machineParams.generateBeamCurrent(n_turns, maxFreq=maxFreq)
machineParams.generateBeamSpectrum()

# Loading scenario at flat top
PS_loader = ImpedanceLoader(MODEL=model,
                            folder=script_path + '/impedance',
                            freq_array=machineParams.freqArray)

# Whole impedance model
PS_loader.importImpedancePS()

# Loading element by element

# # C10
# PS_loader.importC10(n_elements=10,)
# #                             method='rf_cavities/C10/All/Resonators/multi_resonators_h21.txt')

# # C10 1TFB
# PS_loader.importC10_1TFB(n_elements=10,
#                                  main_harmonic_FB=True)

# # C10 ClosedGaps
# PS_loader.importC10_ClosedGapRelay(n_elements=1)

# # C20
# PS_loader.importC20(n_elements=1)

# # C20 MHFB
# PS_loader.importC20_MHFB(n_elements=1)

# # C40
# PS_loader.importC40()
# # method_C40 = 'DFB'  # DFB, Resonators
# # fname = '_DFB'  # _DFB, ''
# # PS_loader.importC40(filename_list='/rf_cavities/C40/Individual/C77/%s/single_resonator%s_b72.txt'%(method_C40, fname),
# #                             method=method_C40)

# # C40 MHFB
# PS_loader.importC40_MHFB()

# # C40 HOMs
# PS_loader.importC40_HOMs()

# # C80
# PS_loader.importC80()
# method_C80 = 'DFB'  # DFB, Resonators
# fname = '_DFB'  # _DFB, ''
# PS_loader.importC80(filename_list='/rf_cavities/C80/Individual/C89/%s/single_resonator%s_b72.txt'%(method_C80, fname),
#                             method=method_C80)

# # C80 MHFB
# PS_loader.importC80_MHFB()

# # C80 HOMs
# PS_loader.importC80_HOMs()

# # C200
# filename = '/rf_cavities/C200/Resonators/from_CST_C200_single_operational.txt'
# # '/rf_cavities/C200/Resonators/from_CST_C200_single_operational.txt'
# # '/rf_cavities/C200/Resonators/from_CST_C200_single_damped.txt'
# PS_loader.importC200(filename=filename)

# # Finemet
# PS_loader.importFinemet()

# # Finemet with MHFB
# PS_loader.importFinemet_MHFB()

# # Kickers
# # filename_Kickers = [
# #             '/kickers/KFA04/Resonators/resonators_broadband_extra_resonances.txt',
# #             '/kickers/KFA13/Resonators/resonators_broadband.txt',
# #             '/kickers/KFA21/Resonators/resonators_broadband.txt',
# #             '/kickers/KFA28/Resonators/resonators_broadband_extra_resonances_best_lowfreq.txt',
# #             '/kickers/KFA45/Resonators/resonators_broadband.txt',
# #             '/kickers/KFA71/Resonators/resonators_broadband.txt',
# #             '/kickers/KFA79/Resonators/resonators_broadband_extra_resonances.txt',
# #             '/kickers/BFA09/Resonators/resonators_broadband_extra_resonances.txt',
# #             '/kickers/BFA21/Resonators/resonators_broadband.txt']
# PS_loader.importKickers()  # filename_list=filename_Kickers

# # Septa
# PS_loader.importSepta()

# # Sector valves
# PS_loader.importSectorValves()

# # MU upstream
# PS_loader.importMU_upstream()

# # MU downstream
# PS_loader.importMU_downstream()

# # Flanges
# PS_loader.importFlanges()

# # Step transitions
# PS_loader.importSteps()

# # Flanges ground loops
# # filename='/vacuum/flanges/ground_loops/Resonators/multi_resonator_1_bypass_low_freq.txt'  # _low_freq
# # # filename='/vacuum/flanges/ground_loops/circuit_model/circuit_model.py'
# # PS_loader.importFlanges_GroundLoops(filename=filename, n_bypasses=200)  # n_bypasses=200
# PS_loader.importFlanges_GroundLoops()

# # Beam instrumentation
# # filename_list = ['/beam_instrumentation/WS/cst_raw_data/PS_SPS_wire_scanner.dat',
# #                  '/beam_instrumentation/Stripline_BPM_SD72/Resonators/multi_resonator.txt',
# #                  '/beam_instrumentation/BGI/Vertical/Resonators/multi_resonator.txt',
# #                  '/beam_instrumentation/WCM_SD03/circuit_model/circuit_model.py']
# # n_elements_list = [4, 1, 1, 2]
# # filename_list = filename_list, n_elements_list = n_elements_list
# PS_loader.importBI()

# # Misc equipment
# PS_loader.importMisc()

# # Resistive wall
# PS_loader.importResistiveWall()

# # Space charge
# PS_loader.importSpaceCharge()

# Using the model (BLonD or imp. toolbox)

imp = PS_loader.export2BLonD()

ResonatorsList = imp.wakeList
ImpedanceTable_list = imp.impedanceList
ImZ_over_f_list = imp.ImZ_over_f_List

# Loading impedance

impedParams = ImpedanceParameters('.', machineParams)

impedParams.importResonatorsList(ResonatorsList)
impedParams.importImpedanceTableList(ImpedanceTable_list)
impedParams.importImZ_over_f_List(ImZ_over_f_list)

impedParams.inducedVoltageGeneration()

# Plotting zone

plt.figure('PS impedance per source')
plt.clf()
plt.plot(impedParams.freqArray,
         np.abs(impedParams.impedance), 'k',
         label='Total', alpha=0.3)
for index_source in range(len(list(PS_loader.table_impedance.keys()))):
    imp.plot_impedance_source(
        list(PS_loader.table_impedance.keys())[index_source],
        figname='PS impedance per source',
        freqArray=impedParams.freqArray)
ax = plt.gca()
ax.xaxis.major.formatter._useMathText = True
plt.ticklabel_format(style='sci', axis='x', scilimits=(-2, 2))
ax.yaxis.major.formatter._useMathText = True
plt.ticklabel_format(style='sci', axis='y', scilimits=(-2, 2))
# plt.savefig(fig_path+'/impedance_per_source.png')

plt.figure('Spectrum and impedance')
plt.clf()
plt.plot(machineParams.freqArray / 1e9, np.abs(machineParams.beamSpectrum))
plt.plot(machineParams.freqArray / 1e9,
         np.abs(machineParams.analyticalSpectrum))
plt.ylim((0, 1.1 * np.max(np.abs(machineParams.beamSpectrum))))
plt.xlabel('Frequency [GHz]')
plt.ylabel('Abs. beam spectrum [$A$]')
plt.twinx()
plt.plot(impedParams.freqArray / 1e9, np.abs(impedParams.impedance) / 1e3, 'r')
plt.ylim((0, 1e-3 * 1.1 * np.max(np.abs(impedParams.impedance))))
plt.ylabel('Abs. impedance [$k\\Omega$]')
plt.tight_layout()
# plt.savefig(fig_path+'/spectrum_and_impedance.png')

plt.figure('Induced Voltage')
plt.clf()
plt.plot(machineParams.timeArray[:len(impedParams.inducedVoltage)] * 1e6,
         machineParams.beamCurrent[:len(impedParams.inducedVoltage)] / np.max(machineParams.beamCurrent) * np.max(np.abs(impedParams.inducedVoltage / 1e3)))
plt.plot(impedParams.timeArray[:len(
    impedParams.inducedVoltage)] * 1e6, impedParams.inducedVoltage / 1e3)
plt.xlim((0, machineParams.timeArray[-1] / n_turns * 2 * 1e6))
plt.xlabel('Time [$\\mu s$]')
plt.ylabel('Induced voltage [kV]')
plt.tight_layout()
# plt.savefig(fig_path+'/induced_voltage.png')
# plt.xlim((0, machineParams.generalParams.t_rev[0]*1e6))
# plt.savefig(fig_path+'/induced_voltage_t_rev.png')


imaginary_Z_over_n = impedParams.impedance.imag / \
    (impedParams.freqArray / machineParams.generalParams.f_rev[0])
imaginary_Z_over_n[0] = imaginary_Z_over_n[1]
plt.figure('ImZ/n')
plt.clf()
plt.plot(machineParams.freqArray / 1e9, machineParams.analyticalSpectrum)
plt.ylim((-1.1 * np.max(np.abs(machineParams.beamSpectrum)),
          1.1 * np.max(np.abs(machineParams.beamSpectrum))))
plt.xlabel('Frequency [GHz]')
plt.ylabel('Single bunch spectrum [$A$]')
plt.twinx()
plt.plot(impedParams.freqArray / 1e9, imaginary_Z_over_n, 'r')
plt.hlines(0, impedParams.freqArray[0] / 1e9,
           impedParams.freqArray[-1] / 1e9, 'k')
plt.ylim((-1.1 * np.max(imaginary_Z_over_n), 1.1 * np.max(imaginary_Z_over_n)))
plt.ylabel('ImZ/n [$\\Omega$]')
plt.tight_layout()
# plt.savefig(fig_path+'/imZ_over_n.png')

plt.show()
