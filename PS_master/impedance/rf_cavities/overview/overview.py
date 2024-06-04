'''
Overview of all the cavities impedance

@author: alasheen
'''

# General import
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# Impedance tools imports
script_path = os.path.dirname(os.path.realpath(__file__))
toolbox_path = script_path+'/../../../impedance_toolbox'
sys.path.insert(0, toolbox_path)
from impedance_toolbox.impedance_params import ImpedanceParameters
from impedance_toolbox.handle_impedance import handleImpedance


freqArray = np.linspace(0, 1.2e9, 10000)
f_rev = 477e3

# Load C10_h7
n_C10_h7 = 10
C10_h7_harmonic = 7

# script_path = os.path.dirname(os.path.realpath(__file__))
# toolbox_path = script_path+'/../../../impedance_toolbox'
sys.path.insert(0, script_path+'/../')
from C10.All.Parametric.All_10MHz_Parametric_model import load_model

C10_h7 = ImpedanceParameters('.')

loaded_data = load_model(C10_h7_harmonic * f_rev)
C10_h7.addResonators(loaded_data[1] * n_C10_h7,
                  loaded_data[0],
                  loaded_data[2],
                  freqArray=freqArray)

# Plotting the impedance
plt.figure('Impedance', figsize=(8,6))
plt.clf()
plt.plot(C10_h7.freqArray/1e6, np.abs(C10_h7.impedance),
         label='C10_h7')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()


# Load C10_h14
n_C10_h14 = 10
C10_h14_harmonic = 14

# script_path = os.path.dirname(os.path.realpath(__file__))
# toolbox_path = script_path+'/../../../impedance_toolbox'
sys.path.insert(0, script_path+'/../')
from C10.All.Parametric.All_10MHz_Parametric_model import load_model

C10_h14 = ImpedanceParameters('.')

loaded_data = load_model(C10_h14_harmonic * f_rev)
C10_h14.addResonators(loaded_data[1] * n_C10_h14,
                  loaded_data[0],
                  loaded_data[2],
                  freqArray=freqArray)

# Plotting the impedance
plt.figure('Impedance')
# plt.clf()
plt.plot(C10_h14.freqArray/1e6, np.abs(C10_h14.impedance),
         label='C10_h14')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()


# Load C10_h21
n_C10_h21 = 10
C10_h21_harmonic = 21

# script_path = os.path.dirname(os.path.realpath(__file__))
# toolbox_path = script_path+'/../../../impedance_toolbox'
sys.path.insert(0, script_path+'/../')
from C10.All.Parametric.All_10MHz_Parametric_model import load_model

C10_h21 = ImpedanceParameters('.')

loaded_data = load_model(C10_h21_harmonic * f_rev)
C10_h21.addResonators(loaded_data[1] * n_C10_h21,
                  loaded_data[0],
                  loaded_data[2],
                  freqArray=freqArray)

# Plotting the impedance
plt.figure('Impedance')
# plt.clf()
plt.plot(C10_h21.freqArray/1e6, np.abs(C10_h21.impedance),
         label='C10_h21')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()


# Load C10_h21_gap_closed
n_C10_h21_gap_closed = 10
C10_h21_gap_closed_harmonic = 21

# script_path = os.path.dirname(os.path.realpath(__file__))
# toolbox_path = script_path+'/../../../impedance_toolbox'
sys.path.insert(0, script_path+'/../')
from C10.ClosedGapRelay.Parametric.All_10MHz_ClosedGapRelay_Parametric_model import load_model

C10_h21_gap_closed = ImpedanceParameters('.')

loaded_data = load_model(C10_h21_gap_closed_harmonic)
C10_h21_gap_closed.addResonators(loaded_data[1] * n_C10_h21_gap_closed,
                  loaded_data[0],
                  loaded_data[2],
                  freqArray=freqArray)

C10_h21_gap_closed.impedance.imag += loaded_data[3] * freqArray * n_C10_h21_gap_closed


# Plotting the impedance
plt.figure('Impedance')
# plt.clf()
plt.plot(C10_h21_gap_closed.freqArray/1e6, np.abs(C10_h21_gap_closed.impedance),
         label='C10_h21_gap_closed')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()


# Load C20
n_C20 = 2
C20_harmonic = 42

filename = '../C20/All/Resonators/single_resonator_20MHz_20kV.txt'
loaded_impedance = np.loadtxt(script_path+'/'+filename)
# im_z_over_f = loaded_impedance[-1, -1] * loaded_impedance.shape[0]

C20_loader = handleImpedance(folder=script_path)
C20_loader.importResonatorFromFile(filename, unitFreq=1, unitRsh=1)

C20 = ImpedanceParameters('.')
C20.addResonators(C20_loader.table_impedance[filename]['Rsh'] * n_C20,
                    C20_loader.table_impedance[filename]['fr'],
                    C20_loader.table_impedance[filename]['Q'],
                    freqArray=freqArray)

# C20.impedance.imag += im_z_over_f * freqArray * n_C20


# Plotting the impedance
plt.figure('Impedance')
# plt.clf()
plt.plot(C20.freqArray/1e6, np.abs(C20.impedance),
         label='C20')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()


# Load C40
C40_harmonic = 84

C40 = ImpedanceParameters('.')

filename = '../C40/C40-77/Resonators/single_resonator_b72.txt'
filename_2 = '../C40/C40-78/Resonators/single_resonator_b72.txt'

loaded_impedance = np.loadtxt(script_path+'/'+filename)
loaded_impedance_2 = np.loadtxt(script_path+'/'+filename_2)
# im_z_over_f = loaded_impedance[-1, -1] * loaded_impedance.shape[0]

C40_loader = handleImpedance(folder=script_path)
C40_loader.importResonatorFromFile(filename, unitFreq=1, unitRsh=1)
C40_loader.importResonatorFromFile(filename_2, unitFreq=1, unitRsh=1)

C40.addResonators(C40_loader.table_impedance[filename]['Rsh'],
                    C40_loader.table_impedance[filename]['fr'],
                    C40_loader.table_impedance[filename]['Q'],
                    freqArray=freqArray)
C40.addResonators(C40_loader.table_impedance[filename_2]['Rsh'],
                    C40_loader.table_impedance[filename_2]['fr'],
                    C40_loader.table_impedance[filename_2]['Q'],
                    freqArray=freqArray)

# C40.impedance.imag += im_z_over_f * freqArray * n_C40


# Plotting the impedance
plt.figure('Impedance')
# plt.clf()
plt.plot(C40.freqArray/1e6, np.abs(C40.impedance),
         label='C40')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()


# Load C80
C80_harmonic = 168

C80 = ImpedanceParameters('.')

filename = '../C80/C80-08/Resonators/single_resonator_b72.txt'
filename_2 = '../C80/C80-88/Resonators/single_resonator_b72.txt'
filename_3 = '../C80/C80-89/Resonators/single_resonator_b72.txt'

loaded_impedance = np.loadtxt(script_path+'/'+filename)
loaded_impedance_2 = np.loadtxt(script_path+'/'+filename_2)
loaded_impedance_3 = np.loadtxt(script_path+'/'+filename_3)
# im_z_over_f = loaded_impedance[-1, -1] * loaded_impedance.shape[0]

C80_loader = handleImpedance(folder=script_path)
C80_loader.importResonatorFromFile(filename, unitFreq=1, unitRsh=1)
C80_loader.importResonatorFromFile(filename_2, unitFreq=1, unitRsh=1)
C80_loader.importResonatorFromFile(filename_3, unitFreq=1, unitRsh=1)

C80.addResonators(C80_loader.table_impedance[filename]['Rsh'],
                    C80_loader.table_impedance[filename]['fr'],
                    C80_loader.table_impedance[filename]['Q'],
                    freqArray=freqArray)
C80.addResonators(C80_loader.table_impedance[filename_2]['Rsh'],
                    C80_loader.table_impedance[filename_2]['fr'],
                    C80_loader.table_impedance[filename_2]['Q'],
                    freqArray=freqArray)
C80.addResonators(C80_loader.table_impedance[filename_3]['Rsh'],
                    C80_loader.table_impedance[filename_3]['fr'],
                    C80_loader.table_impedance[filename_3]['Q'],
                    freqArray=freqArray)

# C80.impedance.imag += im_z_over_f * freqArray * n_C80


# Plotting the impedance
plt.figure('Impedance')
# plt.clf()
plt.plot(C80.freqArray/1e6, np.abs(C80.impedance),
         label='C80')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()


# Load C200_on
n_C200_on = 6
C200_on_harmonic = 420

filename = '../C200/All/Resonators/from_CST_200MHz_single_operational.txt'
loaded_impedance = np.loadtxt(script_path+'/'+filename)
# im_z_over_f = loaded_impedance[-1, -1] * loaded_impedance.shape[0]

C200_on_loader = handleImpedance(folder=script_path)
C200_on_loader.importResonatorFromFile(filename, unitFreq=1, unitRsh=1)

C200_on = ImpedanceParameters('.')
C200_on.addResonators(C200_on_loader.table_impedance[filename]['Rsh'] * n_C200_on,
                    C200_on_loader.table_impedance[filename]['fr'],
                    C200_on_loader.table_impedance[filename]['Q'],
                    freqArray=freqArray)

# C200_on.impedance.imag += im_z_over_f * freqArray * n_C200_on


# Plotting the impedance
plt.figure('Impedance')
# plt.clf()
plt.plot(C200_on.freqArray/1e6, np.abs(C200_on.impedance),
         label='C200_on')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()


# Load C200_off
n_C200_off = 6
C200_off_harmonic = 420

filename = '../C200/All/Resonators/from_CST_200MHz_single_damped.txt'
loaded_impedance = np.loadtxt(script_path+'/'+filename)
# im_z_over_f = loaded_impedance[-1, -1] * loaded_impedance.shape[0]

C200_off_loader = handleImpedance(folder=script_path)
C200_off_loader.importResonatorFromFile(filename, unitFreq=1, unitRsh=1)

C200_off = ImpedanceParameters('.')
C200_off.addResonators(C200_off_loader.table_impedance[filename]['Rsh'] * n_C200_off,
                    C200_off_loader.table_impedance[filename]['fr'],
                    C200_off_loader.table_impedance[filename]['Q'],
                    freqArray=freqArray)

# C200_off.impedance.imag += im_z_over_f * freqArray * n_C200_off


# Plotting the impedance
plt.figure('Impedance')
# plt.clf()
plt.plot(C200_off.freqArray/1e6, np.abs(C200_off.impedance),
         label='C200_off')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()



# Load Finemet
n_Finemet = 1

filename = '../Finemet/Resonators_2018/multi_resonators_impedance.txt'
loaded_impedance = np.loadtxt(script_path+'/'+filename)
# im_z_over_f = loaded_impedance[-1, -1] * loaded_impedance.shape[0]

Finemet_loader = handleImpedance(folder=script_path)
Finemet_loader.importResonatorFromFile(filename, unitFreq=1, unitRsh=1)

Finemet = ImpedanceParameters('.')
Finemet.addResonators(Finemet_loader.table_impedance[filename]['Rsh'] * n_Finemet,
                    Finemet_loader.table_impedance[filename]['fr'],
                    Finemet_loader.table_impedance[filename]['Q'],
                    freqArray=freqArray)

# Finemet.impedance.imag += im_z_over_f * freqArray * n_Finemet


# Plotting the impedance
plt.figure('Impedance')
# plt.clf()
plt.plot(Finemet.freqArray/1e6, np.abs(Finemet.impedance),
         label='Finemet')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()
plt.xlim((0.1, 400))
plt.ylim((10, 1e6))
plt.xscale('log')
plt.yscale('log')
plt.tight_layout()
plt.savefig('overview_cavities.png', dpi=300)



# Load C40_HOM
n_C40_HOM = 2

filename = '../C40/All/Resonators/HOMs_C40.txt'
loaded_impedance = np.loadtxt(script_path+'/'+filename)
# im_z_over_f = loaded_impedance[-1, -1] * loaded_impedance.shape[0]

C40_HOM_loader = handleImpedance(folder=script_path)
C40_HOM_loader.importResonatorFromFile(filename, unitFreq=1, unitRsh=1)

C40_HOM = ImpedanceParameters('.')
C40_HOM.addResonators(C40_HOM_loader.table_impedance[filename]['Rsh'] * n_C40_HOM,
                    C40_HOM_loader.table_impedance[filename]['fr'],
                    C40_HOM_loader.table_impedance[filename]['Q'],
                    freqArray=freqArray)

# C40_HOM.impedance.imag += im_z_over_f * freqArray * n_C40_HOM


# Plotting the impedance
plt.figure('Impedance HOM')
plt.clf()
plt.plot(C40_HOM.freqArray/1e6, np.abs(C40_HOM.impedance),
         label='C40_HOM')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()


# Load C80_HOM
n_C80_HOM = 3

filename = '../C80/All/Resonators/HOMs_C80.txt'
loaded_impedance = np.loadtxt(script_path+'/'+filename)
# im_z_over_f = loaded_impedance[-1, -1] * loaded_impedance.shape[0]

C80_HOM_loader = handleImpedance(folder=script_path)
C80_HOM_loader.importResonatorFromFile(filename, unitFreq=1, unitRsh=1)

C80_HOM = ImpedanceParameters('.')
C80_HOM.addResonators(C80_HOM_loader.table_impedance[filename]['Rsh'] * n_C80_HOM,
                    C80_HOM_loader.table_impedance[filename]['fr'],
                    C80_HOM_loader.table_impedance[filename]['Q'],
                    freqArray=freqArray)

# C80_HOM.impedance.imag += im_z_over_f * freqArray * n_C80_HOM


# Plotting the impedance
plt.figure('Impedance HOM')
# plt.clf()
plt.plot(C80_HOM.freqArray/1e6, np.abs(C80_HOM.impedance),
         label='C80_HOM')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('overview_cavities_HOMs.png', dpi=300)


plt.show()
