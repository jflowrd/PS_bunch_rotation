'''
Overview of all the kickers impedance

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


freqArray = np.linspace(0, 2e9, 10000)
f_rev = 477e3

# Load KFA04
filename = '../KFA04/Resonators/resonators_broadband_extra_resonances.txt'
kfa04_loader = handleImpedance(folder=script_path)
kfa04_loader.importResonatorFromFile(filename, unitFreq=1, unitRsh=1)

kfa04 = ImpedanceParameters('.')
kfa04.addResonators(kfa04_loader.table_impedance[filename]['Rsh'],
                    kfa04_loader.table_impedance[filename]['fr'],
                    kfa04_loader.table_impedance[filename]['Q'],
                    freqArray=freqArray)

# Plotting the impedance
plt.figure('Impedance')
plt.clf()
plt.plot(kfa04.freqArray/1e6, np.abs(kfa04.impedance),
         label='KFA04')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()

impedance_real_summed = kfa04.impedance.real

plt.figure('Impedance real summed')
plt.clf()
plt.plot(kfa04.freqArray/1e6, impedance_real_summed,
         label='KFA04')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance ReZ [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()

impedance_imag_summed = kfa04.impedance.imag

plt.figure('Impedance imag summed')
plt.clf()
plt.plot(kfa04.freqArray/1e6, impedance_imag_summed/(kfa04.freqArray/f_rev),
         label='KFA04')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance ImZ/n [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()

# Load KFA13
filename = '../KFA13/Resonators/resonators_broadband.txt'
kfa13_loader = handleImpedance(folder=script_path)
kfa13_loader.importResonatorFromFile(filename, unitFreq=1, unitRsh=1)

kfa13 = ImpedanceParameters('.')
kfa13.addResonators(kfa13_loader.table_impedance[filename]['Rsh'],
                    kfa13_loader.table_impedance[filename]['fr'],
                    kfa13_loader.table_impedance[filename]['Q'],
                    freqArray=freqArray)

# Plotting the impedance
plt.figure('Impedance')
plt.plot(kfa13.freqArray/1e6, np.abs(kfa13.impedance),
         label='KFA13')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()

impedance_real_summed += kfa13.impedance.real

plt.figure('Impedance real summed')
plt.plot(kfa13.freqArray/1e6, impedance_real_summed,
         label='KFA13')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance ReZ [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()

impedance_imag_summed += kfa13.impedance.imag

plt.figure('Impedance imag summed')
plt.plot(kfa13.freqArray/1e6, impedance_imag_summed/(kfa13.freqArray/f_rev),
         label='KFA13')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance ImZ/n [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()

# Load KFA21
filename = '../KFA21/Resonators/resonators_broadband.txt'
kfa21_loader = handleImpedance(folder=script_path)
kfa21_loader.importResonatorFromFile(filename, unitFreq=1, unitRsh=1)

kfa21 = ImpedanceParameters('.')
kfa21.addResonators(kfa21_loader.table_impedance[filename]['Rsh'],
                    kfa21_loader.table_impedance[filename]['fr'],
                    kfa21_loader.table_impedance[filename]['Q'],
                    freqArray=freqArray)

# Plotting the impedance
plt.figure('Impedance')
plt.plot(kfa21.freqArray/1e6, np.abs(kfa21.impedance),
         label='KFA21')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()

impedance_real_summed += kfa21.impedance.real

plt.figure('Impedance real summed')
plt.plot(kfa21.freqArray/1e6, impedance_real_summed,
         label='KFA21')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance ReZ [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()

impedance_imag_summed += kfa21.impedance.imag

plt.figure('Impedance imag summed')
plt.plot(kfa21.freqArray/1e6, impedance_imag_summed/(kfa21.freqArray/f_rev),
         label='KFA21')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance ImZ/n [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()

# Load KFA28
filename = '../KFA28/Resonators/resonators_broadband_extra_resonances_best_lowfreq.txt'
kfa28_loader = handleImpedance(folder=script_path)
kfa28_loader.importResonatorFromFile(filename, unitFreq=1, unitRsh=1)

kfa28 = ImpedanceParameters('.')
kfa28.addResonators(kfa28_loader.table_impedance[filename]['Rsh'],
                    kfa28_loader.table_impedance[filename]['fr'],
                    kfa28_loader.table_impedance[filename]['Q'],
                    freqArray=freqArray)

# Plotting the impedance
plt.figure('Impedance')
plt.plot(kfa28.freqArray/1e6, np.abs(kfa28.impedance),
         label='KFA28')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()

impedance_real_summed += kfa28.impedance.real

plt.figure('Impedance real summed')
plt.plot(kfa28.freqArray/1e6, impedance_real_summed,
         label='KFA28')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance ReZ [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()

impedance_imag_summed += kfa28.impedance.imag

plt.figure('Impedance imag summed')
plt.plot(kfa28.freqArray/1e6, impedance_imag_summed/(kfa28.freqArray/f_rev),
         label='KFA28')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance ImZ/n [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()

# Load KFA45
filename = '../KFA45/Resonators/resonators_broadband.txt'
kfa45_loader = handleImpedance(folder=script_path)
kfa45_loader.importResonatorFromFile(filename, unitFreq=1, unitRsh=1)

kfa45 = ImpedanceParameters('.')
kfa45.addResonators(kfa45_loader.table_impedance[filename]['Rsh'],
                    kfa45_loader.table_impedance[filename]['fr'],
                    kfa45_loader.table_impedance[filename]['Q'],
                    freqArray=freqArray)

# Plotting the impedance
plt.figure('Impedance')
plt.plot(kfa45.freqArray/1e6, np.abs(kfa45.impedance),
         label='KFA45')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()

impedance_real_summed += kfa45.impedance.real

plt.figure('Impedance real summed')
plt.plot(kfa45.freqArray/1e6, impedance_real_summed,
         label='KFA45')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance ReZ [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()

impedance_imag_summed += kfa45.impedance.imag

plt.figure('Impedance imag summed')
plt.plot(kfa45.freqArray/1e6, impedance_imag_summed/(kfa45.freqArray/f_rev),
         label='KFA45')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance ImZ/n [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()

# Load KFA71
filename = '../KFA71/Resonators/resonators_broadband.txt'
kfa71_loader = handleImpedance(folder=script_path)
kfa71_loader.importResonatorFromFile(filename, unitFreq=1, unitRsh=1)

kfa71 = ImpedanceParameters('.')
kfa71.addResonators(kfa71_loader.table_impedance[filename]['Rsh'],
                    kfa71_loader.table_impedance[filename]['fr'],
                    kfa71_loader.table_impedance[filename]['Q'],
                    freqArray=freqArray)

# Plotting the impedance
plt.figure('Impedance')
plt.plot(kfa71.freqArray/1e6, np.abs(kfa71.impedance),
         label='KFA71')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()

impedance_real_summed += kfa71.impedance.real

plt.figure('Impedance real summed')
plt.plot(kfa71.freqArray/1e6, impedance_real_summed,
         label='KFA71')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance ReZ [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()

impedance_imag_summed += kfa71.impedance.imag

plt.figure('Impedance imag summed')
plt.plot(kfa71.freqArray/1e6, impedance_imag_summed/(kfa71.freqArray/f_rev),
         label='KFA71')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance ImZ/n [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()

# Load KFA79
filename = '../KFA79/Resonators/resonators_broadband_extra_resonances.txt'
kfa79_loader = handleImpedance(folder=script_path)
kfa79_loader.importResonatorFromFile(filename, unitFreq=1, unitRsh=1)

kfa79 = ImpedanceParameters('.')
kfa79.addResonators(kfa79_loader.table_impedance[filename]['Rsh'],
                    kfa79_loader.table_impedance[filename]['fr'],
                    kfa79_loader.table_impedance[filename]['Q'],
                    freqArray=freqArray)

# Plotting the impedance
plt.figure('Impedance')
plt.plot(kfa79.freqArray/1e6, np.abs(kfa79.impedance),
         label='KFA79')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()

impedance_real_summed += kfa79.impedance.real

plt.figure('Impedance real summed')
plt.plot(kfa79.freqArray/1e6, impedance_real_summed,
         label='KFA79')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance ReZ [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()

impedance_imag_summed += kfa79.impedance.imag

plt.figure('Impedance imag summed')
plt.plot(kfa79.freqArray/1e6, impedance_imag_summed/(kfa79.freqArray/f_rev),
         label='KFA79')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance ImZ/n [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()

# Load BFA09
filename = '../BFA09/Resonators/resonators_broadband_extra_resonances.txt'
bfa09_loader = handleImpedance(folder=script_path)
bfa09_loader.importResonatorFromFile(filename, unitFreq=1, unitRsh=1)

bfa09 = ImpedanceParameters('.')
bfa09.addResonators(bfa09_loader.table_impedance[filename]['Rsh'],
                    bfa09_loader.table_impedance[filename]['fr'],
                    bfa09_loader.table_impedance[filename]['Q'],
                    freqArray=freqArray)

# Plotting the impedance
plt.figure('Impedance')
plt.plot(kfa79.freqArray/1e6, np.abs(bfa09.impedance),
         label='BFA09')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()

impedance_real_summed += bfa09.impedance.real

plt.figure('Impedance real summed')
plt.plot(bfa09.freqArray/1e6, impedance_real_summed,
         label='BFA09')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance ReZ [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()

impedance_imag_summed += bfa09.impedance.imag

plt.figure('Impedance imag summed')
plt.plot(bfa09.freqArray/1e6, impedance_imag_summed/(bfa09.freqArray/f_rev),
         label='BFA09')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance ImZ/n [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()

# Load BFA21
filename = '../BFA21/Resonators/resonators_broadband.txt'
bfa21_loader = handleImpedance(folder=script_path)
bfa21_loader.importResonatorFromFile(filename, unitFreq=1, unitRsh=1)

bfa21 = ImpedanceParameters('.')
bfa21.addResonators(bfa21_loader.table_impedance[filename]['Rsh'],
                    bfa21_loader.table_impedance[filename]['fr'],
                    bfa21_loader.table_impedance[filename]['Q'],
                    freqArray=freqArray)

# Plotting the impedance
plt.figure('Impedance')
plt.plot(kfa79.freqArray/1e6, np.abs(bfa21.impedance),
         '--', label='BFA21')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig(script_path+'/overview_indiv_kickers.png')
plt.xlim((0, 300))
plt.ylim((0, 1000))
plt.tight_layout()
plt.savefig(script_path+'/overview_indiv_kickers_zoom.png')

impedance_real_summed += bfa21.impedance.real

plt.figure('Impedance real summed')
plt.plot(bfa21.freqArray/1e6, impedance_real_summed,
         '--', label='BFA21')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance ReZ [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig(script_path+'/overview_real_summed.png')
plt.xlim((0, 300))
plt.ylim((0, 2000))
plt.tight_layout()
plt.savefig(script_path+'/overview_real_summed_zoom.png')

impedance_imag_summed += bfa21.impedance.imag

plt.figure('Impedance imag summed')
plt.plot(bfa21.freqArray/1e6, impedance_imag_summed/(bfa21.freqArray/f_rev),
         '--', label='BFA21')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance ImZ/n [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig(script_path+'/overview_imag_summed.png')
plt.xlim((0, 300))
plt.ylim((-5, 35))
plt.tight_layout()
plt.savefig(script_path+'/overview_imag_summed_zoom.png')

plt.show()
