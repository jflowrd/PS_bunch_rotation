'''
Model of the ground loops on the vacuum flanges with/without the
RF bypasses
'''

import numpy as np
import matplotlib.pyplot as plt
import os


script_path = os.path.dirname(os.path.realpath(__file__))

R = 6  # Ohm
C = 0.65e-11  # Farad
L = 10e-6  # Henri


def load_model(freq_array, R=R, L=L, C=C):

    if freq_array is None:
        freq_array = np.linspace(0, 6e9, 100000)

    omega_array = freq_array * 2*np.pi

    impedance = 1/(1/R + 1/(1j*omega_array*L) + 1j*omega_array*C)

    impedance[freq_array == 0] = 0+0*1j

    return impedance, freq_array


def impedance_to_dB(impedance, reference=None, ratio=False):

    if reference is None:
        return 20*np.log10(np.abs(impedance)/np.nanmax(np.abs(impedance)))
    elif isinstance(reference, float) or (
            isinstance(reference, np.ndarray) and not ratio):
        return 20*np.log10(np.abs(impedance)/np.nanmax(np.abs(reference)))
    elif isinstance(reference, float) or (
            isinstance(reference, np.ndarray) and ratio):
        return 20*np.log10(np.abs(impedance)/np.abs(reference))
    else:
        return None


if __name__ == '__main__':

    f_rev = 477e3

    freq_array = np.linspace(0, 6e9, 100000)
    omega_array = freq_array * 2*np.pi

    n = freq_array/f_rev

    impedance_loop = 1/(1/R + 1/(1j*omega_array*L) + 1j*omega_array*C)
#     impedance_loop[freq_array == 0] = 0+0*1j

    plt.figure('Impedance')
    plt.clf()
    plt.plot(freq_array, impedance_to_dB(impedance_loop),
             label='WCM, R=%.1f $\\Omega$'%(R))
    plt.axhline(-3, alpha=0.5)
    plt.xlabel('Frequency [MHz]')
    plt.ylabel('Abs Impedance [dB]')
    plt.legend(loc='best')
    plt.xscale('log')
    plt.tight_layout()
    plt.savefig(script_path+'/abs_impedance.png')

    plt.show()
