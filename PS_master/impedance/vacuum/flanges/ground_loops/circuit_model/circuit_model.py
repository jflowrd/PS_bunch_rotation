'''
Model of the ground loops on the vacuum flanges with/without the
RF bypasses
'''

import numpy as np
import matplotlib.pyplot as plt
import os


script_path = os.path.dirname(os.path.realpath(__file__))

R = 100  # Ohm
C = 1000e-12  # Farad
L = 10e-6  # Henri

R1 = 1  # Ohm
C1 = 0.4e-6  # Farad
L1 = 12e-9  # Henri


def load_model(freq_array, n_bypasses=1,
               R=R, L=L, C=C, R1=R1, L1=L1, C1=C1):

    if freq_array is None:
        freq_array = np.linspace(0, 500e6, 1000000)

    omega_array = freq_array * 2*np.pi

    impedance_bypass = R1 + 1j*omega_array*L1 + 1/(1j*omega_array*C1)

    impedance = 1/(1/R + 1/(1j*omega_array*L) + 1j*omega_array*C
                   + n_bypasses/impedance_bypass)

    impedance[freq_array == 0] = 0+0*1j

    return impedance, freq_array


if __name__ == '__main__':

    f_rev = 477e3

    freq_array = np.linspace(0, 100e6, 100000)
    omega_array = freq_array * 2*np.pi

    n = freq_array/f_rev

    impedance_loop = 1/(1/R + 1/(1j*omega_array*L) + 1j*omega_array*C)
    impedance_loop[freq_array == 0] = 0+0*1j

    plt.figure('Impedance')
    plt.clf()
    plt.plot(freq_array/1e6, np.abs(impedance_loop), label='Ground loop')

    plt.figure('Impedance Z/n')
    plt.clf()
    plt.plot(freq_array/1e6, np.abs(impedance_loop)/n, label='Ground loop')

    impedance_bypass = R1 + 1j*omega_array*L1 + 1/(1j*omega_array*C1)

#     plt.plot(freq_array, np.abs(impedance_bypass))

    impedance_loop_bypass = 1/(1/R + 1/(1j*omega_array*L) + 1j*omega_array*C
                               + 1/impedance_bypass)
    impedance_loop_bypass[freq_array == 0] = 0+0*1j

    plt.figure('Impedance')
    plt.plot(freq_array/1e6, np.abs(impedance_loop_bypass),
             label='1 RF bypass')
    plt.figure('Impedance Z/n')
    plt.plot(freq_array/1e6, np.abs(impedance_loop_bypass)/n,
             label='1 RF bypass')

    impedance_loop_doublebypass = 1/(1/R + 1/(1j*omega_array*L)
                                     + 1j*omega_array*C
                                     + 1/impedance_bypass
                                     + 1/impedance_bypass)
    impedance_loop_doublebypass[freq_array == 0] = 0+0*1j

    plt.figure('Impedance')
    plt.plot(freq_array/1e6, np.abs(impedance_loop_doublebypass),
             label='2 RF bypasses')
    plt.figure('Impedance Z/n')
    plt.plot(freq_array/1e6, np.abs(impedance_loop_doublebypass)/n,
             label='2 RF bypasses')

    plt.figure('Impedance')
    plt.xlabel('Frequency [MHz]')
    plt.ylabel('Abs Impedance [Ohm]')
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig(script_path+'/abs_impedance.png')
    plt.xlim((0, 5))
    plt.tight_layout()
    plt.savefig(script_path+'/abs_impedance_zoom.png')

    plt.figure('Impedance Z/n')
    plt.xlabel('Frequency [MHz]')
    plt.ylabel('Abs Impedance Z/n [Ohm]')
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig(script_path+'/abs_impedance_z_over_n.png')
    plt.xlim((0, 5))
    plt.tight_layout()
    plt.savefig(script_path+'/abs_impedance_z_over_n_zoom.png')

    plt.show()
