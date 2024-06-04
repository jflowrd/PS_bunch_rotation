'''
Fitting the measured C80-88 cavity impedance with resonators + direct
feedback model

@author: alasheen
'''

# General import
import os
import sys
import numpy as np

# Toolbox import
sys.path.append('../../../impedance_toolbox')
from impedance_toolbox.impedance_toolbox.miniblond.impedances.impedance_sources import Resonators


def feedback(freq, impedance, delay, fb_gain):
    filter_feedback = fb_gain * np.exp(-1j * 2 * np.pi * freq * delay)
    return impedance / (1 + impedance * filter_feedback)


def impedance_calc_FB(freqArray, f_rf, rshunt, Q, delay_n_rf, loop_gain_dB):

    OL_Resonator = Resonators(rshunt, f_rf, Q)
    OL_Resonator.imped_calc(freqArray)

    fb_gain = 10 ** (loop_gain_dB / 20) / rshunt
    delay = delay_n_rf / f_rf

#     fb_gain_max = np.pi/2 * 1/(rshunt/Q) * 1/(2*np.pi*f_rf*delay)

    impedance_FB = feedback(freqArray, OL_Resonator.impedance, delay, fb_gain)

    return impedance_FB


def residue(params, *args):

    (roq, Q, delay_n_rf, loop_gain_dB, f_rf) = params

    roq = np.abs(roq)
    Q = np.abs(Q * 1e4)
    delay_n_rf = np.abs(np.round(delay_n_rf))
    loop_gain_dB = np.abs(loop_gain_dB)
    f_rf = np.abs(f_rf*1e6)

    (freqArray, meas_impedance) = args

    r_shunt = roq * Q

    impedance_FB = impedance_calc_FB(
        freqArray, f_rf, r_shunt, Q, delay_n_rf, loop_gain_dB)

    res = np.sqrt(np.sum((np.abs(meas_impedance) - np.abs(impedance_FB)) ** 2))

    return res


def load_model(filename, folder='./', freqArray=None):

    if freqArray is None:
        freqArray = np.linspace(0, 240e6, 100000)

    loaded_input = np.loadtxt(folder + filename)

    return (impedance_calc_FB(freqArray,
                              loaded_input[0],
                              loaded_input[1],
                              loaded_input[2],
                              loaded_input[3],
                              loaded_input[4]),
            freqArray)


if __name__ == '__main__':

    # General imports

    import matplotlib.pyplot as plt
    # from scipy.fftpack import hilbert
    from scipy.constants import c, e, m_p
    from scipy.optimize import minimize

    # Impedance tools imports
    script_path = os.path.dirname(os.path.realpath(__file__))
    toolbox_path = script_path + '/../../../../../../impedance_toolbox'
    sys.path.insert(0, toolbox_path)
    from impedance_toolbox.impedance_params import ImpedanceParameters
    from impedance_toolbox.miniblond.impedances.impedance_sources import Resonators

    # Loading input impedance and defining the limit in frequency for the
    # multi resonators fit
    n_bunches = np.arange(12, 80, 12)
    # n_bunches = np.array([72])

    cavity_name = 'C80-88'
    data_path = script_path + '/../Measurements/'

    charge = 1  # e
    mass = m_p*c**2/e  # eV
    BField = 1.255e4  # G
    BField = BField*(
                1
                - 1.441767e-6*BField
                + 2.947290e-10*BField**2.
                - 2.357026e-14*BField**3.)
    bending_radius = 70.079  # m
    momentum = BField * bending_radius * charge * c * 1e-4

    t_rev = (2 * np.pi * 100) / (c * momentum /
                                 np.sqrt(momentum ** 2. + mass ** 2.))
    f_rev = 1 / t_rev 
    
    all_cases = np.zeros((len(n_bunches), 5))

    for index_bunch in range(len(n_bunches)):

        loaded_data = np.loadtxt(
            data_path + '/' + cavity_name + '_b' + str(n_bunches[index_bunch]) + '.txt')
        freq_data = loaded_data[:, 0]
        abs_Z_data = loaded_data[:, 1]

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
        plt.plot(impedParams.freqArray / 1e6,
                 np.abs(impedParams.impedance), label='Input')
        plt.xlabel('Frequency [MHz]')
        plt.ylabel('Impedance [$\\Omega$]')
        plt.legend(loc='best')
        plt.tight_layout()

        # Getting a first approximation of the resonances of the impedance
        # table
        f_start = 168*f_rev  # Measured with Xenon beam
        Q_start = 11800
        roq_start = 56
        r_start = roq_start*Q_start
    
        loop_gain_dB = 41
        delay_n_rf = 18

        impedance_FB = impedance_calc_FB(
            impedParams.freqArray, f_start, r_start, Q_start, delay_n_rf, loop_gain_dB)

        x0 = [roq_start, Q_start / 1e4, delay_n_rf, loop_gain_dB, f_start/1e6]

        fitted_params = minimize(residue, x0,
                                 args=(impedParams.freqArray,
                                       impedParams.impedance),
                                 method='Powell', options={'maxiter': 10000, 'disp': True})

        fitted_params = fitted_params['x']

        fitted_params[1] = np.abs(fitted_params[1] * 1e4)
        fitted_params[2] = np.abs(np.round(fitted_params[2]))
        fitted_params[3] = np.abs(fitted_params[3])
        fitted_params[4] = np.abs(fitted_params[4]*1e6)
        print(fitted_params)

        impedance_FB = impedance_calc_FB(impedParams.freqArray,
                                         fitted_params[4],
                                         fitted_params[0] * fitted_params[1],
                                         fitted_params[1],
                                         fitted_params[2],
                                         fitted_params[3])

        # Plotting the result with one resonator
        plt.figure('Impedance')
        plt.plot(impedParams.freqArray / 1e6,
                 np.abs(impedance_FB),
                 label='1 Resonator + DFB')
        plt.xlabel('Frequency [MHz]')
        plt.ylabel('Impedance [$\\Omega$]')
        plt.legend(loc='best')
        plt.tight_layout()
        plt.savefig(script_path + '/fitted_b' +
                    str(n_bunches[index_bunch]) + '.png')

        # Saving result into text file

        R_S_save = np.abs(
            np.array(fitted_params[0] * fitted_params[1], ndmin=2))
        fr_save = np.array(fitted_params[4], ndmin=2)
        Q_save = np.array(fitted_params[1], ndmin=2)
        delay_n_rf_save = np.array(fitted_params[2], ndmin=2)
        gain_save = np.array(fitted_params[3], ndmin=2)

        saved_matrix = np.hstack(
            (fr_save.T, R_S_save.T, Q_save.T, delay_n_rf_save.T, gain_save.T))

        all_cases[index_bunch, :] = np.array(saved_matrix)

        np.savetxt(script_path + '/single_resonator_DFB_b' + str(n_bunches[index_bunch]) + '.txt',
                   saved_matrix,
                   header='Impedance of %s\n' % (cavity_name) +
                   'Single resonator fit for %d bunches\n' % (n_bunches[index_bunch]) +
                   'Impedance damped by Direct Feedback\n' +
                   'Author: A. Lasheen\n' +
                   'f_r [Hz]\tR_s [Ohm]\tQ\tN period delay\tLoop gain [dB]')

        plt.show()

    # Plotting the result for all measurements
    plt.figure('Final R')
    plt.clf()
    plt.plot(n_bunches,
             all_cases[:, 1] / 1e3, 'o--',
             label='Fit')
    plt.axhline(1260, label='Unloaded', color='orange')
    plt.axhline(660, label='Loaded', color='red')
    plt.xlabel('n bunches')
    plt.ylabel('Open loop shunt impedance [kOhm]')
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig(script_path + '/final_rshunt.png')

    plt.figure('Final Q')
    plt.clf()
    plt.plot(n_bunches,
             all_cases[:, 2], 'o--',
             label='Fit')
    plt.axhline(22600, label='Unloaded', color='orange')
    plt.axhline(11800, label='Loaded', color='red')
    plt.xlabel('n bunches')
    plt.ylabel('Open loop Q [kOhm]')
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig(script_path + '/final_Q.png')

    plt.figure('Final RoQ')
    plt.clf()
    plt.plot(n_bunches,
             all_cases[:, 1] / all_cases[:, 2], 'o--',
             label='Fit')
    plt.axhline(56, label='Expected', color='orange')
    plt.xlabel('n bunches')
    plt.ylabel('Open loop RoQ [Ohm]')
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig(script_path + '/final_RoQ.png')

    plt.figure('Final Delay')
    plt.clf()
    plt.plot(n_bunches,
             all_cases[:, 3] / f_start * 1e9, 'o--',
             label='Fit')
    plt.axhline(220, label='Expected', color='orange')
    plt.xlabel('n bunches')
    plt.ylabel('Delay [ns]')
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig(script_path + '/final_delay.png')

    plt.figure('Final Gain')
    plt.clf()
    plt.plot(n_bunches,
             all_cases[:, 4], 'o--',
             label='Fit')
    plt.axhline(41, label='Expected', color='orange')
    plt.xlabel('n bunches')
    plt.ylabel('Loop gain [dB]')
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig(script_path + '/final_gain.png')

    plt.show()
