'''
Created on 24 juin 2019

@author: alasheen
'''

# General imports
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

# Toolbox import
sys.path.append('../../../impedance_toolbox')
from impedance_toolbox.impedance_toolbox.miniblond.impedances.impedance_sources import Resonators

# Main C10 model import
sys.path.append('../../All/Parametric')
# from All_10MHz_Parametric_model import load_model as main_C10
from impedance.rf_cavities.C10.All.Parametric.All_C10_Parametric_model import load_model as main_C10

# Script path
script_path = os.path.dirname(os.path.realpath(__file__))

# Loading LUT
loaded_filter_data = np.load(script_path+'/1TFB_LUT_gains_and_delays.npz')


def load_model(main_harmonic, revolution_frequency,
               impedance_reduction_target, main_harmonic_FB=False,
               frequency_array=None, manual_params=None,
               ZFactor=1.):
    '''
    - Method to load the parametric model
    Requires the main harmonic, the revolution frequency, the impedance
    reduction target on sidebands +-1 f_rev of the central harmonic.
    The feedback can be applied or removed using main_harmonic_FB=True/False
    '''

    if frequency_array is None:
        frequency_array = np.linspace(0, 500e6, int(1e6))

    loaded_main_C10 = main_C10(main_harmonic*revolution_frequency)
    C10_Resonator = Resonators(loaded_main_C10[1],
                               loaded_main_C10[0],
                               loaded_main_C10[2])

    C10_Resonator.imped_calc(frequency_array)

    if impedance_reduction_target < 0:
        if manual_params is None:
            index_harmonic = np.interp(
                main_harmonic,
                loaded_filter_data['main_harmonic_range'],
                np.arange(len(loaded_filter_data['main_harmonic_range'])))
            index_frev = np.interp(
                revolution_frequency,
                loaded_filter_data['revolution_frequency_range'],
                np.arange(len(loaded_filter_data['revolution_frequency_range'])))
            index_red = np.interp(
                -impedance_reduction_target,
                -loaded_filter_data['impedance_reduction_target_range'],
                np.arange(len(loaded_filter_data['impedance_reduction_target_range'])))

            # Bi-linear interpolation to get gain_comb, gain_res and Q_res
            gain_comb = loaded_filter_data['gain_comb_save'][
                int(index_harmonic), int(index_frev), int(index_red)]
            gain_comb += (index_frev-int(index_frev))*(
                loaded_filter_data['gain_comb_save'][
                    int(index_harmonic), int(index_frev)+1, int(index_red)] -
                gain_comb)
            gain_comb += (index_red-int(index_red))*(
                (index_frev-int(index_frev))*(
                    loaded_filter_data['gain_comb_save'][
                        int(index_harmonic), int(index_frev)+1, int(index_red)+1] -
                    loaded_filter_data['gain_comb_save'][
                        int(index_harmonic), int(index_frev), int(index_red)+1]) -
                gain_comb)

            gain_comb *= 10  # The LUT was computed for 10 cavities, should be recomputed!

            if main_harmonic_FB:
                gain_res = 0
                Q_res = 500
            else:
                gain_res = loaded_filter_data['gain_res_save'][
                    int(index_harmonic), int(index_frev), int(index_red)]
                gain_res += (index_frev-int(index_frev))*(
                    loaded_filter_data['gain_res_save'][
                        int(index_harmonic), int(index_frev)+1, int(index_red)] - 
                    gain_res)
                gain_res += (index_red-int(index_red))*(
                    (index_frev-int(index_frev))*(
                        loaded_filter_data['gain_res_save'][
                            int(index_harmonic), int(index_frev)+1, int(index_red)+1] -
                        loaded_filter_data['gain_res_save'][
                            int(index_harmonic), int(index_frev), int(index_red)+1]) -
                    gain_res)

                Q_res = loaded_filter_data['Q_res_save'][
                    int(index_harmonic), int(index_frev), int(index_red)]
                Q_res += (index_frev-int(index_frev))*(
                    loaded_filter_data['Q_res_save'][
                        int(index_harmonic), int(index_frev)+1, int(index_red)] -
                    Q_res)
                Q_res += (index_red-int(index_red))*(
                    (index_frev-int(index_frev))*(
                        loaded_filter_data['Q_res_save'][
                            int(index_harmonic), int(index_frev)+1, int(index_red)+1] -
                        loaded_filter_data['Q_res_save'][
                            int(index_harmonic), int(index_frev), int(index_red)+1]) -
                    Q_res)

                gain_res *= 10  # The LUT was computed for 10 cavities, should be recomputed!

            delay_comb = 1/revolution_frequency
            delay_res = 1/revolution_frequency

            fittedParameters = [delay_comb,
                                gain_comb,
                                delay_res,
                                gain_res,
                                Q_res]
        else:
            fittedParameters = list(manual_params)

        final_filter = final_filter_calc(frequency_array,
                                         revolution_frequency,
                                         main_harmonic,
                                         *fittedParameters)

        final_impedance = final_impedance_calc(C10_Resonator. impedance,
                                               final_filter) * ZFactor

    else:

        final_impedance = C10_Resonator.impedance * ZFactor
        fittedParameters = 0

    return final_impedance, frequency_array, fittedParameters


def final_filter_calc(frequency_array, revolution_frequency, main_harmonic,
                      *params):

    (delay_comb, gain_comb, delay_res, gain_res, Q_res) = params

    # Notch filter for revolution harmonics
    freq_notch = 1/revolution_frequency
    depth = 2**(-4)
    comb_filter = depth/(1-(1-depth)*np.exp(
        -1j*2*np.pi*frequency_array*freq_notch))

    # Resonator to have 0 impedance reduction at central frequency
    main_resonator_filter = Resonators(
        1, main_harmonic*revolution_frequency, Q_res)
    main_resonator_filter.imped_calc(frequency_array)
    resonator_filter = main_resonator_filter.impedance

    final_filter = gain_comb*comb_filter*np.exp(
        -1j*2*np.pi*frequency_array*delay_comb) + \
        gain_res*resonator_filter*np.exp(
            -1j*2*np.pi*frequency_array*delay_res)

    return final_filter


def final_impedance_calc(initial_impedance, final_filter):

    final_impedance = initial_impedance / (1+initial_impedance*final_filter)

    return final_impedance


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

    # Save result
    save_result = False

    main_harmonic_range = [7, 21]
    revolution_frequency_range = [477e3]  # 3e8/(2*np.pi*100)
    impedance_reduction_target_range = [-15]

#     main_harmonic_range = np.arange(7, 21.1)
#     revolution_frequency_range = np.arange(430, 480.1, 10)*1e3
#     impedance_reduction_target_range = -np.arange(1, 21.1, 1)

    delay_comb_save = np.zeros(
        (len(main_harmonic_range),
         len(revolution_frequency_range),
         len(impedance_reduction_target_range)))
    gain_comb_save = np.zeros(
        (len(main_harmonic_range),
         len(revolution_frequency_range),
         len(impedance_reduction_target_range)))
    delay_res_save = np.zeros(
        (len(main_harmonic_range),
         len(revolution_frequency_range),
         len(impedance_reduction_target_range)))
    gain_res_save = np.zeros(
        (len(main_harmonic_range),
         len(revolution_frequency_range),
         len(impedance_reduction_target_range)))
    Q_res_save = np.zeros(
        (len(main_harmonic_range),
         len(revolution_frequency_range),
         len(impedance_reduction_target_range)))
    residue_save = np.zeros(
        (len(main_harmonic_range),
         len(revolution_frequency_range),
         len(impedance_reduction_target_range)))

    for index_h in range(len(main_harmonic_range)):
        for index_freq in range(len(revolution_frequency_range)):
            for index_red in range(len(impedance_reduction_target_range)):

                print('%d/%d-%d/%d-%d/%d' % (
                    index_h+1,
                    len(main_harmonic_range),
                    index_freq+1,
                    len(revolution_frequency_range),
                    index_red+1,
                    len(impedance_reduction_target_range)))

                main_harmonic = main_harmonic_range[index_h]
                revolution_frequency = revolution_frequency_range[index_freq]
                impedance_reduction_target = impedance_reduction_target_range[index_red]

                print(revolution_frequency, main_harmonic,
                      impedance_reduction_target)

                frequency_array = np.linspace(0, 20e6, 100000)
                all_revolution_freq = np.arange(
                    main_harmonic-6, main_harmonic+6+0.1)*revolution_frequency
                adjacent_revolution_freq = np.array(
                    [main_harmonic-1, main_harmonic+1])*revolution_frequency
                main_harmonic_freq = main_harmonic*revolution_frequency

                # Generating main impedance
                loaded_main_C10 = main_C10(main_harmonic*revolution_frequency)
                C10_Resonator = Resonators(loaded_main_C10[1],
                                           loaded_main_C10[0],
                                           loaded_main_C10[2])

                C10_Resonator.imped_calc(frequency_array)

                # Generating filter

                def residue_function(params, *args):

                    (impedance_reduction_target, main_harmonic_freq,
                     adjacent_revolution_freq, C10_impedance, delay_comb_init,
                     delay_res_init, delay_ratio, gain_ratio, Q_ratio) = args

                    (gain_comb, Q_res) = params

                    gain_res = -gain_comb

                    params = (delay_comb_init*delay_ratio,
                              np.abs(gain_comb*gain_ratio),
                              delay_res_init*delay_ratio,
                              gain_res*gain_ratio,
                              np.abs(Q_res*Q_ratio))

                    final_filter = final_filter_calc(frequency_array,
                                                     revolution_frequency,
                                                     main_harmonic,
                                                     *params)

                    final_impedance = final_impedance_calc(C10_impedance,
                                                           final_filter)

                    impedance_ratio_dB = impedance_to_dB(final_impedance,
                                                         C10_impedance,
                                                         True)

                    left_min = np.nanmin(
                        impedance_ratio_dB[
                            (frequency_array >
                             (adjacent_revolution_freq[0]-revolution_frequency/4)) *
                            (frequency_array <
                             (adjacent_revolution_freq[0]+revolution_frequency/4))])

                    right_min = np.nanmin(
                        impedance_ratio_dB[
                            (frequency_array >
                             (adjacent_revolution_freq[1]-revolution_frequency/4)) *
                            (frequency_array <
                             (adjacent_revolution_freq[1]+revolution_frequency/4))])

                    min_center = np.nanmin(
                        impedance_ratio_dB[
                            (frequency_array >
                             (main_harmonic_freq-revolution_frequency/4)) *
                            (frequency_array <
                             (main_harmonic_freq+revolution_frequency/4))])

                    max_center = np.nanmax(
                        impedance_ratio_dB[
                            (frequency_array >
                             (main_harmonic_freq-revolution_frequency/4)) *
                            (frequency_array <
                             (main_harmonic_freq+revolution_frequency/4))])

                    residue_left = (left_min-impedance_reduction_target)**2.
                    residue_right = (right_min-impedance_reduction_target)**2.
                    residue_center = (max_center-min_center)**2.

                    residue = np.sqrt(residue_left + residue_right +
                                      residue_center)

                    return residue

                delay_ratio = 1e5
                gain_ratio = 1e-4
                Q_ratio = 1e3

                delay_comb_init = 1/revolution_frequency / delay_ratio
                gain_comb_init = 1e-3 / gain_ratio * 10**(
                    impedance_reduction_target/20)
                delay_res_init = 1/revolution_frequency / delay_ratio
                Q_res_init = (350+(990-350)/(21-7)*main_harmonic) / Q_ratio

                startingParams = [gain_comb_init,
                                  Q_res_init]

                fittedParameters = minimize(
                            residue_function, startingParams,
                            args=(impedance_reduction_target,
                                  main_harmonic_freq,
                                  adjacent_revolution_freq,
                                  C10_Resonator.impedance,
                                  delay_comb_init,
                                  delay_res_init,
                                  delay_ratio,
                                  gain_ratio,
                                  Q_ratio),
                            options={'disp': True,
                                     'maxiter': 1000,
                                     'maxfev': 10000},
                            method='Nelder-Mead')['x']

                (gain_comb, Q_res) = fittedParameters
                residue = residue_function(fittedParameters,
                                           *(impedance_reduction_target,
                                             main_harmonic_freq,
                                             adjacent_revolution_freq,
                                             C10_Resonator.impedance,
                                             delay_comb_init,
                                             delay_res_init,
                                             delay_ratio,
                                             gain_ratio,
                                             Q_ratio))
                fittedParameters = (delay_comb_init*delay_ratio,
                                    np.abs(gain_comb*gain_ratio),
                                    delay_res_init*delay_ratio,
                                    -gain_comb*gain_ratio,
                                    np.abs(Q_res*Q_ratio))
                print(fittedParameters)

                delay_comb_save[index_h,
                                index_freq,
                                index_red] = fittedParameters[0]
                gain_comb_save[index_h,
                               index_freq,
                               index_red] = fittedParameters[1]
                delay_res_save[index_h,
                               index_freq,
                               index_red] = fittedParameters[2]
                gain_res_save[index_h,
                              index_freq,
                              index_red] = fittedParameters[3]
                Q_res_save[index_h,
                           index_freq,
                           index_red] = fittedParameters[4]
                residue_save[index_h,
                             index_freq,
                             index_red] = residue

                if save_result:
                    np.savez(script_path+'/1TFB_LUT_gains_and_delays.npz',
                             main_harmonic_range=main_harmonic_range,
                             revolution_frequency_range=revolution_frequency_range,
                             impedance_reduction_target_range=impedance_reduction_target_range,
                             delay_comb_save=delay_comb_save,
                             gain_comb_save=gain_comb_save,
                             delay_res_save=delay_res_save,
                             gain_res_save=gain_res_save,
                             Q_res_save=Q_res_save)

                if np.isclose(revolution_frequency, 477e3) and \
                        np.isclose(impedance_reduction_target, -15):
#                         np.isclose(main_harmonic, 21) and \

                    # Final impedance
                    final_filter = final_filter_calc(frequency_array,
                                                     revolution_frequency,
                                                     main_harmonic,
                                                     *fittedParameters)

                    final_impedance = final_impedance_calc(
                        C10_Resonator.impedance,
                        final_filter)

                    C10_impedance_dB = impedance_to_dB(C10_Resonator.impedance)

                    final_impedance_dB = impedance_to_dB(
                        final_impedance,
                        np.abs(C10_Resonator.impedance))

                    impedance_ratio_dB = impedance_to_dB(
                        final_impedance,
                        np.abs(C10_Resonator.impedance),
                                                         True)

#                     plt.figure('Filter')
#                     plt.clf()
#                     plt.plot(frequency_array, np.abs(final_filter))
#
#                     plt.figure('Filter phase')
#                     plt.clf()
#                     plt.plot(frequency_array, np.angle(final_filter))

                    plt.figure('Impedance')
                    plt.clf()
                    plt.plot(frequency_array/1e6,
                             np.abs(C10_Resonator.impedance)/1e3)
                    plt.plot(frequency_array/1e6, np.abs(final_impedance)/1e3)
                    plt.xlim(
                        (main_harmonic*revolution_frequency/1e6-5/2,
                         main_harmonic*revolution_frequency/1e6+5/2))
                    plt.hlines(0,
                               main_harmonic*revolution_frequency/1e6-5/2,
                               main_harmonic*revolution_frequency/1e6+5/2,
                               linewidth=0.5, alpha=0.5)
                    plt.vlines(all_revolution_freq/1e6, 0,
                               np.nanmax(np.abs(C10_Resonator.impedance)/1e3),
                               linewidth=0.5, alpha=0.5)
                    plt.xlabel('Frequency [MHz]')
                    plt.ylabel('Magnitude [$k\\Omega$]')
                    plt.tight_layout()
                    plt.savefig(
                        script_path+'/impedance_with_1TFB_h%d_-15dB.png' % (main_harmonic))

                    min_final_impedance_dB = np.nanmin(final_impedance_dB[
                            (frequency_array >
                             (main_harmonic*revolution_frequency-5e6/2)) *
                            (frequency_array <
                             (main_harmonic*revolution_frequency+5e6/2))])
                    max_final_impedance_dB = np.nanmax(final_impedance_dB[
                            (frequency_array >
                             (main_harmonic*revolution_frequency-5e6/2)) *
                            (frequency_array <
                             (main_harmonic*revolution_frequency+5e6/2))])
                    plt.figure('Impedance dB')
                    plt.clf()
                    plt.plot(frequency_array/1e6, C10_impedance_dB)
                    plt.plot(frequency_array/1e6, final_impedance_dB)
                    plt.xlim(
                        (main_harmonic*revolution_frequency/1e6-5/2,
                         main_harmonic*revolution_frequency/1e6+5/2))
                    plt.ylim((1.1*min_final_impedance_dB, 1.1*max_final_impedance_dB))
                    plt.vlines(all_revolution_freq/1e6, min_final_impedance_dB,
                               max_final_impedance_dB, linewidth=0.5, alpha=0.5)
                    plt.xlabel('Frequency [MHz]')
                    plt.ylabel('Magnitude [dB]')
                    plt.tight_layout()
                    plt.savefig(
                        script_path+'/impedance_dB_with_1TFB_h%d_-15dB.png' % (main_harmonic))

                    min_impedance_ratio_dB = np.nanmin(impedance_ratio_dB[
                            (frequency_array >
                             (main_harmonic*revolution_frequency-5e6/2)) *
                            (frequency_array <
                             (main_harmonic*revolution_frequency+5e6/2))])
                    max_impedance_ratio_dB = np.nanmax(impedance_ratio_dB[
                            (frequency_array >
                             (main_harmonic*revolution_frequency-5e6/2)) *
                            (frequency_array <
                             (main_harmonic*revolution_frequency+5e6/2))])
                    plt.figure('Impedance ratio')
                    plt.clf()
                    plt.plot(frequency_array/1e6, impedance_ratio_dB)
                    plt.xlim(
                        (main_harmonic*revolution_frequency/1e6-5/2,
                         main_harmonic*revolution_frequency/1e6+5/2))
                    plt.vlines(all_revolution_freq/1e6, min_impedance_ratio_dB,
                               max_impedance_ratio_dB, linewidth=0.5,
                               alpha=0.5)
                    plt.hlines(0,
                               main_harmonic*revolution_frequency/1e6-5/2,
                               main_harmonic*revolution_frequency/1e6+5/2,
                               linewidth=0.5, alpha=0.5)
                    plt.xlabel('Frequency [MHz]')
                    plt.ylabel('Magnitude [dB]')
                    plt.tight_layout()
                    plt.savefig(
                        script_path+'/impedance_ratio_dB_with_1TFB_h%d_-15dB.png' % (main_harmonic))

#                     plt.figure('Impedance phase')
#                     plt.clf()
#                     plt.plot(frequency_array/1e6, np.angle(C10_Resonator.impedance))
#                     plt.plot(frequency_array/1e6, np.angle(final_impedance))
#                     plt.xlabel('Frequency [MHz]')
#                     plt.ylabel('Phase')
#                     plt.tight_layout()
#                     plt.show()
