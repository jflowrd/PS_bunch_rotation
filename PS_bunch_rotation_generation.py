'''
PS bunch rotation RF generation
'''


# General imports
from __future__ import division
import numpy as np

def C40_77_prog(time_array, start_bunch_rotation_h84, normalizing_factor = 1.):


    time_linear_fit_40_78 = np.array([-247], dtype=float)
    time_linear_fit_40_77 = np.array([-231, -203, -192, 0], dtype=float)
    voltage_linear_fit_40_77 = np.array([0, 279, 296, 299], dtype=float)

    return normalizing_factor*np.interp(time_array,
                                        time_linear_fit_40_77-time_linear_fit_40_78[0]+start_bunch_rotation_h84,
                                        voltage_linear_fit_40_77)


def C40_78_prog(time_array, start_bunch_rotation_h84, normalizing_factor = 1.):

    time_linear_fit_40_78 = np.array([-247, -210, -185, -145, 0], dtype=float)
    voltage_linear_fit_40_78 = np.array([95, 280, 293, 298, 298], dtype=float)

    return normalizing_factor*np.interp(time_array,
                                        time_linear_fit_40_78-time_linear_fit_40_78[0]+start_bunch_rotation_h84,
                                        voltage_linear_fit_40_78)



def C80_08_prog(time_array, start_bunch_rotation_h168, length_step, step_amplitude, normalizing_factor = 1.):

#    time_linear_fit_80_08 = np.array([-249, -237, -96, -86, -70, 0], dtype=float)
#    voltage_linear_fit_80_08 = np.array([0, 38, 38, 25, 286, 308], dtype=float)

    time_linear_fit_80_08 = np.array([-249, -237, -96, -86, -55, -40, 0], dtype=float)
    voltage_linear_fit_80_08 = np.array([0, 38, 38, 38, 283, 297, 300], dtype=float)

    time_linear_fit_80_08[2:] += time_linear_fit_80_08[1]-time_linear_fit_80_08[2]+length_step
    voltage_linear_fit_80_08[1:4] = voltage_linear_fit_80_08[1:4]/np.mean(voltage_linear_fit_80_08[1:4])
    voltage_linear_fit_80_08[1:4] = voltage_linear_fit_80_08[1:4]*step_amplitude

    return normalizing_factor*np.interp(time_array,
                                        time_linear_fit_80_08-time_linear_fit_80_08[0]+start_bunch_rotation_h168,
                                        voltage_linear_fit_80_08)


def C80_88_prog(time_array, start_bunch_rotation_h168, length_step, step_amplitude, normalizing_factor = 1.):

#    time_linear_fit_80_08 = np.array([-249, -237, -96, -86, -70, 0], dtype=float)
#    voltage_linear_fit_80_08 = np.array([0, 38, 38, 25, 286, 308], dtype=float)
#
#    time_linear_fit_80_88 = np.array([-224, -215, -96, -86, -65, 0], dtype=float)
#    voltage_linear_fit_80_88 = np.array([0, 33, 33, 14, 302, 296], dtype=float)


    time_linear_fit_80_08 = np.array([-249, -237, -96, -86, -55, -40, 0], dtype=float)
    voltage_linear_fit_80_08 = np.array([0, 38, 38, 38, 283, 297, 300], dtype=float)

    time_linear_fit_80_88 = np.array([-224, -237, -96, -86, -55, -40, 0], dtype=float)
    voltage_linear_fit_80_88 = np.array([0, 38, 38, 38, 283, 297, 300], dtype=float)

    time_linear_fit_80_88[2:] += time_linear_fit_80_08[1]-time_linear_fit_80_08[2]+length_step

    voltage_linear_fit_80_88[1:4] = voltage_linear_fit_80_88[1:4]/np.mean(voltage_linear_fit_80_08[1:4])
    voltage_linear_fit_80_88[1:4] = voltage_linear_fit_80_88[1:4]*step_amplitude

    return normalizing_factor*np.interp(time_array,
                                        time_linear_fit_80_88-time_linear_fit_80_08[0]+start_bunch_rotation_h168,
                                        voltage_linear_fit_80_88)


def C80_phase_prog(time_array, start_bunch_rotation_h168, length_step, phase_0, phase_1):

#    time_linear_fit_80_08 = np.array([-249, -237, -96, -86, -70, 0], dtype=float)
#    voltage_linear_fit_80_08 = np.array([0, 38, 38, 25, 286, 308], dtype=float)
#
#    time_linear_fit_80_88 = np.array([-224, -215, -96, -86, -65, 0], dtype=float)
#    voltage_linear_fit_80_88 = np.array([0, 33, 33, 14, 302, 296], dtype=float)


    length_phase_jump = 10

    time_linear_fit_80_08 = np.array([-249, -237, -86-length_phase_jump, -86, -70, 0], dtype=float)

    time_linear_fit_80_08[2:] += time_linear_fit_80_08[1]-time_linear_fit_80_08[2]+length_step

    phase_linear_fit_80 = np.array([phase_0, phase_0, phase_0, phase_1, phase_1, phase_1], dtype=float)#np.array([phase_0, phase_1], dtype=float)


    return np.interp(time_array,
                     time_linear_fit_80_08-time_linear_fit_80_08[0]+start_bunch_rotation_h168,
                     phase_linear_fit_80)


def C40_phase_prog(time_array, start_bunch_rotation_h84, phase_0, phase_1):

    length_phase_jump = 10

    time_phase_jump = np.array([-length_phase_jump, 0], dtype=float)

    phase_jump = np.array([phase_0, phase_1], dtype=float)


    return np.interp(time_array,
                     time_phase_jump+start_bunch_rotation_h84,
                     phase_jump)


