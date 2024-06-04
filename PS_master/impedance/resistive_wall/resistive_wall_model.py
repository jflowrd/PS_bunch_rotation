'''
Created on 14 feb. 2018

@author: alasheen
'''

# General imports
import os
import numpy as np
from scipy.constants import c, epsilon_0, mu_0
Z0 = 1 / (epsilon_0 * c)

# Toolbox import
script_path = os.path.dirname(os.path.realpath(__file__))


def load_model(frequency):
    '''
    Method to compute resistive wall impedance for a program of beta
    '''

    ring_length = 2 * np.pi * 100
    chamber_radius = 35.e-3
    form_factor = 0.96

    # 70% of the ring in Stainless steel 316
    length_percent_1 = 0.7
    conductivity_1 = 1.3e6  # S/m
#     thickness_1 = 2e-3  # m

    impedance_1 = length_percent_1 * _resistive_wall_impedance(
        frequency, ring_length, conductivity_1, form_factor, chamber_radius)

    # 20% of the ring in Inconel X750
    length_percent_2 = 0.2
    conductivity_2 = 8.3e5  # S/m
#     thickness_2 = 1.5e-3  # m

    impedance_2 = length_percent_2 * _resistive_wall_impedance(
        frequency, ring_length, conductivity_2, form_factor, chamber_radius)

    total_impedance = impedance_1 + impedance_2

    return total_impedance


def _resistive_wall_impedance(frequency, length, conductivity, form_factor,
                              chamber_radius):
    '''
    Resistive wall impedance based on ...
    '''

    skin_depth = np.sqrt(2 / (conductivity * mu_0 * 2 * np.pi * frequency))

    Re_RW_1 = length / (2 * np.pi * conductivity *
                        skin_depth * form_factor * chamber_radius)

    Im_RW_1 = length / (2 * np.pi * conductivity *
                        skin_depth * form_factor * chamber_radius)

    return Re_RW_1 + 1j * Im_RW_1
