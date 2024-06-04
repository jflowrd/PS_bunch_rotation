'''
Loading the PS impedance model

@author: alasheen
'''

# General imports
# import sys
import os
import numpy as np
import warnings
from scipy.constants import m_p, c, e


# Toolbox import
from impedance_toolbox.impedance_toolbox.handle_impedance import handleImpedance, impedance2blond


class ImpedanceLoader(handleImpedance):
    """Class containing the functions to load the PS impedance

    This class inherits of the handleImpedance class which contains the methods
    used to load the data. This class which files are used, depending on the
    scenario and apply various modifications (number of
    element in the ring, damping wanted etc...)

    Attributes
    ----------
    MODEL : str
        scenario, can be 'top'. In case you want a user defined
        scenario, you can overwrite all the various options (see __init__).
        You can also use the method in handleImpedance to modify the impedance
        before converting it in BLonD variable with the impedance2blond class

    """

    def __init__(self, MODEL=None, particle_type='proton', **kwargs):

        folder = kwargs.pop(
            'folder',
            os.path.dirname(os.path.realpath(__file__)) + '/impedance/')

        handleImpedance.__init__(self, folder=folder)

        # Default parameters
        if particle_type == 'proton':

            self.mass = (m_p * c**2. / e)
            self.bending_radius = 70.07887
            self.circumference = 100

            self.momentum_top = 26.366945556921127e9
            self.b_field_top = self.momentum_top / \
                (self.bending_radius * c) * 1e4  # G
            self.b_field_top = self.b_field_top * (
                1
                - 1.441767e-6 * self.b_field_top
                + 2.947290e-10 * self.b_field_top**2.
                - 2.357026e-14 * self.b_field_top**3.)
            self.momentum_top = self.b_field_top * self.bending_radius * c * 1e-4
            self.total_energy_top = np.sqrt(
                self.momentum_top**2. + self.mass**2.)
            self.t_rev_top = (2 * np.pi * self.circumference) / (
                c * self.momentum_top / self.total_energy_top)
            self.f_rev_top = 1 / self.t_rev_top

            self.momentum_bottom = 2.139992595953578e9
            self.b_field_bottom = self.momentum_bottom / \
                (self.bending_radius * c) * 1e4  # G
            self.b_field_bottom = self.b_field_bottom * (
                1
                - 1.441767e-6 * self.b_field_bottom
                + 2.947290e-10 * self.b_field_bottom**2.
                - 2.357026e-14 * self.b_field_bottom**3.)
            self.momentum_bottom = self.b_field_bottom * self.bending_radius * c * 1e-4
            self.total_energy_bottom = np.sqrt(
                self.momentum_bottom**2. + self.mass**2.)
            self.t_rev_bottom = (2 * np.pi * self.circumference) / (
                c * self.momentum_bottom / self.total_energy_bottom)
            self.f_rev_bottom = 1 / self.t_rev_bottom

            self.momentum_plateau = 3.3011686295129175e9
            self.b_field_plateau = self.momentum_plateau / \
                (self.bending_radius * c) * 1e4  # G
            self.b_field_plateau = self.b_field_plateau * (
                1
                - 1.441767e-6 * self.b_field_plateau
                + 2.947290e-10 * self.b_field_plateau**2.
                - 2.357026e-14 * self.b_field_plateau**3.)
            self.momentum_plateau = self.b_field_plateau * self.bending_radius * c * 1e-4
            self.total_energy_plateau = np.sqrt(
                self.momentum_plateau**2. + self.mass**2.)
            self.t_rev_plateau = (2 * np.pi * self.circumference) / (
                c * self.momentum_plateau / self.total_energy_plateau)
            self.f_rev_plateau = 1 / self.t_rev_plateau

#             # Lead ions parameters
#             charge = 54  # e
#             mass = 193.702e9  # eV
#             BField = 1.25475  # T

#             # Xenon ions parameters
#             charge = 39  # e
#             mass = 120.054e9  # eV
#             BField = 1.13615  # T

        # Getting base model configuration (mainly for cavities)
        self.MODEL = MODEL

        # Getting kwargs based on model
        model_kwargs = self._MODEL_kwargs()
        for key in model_kwargs.keys():
            if key not in kwargs:
                kwargs[key] = model_kwargs.get(key)

        # Getting f_rev and momentum
        self.f_rev = kwargs.pop('f_rev')
        self.momentum = kwargs.pop('momentum')

        # Getting freq_array for generation and interpolation
        self.freq_array = kwargs.pop('freq_array', None)

        # Getting all settings
        # C10 cavities
        self.method_C10 = kwargs.pop('method_C10', 'parametric')
        self.n_elements_C10 = kwargs.pop('n_elements_C10', 10)
        self.harmonic_C10 = kwargs.pop('harmonic_C10', 21)
        self.freq_C10 = kwargs.pop('freq_C10',
                                   self.harmonic_C10 * self.f_rev)
        self.RshFactor_C10 = kwargs.pop('RshFactor_C10', 1.)
        self.QFactor_C10 = kwargs.pop('QFactor_C10', 1.)

        # C10 with 1TFB
        self.enable_C10_1TFB = kwargs.pop('enable_C10_1TFB', False)
        self.impedance_reduction_C10_1TFB = kwargs.pop(
            'impedance_reduction_C10_1TFB', -15.)
        self.main_harmonic_enableFB_C10_1TFB = kwargs.pop(
            'main_harmonic_enableFB_C10_1TFB', False)
        self.ZFactor_C10_1TFB = kwargs.pop('ZFactor_C10_1TFB', 1.)
        self.ZFactor_C10_1TFB_cav = kwargs.pop('ZFactor_C10_1TFB_cav', 1.)
        self.manual_params_C10_1TFB = kwargs.pop('manual_params_C10_1TFB',
                                                 None)

        # C10 cavities with closed gap relay
        self.method_C10_ClosedGap = kwargs.pop('method_C10_ClosedGap',
                                               'parametric')
        self.n_elements_C10_ClosedGap = kwargs.pop(
            'n_elements_C10_ClosedGap',
            int(11 - np.sum(self.n_elements_C10)))
        self.harmonic_C10_ClosedGap = kwargs.pop('harmonic_C10_ClosedGap',
                                                 21)
        self.RshFactor_C10_ClosedGap = kwargs.pop('RshFactor_C10_ClosedGap',
                                                  1.)
        self.QFactor_C10_ClosedGap = kwargs.pop('QFactor_C10_ClosedGap', 1.)

        # C20 cavities
        self.n_elements_C20 = kwargs.pop('n_elements_C20', 0)
        self.freq_C20 = kwargs.pop('freq_C20', 42 * self.f_rev_top)
        self.filename_C20 = kwargs.pop(
            'filename_C20',
            '/rf_cavities/C20/All/Resonators/single_resonator_20MHz_100V.txt')
        self.RshFactor_C20 = kwargs.pop('RshFactor_C20', 1.)
        self.QFactor_C20 = kwargs.pop('QFactor_C20', 1.)

        # C20 cavities with MHFB
        self.enable_C20_MHFB = kwargs.pop('enable_C20_MHFB', False)
        self.impedance_reduction_C20_MHFB = kwargs.pop(
            'impedance_reduction_C20_MHFB', -20)
        self.bandwidth_target_C20_MHFB = kwargs.pop(
            'bandwidth_target_C20_MHFB', 5e3)
        self.main_harmonic_enableFB_C20_MHFB = kwargs.pop(
            'main_harmonic_enableFB_C20_MHFB', False)
        self.manual_params_C20_LV = kwargs.pop(
            'manual_params_C20_LV', None)
        self.manual_params_C20_HV = kwargs.pop(
            'manual_params_C20_HV', None)
        self.ZFactor_C20_MHFB = kwargs.pop('ZFactor_C20_MHFB', 1.)

        # C40 cavities
        self.freq_C40 = kwargs.pop('freq_C40', 84 * self.f_rev_top)
        self.method_C40 = 'DFB'
        self.filename_C40 = kwargs.pop(
            'filename_C40',
            ['/rf_cavities/C40/Individual/C77/DFB/single_resonator_DFB_b72.txt',
             '/rf_cavities/C40/Individual/C78/DFB/single_resonator_DFB_b72.txt'])
        self.RshFactor_C40 = kwargs.pop('RshFactor_C40', 1.)
        self.QFactor_C40 = kwargs.pop('QFactor_C40', 1.)
        self.ZFactor_C40 = kwargs.pop('ZFactor_C40', 1.)

        # C40 cavities with MHFB
        self.enable_C40_MHFB = kwargs.pop('enable_C40_MHFB', False)
        self.filename_C40_MHFB = kwargs.pop(
            'filename_C40_MHFB',
            ['/rf_cavities/C40/Individual/C77/Resonators/single_resonator_b72.txt',
             '/rf_cavities/C40/Individual/C78/Resonators/single_resonator_b72.txt'])
        self.main_harmonic_enableFB_C40_MHFB = kwargs.pop(
            'main_harmonic_enableFB_C40_MHFB',
            ['C77', 'C78'])
        self.impedance_reduction_target_C40_MHFB = kwargs.pop(
            'impedance_reduction_target_C40_MHFB', -20)
        self.bandwidth_target_C40_MHFB = kwargs.pop(
            'bandwidth_target_C40_MHFB', 5e3)
        self.manual_params_C40_77_MHFB = kwargs.pop(
            'manual_params_C40_77_MHFB', None)
        self.manual_params_C40_78_MHFB = kwargs.pop(
            'manual_params_C40_78_MHFB', None)
        self.RshFactor_C40_MHFB = kwargs.pop('RshFactor_C40_MHFB', 1.)
        self.QFactor_C40_MHFB = kwargs.pop('QFactor_C40_MHFB', 1.)
        self.ZFactor_C40_MHFB = kwargs.pop('ZFactor_C40_MHFB', 1.)

        # C40 cavities HOMs
        self.n_elements_C40_HOMs = kwargs.pop('n_elements_C40_HOMs',
                                              len(self.filename_C40))
        self.filename_C40_HOMs = kwargs.pop(
            'filename_C40_HOMs',
            '/rf_cavities/C40/High_Order_Modes/Resonators/HOMs_C40.txt')
        self.RshFactor_C40_HOMs = kwargs.pop('RshFactor_C40_HOMs', 1.)
        self.QFactor_C40_HOMs = kwargs.pop('QFactor_C40_HOMs', 1.)

        # C80 cavities
        self.freq_C80 = kwargs.pop('freq_C80', 168 * self.f_rev_top)
        self.method_C80 = 'DFB'
        self.filename_C80 = kwargs.pop(
            'filename_C80',
            ['/rf_cavities/C80/Individual/C08/DFB/single_resonator_DFB_b72.txt',
             '/rf_cavities/C80/Individual/C88/DFB/single_resonator_DFB_b72.txt'])
#              '/rf_cavities/C80/Individual/C89/DFB/single_resonator_DFB_b72.txt'
        self.RshFactor_C80 = kwargs.pop('RshFactor_C80', 1.)
        self.QFactor_C80 = kwargs.pop('QFactor_C80', 1.)
        self.ZFactor_C80 = kwargs.pop('ZFactor_C80', 1.)

        # C80 cavities with MHFB
        self.enable_C80_MHFB = kwargs.pop('enable_C80_MHFB', False)
        self.filename_C80_MHFB = kwargs.pop(
            'filename_C80_MHFB',
            ['/rf_cavities/C80/Individual/C08/Resonators/single_resonator_b72.txt',
             '/rf_cavities/C80/Individual/C88/Resonators/single_resonator_b72.txt'])
#              '/rf_cavities/C80/Individual/C89/Resonators/single_resonator_b72.txt'
        self.main_harmonic_enableFB_C80_MHFB = kwargs.pop(
            'main_harmonic_enableFB_C80_MHFB',
            ['C08', 'C88', 'C89'])
        self.impedance_reduction_target_C80_MHFB = kwargs.pop(
            'impedance_reduction_target_C80_MHFB', -20)
        self.bandwidth_target_C80_MHFB = kwargs.pop(
            'bandwidth_target_C80_MHFB', 5e3)
        self.manual_params_C80_08_MHFB = kwargs.pop(
            'manual_params_C80_08_MHFB', None)
        self.manual_params_C80_88_MHFB = kwargs.pop(
            'manual_params_C80_88_MHFB', None)
        self.manual_params_C80_89_MHFB = kwargs.pop(
            'manual_params_C80_89_MHFB', None)
        self.RshFactor_C80_MHFB = kwargs.pop('RshFactor_C80_MHFB', 1.)
        self.QFactor_C80_MHFB = kwargs.pop('QFactor_C80_MHFB', 1.)
        self.ZFactor_C80_MHFB = kwargs.pop('ZFactor_C80_MHFB', 1.)

        # C80 cavities HOMs
        self.n_elements_C80_HOMs = kwargs.pop('n_elements_C80_HOMs',
                                              len(self.filename_C80))
        self.filename_C80_HOMs = kwargs.pop(
            'filename_C80_HOMs',
            '/rf_cavities/C80/High_Order_Modes/Resonators/HOMs_C80.txt')
        self.RshFactor_C80_HOMs = kwargs.pop('RshFactor_C80_HOMs', 1.)
        self.QFactor_C80_HOMs = kwargs.pop('QFactor_C80_HOMs', 1.)

        # C200 cavities
        self.n_elements_C200 = kwargs.pop('n_elements_C200', 6)
        self.filename_C200 = kwargs.pop(
            'filename_C200',
            '/rf_cavities/C200/Resonators/from_CST_C200_single_damped.txt')
        self.RshFactor_C200 = kwargs.pop('RshFactor_C200', 1.)
        self.QFactor_C200 = kwargs.pop('QFactor_C200', 1.)

        # Finemet cavity
        self.filename_Finemet = '/rf_cavities/Finemet/Resonators_2018/multi_resonators_impedance.txt'
        self.RshFactor_Finemet = 1.
        self.QFactor_Finemet = 1.

        # Finemet with MHFB
        self.enable_Finemet_MHFB = kwargs.pop('enable_Finemet_MHFB', False)
        self.filename_Finemet_MHFB = kwargs.pop(
            'filename_Finemet_MHFB',
            '/rf_cavities/Finemet/Resonators_2018/multi_resonators_impedance.txt')
        self.impedance_reduction_target_Finemet_MHFB = kwargs.pop(
            'impedance_reduction_target_Finemet_MHFB', -15)
        self.bandwidth_target_Finemet_MHFB = kwargs.pop(
            'bandwidth_target_Finemet_MHFB', 8e3)
        self.manual_params_Finemet_MHFB = kwargs.pop(
            'manual_params_Finemet_MHFB', None)
        self.RshFactor_Finemet_MHFB = kwargs.pop('RshFactor_Finemet_MHFB', 1.)
        self.QFactor_Finemet_MHFB = kwargs.pop('QFactor_Finemet_MHFB', 1.)
        self.ZFactor_Finemet_MHFB = kwargs.pop('ZFactor_Finemet_MHFB', 1.)

        # Kickers
        self.filename_Kickers = kwargs.pop('filename_Kickers', [
            '/kickers/KFA04/Resonators/resonators_broadband_extra_resonances.txt',
            '/kickers/KFA13/Resonators/resonators_broadband.txt',
            '/kickers/KFA21/Resonators/resonators_broadband.txt',
            '/kickers/KFA28/Resonators/resonators_broadband_extra_resonances_best_lowfreq.txt',
            '/kickers/KFA45/Resonators/resonators_broadband.txt',
            '/kickers/KFA71/Resonators/resonators_broadband.txt',
            '/kickers/KFA79/Resonators/resonators_broadband_extra_resonances.txt',
            '/kickers/BFA09/Resonators/resonators_broadband_extra_resonances.txt',
            '/kickers/BFA21/Resonators/resonators_broadband.txt'])

        self.RshFactor_Kickers = kwargs.pop('RshFactor_Kickers', 1.)
        self.QFactor_Kickers = kwargs.pop('QFactor_Kickers', 1.)

        # Septa
        self.filename_Septa = kwargs.pop('filename_Septa', [
            '/septa/SMH42/Resonators/multi_resonator_and_ImZ_over_f.txt'])

        self.RshFactor_Septa = kwargs.pop('RshFactor_Septa', 1.)
        self.QFactor_Septa = kwargs.pop('QFactor_Septa', 1.)

        # Vacuum
        # Sector valve
        self.filename_SectorValves = kwargs.pop(
            'filename_SectorValves',
            ['/vacuum/sector_valve/cst_raw_data/eigen'])
        self.n_elements_SectorValves = kwargs.pop('n_elements_SectorValves',
                                                  10.)
        self.RshFactor_SectorValves = kwargs.pop('RshFactor_SectorValves', 1.)
        self.QFactor_SectorValves = kwargs.pop('QFactor_SectorValves', 1.)

        # MU sections upstream
        self.filename_MU_upstream = kwargs.pop(
            'filename_MU_upstream',
            ['/vacuum/bellows/Inductive/ImZ_over_f.txt'])
        self.n_elements_MU_upstream = kwargs.pop('n_elements_MU_upstream', 99.)
        self.RshFactor_MU_upstream = kwargs.pop('RshFactor_MU_upstream', 1.)
        self.QFactor_MU_upstream = kwargs.pop('QFactor_MU_upstream', 1.)

        # MU sections downstream
        self.filename_MU_downstream = kwargs.pop(
            'filename_MU_downstream',
            ['/vacuum/pumping_manifold_empty/Resonators/multi_resonator_and_ImZ_over_f.txt',
             '/vacuum/pumping_manifold_CODD/Resonators/multi_resonator_and_ImZ_over_f.txt'])
        self.n_elements_MU_downstream = kwargs.pop('n_elements_MU_downstream',
                                                   [36, 100 - 36])
        self.RshFactor_MU_downstream = kwargs.pop(
            'RshFactor_MU_downstream', 1.)
        self.QFactor_MU_downstream = kwargs.pop('QFactor_MU_downstream', 1.)

        # Vacuum Flanges
        self.filename_Flanges = kwargs.pop(
            'filename_Flanges',
            '/vacuum/flanges/gap_ps195/cst_raw_data/eigen')
        self.n_elements_Flanges = kwargs.pop(
            'n_elements_Flanges',
            int(259 - np.sum(self.n_elements_MU_downstream)))
        self.RshFactor_Flanges = kwargs.pop(
            'RshFactor_Flanges', 1.)
        self.QFactor_Flanges = kwargs.pop('QFactor_Flanges', 1.)

        # Step transition
        self.filename_Steps = kwargs.pop(
            'filename_Steps',
            '/vacuum/flanges/steps/Inductive/ImZ_over_f.txt')
        self.ZFactor_Steps = kwargs.pop(
            'ZFactor_Steps', 1.)

        # Flanges ground loops and rf bypasses
        self.filename_Flanges_GroundLoops = kwargs.pop(
            'filename_Flanges_GroundLoops',
            '/vacuum/flanges/ground_loops/Resonators/multi_resonator_1_bypass_low_freq.txt')
        self.n_elements_Flanges_GroundLoops = kwargs.pop(
            'n_elements_Flanges_GroundLoops', 200)
        self.n_bypasses_Flanges_GroundLoops = kwargs.pop(
            'n_bypasses_Flanges_GroundLoops', 200)
        self.ZFactor_Flanges_GroundLoops = kwargs.pop(
            'ZFactor_Flanges_GroundLoops', 1.)

        # Beam instrumentation
        self.filename_BI = kwargs.pop(
            'filename_BI',
            ['/beam_instrumentation/WS/cst_raw_data/PS_SPS_wire_scanner.dat',
             '/beam_instrumentation/Stripline_BPM_SD72/Resonators/multi_resonator_and_ImZ_over_f.txt',
             '/beam_instrumentation/BGI/Vertical/Resonators/multi_resonator.txt',
             '/beam_instrumentation/WCM_SD03/circuit_model/circuit_model.py'])
        self.n_elements_BI = kwargs.pop('n_elements_BI', [4, 1, 1, 2])
        self.RshFactor_BI = kwargs.pop(
            'RshFactor_BI', 1.)
        self.QFactor_BI = kwargs.pop('QFactor_BI', 1.)

        # Misc equipment
        self.filename_Misc = kwargs.pop(
            'filename_Misc',
            ['/miscellaneous/internal_dump/PreLS2/Resonators/multi_resonator_and_ImZ_over_f.txt',
             '/miscellaneous/TFB_Kicker_SD97/Resonators/multi_resonator_and_ImZ_over_f.txt',
             '/miscellaneous/ralentisseur/Resonators/cst_eig.txt'])
        self.n_elements_Misc = kwargs.pop('n_elements_Misc', [2, 1, 1])
        self.RshFactor_Misc = kwargs.pop(
            'RshFactor_Misc', 1.)
        self.QFactor_Misc = kwargs.pop('QFactor_Misc', 1.)

        # Resistive wall
        self.ZFactor_ResistiveWall = kwargs.pop('ZFactor_ResistiveWall', 1.)

        # Space charge
        self.load_SC = kwargs.pop('load_SC', True)
        self.emittance_x_norm_SC = kwargs.pop('emittance_x_norm_SC', 2e-6)
        self.emittance_y_norm_SC = kwargs.pop('emittance_y_norm_SC', 2e-6)
        self.momentum_spread_SC = kwargs.pop('momentum_spread_SC', 1e-3)
        self.method_SC = kwargs.pop('method_SC', 'rectangle')
        self.aperture_X_SC = kwargs.pop('aperture_X_SC', None)
        self.aperture_Y_SC = kwargs.pop('aperture_Y_SC', None)
        self.beta_X_SC = kwargs.pop('beta_X_SC', None)
        self.beta_Y_SC = kwargs.pop('beta_Y_SC', None)
        self.disp_X_SC = kwargs.pop('disp_X_SC', None)
        self.optics_S_SC = kwargs.pop('optics_S_SC', None)
        self.reinterpolate_SC = kwargs.pop('reinterpolate_SC', True)
        self.return_diagnose_SC = kwargs.pop('return_diagnose_SC', False)
        self.BI_device_x_SC = kwargs.pop('BI_device_x_SC', 'BWSH65')
        self.BI_device_y_SC = kwargs.pop('BI_device_y_SC', 'BWSV85')
        self.ZFactor_SC = kwargs.pop('ZFactor_SC', 1.)

    def _MODEL_kwargs(self):
        '''
        Getting the cavities settings and setting them in the init
        if not defined by the user
        '''

        model_kwargs = {}

        # Model based cavity settings
        if self.MODEL is None:

            pass

        elif self.MODEL == 'top':

            model_kwargs['f_rev'] = self.f_rev_top
            model_kwargs['momentum'] = self.momentum_top

            # C10 cavities
            model_kwargs['harmonic_C10'] = 21
            model_kwargs['freq_C10'] = model_kwargs['harmonic_C10'] * \
                model_kwargs['f_rev']
            model_kwargs['n_elements_C10'] = 4

        elif self.MODEL == 'ramp':

            # C10 cavities
            model_kwargs['harmonic_C10'] = 21
            model_kwargs['freq_C10'] = model_kwargs['harmonic_C10'] * \
                model_kwargs['f_rev']
            model_kwargs['n_elements_C10'] = 10

        elif self.MODEL == 'foursplit':

            model_kwargs['f_rev'] = self.f_rev_top
            model_kwargs['momentum'] = self.momentum_top

            # C10 cavities
            model_kwargs['harmonic_C10'] = 21
            model_kwargs['freq_C10'] = model_kwargs['harmonic_C10'] * \
                model_kwargs['f_rev']
            model_kwargs['n_elements_C10'] = 1

            # C20 cavities
            model_kwargs['n_elements_C20'] = 1

            # C40 MHFB
            model_kwargs['main_harmonic_FB_C40'] = ['C40-77']

        elif self.MODEL == 'h84':

            model_kwargs['f_rev'] = self.f_rev_top
            model_kwargs['momentum'] = self.momentum_top

            # C10 cavities
            model_kwargs['harmonic_C10'] = 21
            model_kwargs['freq_C10'] = model_kwargs['harmonic_C10'] * \
                model_kwargs['f_rev']
            model_kwargs['n_elements_C10'] = 0

            # C40 MHFB
            model_kwargs['main_harmonic_FB_C40'] = ['C40-77']

        elif self.MODEL == 'bunch_rotation':

            model_kwargs['f_rev'] = self.f_rev_top
            model_kwargs['momentum'] = self.momentum_top

            # C10 cavities
            model_kwargs['harmonic_C10'] = 21
            model_kwargs['freq_C10'] = model_kwargs['harmonic_C10'] * \
                model_kwargs['f_rev']
            model_kwargs['n_elements_C10'] = 0

            # C40 MHFB
            model_kwargs['main_harmonic_FB_C40'] = ['C40-77', 'C40-78']

            # C80 MHFB
            model_kwargs['main_harmonic_FB_C80'] = ['C80-08', 'C80-88']

        elif self.MODEL == 'bottom_h7':

            model_kwargs['f_rev'] = self.f_rev_bottom
            model_kwargs['momentum'] = self.momentum_bottom

            # C10 cavities
            model_kwargs['harmonic_C10'] = 7
            model_kwargs['freq_C10'] = model_kwargs['harmonic_C10'] * \
                model_kwargs['f_rev']
            model_kwargs['n_elements_C10'] = 10

        elif self.MODEL == 'plateau_h21':

            model_kwargs['f_rev'] = self.f_rev_plateau
            model_kwargs['momentum'] = self.momentum_plateau

            # C10 cavities
            model_kwargs['harmonic_C10'] = 21
            model_kwargs['freq_C10'] = model_kwargs['harmonic_C10'] * \
                model_kwargs['f_rev']
            model_kwargs['n_elements_C10'] = 0

        elif self.MODEL == 'trisplit':

            model_kwargs['f_rev'] = self.f_rev_plateau
            model_kwargs['momentum'] = self.momentum_plateau

            # C10 cavities
            model_kwargs['harmonic_C10'] = np.array([7, 14, 21])
            model_kwargs['freq_C10'] = model_kwargs['harmonic_C10'] * \
                model_kwargs['f_rev']
            model_kwargs['n_elements_C10'] = np.array([3, 4, 3])

        else:

            warnings.warn(
                'The model %s does not belong to the list !!'
                % (self.MODEL))

        return model_kwargs

    def importC10(self, n_elements=None, freq=None, method=None,
                  RshFactor=None, QFactor=None):
        '''
        Determines the impedance model used to represent the C10 cavities
        '''

        if n_elements is None:
            n_elements = self.n_elements_C10
        if isinstance(n_elements, (float, int)) and (n_elements == 0):
            return
        if freq is None:
            freq = self.freq_C10
            if isinstance(freq, float):
                freq = [freq]
            if isinstance(n_elements, (float, int)):
                n_elements = [n_elements]
            if (len(n_elements) > 1) \
                    and (len(freq) == 1):
                freq = [freq] * len(n_elements)
        if method is None:
            method = self.method_C10
        if RshFactor is None:
            RshFactor = self.RshFactor_C10
        if QFactor is None:
            QFactor = self.QFactor_C10

        if method == 'parametric':
            from impedance.rf_cavities.C10.All.Parametric.All_C10_Parametric_model import load_model as load_C10

            for index_freq, single_freq in enumerate(freq):
                loaded_resonator = load_C10(single_freq)

                self.importResonatorsFromList(
                    loaded_resonator,
                    'rf_C10_%.0fMHz' % (single_freq / 1e6),
                    unitFreq=1., unitRsh=1.,
                    RshFactor=RshFactor * n_elements[index_freq],
                    QFactor=QFactor)

        else:

            filename = method

            self.importResonatorFromFile(filename, unitFreq=1., unitRsh=1.,
                                         RshFactor=RshFactor * n_elements,
                                         QFactor=QFactor)

    def importC10_1TFB(self, n_elements=None,
                       main_harmonic=None, f_rev=None,
                       impedance_reduction_target=None,
                       main_harmonic_enableFB=None,
                       frequency_array=None, manual_params=None,
                       ZFactor_cav=None,  ZFactor=None):
        '''
        Determines the impedance model used to represent the C10 cavities
        with the 1TFB
        '''

        if n_elements is None:
            n_elements = self.n_elements_C10
        if isinstance(n_elements, (float, int)) and (n_elements == 0):
            return
        if main_harmonic is None:
            main_harmonic = self.harmonic_C10
            if (isinstance(n_elements, list)
                or isinstance(n_elements, np.ndarray)) \
                    and isinstance(main_harmonic, float):
                main_harmonic = [main_harmonic] * len(n_elements)
        if f_rev is None:
            f_rev = self.f_rev
        if impedance_reduction_target is None:
            impedance_reduction_target = self.impedance_reduction_C10_1TFB
            if (isinstance(n_elements, list)
                or isinstance(n_elements, np.ndarray)) \
                    and isinstance(impedance_reduction_target, float):
                impedance_reduction_target = [impedance_reduction_target] * len(n_elements)
        if main_harmonic_enableFB is None:
            main_harmonic_enableFB = self.main_harmonic_enableFB_C10_1TFB
            if (isinstance(n_elements, list)
                or isinstance(n_elements, np.ndarray)) \
                    and isinstance(main_harmonic_enableFB, bool):
                main_harmonic_enableFB = [main_harmonic_enableFB] * len(n_elements)
        if frequency_array is None:
            frequency_array = self.freq_array
        if manual_params is None:
            manual_params = self.manual_params_C10_1TFB
        if ZFactor_cav is None:
            ZFactor_cav = self.ZFactor_C10_1TFB_cav
        if ZFactor is None:
            ZFactor = self.ZFactor_C10_1TFB

        from impedance.rf_cavities.C10.OTFB.Parametric.All_C10_1TFB_Parametric_model import load_model as load_C10_1TFB

        for index_h, single_h in enumerate(main_harmonic):
            final_impedance, frequency_array, fittedParameters = load_C10_1TFB(
                single_h, f_rev, impedance_reduction_target[index_h],
                main_harmonic_FB=main_harmonic_enableFB[index_h],
                frequency_array=frequency_array, manual_params=manual_params,
                ZFactor=ZFactor_cav)

            self.importInputTableFromList(
                [frequency_array, final_impedance.real, final_impedance.imag],
                'rf_C10_%.0fMHz_%d_1TFB' % ((single_h * f_rev) / 1e6, index_h),
                unitFreq=1., ZFactor=ZFactor * n_elements[index_h])

        return fittedParameters

    def importC10_ClosedGapRelay(self, n_elements=None,
                                 method=None,
                                 harmonic=None,
                                 RshFactor=None, QFactor=None):
        '''
        Determines the impedance model used to represent the C10 cavities
        '''

        if n_elements is None:
            n_elements = self.n_elements_C10_ClosedGap
        if isinstance(n_elements, (float, int)) and (n_elements == 0):
            return
        if method is None:
            method = self.method_C10
        if harmonic is None:
            harmonic = self.harmonic_C10_ClosedGap
        if RshFactor is None:
            RshFactor = self.RshFactor_C10
        if QFactor is None:
            QFactor = self.QFactor_C10

        if method == 'parametric':
            from impedance.rf_cavities.C10.Closed_Gap_Relay.Parametric.All_C10_ClosedGapRelay_Parametric_model import load_model as load_C10_ClosedGap

            loaded_resonator = load_C10_ClosedGap(harmonic)
            im_z_over_n = loaded_resonator[-1]

            self.importResonatorsFromList(
                list(loaded_resonator[0:3]), 'rf_C10_ClosedGapRelay',
                unitFreq=1., unitRsh=1.,
                RshFactor=RshFactor * n_elements,
                QFactor=QFactor)

            self.importImZ_over_f('rf_C10_ClosedGapRelay_ImZ/f',
                                  im_z_over_n * RshFactor * n_elements)

        else:

            warnings.warn("Need to run the parametric in " +
                          "importCavitiesC10_ClosedGapRelay at the moment")

    def importC20(self, n_elements=None, freq=None, filename=None,
                  RshFactor=None, QFactor=None):
        '''
        Determines the impedance model used to represent the C20 cavities
        '''

        if n_elements is None:
            n_elements = self.n_elements_C20
        if isinstance(n_elements, (float, int)) and (n_elements == 0):
            return
        if freq is None:
            freq = self.freq_C20
        if filename is None:
            filename = self.filename_C20
        if RshFactor is None:
            RshFactor = self.RshFactor_C20
        if QFactor is None:
            QFactor = self.QFactor_C20

        self.importResonatorFromFile(filename, unitFreq=1., unitRsh=1.,
                                     RshFactor=RshFactor, QFactor=QFactor)

        self.set_fr_resonatorOrTWC(filename, freq)

    def importC20_MHFB(self, n_elements=None, filename_C20=None,
                       f_rev=None,
                       impedance_reduction_target=None,
                       bandwidth_target=None,
                       main_harmonic_enableFB=False,
                       frequency_array=None,
                       manual_params_C20_LV=None,
                       manual_params_C20_HV=None,
                       RshFactor=None, QFactor=None, ZFactor=None):
        '''
        Determines the impedance model used to represent the C20 cavity
        with MHFB
        '''

        if n_elements is None:
            n_elements = self.n_elements_C20
        if isinstance(n_elements, (float, int)) and (n_elements == 0):
            return
        if filename_C20 is None:
            filename_C20 = self.filename_C20
        if f_rev is None:
            f_rev = self.f_rev
        if impedance_reduction_target is None:
            impedance_reduction_target = self.impedance_reduction_C20_MHFB
        if bandwidth_target is None:
            bandwidth_target = self.bandwidth_target_C20_MHFB
        if main_harmonic_enableFB is None:
            main_harmonic_enableFB = self.main_harmonic_enableFB_C20_MHFB
        if frequency_array is None:
            frequency_array = self.freq_array
        if manual_params_C20_LV is None:
            manual_params_C20_LV = self.manual_params_C20_LV
        if manual_params_C20_HV is None:
            manual_params_C20_HV = self.manual_params_C20_HV
        if RshFactor is None:
            RshFactor = self.RshFactor_C20
        if QFactor is None:
            QFactor = self.QFactor_C20
        if ZFactor is None:
            ZFactor = self.ZFactor_C20_MHFB

        from impedance.rf_cavities.C20.MHFB.Resonators.All_C20_MHFB_Resonators import load_model as load_C20_MHFB

        final_impedance, frequency_array, fittedParameters = load_C20_MHFB(
            filename_C20, f_rev,
            impedance_reduction_target, bandwidth_target,
            main_harmonic_FB=main_harmonic_enableFB,
            frequency_array=frequency_array,
            manual_params_C20_LV=manual_params_C20_LV,
            manual_params_C20_HV=manual_params_C20_HV,
            RshFactor=RshFactor, QFactor=QFactor)

        self.importInputTableFromList(
            [frequency_array, final_impedance.real, final_impedance.imag],
            'rf_C20_MHFB', unitFreq=1., ZFactor=ZFactor)

        return fittedParameters

    def importC40(self, method=None, filename_list=None,
                  frequency_array=None,
                  RshFactor=None, QFactor=None, ZFactor=None):
        '''
        Determines the impedance model used to represent the C40 cavities
        without MHFB
        '''

        if method is None:
            method = self.method_C40
        if filename_list is None:
            filename_list = self.filename_C40
        if frequency_array is None:
            frequency_array = self.freq_array
        if RshFactor is None:
            RshFactor = self.RshFactor_C40
        if QFactor is None:
            QFactor = self.QFactor_C40
        if ZFactor is None:
            ZFactor = self.ZFactor_C40

        if isinstance(filename_list, str):
            filename_list = [filename_list]

        if method == 'Resonators':
            for filename in filename_list:
                self.importResonatorFromFile(filename, unitFreq=1., unitRsh=1.,
                                             RshFactor=RshFactor,
                                             QFactor=QFactor)
        elif method == 'DFB':
            for filename in filename_list:

                if 'C77' in filename:
                    from impedance.rf_cavities.C40.Individual.C77.DFB.C77_fit_DFB import load_model as load_C77
                    load_C40 = load_C77
                elif 'C78' in filename:
                    from impedance.rf_cavities.C40.Individual.C78.DFB.C78_fit_DFB import load_model as load_C78
                    load_C40 = load_C78
                else:
                    continue

                loaded_impedance, freqArray = load_C40(
                    filename, folder=self.impedanceFolder,
                    freqArray=frequency_array)
                self.importInputTableFromList([freqArray,
                                               loaded_impedance.real,
                                               loaded_impedance.imag],
                                              filename,
                                              unitFreq=1, ZFactor=ZFactor)

        else:
            raise RuntimeError('Method for importC40 not recognized')

    def importC40_MHFB(self, filename_C40=None, f_rev=None,
                       impedance_reduction_target=None,
                       bandwidth_target=None,
                       main_harmonic_enableFB=None,
                       frequency_array=None,
                       manual_params_C40_77=None,
                       manual_params_C40_78=None,
                       RshFactor=None, QFactor=None, ZFactor=None):
        '''
        Determines the impedance model used to represent the C40 cavities
        with MHFB
        '''

        if filename_C40 is None:
            filename_C40 = self.filename_C40_MHFB
        if f_rev is None:
            f_rev = self.f_rev
        if impedance_reduction_target is None:
            impedance_reduction_target = self.impedance_reduction_target_C40_MHFB
        if bandwidth_target is None:
            bandwidth_target = self.bandwidth_target_C40_MHFB
        if main_harmonic_enableFB is None:
            main_harmonic_enableFB = self.main_harmonic_enableFB_C40_MHFB
        if frequency_array is None:
            frequency_array = self.freq_array
        if manual_params_C40_77 is None:
            manual_params_C40_77 = self.manual_params_C40_77_MHFB
        if manual_params_C40_78 is None:
            manual_params_C40_78 = self.manual_params_C40_78_MHFB
        if RshFactor is None:
            RshFactor = self.RshFactor_C40_MHFB
        if QFactor is None:
            QFactor = self.QFactor_C40_MHFB
        if ZFactor is None:
            ZFactor = self.ZFactor_C40_MHFB

        from impedance.rf_cavities.C40.MHFB.Resonators.All_C40_MHFB_Resonators import load_model

        final_impedance, frequency_array, fittedParameters = load_model(
            filename_C40, f_rev,
            impedance_reduction_target, bandwidth_target,
            main_harmonic_FB=main_harmonic_enableFB,
            frequency_array=frequency_array,
            manual_params_C40_77=manual_params_C40_77,
            manual_params_C40_78=manual_params_C40_78,
            RshFactor=RshFactor, QFactor=QFactor)

        self.importInputTableFromList(
            [frequency_array, final_impedance.real, final_impedance.imag],
            'rf_C40_MHFB', unitFreq=1., ZFactor=ZFactor)

        return fittedParameters

    def importC40_HOMs(self, n_elements=None, filename=None,
                       RshFactor=None, QFactor=None):
        '''
        Determines the impedance model used to represent the C40 cavities
        HOMs.
        '''

        if n_elements is None:
            n_elements = self.n_elements_C40_HOMs
        if filename is None:
            filename = self.filename_C40_HOMs
        if RshFactor is None:
            RshFactor = self.RshFactor_C40_HOMs
        if QFactor is None:
            QFactor = self.QFactor_C40_HOMs

        self.importResonatorFromFile(filename, unitFreq=1., unitRsh=1.,
                                     RshFactor=RshFactor * n_elements,
                                     QFactor=QFactor)

    def importC80(self, method=None, filename_list=None,
                  frequency_array=None,
                  RshFactor=None, QFactor=None, ZFactor=None):
        '''
        Determines the impedance model used to represent the C80 cavities
        without MHFB
        '''

        if method is None:
            method = self.method_C80
        if filename_list is None:
            filename_list = self.filename_C80
        if frequency_array is None:
            frequency_array = self.freq_array
        if RshFactor is None:
            RshFactor = self.RshFactor_C80
        if QFactor is None:
            QFactor = self.QFactor_C80
        if ZFactor is None:
            ZFactor = self.ZFactor_C80

        if isinstance(filename_list, str):
            filename_list = [filename_list]

        if method == 'Resonators':
            for filename in filename_list:
                self.importResonatorFromFile(filename, unitFreq=1., unitRsh=1.,
                                             RshFactor=RshFactor,
                                             QFactor=QFactor)
        elif method == 'DFB':
            for filename in filename_list:

                if 'C08' in filename:
                    from impedance.rf_cavities.C80.Individual.C08.DFB.C08_fit_DFB import load_model as load_C08
                    load_C80 = load_C08
                elif 'C88' in filename:
                    from impedance.rf_cavities.C80.Individual.C88.DFB.C88_fit_DFB import load_model as load_C88
                    load_C80 = load_C88
                elif 'C89' in filename:
                    from impedance.rf_cavities.C80.Individual.C89.DFB.C89_fit_DFB import load_model as load_C89
                    load_C80 = load_C89
                else:
                    continue

                loaded_impedance, freqArray = load_C80(
                    filename, folder=self.impedanceFolder,
                    freqArray=frequency_array)
                self.importInputTableFromList([freqArray,
                                               loaded_impedance.real,
                                               loaded_impedance.imag],
                                              filename,
                                              unitFreq=1, ZFactor=ZFactor)

        else:
            raise RuntimeError('Method for importC80 not recognized')

    def importC80_MHFB(self, filename_C80=None, f_rev=None,
                       impedance_reduction_target=None,
                       bandwidth_target=None,
                       main_harmonic_enableFB=None,
                       frequency_array=None,
                       manual_params_C80_08=None,
                       manual_params_C80_88=None,
                       manual_params_C80_89=None,
                       RshFactor=None, QFactor=None, ZFactor=None):
        '''
        Determines the impedance model used to represent the C80 cavities
        with MHFB
        '''

        if filename_C80 is None:
            filename_C80 = self.filename_C80_MHFB
        if f_rev is None:
            f_rev = self.f_rev
        if impedance_reduction_target is None:
            impedance_reduction_target = self.impedance_reduction_target_C80_MHFB
        if bandwidth_target is None:
            bandwidth_target = self.bandwidth_target_C80_MHFB
        if main_harmonic_enableFB is None:
            main_harmonic_enableFB = self.main_harmonic_enableFB_C80_MHFB
        if frequency_array is None:
            frequency_array = self.freq_array
        if manual_params_C80_08 is None:
            manual_params_C80_08 = self.manual_params_C80_08_MHFB
        if manual_params_C80_88 is None:
            manual_params_C80_88 = self.manual_params_C80_88_MHFB
        if manual_params_C80_89 is None:
            manual_params_C80_89 = self.manual_params_C80_89_MHFB
        if RshFactor is None:
            RshFactor = self.RshFactor_C80_MHFB
        if QFactor is None:
            QFactor = self.QFactor_C80_MHFB
        if ZFactor is None:
            ZFactor = self.ZFactor_C80_MHFB

        from impedance.rf_cavities.C80.MHFB.Resonators.All_80MHz_MHFB_Resonators import load_model

        final_impedance, frequency_array, fittedParameters = load_model(
            filename_C80, f_rev,
            impedance_reduction_target, bandwidth_target,
            main_harmonic_FB=main_harmonic_enableFB,
            frequency_array=frequency_array,
            manual_params_C80_08=manual_params_C80_08,
            manual_params_C80_88=manual_params_C80_88,
            manual_params_C80_89=manual_params_C80_89,
            RshFactor=RshFactor, QFactor=QFactor)

        self.importInputTableFromList(
            [frequency_array, final_impedance.real, final_impedance.imag],
            'rf_C80_MHFB', unitFreq=1., ZFactor=ZFactor)

        return fittedParameters

    def importC80_HOMs(self, n_elements=None, filename=None,
                       RshFactor=None, QFactor=None):
        '''
        Determines the impedance model used to represent the C80 cavities
        HOMs.
        '''

        if n_elements is None:
            n_elements = self.n_elements_C80_HOMs
        if filename is None:
            filename = self.filename_C80_HOMs
        if RshFactor is None:
            RshFactor = self.RshFactor_C80_HOMs
        if QFactor is None:
            QFactor = self.QFactor_C80_HOMs

        self.importResonatorFromFile(filename, unitFreq=1., unitRsh=1.,
                                     RshFactor=RshFactor * n_elements,
                                     QFactor=QFactor)

    def importC200(self, n_elements=None, filename=None,
                   RshFactor=None, QFactor=None):
        '''
        Determines the impedance model used to represent the C200 cavities
        main harmonic
        '''

        if n_elements is None:
            n_elements = self.n_elements_C200
        if filename is None:
            filename = self.filename_C200
        if RshFactor is None:
            RshFactor = self.RshFactor_C200
        if QFactor is None:
            QFactor = self.QFactor_C200

        self.importResonatorFromFile(filename, unitFreq=1., unitRsh=1.,
                                     RshFactor=RshFactor * n_elements,
                                     QFactor=QFactor)

    def importFinemet(self, filename=None, RshFactor=None, QFactor=None):
        '''
        Determines the impedance model used to represent the Finemet cavity.
        '''

        if filename is None:
            filename = self.filename_Finemet
        if RshFactor is None:
            RshFactor = self.RshFactor_Finemet
        if QFactor is None:
            QFactor = self.QFactor_Finemet

        self.importResonatorFromFile(filename, unitFreq=1., unitRsh=1.,
                                     RshFactor=RshFactor, QFactor=QFactor)

    def importFinemet_MHFB(self, filename_Finemet=None, f_rev=None,
                           impedance_reduction_target=None,
                           bandwidth_target=None,
                           frequency_array=None,
                           manual_params_finemet=None,
                           RshFactor=None, QFactor=None, ZFactor=None):
        '''
        Determines the impedance model used to represent the Finemet cavity
        with MHFB.
        '''

        if filename_Finemet is None:
            filename_Finemet = self.filename_Finemet
        if f_rev is None:
            f_rev = self.f_rev
        if impedance_reduction_target is None:
            impedance_reduction_target = self.impedance_reduction_target_Finemet_MHFB
        if bandwidth_target is None:
            bandwidth_target = self.bandwidth_target_Finemet_MHFB
        if frequency_array is None:
            frequency_array = self.freq_array
        if manual_params_finemet is None:
            manual_params_finemet = self.manual_params_Finemet_MHFB
        if RshFactor is None:
            RshFactor = self.RshFactor_Finemet_MHFB
        if QFactor is None:
            QFactor = self.QFactor_Finemet_MHFB
        if ZFactor is None:
            ZFactor = self.ZFactor_Finemet_MHFB

        from impedance.rf_cavities.Finemet.MHFB.Resonators.Finemet_MHFB_Resonators import load_model

        final_impedance, frequency_array, fittedParameters = load_model(
            filename_Finemet, f_rev,
            impedance_reduction_target, bandwidth_target,
            frequency_array=frequency_array,
            manual_params_finemet=manual_params_finemet,
            RshFactor=RshFactor, QFactor=QFactor)

        self.importInputTableFromList(
            [frequency_array, final_impedance.real, final_impedance.imag],
            'rf_Finemet_MHFB', unitFreq=1., ZFactor=ZFactor)

        return fittedParameters

    def importKickers(self, filename_list=None, RshFactor=None, QFactor=None):
        '''
        Determines the impedance model used to represent the Kickers
        '''

        if filename_list is None:
            filename_list = self.filename_Kickers
        if RshFactor is None:
            RshFactor = self.RshFactor_Kickers
        if QFactor is None:
            QFactor = self.QFactor_Kickers

        if isinstance(filename_list, str):
            filename_list = [filename_list]

        for filename in filename_list:
            self.importResonatorFromFile(filename, unitFreq=1., unitRsh=1.,
                                         RshFactor=RshFactor,
                                         QFactor=QFactor)

    def importSepta(self, filename_list=None, RshFactor=None, QFactor=None):
        '''
        Determines the impedance model used to represent the Septa
        '''

        if filename_list is None:
            filename_list = self.filename_Septa
        if RshFactor is None:
            RshFactor = self.RshFactor_Septa
        if QFactor is None:
            QFactor = self.QFactor_Septa

        if isinstance(filename_list, str):
            filename_list = [filename_list]

        for filename in filename_list:
            loaded_data = np.loadtxt(self.impedanceFolder + '/' + filename)

            im_z_over_n = np.sum(loaded_data[:, -1])

            self.importResonatorsFromList(
                [loaded_data[:, 0],
                 loaded_data[:, 1],
                 loaded_data[:, 2]], filename,
                unitFreq=1., unitRsh=1.,
                RshFactor=RshFactor,
                QFactor=QFactor)

            self.importImZ_over_f(filename + '_ImZ/f',
                                  im_z_over_n * RshFactor)

    def importSectorValves(self, n_elements=None, filename_list=None,
                           RshFactor=None, QFactor=None):
        '''
        Determines the impedance model used to represent the sector valve
        '''

        if n_elements is None:
            n_elements = self.n_elements_SectorValves
        if filename_list is None:
            filename_list = self.filename_SectorValves
        if RshFactor is None:
            RshFactor = self.RshFactor_SectorValves
        if QFactor is None:
            QFactor = self.QFactor_SectorValves

        for filename in filename_list:
            self.importEigenFromCST(filename,
                                    RshFactor=RshFactor * n_elements,
                                    QFactor=QFactor)

    def importMU_upstream(self, n_elements=None, filename_list=None,
                          ZFactor=None):
        '''
        Determines the impedance model used to represent the MU section
        (upstream), now represented as bellows only
        '''

        if n_elements is None:
            n_elements = self.n_elements_MU_upstream
        if filename_list is None:
            filename_list = self.filename_MU_upstream
        if ZFactor is None:
            ZFactor = self.RshFactor_MU_upstream

        for filename in filename_list:
            loaded_data = np.loadtxt(self.impedanceFolder + '/' + filename)

            im_z_over_n = loaded_data

            self.importImZ_over_f(filename + '_ImZ/f',
                                  im_z_over_n * ZFactor * n_elements)

    def importMU_downstream(self, n_elements_list=None, filename_list=None,
                            RshFactor=None, QFactor=None):
        '''
        Determines the impedance model used to represent the MU section
        (downstream), now represented as pumping manifolds, empty or with
        CODD
        '''

        if n_elements_list is None:
            n_elements_list = self.n_elements_MU_downstream
        if filename_list is None:
            filename_list = self.filename_MU_downstream
        if RshFactor is None:
            RshFactor = self.RshFactor_MU_downstream
        if QFactor is None:
            QFactor = self.QFactor_MU_downstream

        if isinstance(n_elements_list, (int, float)):
            n_elements_list = [n_elements_list]
        if isinstance(filename_list, str):
            filename_list = [filename_list]
        if len(n_elements_list) != len(filename_list):
            raise RuntimeError(
                'n_elements_list should be same length as filename_list')

        for index, filename in enumerate(filename_list):
            loaded_data = np.loadtxt(self.impedanceFolder + '/' + filename)

            im_z_over_n = np.sum(loaded_data[:, -1])

            self.importResonatorsFromList(
                [loaded_data[:, 0],
                 loaded_data[:, 1],
                 loaded_data[:, 2]], filename,
                unitFreq=1., unitRsh=1.,
                RshFactor=RshFactor * n_elements_list[index],
                QFactor=QFactor)

            self.importImZ_over_f(filename + '_ImZ/f',
                                  im_z_over_n * RshFactor * n_elements_list[index])

    def importFlanges(self, n_elements_list=None, filename_list=None,
                      RshFactor=None, QFactor=None):
        '''
        Model of vacuum flanges (high frequency resonance from the gap)
        '''

        if n_elements_list is None:
            n_elements_list = self.n_elements_Flanges
        if filename_list is None:
            filename_list = self.filename_Flanges
        if RshFactor is None:
            RshFactor = self.RshFactor_Flanges
        if QFactor is None:
            QFactor = self.QFactor_Flanges

        if isinstance(n_elements_list, (int, float)):
            n_elements_list = [n_elements_list]
        if isinstance(filename_list, str):
            filename_list = [filename_list]
        if len(n_elements_list) != len(filename_list):
            raise RuntimeError(
                'n_elements_list should be same length as filename_list')

        for index, filename in enumerate(filename_list):

            self.importEigenFromCST(filename,
                                    RshFactor=RshFactor *
                                    n_elements_list[index],
                                    QFactor=QFactor,
                                    unitFreq=1e9)

    def importSteps(self, filename=None, ZFactor=None):
        '''
        Model of step transitions
        '''

        if filename is None:
            filename = self.filename_Steps
        if ZFactor is None:
            ZFactor = self.ZFactor_Steps

        loaded_data = np.loadtxt(self.impedanceFolder + '/' + filename)

        im_z_over_n = loaded_data

        self.importImZ_over_f(
            filename, im_z_over_n * ZFactor)

    def importFlanges_GroundLoops(self, filename=None,
                                  n_elements=None,
                                  n_bypasses=None,
                                  frequency_array=None,
                                  ZFactor=None):
        '''
        Determines the impedance model used to represent the ground
        loops around the vacuum flanges
        '''

        if filename is None:
            filename = self.filename_Flanges_GroundLoops
        if n_elements is None:
            n_elements = self.n_elements_Flanges_GroundLoops
        if n_bypasses is None:
            n_bypasses = self.n_bypasses_Flanges_GroundLoops
        if frequency_array is None:
            frequency_array = self.freq_array
        if ZFactor is None:
            ZFactor = self.ZFactor_Flanges_GroundLoops

        if filename.split('.')[-1] == 'py':
            from impedance.vacuum.flanges.ground_loops.circuit_model.circuit_model import load_model

            if n_bypasses > n_elements:
                raise RuntimeError('More bypasses than number of flanges')

            impedance_bypass, freq_array = load_model(
                frequency_array, n_bypasses=1)

            self.importInputTableFromList(
                [freq_array, impedance_bypass.real, impedance_bypass.imag],
                'rf_Flanges_GroundLoops_Withbypass', unitFreq=1.,
                ZFactor=ZFactor * n_bypasses)

            if (n_elements - n_bypasses) > 0:
                impedance_nobypass, _ = load_model(
                    frequency_array, n_bypasses=0)

                self.importInputTableFromList(
                    [freq_array, impedance_nobypass.real, impedance_nobypass.imag],
                    'rf_Flanges_GroundLoops_Nobypass', unitFreq=1.,
                    ZFactor=ZFactor * (n_elements - n_bypasses))

        else:

            self.importResonatorFromFile(filename, unitFreq=1., unitRsh=1.,
                                         RshFactor=ZFactor * n_bypasses)

            if (n_elements - n_bypasses) > 0:

                self.importResonatorFromFile(
                    '/vacuum/flanges/ground_loops/Resonators/multi_resonator_0_bypass.txt',
                    unitFreq=1., unitRsh=1.,
                    RshFactor=ZFactor * (n_elements - n_bypasses))

    def importBI(self, n_elements_list=None, filename_list=None,
                 frequency_array=None,
                 RshFactor=None, QFactor=None):
        '''
        Determines the impedance model used to represent the beam
        instrumentation devices (except CODD which is handled separately)
        '''

        if n_elements_list is None:
            n_elements_list = self.n_elements_BI
        if filename_list is None:
            filename_list = self.filename_BI
        if frequency_array is None:
            frequency_array = self.freq_array
        if RshFactor is None:
            RshFactor = self.RshFactor_BI
        if QFactor is None:
            QFactor = self.QFactor_BI

        if isinstance(n_elements_list, (int, float)):
            n_elements_list = [n_elements_list]
        if isinstance(filename_list, str):
            filename_list = [filename_list]
        if len(n_elements_list) != len(filename_list):
            raise RuntimeError(
                'n_elements_list should be same length as filename_list')

        for index, filename in enumerate(filename_list):

            if filename.split('.')[-1] == 'py':

                load_model = __import__(
                    'impedance' + filename.replace('/', '.')[:-3], fromlist=['']).load_model

                loaded_impedance, freq_array = load_model(frequency_array)

                self.importInputTableFromList(
                    [freq_array, loaded_impedance.real, loaded_impedance.imag],
                    filename, unitFreq=1.,
                    ZFactor=RshFactor * n_elements_list[index])

            else:
                if filename.split('.')[-1] == 'dat':
                    delimiter = '\t'
                    unitFreq = 1e9
                    unitRsh = 1e3
                else:
                    delimiter = ' '
                    unitFreq = 1.
                    unitRsh = 1.
                loaded_data = np.array(
                    np.loadtxt(self.impedanceFolder + '/' + filename,
                               delimiter=delimiter,
                               comments=['#', '!']), ndmin=2)

                self.importResonatorsFromList(
                    [loaded_data[:, 0],
                     loaded_data[:, 1],
                     loaded_data[:, 2]], filename,
                    unitFreq=unitFreq, unitRsh=unitRsh,
                    RshFactor=RshFactor * n_elements_list[index],
                    QFactor=QFactor)

                if loaded_data.shape[1] == 4:
                    im_z_over_n = np.sum(loaded_data[:, -1])

                    self.importImZ_over_f(filename + '_ImZ/f',
                                          im_z_over_n * RshFactor * n_elements_list[index])

    def importMisc(self, n_elements_list=None, filename_list=None,
                   RshFactor=None, QFactor=None):
        '''
        Determines the impedance model used to represent the miscellaneous
        equipment (dump, TFB, ralentisseur...)
        '''

        if n_elements_list is None:
            n_elements_list = self.n_elements_Misc
        if filename_list is None:
            filename_list = self.filename_Misc
        if RshFactor is None:
            RshFactor = self.RshFactor_Misc
        if QFactor is None:
            QFactor = self.QFactor_Misc

        if isinstance(n_elements_list, (int, float)):
            n_elements_list = [n_elements_list]
        if isinstance(filename_list, str):
            filename_list = [filename_list]
        if len(n_elements_list) != len(filename_list):
            raise RuntimeError(
                'n_elements_list should be same length as filename_list')

        for index, filename in enumerate(filename_list):
            loaded_data = np.array(
                np.loadtxt(self.impedanceFolder + '/' + filename), ndmin=2)

            self.importResonatorsFromList(
                [loaded_data[:, 0],
                 loaded_data[:, 1],
                 loaded_data[:, 2]], filename,
                unitFreq=1., unitRsh=1.,
                RshFactor=RshFactor * n_elements_list[index],
                QFactor=QFactor)

            if loaded_data.shape[1] == 4:
                im_z_over_n = np.sum(loaded_data[:, -1])

                self.importImZ_over_f(filename + '_ImZ/f',
                                      im_z_over_n * RshFactor * n_elements_list[index])

    def importResistiveWall(self, frequency_array=None, ZFactor=None):
        '''
        Load the resistive wall impedance
        '''

        if frequency_array is None:
            frequency_array = self.freq_array
        if ZFactor is None:
            ZFactor = self.ZFactor_ResistiveWall

        from impedance.resistive_wall.resistive_wall_model import load_model as load_RW

        resistive_wall_impedance = load_RW(frequency_array)

        self.importInputTableFromList([frequency_array,
                                       resistive_wall_impedance.real,
                                       resistive_wall_impedance.imag],
                                      'resistive_wall',
                                      unitFreq=1.,
                                      ZFactor=ZFactor)

    def importSpaceCharge(self, emittance_x_norm=None, emittance_y_norm=None,
                          mass=None, momentum=None, momentum_spread=None,
                          method=None,
                          aperture_X=None, aperture_Y=None,
                          optics_S=None, beta_X=None,
                          beta_Y=None, disp_X=None,
                          reinterpolate=None,
                          return_diagnose=None,
                          BI_device_x=None,
                          BI_device_y=None,
                          ZFactor=None):
        '''
        Returns ImZ/n for space charge
        '''

        if emittance_x_norm is None:
            emittance_x_norm = self.emittance_x_norm_SC
        if emittance_y_norm is None:
            emittance_y_norm = self.emittance_y_norm_SC
        if mass is None:
            mass = self.mass
        if momentum is None:
            momentum = self.momentum
        if momentum_spread is None:
            momentum_spread = self.momentum_spread_SC
        if method is None:
            method = self.method_SC
        if aperture_X is None:
            aperture_X = self.aperture_X_SC
        if aperture_Y is None:
            aperture_Y = self.aperture_Y_SC
        if optics_S is None:
            optics_S = self.optics_S_SC
        if beta_X is None:
            beta_X = self.beta_X_SC
        if beta_Y is None:
            beta_Y = self.beta_Y_SC
        if disp_X is None:
            disp_X = self.disp_X_SC
        if reinterpolate is None:
            reinterpolate = self.reinterpolate_SC
        if return_diagnose is None:
            return_diagnose = self.return_diagnose_SC
        if BI_device_x is None:
            BI_device_x = self.BI_device_x_SC
        if return_diagnose is None:
            BI_device_y = self.BI_device_y_SC
        if ZFactor is None:
            ZFactor = self.ZFactor_SC

        from impedance.space_charge.space_charge_model import load_model

        space_charge_model = load_model(
            emittance_x_norm, emittance_y_norm, mass, momentum,
            momentum_spread, method,
            aperture_X=aperture_X, aperture_Y=aperture_Y,
            optics_S=optics_S, beta_X=beta_X,
            beta_Y=beta_Y, disp_X=disp_X,
            reinterpolate=reinterpolate,
            return_diagnose=return_diagnose,
            BI_device_x=BI_device_x,
            BI_device_y=BI_device_y)

        if return_diagnose:
            im_z_over_n = space_charge_model[0]
        else:
            im_z_over_n = space_charge_model

        self.importImZ_over_f('space_charge',
                              im_z_over_n / self. f_rev * ZFactor)

        return space_charge_model

    def importImpedancePS(self):

        # RF systems
        if self.enable_C10_1TFB:
            self.importC10_1TFB()
        else:
            self.importC10()

        self.importC10_ClosedGapRelay()

        if self.enable_C20_MHFB:
            self.importC20_MHFB()
        else:
            self.importC20()

        if self.enable_C40_MHFB:
            self.importC40_MHFB()
        else:
            self.importC40()

        self.importC40_HOMs()

        if self.enable_C80_MHFB:
            self.importC80_MHFB()
        else:
            self.importC80()

        self.importC80_HOMs()

        self.importC200()

        self.importFinemet()

        # Kickers and Setpa
        self.importKickers()

        self.importSepta()

        # Vacuum equipment
        self.importSectorValves()

        self.importMU_upstream()

        self.importMU_downstream()

        self.importFlanges()

        self.importSteps()

        self.importFlanges_GroundLoops()

        # Beam instrumentation
        self.importBI()

        # Misc equipment
        self.importMisc()

        # Resistive wall
        self.importResistiveWall()

        # Space charge
        if self.load_SC:
            self.importSpaceCharge()

    def export2BLonD(self):

        if hasattr(self, 'table_impedance'):
            return impedance2blond(self.table_impedance)
        else:
            return 0
