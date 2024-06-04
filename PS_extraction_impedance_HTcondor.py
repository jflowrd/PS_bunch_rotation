# coding: utf-8
'''
Author: A. Lasheen
PS non adiabatic bunch shortening
Single bunch
With intensity effects
'''
'''
Edited by: J.Flowerdew
'''

# %% Imports
#Use env for numba backend

# General imports
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.stats import gaussian_kde

# BLonD imports
# sys.path.append('../../pymodules/blond')
from blond.input_parameters.ring import Ring
from blond.input_parameters.rf_parameters import RFStation
from blond.beam.beam import Proton, Beam
from blond.beam.profile import Profile, CutOptions
from blond.trackers.tracker import RingAndRFTracker, FullRingAndRF
from blond.beam.distributions_multibunch import match_beam_from_distribution
from blond.plots.plot import Plot
from blond.impedances.impedance import InductiveImpedance, InducedVoltageTime, InducedVoltageFreq, TotalInducedVoltage
from blond.impedances.impedance_sources import Resonators
from blond.impedances import induced_voltage_analytical
from PS_bunch_rotation_generation import C40_77_prog, C40_78_prog, C80_08_prog, C80_88_prog, C80_phase_prog, C40_phase_prog

from blond_common.fitting.profile import gaussian_fit

# Other imports
from PS_master.plots.mpl import uber_plot, uber_imshow

# Toolbox import
from PS_master.impedance_toolbox.impedance_toolbox.machine_params import MachineParameters
from PS_master.impedance_toolbox.impedance_toolbox.impedance_params import ImpedanceParameters

# Scenario import
from PS_master.impedance_loader import ImpedanceLoader

plot_turns = False
plot_animation = False
plot_save_turns = False
plot_final_turn = True
plot_voltage_ramp = True
plot_impeadance = False
plot_induced_voltage = True

multi_turn_wake = False

def run(output_dir='/eos/user/j/jflowerd/Work/', 
        n_bunches = 1, n_macroparticles_per_bunch = 1e6, length_step_h168 = 133, start_bunch_rotation_h84 = 40 ):

    case_name = 'PS_impedance_HTcondor'
    results_folder = output_dir + 'results/' + case_name

    ''' Plotting parameters '''
    dt_plt = 100 # <----- Number of turns between phase space plots
    E_min = -100e6
    E_max = 100e6

    start_bunch_rotation_h84 = start_bunch_rotation_h84*1e-6

    # %% BLonD setting up

    # Defining PS parameters
    # Ring
    particle_type = Proton()
    momentum = 25.92e9                                      # Momentum [eV/c]
    circumference_ps = 2. * np.pi * 100.                    # Circumference [m]
    gamma_transition_ps = 6.1                               # Transition gamma
    momentum_compaction_ps = 1. / gamma_transition_ps**2    # Momentum compaction

    # RF parameters PS
    n_rf_systems_ps = 2                     # Number of rf systems
    harmonic_numbers_ps = [84, 168]         # Harmonic numbers
    #length_step_h168 = 133                  # First Step time [us] (old value: 133)
    amplitude_step_h168 = 0.                # Lin. step amplitude [kV]
    C40_77_prog_factor = 271.34/300         # Calibration factor C40-77 (old value: 284.700 / 300)
    C40_78_prog_factor = 231.46/300         # Calibration factor C40-78 (old value: 1.)
    C80_08_prog_factor = 0                  # Calibration factor C80-08 (old value: 1.) # 279.14/300
    C80_88_prog_factor = 311.56/300         # Calibration factor C80-88 (old value: 1.)
    C80_89_prog_factor = 277.58/300         # Calibration factor C80-89 (old value: 0.) # 277.58/300
    intial_phase_h84 = 0
    final_phase_h84 = 0
    initial_phase_h168 = 0
    step_adjust_h168 = 0
    final_phase_h168 = initial_phase_h168 - np.pi + step_adjust_h168

    # Beam parameters
    n_bunches = 1
    n_macroparticles_per_bunch = 1e6 #1e6
    intensity_per_bunch = 2.6e11
    intensity = intensity_per_bunch * n_bunches
    n_macroparticles = int(n_macroparticles_per_bunch * n_bunches)
    adjusted_4sig = 11.4e-9
    adjusted_full_bl = 14.6e-9
    adjusted_exponent = ((2 * adjusted_full_bl / adjusted_4sig)**2. - 3) / 2. - 0.5

    # Profile parameters
    n_bins_per_bunch_ps = 2**8
    bunch_spacing_buckets_ps = 1
    n_bins_ps = int(n_bins_per_bunch_ps *
                    (bunch_spacing_buckets_ps * (harmonic_numbers_ps[0] - 1) + 1))

    # Impedance
    filter_front_wake = 0.5
    n_turns_memory = 100
    gene_iterations = 8

    # %% Building BLonD objects and running the simulation

    # Creating output folders
    try:
        os.mkdir('results')
    except FileExistsError:
        pass

    try:
        os.mkdir(results_folder)
    except FileExistsError:
        pass

    # Iterate the simulation twice to extract when the beam is the shortest
    number_iterations_for_extraction = 2
    extraction_turn = -1
    idx_sim = 0

    data_dt = []
    data_dE = []
    for idx_sim in range(number_iterations_for_extraction):

        # Ring
        ring_ps = Ring(circumference_ps, momentum_compaction_ps,
                    momentum, particle_type)
        #print('extraction_turn:', extraction_turn)

        if extraction_turn == -1:
            time_shortening = 300e-6 #1500e-6 for filamentation
            n_turns_ps = int(round(time_shortening / ring_ps.t_rev[0]))
        else:
            n_turns_ps = extraction_turn[0] + 1

        ring_ps = Ring(circumference_ps, momentum_compaction_ps, momentum,
                    particle_type, n_turns_ps)

        # RF parameters
        # start_bunch_rotation_h84 = 40e-6

        C40_77_RF_program = 1e3 * C40_77_prog(ring_ps.cycle_time * 1e6,
                                            start_bunch_rotation_h84 * 1e6,
                                            C40_77_prog_factor)

        C40_78_RF_program = 1e3 * C40_78_prog(ring_ps.cycle_time * 1e6,
                                            start_bunch_rotation_h84 * 1e6,
                                            C40_78_prog_factor)

        rf_prog_h84 = C40_77_RF_program + C40_78_RF_program

        phi_prog_h84 = C40_phase_prog(ring_ps.cycle_time * 1e6,
                                    start_bunch_rotation_h84 * 1e6,
                                    intial_phase_h84, final_phase_h84)

        C80_08_RF_program = 1e3 * C80_08_prog(ring_ps.cycle_time * 1e6,
                                            start_bunch_rotation_h84 * 1e6,
                                            length_step_h168,
                                            amplitude_step_h168,
                                            C80_08_prog_factor)

        C80_88_RF_program = 1e3 * C80_88_prog(ring_ps.cycle_time * 1e6,
                                            start_bunch_rotation_h84 * 1e6,
                                            length_step_h168,
                                            amplitude_step_h168,
                                            C80_88_prog_factor)

        C80_89_RF_program = 1e3 * C80_88_prog(ring_ps.cycle_time * 1e6,
                                            start_bunch_rotation_h84 * 1e6,
                                            length_step_h168,
                                            amplitude_step_h168,
                                            C80_89_prog_factor)

        rf_prog_h168 = C80_08_RF_program + C80_88_RF_program + C80_89_RF_program

        phi_prog_h168 = C80_phase_prog(ring_ps.cycle_time * 1e6,
                                    start_bunch_rotation_h84 * 1e6,
                                    length_step_h168,
                                    initial_phase_h168,
                                    final_phase_h168)

        rf_params_ps = RFStation(ring_ps, harmonic_numbers_ps,
                                [rf_prog_h84, rf_prog_h168],
                                [phi_prog_h84, phi_prog_h168],
                                n_rf_systems_ps)

        omega_prog_h168 = np.diff(phi_prog_h168) / \
            rf_params_ps.t_rev[0] + rf_params_ps.omega_rf[1, 0]
        omega_prog_h168 = np.append(omega_prog_h168, omega_prog_h168[-1])

        rf_params_ps.omega_rf[1, :] = omega_prog_h168
        rf_params_ps.omega_rf_d[1, :] = omega_prog_h168

        # Beam
        beam_ps = Beam(ring_ps, n_macroparticles, intensity)


        # Profile
        # cut_opt = CutOptions(n_slices=n_bins_per_bunch_ps,
        #                      cut_left=0, cut_right=rf_params_ps.t_rf[0, 0],
        #                      cuts_unit='s')
        cut_opt = CutOptions(n_slices=n_bins_ps,
                            cut_left=0, cut_right=ring_ps.t_rev[0],
                            cuts_unit='s')
        profile_ps = Profile(beam_ps, CutOptions=cut_opt)

        # Impedance
        script_path = os.path.dirname(os.path.realpath(__file__))
        toolbox_path = script_path + '/impedance_toolbox'

        # Case selector
        model = 'top'  # top, trisplit
        beam = 'LHC25ns_flat_top_after_four_split'  # LHC25ns_C2595_flat_top, LHC25ns_flat_top_after_four_split, LHC25ns_C1940_end_flat_plateau
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
        PS_loader = ImpedanceLoader(MODEL=model, f_rev=ring_ps.f_rev[0],
                                    momentum=ring_ps.momentum[0,0],
                                    folder=script_path + '/impedance',
                                    freq_array=np.linspace(0, 5.2e9, int(1e7))) # machineParams.freqArray


        PS_loader.importImpedancePS() 

        imp = PS_loader.export2BLonD()

        ResonatorsList = imp.wakeList
        ImpedanceTable_list = imp.impedanceList
        ImZ_over_f_list = imp.ImZ_over_f_List

        #Z_over_n = np.sum(np.array(ImZ_over_f_list[:-1])*ring_ps.f_rev[0]) * np.ones(ring_ps.n_turns+1)
        Z_over_n = np.zeros(ring_ps.n_turns + 1)

        if multi_turn_wake:
            frequency_step = 1/(ring_ps.t_rev[0]*n_turns_memory) # [Hz]
            front_wake_length = filter_front_wake * ring_ps.t_rev[0]*n_turns_memory
            impedance_sources = ResonatorsList+ImpedanceTable_list
        else:
            frequency_step = None
            front_wake_length = None
            impedance_sources = ResonatorsList+ImpedanceTable_list


        PS_intensity_freq = InducedVoltageFreq(beam_ps,
                                            profile_ps,
                                            impedance_sources, #[my_res_10, my_res_20, my_res_40, my_res_80], #[my_res], #ResonatorsList+ImpedanceTable_list
                                            RFParams=rf_params_ps,
                                                # use_regular_fft=False,
                                                frequency_resolution=frequency_step,
                                            multi_turn_wake=multi_turn_wake,
                                            front_wake_length=front_wake_length)


        PS_inductive = InductiveImpedance(beam_ps, profile_ps, Z_over_n, rf_params_ps, deriv_mode='gradient') # Z_over_n

        PS_longitudinal_intensity = TotalInducedVoltage(beam_ps, profile_ps, [PS_intensity_freq, PS_inductive])

        uber_plot(PS_intensity_freq.freq/1e6,
            np.abs(PS_intensity_freq.total_impedance*profile_ps.bin_size)/1e3,
            figname='Impedance',
            xlabel='Frequency [MHz]',
            ylabel='Impedance $\\|Z\\|$ [$\\mathrm{k\\Omega}$]',
            savefig=results_folder+'/impedance.png')
        plt.xlim((0.4, 1.02*np.max(PS_intensity_freq.freq)/1e6))
        plt.ylim((0.5, 1.1*np.max(np.abs(PS_intensity_freq.total_impedance*profile_ps.bin_size))/1e3))
        plt.xscale('log')
        plt.yscale('log')
        plt.tight_layout()
        plt.savefig(results_folder+'/impedance_log.png')

        # RF tracker
        longitudinal_tracker_ps = RingAndRFTracker(rf_params_ps, beam_ps, interpolation=True,
                                                Profile=profile_ps, TotalInducedVoltage=PS_longitudinal_intensity)
        full_tracker_ps = FullRingAndRF([longitudinal_tracker_ps])

        if idx_sim == 0:
            distribution_options = {'type': 'binomial', 'exponent': adjusted_exponent,
                                    'bunch_length': adjusted_full_bl,
                                    'bunch_length_fit': 'full',
                                    'density_variable': 'Hamiltonian'}

            match_beam_from_distribution(beam_ps, full_tracker_ps, ring_ps,
                                        distribution_options, n_bunches,
                                        bunch_spacing_buckets_ps,
                                        TotalInducedVoltage=None, #PS_longitudinal_intensity,
                                        n_iterations=gene_iterations,
                                        n_points_potential=int(1e3))

            saved_dt = np.array(beam_ps.dt)
            saved_dE = np.array(beam_ps.dE)

            ## PLotting initial distribution
            # plt.figure('Initial bunch distribution')
            # plt.clf()
            # plt.plot(saved_dt*1e9, saved_dE/1e6, '.')
            # plt.xlabel('$\\Delta t$ [ns]')
            # plt.ylabel('$\\Delta E$ [MeV]')
            # plt.tight_layout()
            # plt.savefig(results_folder + '/bunch_init_distribution.png')
        else:
            beam_ps.dt[:] = np.array(saved_dt)
            beam_ps.dE[:] = np.array(saved_dE)

        profile_ps.track()
        if multi_turn_wake:
            for index_gene in range(n_turns_memory):
                PS_longitudinal_intensity.induced_voltage_sum() # Update induced voltage
        else:
            PS_longitudinal_intensity.induced_voltage_sum() # Update induced voltage

        if plot_induced_voltage:
            # profile_ps.rms()
            time_array = np.linspace(0, ring_ps.t_rev[0], len(PS_longitudinal_intensity.induced_voltage))
            # analytical_induced_voltage = induced_voltage_analytical.analytical_gaussian_resonator(profile_ps.bunchLength/4, my_res_40.Q, my_res_40.R_S, my_res_40.omega_R, time_array-0.5*rf_params_ps.t_rf[0, 0], intensity) \
            # + induced_voltage_analytical.analytical_gaussian_resonator(profile_ps.bunchLength/4, my_res_80.Q, my_res_80.R_S, my_res_80.omega_R, time_array-0.5*rf_params_ps.t_rf[0, 0], intensity)
            # + induced_voltage_analytical.analytical_gaussian_resonator(profile_ps.bunchLength/4, my_res_10.Q, my_res_10.R_S, my_res_10.omega_R, time_array-0.5*rf_params_ps.t_rf[0, 0], intensity)
            # + induced_voltage_analytical.analytical_gaussian_resonator(profile_ps.bunchLength/4, my_res_20.Q, my_res_20.R_S, my_res_20.omega_R, time_array-0.5*rf_params_ps.t_rf[0, 0], intensity)
            # print('bunch length: ', profile_ps.bunchLength)
            plt.figure()
            plt.plot(time_array , PS_longitudinal_intensity.induced_voltage, label='Induced Voltage')
            # plt.plot(time_array , analytical_induced_voltage, '--', label= 'Analytical Induced Voltage') # Analyticaly calculate the induced voltage
            plt.ylabel('Induced voltage [V]')
            plt.xlabel('Time [s]')
            plt.legend(loc='best')
            plt.savefig(results_folder + '/induced_voltage.png')

            if multi_turn_wake:
                plt.figure()
                plt.plot(PS_intensity_freq.time_mtw, PS_intensity_freq.induced_voltage)
                plt.ylabel('Induced voltage [V]')
                plt.xlabel('Time [s]')
                plt.legend(loc='best')
                plt.savefig(results_folder + '/induced_voltage_mtw.png')
        
        # Saving beam parameters
        fwhm_save = np.zeros((n_bunches, ring_ps.n_turns + 1))
        profile_ps.fwhm()
        fwhm_save[0, 0] = profile_ps.bunchLength

        gaussian_4sig_save = np.zeros((n_bunches, ring_ps.n_turns + 1))
        gaussian_4sig_save[0, 0] = gaussian_fit(
            profile_ps.bin_centers,
            profile_ps.n_macroparticles)[-1] * 4

        # Plots
        # plots = Plot(ring_ps, rf_params_ps, beam_ps, dt_plt, ring_ps.n_turns,
        #         0, (n_bunches)*rf_params_ps.t_rf[0, 0], # <----- dt limits of the phase space plot
        #         E_min, E_max, # <---------------- dE limits of the phase space plot
        #         show_plots=False,
        #         separatrix_plot=True)
        

        # Tracking
        for turn in range(ring_ps.n_turns):
            # plt.plot(abs(PS_longitudinal_intensity.induced_voltage[:1000]))
            # plt.show()
            full_tracker_ps.track()
            profile_ps.track()
            PS_longitudinal_intensity.induced_voltage_sum()

            profile_ps.fwhm()
            fwhm_save[0, turn + 1] = profile_ps.bunchLength

            gaussian_4sig_save[0, turn + 1] = gaussian_fit(
                profile_ps.bin_centers,
                profile_ps.n_macroparticles)[-1] * 4

            if plot_turns and idx_sim == 0:
                plots.track() 

            if plot_animation:
                if turn%dt_plt == 0: 
                    data_dt.append(beam_ps.dt*1e9)
                    data_dE.append(beam_ps.dE/1e6)

            if plot_save_turns:
                if turn%dt_plt == 0: #or (turn>120 and turn<130):
                    plt.figure(f'Bunch distribution turn {turn}')
                    plt.clf()
                    plt.plot(beam_ps.dt*1e9, beam_ps.dE/1e6, '.')
                    plt.xlabel('$\\Delta t$ [ns]')
                    plt.ylabel('$\\Delta E$ [MeV]')
                    plt.tight_layout()
                    plt.savefig(results_folder + '/bunch_distribution_turn_' + str(turn) + '.png')

        plt.figure(f'Bunch distribution final')
        plots = Plot(ring_ps, rf_params_ps, beam_ps, 1, ring_ps.n_turns,
                0, (n_bunches)*rf_params_ps.t_rf[0, 0], # <----- dt limits of the phase space plot
                E_min, E_max, # <---------------- dE limits of the phase space plot
                show_plots=plot_final_turn,
                separatrix_plot=True)
        plt.savefig(results_folder + '/bunch_distribution_final.png')
        
        extraction_turn = np.where(
            np.nanmean(gaussian_4sig_save, axis=0) ==
            np.nanmin(np.nanmean(gaussian_4sig_save, axis=0)))[0]
        #print("Minimum bunch length: ",  np.nanmin(np.nanmean(gaussian_4sig_save, axis=0)))
        
    # profile_ps.rms()
    # time_array = np.linspace(0, ring_ps.t_rev[0], len(PS_longitudinal_intensity.induced_voltage))
    # print('bunch length: ', profile_ps.bunchLength)
    # plt.clf()
    # plt.plot(time_array , PS_longitudinal_intensity.induced_voltage, label='Induced Voltage')
    # #plt.plot(time_array , analytical_induced_voltage, '--', label= 'Analytical Induced Voltage') # Analyticaly calculate the induced voltage
    # plt.ylabel('Induced voltage [V]')
    # plt.xlabel('Time [s]')
    # plt.legend(loc='best')
    # plt.show()
        
    if plot_voltage_ramp:
        plt.figure('Voltage ramp')
        plt.clf()
        plt.plot(ring_ps.cycle_time * 1e6, C40_77_RF_program*1e-3, label='C40_77_RF_program')
        plt.plot(ring_ps.cycle_time * 1e6, C40_78_RF_program*1e-3, label='C40_78_RF_program')
        plt.plot(ring_ps.cycle_time * 1e6, C80_08_RF_program*1e-3, label='C80_08_RF_program')
        plt.plot(ring_ps.cycle_time * 1e6, C80_88_RF_program*1e-3, label='C80_88_RF_program')
        plt.plot(ring_ps.cycle_time * 1e6, C80_89_RF_program*1e-3, label='C80_89_RF_program')
        # plt.plot(ring_ps.cycle_time * 1e6, rf_prog_h84*1e-3, 'k--', label='rf_prog_h84')
        # plt.plot(ring_ps.cycle_time * 1e6, rf_prog_h168*1e-3, 'r--', label='rf_prog_h168')
        plt.xlabel('Time [us]')
        plt.ylabel('Voltage [kV]')
        plt.legend(loc='best')
        plt.tight_layout()
        plt.savefig(results_folder + '/voltage_ramp.png')


    # PLotting results
    plt.figure('Bunch lengths')
    plt.clf()
    plt.plot(fwhm_save[0, :] * 1e9, label='FWHM->4$\\sigma$')
    plt.plot(gaussian_4sig_save[0, :] * 1e9, label='Gaussian fit 4$\\sigma$')
    plt.xlabel('Turn')
    plt.ylabel('Bunch length [ns]')
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig(results_folder + '/bunch_lengths.png')


    # PLotting results
    plt.figure('Bunch distribution')
    plt.clf()
    plt.plot(beam_ps.dt*1e9, beam_ps.dE/1e6, '.')
    plt.xlabel('$\\Delta t$ [ns]')
    plt.ylabel('$\\Delta E$ [MeV]')
    plt.tight_layout()
    plt.savefig(results_folder + '/bunch_distribution.png')


    # PLotting results heat map
    x = beam_ps.dt*1e9
    y = beam_ps.dE/1e6
    #xy = np.vstack([x, y])
    #kde = gaussian_kde(xy)
    sample_size = 10000  # Adjust the sample size as needed as heat map takes a while to plot
    idx = np.random.choice(len(x), size=sample_size, replace=False)
    x_sample = x[idx]
    y_sample = y[idx]
    xy_sample = np.vstack([x_sample, y_sample])
    kde = gaussian_kde(xy_sample)
    x_grid, y_grid = np.mgrid[0:(n_bunches)*rf_params_ps.t_rf[0, 0]/1e-9:100j, y.min():y.max():100j]
    positions = np.vstack([x_grid.ravel(), y_grid.ravel()])
    density = np.reshape(kde(positions).T, x_grid.shape)
    plt.figure('Bunch distribution: heat map')
    plt.clf()
    #plt.pcolormesh(x_grid, y_grid, density, shading='auto', cmap='gist_heat_r')
    plt.contourf(x_grid, y_grid, density, levels=100, cmap='turbo')
    plt.colorbar(label='Density')
    plt.xlim(10,15)
    plt.xlabel('$\\Delta t$ [ns]')
    plt.ylabel('$\\Delta E$ [MeV]')
    plt.tight_layout()
    plt.savefig(results_folder + '/bunch_distribution_density.png')


    # Saving results
    np.savez(results_folder + '/bunch_lengths.npz',
            fwhm_save=fwhm_save,
            gaussian_4sig_save=gaussian_4sig_save,
            cycle_time=ring_ps.cycle_time)

    np.savez(results_folder + '/bunch_distribution.npz',
            beam_dt=beam_ps.dt,
            beam_dE=beam_ps.dE)

    #Plotting animation
    if plot_animation:
        fig = plt.figure('Animation') 
        
        # marking the x-axis and y-axis 
        axis = plt.axes(xlim=(5, 20),  ylim=(-80, 80), xlabel = '$\\Delta t$ [ns]', ylabel ='$\\Delta E$ [MeV]' ) 
        
        # initializing a line variable 
        scatter_plot, = axis.plot([], [], linestyle="",  marker='.') 
        
        # data which the line will 
        # contain (x, y) 
        def init(): 
            scatter_plot.set_data([], []) 
            return scatter_plot, 
        
        def animate(i): 
            x = data_dt[i]
            y = data_dE[i]
            scatter_plot.set_data(x, y) 
            return scatter_plot, 
        
        anim = FuncAnimation(fig, animate, 
                            init_func = init, 
                            frames = len(data_dt), 
                            interval = dt_plt,
                            blit = True) 
        
        anim.save(results_folder +'/bunch_rotation.gif', 
                writer = 'pillow', fps = 3)
        

        
    #Plotting impedance
    if plot_impeadance:

        # # Loading impedance
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
        plt.savefig(fig_path+'/impedance_per_source.png')

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
        plt.savefig(fig_path+'/spectrum_and_impedance.png')

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
        plt.savefig(fig_path+'/induced_voltage.png')
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
        plt.savefig(fig_path+'/imZ_over_n.png')

