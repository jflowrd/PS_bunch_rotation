# coding: utf-8
'''
Author: A. Lasheen
PS non adiabatic bunch shortening
Single bunch
No intensity effects
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

from PS_bunch_rotation_generation import C40_77_prog, C40_78_prog, C80_08_prog, C80_88_prog, C80_phase_prog, C40_phase_prog

def run(output_dir='/eos/user/j/jflowerd/Work/', n_bunches = 1, n_macroparticles_per_bunch = 1e5):
    # BLonD common imports
    #sys.path.append('../../pymodules/')
    from blond_common.fitting.profile import gaussian_fit

    # PS RF program imports
    #sys.path.append('../pymodules/')


    # %% Simulation settings

    #case_name = 'first_test'
    case_name = 'HT_condor_test'
    results_folder = output_dir + 'results/' + case_name

    plot_turns = False
    plot_animation = False
    plot_save_turns = False

    ''' Plotting parameters '''
    dt_plt = 10 # <----- Number of turns between phase space plots
    t_offset = 1e-9
    E_min = -100e6
    E_max = 100e6

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
    length_step_h168 = 133                  # First Step time [us]
    amplitude_step_h168 = 0.                # Lin. step amplitude [kV]
    C40_77_prog_factor = 284.700 / 300      # Calibration factor C40-77
    C40_78_prog_factor = 1.                 # Calibration factor C40-78
    C80_08_prog_factor = 1.                 # Calibration factor C80-08
    C80_88_prog_factor = 1.                 # Calibration factor C80-88
    C80_89_prog_factor = 0.                 # Calibration factor C80-89
    intial_phase_h84 = 0
    final_phase_h84 = 0
    initial_phase_h168 = 0
    step_adjust_h168 = 0
    final_phase_h168 = initial_phase_h168 - np.pi + step_adjust_h168

    # Beam parameters
    #n_bunches = 1
    #n_macroparticles_per_bunch = 1e5 #1e6
    intensity_per_bunch = 1.3e11
    intensity = intensity_per_bunch * n_bunches
    n_macroparticles = int(n_macroparticles_per_bunch * n_bunches)
    adjusted_4sig = 11.4e-9
    adjusted_full_bl = 14.6e-9
    adjusted_exponent = ((2 * adjusted_full_bl / adjusted_4sig)**2. - 3) / 2. - 0.5

    # Profile parameters
    n_bins_per_bunch_ps = 2**8
    bunch_spacing_buckets_ps = 1

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
    number_iterations_for_extraction = 1
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
            time_shortening = 350e-6
            n_turns_ps = int(round(time_shortening / ring_ps.t_rev[0]))
        else:
            n_turns_ps = extraction_turn[0] + 1

        ring_ps = Ring(circumference_ps, momentum_compaction_ps, momentum,
                    particle_type, n_turns_ps)

        # RF parameters
        start_bunch_rotation_h84 = 40e-6

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
                                    start_bunch_rotation_h84 *
                                    1e6,
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
        cut_opt = CutOptions(n_slices=n_bins_per_bunch_ps,
                            cut_left=0, cut_right=rf_params_ps.t_rf[0, 0],
                            cuts_unit='s')
        profile_ps = Profile(beam_ps, CutOptions=cut_opt)

        # RF tracker
        longitudinal_tracker_ps = RingAndRFTracker(rf_params_ps, beam_ps, interpolation=False,
                                                Profile=profile_ps, TotalInducedVoltage=None)
        full_tracker_ps = FullRingAndRF([longitudinal_tracker_ps])

        if idx_sim == 0:
            distribution_options = {'type': 'binomial', 'exponent': adjusted_exponent,
                                    'bunch_length': adjusted_full_bl,
                                    'bunch_length_fit': 'full',
                                    'density_variable': 'Hamiltonian'}

            match_beam_from_distribution(beam_ps, full_tracker_ps, ring_ps,
                                        distribution_options, n_bunches,
                                        bunch_spacing_buckets_ps,
                                        TotalInducedVoltage=None,
                                        n_iterations=1,
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

        # Saving beam parameters
        profile_ps.track()
        fwhm_save = np.zeros((n_bunches, ring_ps.n_turns + 1))
        profile_ps.fwhm()
        fwhm_save[0, 0] = profile_ps.bunchLength

        gaussian_4sig_save = np.zeros((n_bunches, ring_ps.n_turns + 1))
        gaussian_4sig_save[0, 0] = gaussian_fit(
            profile_ps.bin_centers,
            profile_ps.n_macroparticles)[-1] * 4

        # Plots
        plots = Plot(ring_ps, rf_params_ps, beam_ps, dt_plt, ring_ps.n_turns,
                0, (n_bunches)*rf_params_ps.t_rf[0, 0], # <----- dt limits of the phase space plot
                E_min, E_max, # <---------------- dE limits of the phase space plot
                show_plots=True,
                separatrix_plot=True)

        # Tracking
        for turn in range(ring_ps.n_turns):

            full_tracker_ps.track()
            profile_ps.track()

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

        extraction_turn = np.where(
            np.nanmean(gaussian_4sig_save, axis=0) ==
            np.nanmin(np.nanmean(gaussian_4sig_save, axis=0)))[0]
        

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
    sample_size = 10000  # Adjust the sample size as needed
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
    plt.contourf(x_grid, y_grid, density, levels=100, cmap='gist_heat_r')
    plt.colorbar(label='Density')
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
        fig = plt.figure() 
        
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
                            interval = 1, #dt_plt
                            blit = True) 
        
        anim.save(results_folder +'/bunch_rotation.gif', 
                writer = 'pillow', fps = 2)
        
