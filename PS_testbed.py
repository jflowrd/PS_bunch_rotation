import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import m_p, c, e
import blond
#blond.test()
from blond.input_parameters.ring import Ring
from blond.input_parameters.rf_parameters import RFStation
from blond.beam.beam import Proton, Beam
from blond.beam.distributions import parabolic, bigaussian
from blond.beam.profile import Profile, CutOptions, FitOptions
from blond.trackers.tracker import RingAndRFTracker, FullRingAndRF
from blond.plots.plot import Plot
from blond.monitors.monitors import BunchMonitor
from support_functions import plot_phase_space_trajectory, plot_phase_space_distribution

single_particle = False               #Track a single particle
single_particle_track = True          #Track multiple single particle
bunch_plot = False                    #Plot a bunch after a certain number of turns (n_turns)
bunch_tracking = True                 #Track a bunch for n_turns plotting every dt_plt turns
bunch_tracking_acc = False            #Track an accelerating bunch for n_turns plotting every dt_plt turns
bunch_parameters = False               #Plot bunch parameters
bunch_monitor = False                #Plot bunch monitor data

''' Tracking parameters '''
n_particles = 30            #Single particle tracking 
n_macroparticles = 1e6    #Bunch tracking 
n_turns = 1000              #Number of turns tracked for
intensity = 2.5e11          #Intensity of bunch

energy_offset = 1e7        #Energy offset in bucket [eV]
energy_spread = 0        #Energy spread in bucket [eV]
bunch_len_factor = 1     #Factor to change length of bunch
bunch_offset = 0e-9         #Timing offset to bunch [s]

''' Acceleration parameters '''
acc_time = 50e-3           #Acceleration time [s]
momentum_min = 26e9         #Injection mometum [eV/c]
momentum_max = 27e9         #Extraction momentum [eV/c]
momentum_time = (0, acc_time) # <----------------------------------- Time values of the momentum program
momentum_program = (momentum_min, momentum_max) # <------------------------------- Momentum values 
momentum_program_input = (momentum_time, momentum_program) # <---- Declare momentum program as a tuple

''' Plotting parameters '''
dt_plt = 500 # <----- Number of turns between phase space plots
t_offset = 1e-9
E_min = -100e6
E_max = 100e6

''' Machine parameters '''
Ekin = 26e9 # [eV]

charge = 1
E0 = m_p * c**2. / e
circumference_PS = 2. * np.pi * 100.  # Machine circumference [m]
circumference = circumference_PS
energy = Ekin + E0 # [eV]
momentum = np.sqrt(energy**2. - E0**2.) # [eV/c]
momentum = 25.92e9      # Momentum [GeV/c]
beta = momentum / energy
gamma = energy / E0

t_rev = circumference / (beta*c)
f_rev = 1 / t_rev
w_0 = 2*np.pi*f_rev #revolution frequency
print(f_rev)

# RF parameters PS
n_rf_systems_PS = 2                        # Number of rf systems
harmonic_numbers_PS = [84, 168]                 # Harmonic numbers
phi_prog_h84 = 0

h = 84 #####################################################################################
voltage = 4.5e6 # [V]
f_rf = h*f_rev
t_rf = 1/ f_rf
phi_s = 0
phase = phi_s
w_r = h*w_0

gamma_t = 6.1
alpha_c = 1 / gamma_t**2.
eta = alpha_c - 1 / gamma**2. 

Qs = np.sqrt((charge*voltage*h*eta*np.cos(phi_s))/(2*np.pi*beta**2 * energy)) #Note this should have a minus sign withing square root
w_s0 = Qs*w_0 # synchrotron frequency

particle = Proton() # <------- Predefined particle with definition of mass and charge
alpha_0 = 1 / gamma_t**2 # <-- Momentum compaction factor

ring_PS = Ring(circumference, alpha_0, momentum, particle, n_turns=n_turns)

from PS_bunch_rotation_generation import C40_77_prog, C40_78_prog, C80_08_prog, C80_88_prog, C80_phase_prog, C40_phase_prog

start_bunch_rotation_h84 = 40e-6

# Linearized bunch rotation
length_step_h168 = 141  # us
amplitude_step_h168 = 40.  # kV
C40_77_prog_factor = 1. # 284.700/300

# Other progs
C40_78_prog_factor = 1. # 256.607/300
C80_08_prog_factor = 1.
C80_88_prog_factor = 1.
C80_89_prog_factor = 0.

# Phase
final_phase_h84 = 0
initial_phase_h168 = 0
step_adjust_h168 = 0
final_phase_h168 = initial_phase_h168 - np.pi + step_adjust_h168

C40_77_RF_program = 1e3 * \
    C40_77_prog(ring_PS.cycle_time * 1e6,
                start_bunch_rotation_h84 * 1e6, C40_77_prog_factor)

C40_78_RF_program = 1e3 * \
    C40_78_prog(ring_PS.cycle_time * 1e6, start_bunch_rotation_h84 * 1e6,
                C40_78_prog_factor)
# C40_78_RF_program[ring_PS.cycle_time<=start_bunch_rotation_h84] = 0

rf_prog_h84 = C40_77_RF_program + C40_78_RF_program

phi_prog_h84 = C40_phase_prog(ring_PS.cycle_time * 1e6,
                                start_bunch_rotation_h84 * 1e6,
                                0, final_phase_h84)

C80_08_RF_program = 1e3 * C80_08_prog(ring_PS.cycle_time * 1e6,
                                        start_bunch_rotation_h84 * 1e6,
                                        length_step_h168,
                                        amplitude_step_h168,
                                        C80_08_prog_factor)

C80_88_RF_program = 1e3 * C80_88_prog(ring_PS.cycle_time * 1e6,
                                        start_bunch_rotation_h84 * 1e6,
                                        length_step_h168,
                                        amplitude_step_h168,
                                        C80_88_prog_factor)

C80_89_RF_program = 1e3 * C80_88_prog(ring_PS.cycle_time * 1e6,
                                        start_bunch_rotation_h84 * 1e6,
                                        length_step_h168,
                                        amplitude_step_h168,
                                        C80_89_prog_factor)

rf_prog_h168 = C80_08_RF_program + C80_88_RF_program + C80_89_RF_program

phi_prog_h168 = C80_phase_prog(ring_PS.cycle_time * 1e6, start_bunch_rotation_h84 *
                                1e6, length_step_h168, initial_phase_h168, final_phase_h168)

#RF station
rf_params_PS = RFStation(ring_PS, harmonic_numbers_PS,
                             [rf_prog_h84, rf_prog_h168],
                             [phi_prog_h84, phi_prog_h168],
                             n_rf_systems_PS)

# Get the calculated parameters at a given turn, for a given rf harmonic
# The RF station is located after a ring section
# Some are done for you, try printing the RF voltage, and frequency
# Hint: type "rf_station." and then press tab to check the content

idx_section = 0
idx_turn = 0
idx_rf = 0
idx_turn = 0


print("-------Ring parameters-------")
print('Total energy [GeV/c]: ', ring_PS.energy[idx_section, idx_turn]/1e9)
print('Kinetic energy [GeV/c]: ', ring_PS.kin_energy[idx_section, idx_turn]/1e9)
print('Revolution period [us]: ', ring_PS.t_rev[idx_turn]*1e6)
print('Relativistic velocity: ', ring_PS.beta[idx_section, idx_turn])
print('Lorentz factor: ', ring_PS.gamma[idx_section, idx_turn])
print('Linear slippage factor: ', ring_PS.eta_0[idx_section, idx_turn])
print("------------------------------")
print('   ')
print("-------RF parameters-------")
print('RF harmonic number: ', rf_params_PS.harmonic[idx_rf, idx_turn])
print('RF period [s]: ', rf_params_PS.t_rf[idx_rf, idx_turn])
print('Synchrotron frequency [Hz]: ', rf_params_PS.omega_s0[idx_turn] / (2 * np.pi))
print('Synchrotron tune: ', rf_params_PS.Q_s[idx_turn])
print('RF Voltage [MV]: ', rf_params_PS.voltage[idx_rf, idx_turn]/1e6)
print('Omega RF: ', rf_params_PS.omega_rf[idx_rf, idx_turn])
print('Revolution Frequency [Hz]: ', 1/rf_params_PS.t_rev[idx_rf])
print('RF Frequency [MHz]: ', 1e-6/rf_params_PS.t_rf[idx_rf, idx_turn])
print("------------------------------")
print('   ')
print("-------Particle parameters-------")
print('Number of particles (Single particle only): ', n_particles)
print('Number of macroparticles: ', n_macroparticles)
print('Number of turns: ' , n_turns)
print('Energy spread: ', energy_spread)
print('Energy offset: ', energy_offset)
print('Bunch length: ', (rf_params_PS.t_rf[0, 0] / 4)*bunch_len_factor)
print('Bunch length factor: ', bunch_len_factor)
print('Bunch offset (time delay): ', bunch_offset)
print("------------------------------")

#Track a single particle
if single_particle:
    n_particles = 1
    intensity = 2.5e11

    ring = ring_PS
    rf_station = rf_params_PS

    beam = Beam(ring, n_particles, intensity)
    beam.dt[0] = 2.965E-9
    beam.dE[0] = 1E6

    rf_tracker = RingAndRFTracker(rf_station, beam)
    full_tracker = FullRingAndRF([rf_tracker])

    saved_dt = np.zeros(ring.n_turns) # <----- "Empty" array to accumulate particle time coordinates
    saved_dE = np.zeros(ring.n_turns) # <----- "Empty" array to accumulate particle energy coordinates

    for turn in range(ring.n_turns):
        full_tracker.track() # <-------------- Apply BLonD kick & drift functions
        saved_dt[turn] = beam.dt # <---------- Saving the time coordinate
        saved_dE[turn] = beam.dE # <---------- Saving the energy coordinate
    
    plot_phase_space_trajectory(saved_dt, saved_dE)
    plt.show()

#Track multiple single particle
if single_particle_track:
    ring = ring_PS
    rf_station = rf_params_PS

    beam = Beam(ring, n_particles, intensity)

    # Linearly distribute particles in time
    beam.dt[:] = np.linspace(0, rf_station.t_rf[0, 0] / 2, n_particles)
    beam.dE[:] = np.zeros(n_particles) 

    rf_tracker = RingAndRFTracker(rf_station, beam)
    full_tracker = FullRingAndRF([rf_tracker])

    saved_dt = np.zeros((n_turns, n_particles))
    saved_dE = np.zeros((n_turns, n_particles))

    for turn in range(ring.n_turns):
        full_tracker.track()
        saved_dt[turn] = beam.dt
        saved_dE[turn] = beam.dE

    plot_phase_space_trajectory(saved_dt, saved_dE)
    plt.show()

#Plot a bunch after a certain number of turns (n_turns)
if bunch_plot:
    ring = ring_PS
    rf_station = rf_params_PS

    beam = Beam(ring, n_macroparticles, intensity)

    bunch_position = rf_station.t_rf[0, 0] / 2 + bunch_offset # <--------- Optional position offset
    bunch_length = (rf_station.t_rf[0, 0] / 4)*bunch_len_factor

    energy_position = energy_offset
    energy_spread = energy_spread # <------- Optional change of energy spread

    if energy_spread == 0 or energy_spread == 0.0:
        parabolic(
            ring,
            rf_station,
            beam,
            bunch_length,
            bunch_position=bunch_position, # <------ Uncomment to add position offset
            bunch_energy=energy_position
        )
    else:
        parabolic(
            ring,
            rf_station,
            beam,
            bunch_length,
            bunch_position=bunch_position, # <------ Uncomment to add position offset
            bunch_energy=energy_position,
            energy_spread=energy_spread # <-------- Uncomment to change of energy spread
        )

    rf_tracker = RingAndRFTracker(rf_station, beam)
    full_tracker = FullRingAndRF([rf_tracker])

    for turn in range(ring.n_turns):
        full_tracker.track()

    plot_phase_space_distribution(beam.dt, beam.dE)
    plt.show()

#Track a bunch for n_turns plotting every dt_plt turns
if bunch_tracking:
    ring = ring_PS
    rf_station = rf_params_PS

    beam = Beam(ring, n_macroparticles, intensity)

    bunch_position = rf_station.t_rf[0, 0] / 2 + bunch_offset # <--------- Optional position offset
    bunch_length = (rf_station.t_rf[0, 0] / 4)*bunch_len_factor

    energy_position = energy_offset
    energy_spread = energy_spread # <------- Optional change of energy spread

    if energy_spread == 0 or energy_spread == 0.0:
        parabolic(
            ring,
            rf_station,
            beam,
            bunch_length,
            bunch_position=bunch_position, # <------ Uncomment to add position offset
            bunch_energy=energy_position
        )
    else:
        parabolic(
            ring,
            rf_station,
            beam,
            bunch_length,
            bunch_position=bunch_position, # <------ Uncomment to add position offset
            bunch_energy=energy_position,
            energy_spread=energy_spread # <-------- Uncomment to change of energy spread
        )
    
    cut_opt = CutOptions(cut_left=0, cut_right=rf_station.t_rf[0, 0], n_slices=60)
    fit_opt = FitOptions(fit_option='gaussian') # <----- Specify the type of fit to be applied
    profile = Profile(beam, CutOptions=cut_opt, FitOptions=fit_opt)

    rf_tracker = RingAndRFTracker(rf_station, beam)
    full_tracker = FullRingAndRF([rf_tracker])

    plots = Plot(ring, rf_station, beam, dt_plt, ring.n_turns,
             0-t_offset, rf_station.t_rf[0, 0]+t_offset, # <----- dt limits of the phase space plot
             E_min, E_max, # <---------------- dE limits of the phase space plot
             show_plots=True,
             separatrix_plot=True)

    if bunch_parameters:
        for turn in range(ring.n_turns):
            full_tracker.track()
            profile.track()
            plots.track()
        plt.figure('profile')
        plt.clf()
        plt.plot(profile.bin_centers, profile.n_macroparticles)
        plt.show()

    elif bunch_monitor:
        h5file = './data/saved_data'
        bunchmonitor = BunchMonitor( # <------- Declare bunch monitor, used to compute bunch parameters
            ring, rf_station, beam,
            h5file,
            Profile=profile)
        bunchmonitor.track()
        rf_tracker = RingAndRFTracker(rf_station, beam)
        full_tracker = FullRingAndRF([rf_tracker])
        for turn in range(ring.n_turns):
            full_tracker.track()
            profile.track()
            plots.track()
            bunchmonitor.track()
        plt.figure('profile')
        plt.clf()
        plt.plot(profile.bin_centers, profile.n_macroparticles)
        plt.show()

    else:
         for turn in range(ring.n_turns):
            full_tracker.track()
            plots.track()        

    
#Track an accelerating bunch for n_turns plotting every dt_plt turns
if bunch_tracking_acc:
    ring = ring_PS
    rf_station = rf_params_PS

    bunch_position = rf_station.t_rf[0, 0] / 2 + bunch_offset # <--------- Optional position offset
    bunch_length = (rf_station.t_rf[0, 0] / 4)*bunch_len_factor

    energy_position = energy_offset
    energy_spread = energy_spread # <------- Optional change of energy spread

    if energy_spread == 0 or energy_spread == 0.0:
        parabolic(
            ring,
            rf_station,
            beam,
            bunch_length,
            bunch_position=bunch_position, # <------ Uncomment to add position offset
            bunch_energy=energy_position
        )
    else:
        parabolic(
            ring,
            rf_station,
            beam,
            bunch_length,
            bunch_position=bunch_position, # <------ Uncomment to add position offset
            bunch_energy=energy_position,
            energy_spread=energy_spread # <-------- Uncomment to change of energy spread
        )

    rf_tracker = RingAndRFTracker(rf_station, beam)
    full_tracker = FullRingAndRF([rf_tracker])


    plots = Plot(ring, rf_station, beam, dt_plt, ring.n_turns,
             0-t_offset, rf_station.t_rf[0, 0]+t_offset, # <----- dt limits of the phase space plot
             E_min, E_max, # <---------------- dE limits of the phase space plot
             show_plots=True,
             separatrix_plot=True)

    for turn in range(ring.n_turns):
        full_tracker.track()
        plots.track() 


