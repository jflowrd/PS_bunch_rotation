# BLonD PS bunch rotation
A script to simulate the non-adiabatic bunch shortening in the PS, with intensity effects.

## How to run
Main file is 'PS_master/PS_extraction_impedance.py'
- The script uses the following modules:
  - SPS_parameters.py for plotting SPS separatrix and calculating losses
  - PS_bunch_rotation_generation.py to define the 40 MHz and 80 MHz voltage and phase programs
  - impedance_loader.py to load the PS impedance model
  - /plots and /blond_common for plotting and fitting
- Full list of packages are found in requirements.txt (although many are not needed and are for other experimentations)
