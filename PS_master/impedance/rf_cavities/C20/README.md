# Impedance of 13.3-20 MHz cavities

Those cavities are used for rf manipulations (ferrite loaded coaxial resonator)

Each cavity consists in two rf gaps providing for a max total of 20 kV. 
The impedance files are given for one full cavity, including both gaps.

There are two cavities in total, named after there position in the ring (C20-80 and C20-92).

The cavities can be tuned around 13.3 MHz or 20 MHz. The characteristics of
the cavity are (see reference):

| Parameter \ Frequency [MHz] 			  | 13.3	| 20	  |
| -----------------------------------	|:-----:|:-----:|
| Maximum VRF [kVp]       				    | 20	  | 20 	  |
| Shunt resistance R at 20 kVp [kOhm]	| 1.7  	| 2.0 	|
| Quality factor Q at 20 kVp  			  | 82 	  | 63 	  |
| R/Q at 20 kVp [Ohm] 					      | 20.73	| 31.75	|
| Quality factor Q at 100 Vp  			  | 163	  | 100 	|

When not in use, the cavities are shorted by the gap relay (always the case
with the exception of when one cavity is used for rf manipulation).

The impedance model in this database includes the effect of the fast feedback 
loop:

| Parameter \ Frequency [MHz] 			    | 13.3		  | 20  		|
| -------------------------------------	|:---------:|:-------:|
| Reduced quality factor Q at 20 kVp  	| 4.85 		  | 4.2239 	|
| Reduced quality factor Q at 100 Vp  	| 5.0222	  | 4.5801	|

A test model of a multi-Harmonic feedback was included for studies purpose,
but is not implemented in the machine or planned. The implementation is
based on the same idea as [C40-MHFB](../C40/MHFB/README.md) and [C80-MHFB](../C80/MHFB/README.md)

<details>
  <summary>C10 in the ring</summary>
  <img src="http://cern.ch/psring/psring/pictures/fullsize/ss80.jpg">
  <img src="http://cern.ch/psring/psring/pictures/fullsize/ss92.jpg">
</details>

## Subfolders and usage

- Resonators

The resonators are defined as in the general description.
The central frequency of the cavity can be modified but should be adjusted
to match the harmonic $`h=42`$ with a proton beam at top energy for the
20 MHz case.

- MHFB

The model is based on the different scenarios defined by resonators,
on which extra parameters were included to model the MHFB.

The user specifies the scenario to load the base impedance , the revolution period, 
the impedance reduction target and the bandwidth at the central frequency of the cavity.
The relevant feedback parameters (gain, delay, bandwidth) are interpolated from a lookup table.
The user can specify whether the notch at the central frequency should be removed
(e.g. in the eventuality where the cavity is pulsing).

## Plots

<details>
  <summary>Impedance with direct feedback at 13.3 MHz high voltage</summary>
  <img src="All/Resonators/single_resonator_13.3MHz_20kV.png">
</details>

<details>
  <summary>Impedance with direct feedback at 20 MHz high voltage</summary>
  <img src="All/Resonators/single_resonator_20MHz_20kV.png">
</details>

<details>
  <summary>Feedback with -20dB impedance reduction with central notch</summary>
  <img src="MHFB/Resonators/C20-HV_impedance_dB_with_MHFB_-20dB_477kHz_central_False.png">
  <img src="MHFB/Resonators/C20-HV_impedance_with_MHFB_-20dB_477kHz_central_False.png">
</details>

<details>
  <summary>Feedback with -20dB impedance reduction without central notch</summary>
  <img src="MHFB/Resonators/C20-HV_impedance_dB_with_MHFB_-20dB_477kHz_central_True.png">
  <img src="MHFB/Resonators/C20-HV_impedance_with_MHFB_-20dB_477kHz_central_True.png">
</details>

## Further work

- The shunt impedance is likely to be underestimated like the 10 MHz cavities were.
This is to be confirmed.
- Include parametrized model of open-loop impedance with fast feedback

To be determined (beam based measurement of the impedance ? 
residual impedance with gap relay ? impedance for low vs. high voltage in the gap ? 
simulation based modelling ?)

## References

### Measurements
- [M. Morvillo et al., The PS 13.3-20 MHZ RF Systems for LHC](http://cds.cern.ch/record/620368)

### Other
- [Fifty years of the CERN Proton Synchrotron : Volume 1](https://cds.cern.ch/record/1359959)

