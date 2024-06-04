# Impedance of 10 MHz cavities

The main cavities for acceleration in PS (ferrite loaded coaxial resonator).

Each cavity consists in two accelerating gaps providing for a max total of 20 kV. 
The impedance files are given for one full cavity, including both gaps.

There are eleven cavities in total, named after their position in the ring. The
cavities in the same tuning groups share the same harmonic number:

- Group A: C10-36, C10-46, C10-51
- Group B: C10-56, C10-66, C10-76, C10-81
- Group C: C10-86, C10-91, C10-96
- Group D: C10-11 (this cavity can replace the cavity from any other group, in
nominal operation the cavity is shorted by a gap relay and its impedance is
assumed to be reduced to zero)

When not in use, the cavities from a same group can be detuned to a frequency
out from the bandwidth of the beam spectrum or shorted by the gap relay. The
cavity impedance with the gap relay is reduced, but not completely negligible. 

Two feedback system are implemented around each cavity. A direct feedback loop
with small delay and a 1-turn delay feedback. By default, the impedance
model includes the effect of the fast feedback loop. The Amplitude Voltage
Control (AVC) loop is not included.

The power amplifiers are being updated during the Long Shutdown 2 (2019-2020).
The shunt impedance is expected to be reduced by a factor 2. The upgraded
power amplifier was tested on C10-11, its impedance model includes the effect
of the upgraded amplifier.

<details>
  <summary>C10 in the ring</summary>
  <img src="http://cern.ch/psring/psring/pictures/fullsize/ss11.jpg">
  <img src="http://cern.ch/psring/psring/pictures/fullsize/ss36.jpg">
  <img src="http://cern.ch/psring/psring/pictures/fullsize/ss46.jpg">
  <img src="http://cern.ch/psring/psring/pictures/fullsize/ss51.jpg">
  <img src="http://cern.ch/psring/psring/pictures/fullsize/ss56.jpg">
  <img src="http://cern.ch/psring/psring/pictures/fullsize/ss66.jpg">
  <img src="http://cern.ch/psring/psring/pictures/fullsize/ss76.jpg">
  <img src="http://cern.ch/psring/psring/pictures/fullsize/ss81.jpg">
  <img src="http://cern.ch/psring/psring/pictures/fullsize/ss86.jpg">
  <img src="http://cern.ch/psring/psring/pictures/fullsize/ss91.jpg">
  <img src="http://cern.ch/psring/psring/pictures/fullsize/ss96.jpg">
</details>

## Subfolders and usage

Various possibilities are available to load the impedance for calculations and
simulations.

- Individual

The impedance of each single cavity can be loaded individually (including
direct feedback).

The impedance is from beam-based measurements and all contributions were fitted
individually with multiple resonators.

More details [here](Individual/README.md)

- All

The impedance is averaged over all 10 cavities (excluding C10-11).

The resulting impedance is fitted with resonators and a parametric model is generated
to get the impedance for all possible tuning parameter.

More details [here](All/README.md)

- OTFB

The impedance including the contribution of the 1-turn delay feedback.

The averaged impedance over the 10 cavities, fitted as a single resonator, is used as a
basis and filtered using a comb-filter.

More details [here](OTFB/README.md)
 
- ClosedGapRelay

The impedance of the cavity with closed gap relays.

The impedance is from beam-based measurements of a single gap. The impedance of both
gaps both shorted with gap relays is assumed identical. The resulting impedance is fitted
with resonators.

More details [here](Closed_Gap_Relay/README.md)

## History

To be completed

## Further work

To be completed (impedance for low vs. high voltage in the gap ? simulation based modelling ?)

## References

To be completed

### Measurements
- [G. Favia, Study of the beam-cavity interaction in the CERN PS 10 MHz 
cavities and investigation of hardware solutions to reduce beam loading, 
PhD Thesis](https://cds.cern.ch/record/2286835) 
(NB: the data available in this database was cut to cover the frequency 
range 0-25 MHz, the frequency range of the original data after analysis and 
interpolation was 62.5 GHz)

### Other
- [G. Favia et al., Study of the beam-cavity interaction in the CERN PS 10 MHz 
rf system](https://cds.cern.ch/record/2207324/files/mopor012.pdf)
- [Fifty years of the CERN Proton Synchrotron : Volume 1](https://cds.cern.ch/record/1359959)
- [H. Damerau et al., Longitudinal Coupled-Bunch Instabilities in the CERN PS](https://cds.cern.ch/record/1055555/)
- [R. Garoby et al., RF system for high beam intensity acceleration in the CERN PS](https://cds.cern.ch/record/196412)
- [D. Grier, The PS 10 MHz cavity and power amplifier](https://cds.cern.ch/record/960421/)
- [L. Ventura, Studies of Longitudinal Coupled-Bunch Instabilities in the LHC Injectors Chain](https://cds.cern.ch/record/2253783)

## Credits

To be completed

- Measured data: G. Favia
- Resonator fit: A. Lasheen

