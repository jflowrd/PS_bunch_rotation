# Impedance of Finemet cavity

The Finement cavity in the PS is used as a longitudinal damper for dipole mode of oscillation.

One cavity is installed and consists of six rf cells.

The impedance is based on measurements and should include the contribution of the amplifiers and the transmission line.
Presently, two fits are proposed based on measured data from 2013 and 2018. Both
were kept in the present model since the flatness of the impedance is different
and the effect on the beam can be assessed.

<details>
  <summary>Finemet cavity in the ring</summary>
  <img src="http://cern.ch/psring/psring/pictures/fullsize/ss02.jpg">
</details>

## Subfolders and usage

- Resonators

The measured data is fitted with resonators to catch the more global features
of the impedance and allow to compute the impedance on a custom frequency range.

The measured data from 2013 contains both amplitude and phase information.
The measured data from 2018 consists of only amplitude information.
The fit was adjusted to have a strictly positive resistive impedance.

The fit aims at having a good fit up to ~30MHz (the resonance above comes from the measurement setup). 

- MHFB

The impedance of the Finemet cavity reduced using a Multi-Harmonic Feedback. This
is not implemented presently in the machine. This was included in the present
model for testing purposes.

The MHFB implementation was done using the same principles as for the high frequency cavities.

## Plots

<details>
  <summary>Resonator fit on 2013 measurements</summary>
  <img src="Resonators_2013/fitted_multi_resonators_abs_log.png">
  <img src="Resonators_2013/fitted_multi_resonators_real_imag_log.png">
</details>

<details>
  <summary>Resonator fit on 2018 measurements</summary>
  <img src="Resonators_2018/fitted_multi_resonators_abs_log.png">
  <img src="Resonators_2018/fitted_multi_resonators_real_imag_log.png">
</details>

<details>
  <summary>MHFB on Resonators_2018 model</summary>
  <img src="MHFB/Resonators/2018_impedance_dB_with_MHFB_-15dB_458kHz.png">
</details>


## History

To be completed

## Further work

- Results with CST Simulation model to be included

## References

- [S. Persichelli et al., Impedance studies for the PS Finemet loaded
longitudinal damper,
IPAC'14](https://cds.cern.ch/record/1742165/files/CERN-ACC-2014-0113.pdf) 

## Credits

- M. Paoluzzi for the measured data
