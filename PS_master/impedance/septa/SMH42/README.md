# Impedance of SMH42

The SMH42 is the injection septum from PSB to PS for the proton beam.

The impedance is difficult to obtain from CST simulations due
to the lamniation material which characteristics are difficult to
get and model accurately.

The impedance was evaluated using wire measurements.

<details>
  <summary>SMH42 in the ring</summary>
  <img src="http://cern.ch/psring/psring/pictures/fullsize/ss42.jpg">
</details>

## Subfolders and usage

- Measurements

The impedance source can be loaded directly from the measured data (impedance
table).

This model may not be the most suitable for particle simulations since the
measured data is noisy and cut at 1 GHz.

- Resonators

The impedance is fitted with resonators to get the parameters of the resonances
at low frequency, and a contribution in ImZ/f is added to get the reactive
impedance up to 1 GHz.

Note that the resistive impedance was offset by the value at 1 MHz in order
to have 0 contribution at DC (appearent systematic offset).

## Plots

<details>
  <summary>Resonator and ImZ/f fit</summary>
  <img src="Resonators/fitted_resonances_realimag.png">
  <img src="Resonators/fitted_broadband_realimag.png">
</details>

## History

To be completed

## Further work

To be completed

## References

To be completed

- [B. Popovic, 36th Impedance Working Group meeting](https://indico.cern.ch/event/860185/contributions/3623056/attachments/1938347/3212946/IWG_SMH42_041119.pdf)

## Credits

- Measurements: B. Popovic
