# Impedance of 40 MHz cavities

Those cavities are used for rf manipulations (narrow-band resonator)

Each cavity consists in one rf gap providing a max voltage of 300 kV (the effective
rf voltage measured from quadrupole oscillations indicates a lower rf voltage, ref necessary). 
The impedance files are given for one full cavity.

There are two cavities in total, named after there position in the ring (C40-77 and C40-78).

When not in use, the cavity gap can be closed mecanically. Presently both
cavities are usually open during operation.

<details>
  <summary>C40 in the ring</summary>
  <img src="http://cern.ch/psring/psring/pictures/fullsize/ss77.jpg">
  <img src="http://cern.ch/psring/psring/pictures/fullsize/ss78.jpg">
</details>

## Subfolders and usage

- Individual

This is a model based on beam-based measurements (main harmonic). The 
beam based measurements display important difference with the expected design
value and is the prefered option till an explanation to the discrepancy is
found.

More details [here](Individual/README.md)

- MHFB

The model including the Multi Harmonic Feedback. The model is based on the
beam based impedance modified with the MHFB transfer function.

More details [here](MHFB/README.md)

- High Order Modes

The impedance of the HOMs is a list of resonators and should be imported on top
of the main harmonic impedance above

More details [here](High_Order_Modes/README.md)

## Further work

- Understand the discrepancy between the expectation in terms of Rshunt with
respect to the measured one.
- Include the impedance with mechanical short
- Other tasks to be determined (continue beam based measurement of the impedance ? 
modelling of the HOMs with couplers ? impedance for low vs. high voltage in the gap ? 
simulation based modelling ?)

## History

## References

### Design
- [The PS complex as proton pre-injector for the LHC : design and implementation report](https://cds.cern.ch/record/449242)

### Other
- [Fifty years of the CERN Proton Synchrotron : Volume 1](https://cds.cern.ch/record/1359959)

## Credits

To be completed
