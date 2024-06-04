# Impedance of Bellows

The impedance of the standard bellows upstream of the main Magnet Units.

There are 99 bellows upstream of the MUs.
The CST model is based on the **Standard Upstream Assembly**
which is the design for 67 of the 99 elements.

The remaining elements were not yet evaluated and the **Standard Upstream Assembly**
model is used to represent all the bellows upstream of the MUs.

## Subfolders and usage

- CST Raw Data

The wakefield is obtained using CST. The impedance is mostly a pure inductive
component up to high frequencies (>1.6GHz).

<details>
  <summary>CST model view</summary>
  <img src="cst_raw_data/cst_model_images/std_upstream_assy_Model_View1.png">
</details>

- Inductive

The result of the CST simulation is fitted with constant ImZ/f.

## Plots

<details>
  <summary>ImZ/f fit</summary>
  <img src="Inductive/fitted_impedance.png">
</details>

## History

To be completed

## Further work

To be completed

- Modeling the different bellows
- Extend the frequency range of the simulation

## References

To be completed

- [B. Popovic, Longitudinal limitations with LIU-PS RF upgrades and mitigation strategy](https://indico.cern.ch/event/750790/contributions/3108016/attachments/1719965/2776247/RG_PS_Impedance_Model_Meeting.pdf)
- [B. Popovic, SPS Injection Losses Review](https://indico.cern.ch/event/672967/contributions/2753563/attachments/1567554/2471183/PS_Impedance_Status_Injection_Losses_Meeting_301117.pdf)

## Credits

To be completed

