# Impedance of Sector Valves

The impedance of the sector valves.

There are 10 sector valves located after straight sections 00,10,20,30,40,50,60,70,80,90.

The valve was modelled in CST and is available on
[EDMS](https://edms.cern.ch/ui/#!master/navigator/document?P:100327819:100512607:subDocs).
with the layout Drawings:

- ps_lm___0011-vAC
- ps_lm___0030-vAD
- ps_lm___0050-vAE
- ps_lm___0070-vAC
- ps_lm___0090
- ps_lm___0110vAB
- ps_lm___0130vAD
- ps_lm___0150vAC

## Subfolders and usage

- CST Raw Data

The wakefield and eigenmodes are obtained using CST.

The internal structure is proprietary design and was drawn in the CST model using datasheets.

The CST model was benchmarked with measurements to ensure a good modelling of the more
important features of the internal components.

The output result is the impedance in frequency domain and
can directly be used in particle simulations. The impedance
is noisy, to obtain a smoother response the fit with resonators
can be preferred.

The eigenmode ouput can be used directly as resonators.

<details>
  <summary>CST model view</summary>
  <img src="cst_raw_data/cst_model_images/SectorValve_CST_Model_View1.png">
  <img src="cst_raw_data/cst_model_images/SectorValve_CST_Model_View2.png">
</details>

- Resonators

The result of the CST simulation is fitted with resonators to get
the more global features of the impedance.

The impedance is fitted with 5 resonators (a resonator added at 2.28 GHz to have a better fit below <2GHz).

The eigenmode calculation in CST ranges to larger frequencies than the wake calculation and is the recommended usage.

## Plots

<details>
  <summary>Comparison eigen/wake</summary>
  <img src="cst_raw_data/plot_comparison/cst_comparison.png">
</details>

<details>
  <summary>Resonator fit</summary>
  <img src="Resonators/fitted_impedance_wake.png">
</details>

## History

To be completed

## Further work

To be completed

- Variations with neighboring bellows and flanges in the CST model

## References

To be completed

- [B. Popovic, Longitudinal limitations with LIU-PS RF upgrades and mitigation strategy](https://indico.cern.ch/event/750790/contributions/3108016/attachments/1719965/2776247/RG_PS_Impedance_Model_Meeting.pdf)
- [B. Popovic, SPS injection losses review](https://indico.cern.ch/event/672967/contributions/2753563/attachments/1567554/2471183/PS_Impedance_Status_Injection_Losses_Meeting_301117.pdf)
- [B. Popovic, LIU-PS Beam Dynamics WG meeting #1](https://indico.cern.ch/event/662292/contributions/2704441/attachments/1516648/2367108/BP_PS_Impedance_Model_Update_310817_Final.pdf)

## Credits

To be completed

