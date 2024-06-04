# Impedance of Transverse Feedback Kicker

The impedance of the Transverse Feedback Kicker in straight section 97.

The TFB Kicker was modeled in CST and is stored on
[EDMS](https://edms.cern.ch/ui/#!master/navigator/project?P:100327819:100578987:subDocs).

<details>
  <summary>TFB Kicker in the ring</summary>
  <img src="http://cern.ch/psring/psring/pictures/fullsize/ss97.jpg">
</details>


## Subfolders and usage

- CST Raw Data

The wakefield and eigenmodes are obtained using CST.

The eigenmode and wakefield give comparable results. However, a resonance
at higher frequency than the maximum frequency considered in simulations
seem to contribute with broandband inductive impedance.

For eigenmodes, the contributions from all the folders should be loaded.

<details>
  <summary>CST model view</summary>
  <img src="cst_raw_data/cst_model_images/TFBKicker_CST_Model_View1.png">
</details>

- Resonators

The result of the CST simulation is fitted with resonators to get
the more global features of the impedance.

The impedance is fitted with 10 resonators and an extra contribution of 
constance ImZ/f.

The eigenmode calculation does not provide with the low frequency inductive
contribution, the resonator fit is recommended as default usage.

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

# References

- [B. Popovic, LIU-PS Beam Dynamics WG meeting #35](https://indico.cern.ch/event/851693/contributions/3581115/attachments/1922953/3181642/LIU_PS_BD_WG_Trans_Damper_091019_FINAL.pdf)

