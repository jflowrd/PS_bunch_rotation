# Impedance of Stripline BPM

The impedance of the Stripline BPM in straight section 72.

The BPM was modeled in CST and is stored on
[EDMS](https://edms.cern.ch/ui/#!master/navigator/project?P:100327819:100512463:subDocs).
The original model is from the impedance website \[1\].

Wakefield and eigenmode simulations were run using the model shown in \[1\].

<details>
  <summary>Stripline BPM in the ring</summary>
  <img src="http://cern.ch/psring/psring/pictures/fullsize/ss72.jpg">
</details>


## Subfolders and usage

- CST Raw Data

The wakefield and eigenmodes are obtained using CST.

The main modes are well identified, but the R and Q parameters of the resonances
are not identical. Some low frequency contribution is visible with wakefield
simulations, hitting towards the presence of a high frequency resonance
above the highest frequency considered in the simulation.

<details>
  <summary>CST model view</summary>
  <img src="cst_raw_data/cst_model_images/StriplineBPM_CST_Model_View1.png">
</details>

- Resonators

The result of the CST simulation is fitted with resonators to get
the more global features of the impedance.

The impedance is fitted with 5 resonators and an extra contribution of 
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

- \[1\] http://impedance.web.cern.ch/impedance/PS.htm
- \[2\] http://psring.web.cern.ch/psring/psring/showpicture.php?section=72
- \[3\] ECR: 1234117
- [B. Popovic, Longitudinal limitations with LIU-PS RF upgrades and mitigation strategy](https://indico.cern.ch/event/750790/contributions/3108016/attachments/1719965/2776247/RG_PS_Impedance_Model_Meeting.pdf)

