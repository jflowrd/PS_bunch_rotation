# Impedance of Empty Pumping Manifolds

The impedance of the empty pumping manifolds downstream of the main Magnet Units.

Out of the 100 elements, 36 are empty.
 
The CST model is based on the **Standard Downstream Assembly** as representative for
the overall empty pumping manifolds. 

The remaining elements were not yet evaluated and the **Standard Downstream Assembly**
model is used to represent all the pumping manifolds downstream of the MUs.

Uses following Stepfiles:
- ST0371891_01.1
- ST0397134_01.1

Damping resistors are inserted in the arm of the pumping manifold \[1-2\].
Their contribution is not included in the model.

## Subfolders and usage

- CST Raw Data

The model is based on the Standard Downstream Assembly. The model includes a flange,
a bellow, and the armature to the pump. The whole armature to the pump is included
in the model.

The damping resistors in the arm of the pumping manifold are not included in the model.

The output result is the impedance in frequency domain and
can directly be used in particle simulations. The impedance
is noisy, to obtain a smoother response the fit with resonators
can be preferred.

<details>
  <summary>CST model view</summary>
  <img src="cst_raw_data/cst_model_images/std_int_downstream_assy_Model_View1.png">
  <img src="cst_raw_data/cst_model_images/std_int_downstream_assy_Model_View2.png">
  <img src="cst_raw_data/cst_model_images/std_int_downstream_assy_Model_View3.png">
  <img src="cst_raw_data/cst_model_images/std_int_downstream_assy_Model_View4.png">
</details>

- Resonators

The output CST impedance is fitted with 8 resonators and an extra contribution in
ImZ/f.

Some smaller resonances were not included in the fit to ensure convergence.

## Plots

<details>
  <summary>Resonator fit</summary>
  <img src="Resonators/fitted_impedance.png">
</details>

## History

To be completed

## Further work

To be completed

- Modelling the variations of the pumping manifolds
- Draw a list of the damping resistors insertion and evaluate their effect

## References

To be completed

- \[1\] [D. Boussard, Observation of microwave longitudinal instabilities in the CPS](https://cds.cern.ch/record/872559/files/cer-002556319.pdf)
- \[2\] [J. Belleman pages](http://jeroen.web.cern.ch/jeroen/codd/damperresistors.html)
- [B. Popovic, Longitudinal limitations with LIU-PS RF upgrades and mitigation strategy](https://indico.cern.ch/event/750790/contributions/3108016/attachments/1719965/2776247/RG_PS_Impedance_Model_Meeting.pdf)
- [B. Popovic, SPS Injection Losses Review](https://indico.cern.ch/event/672967/contributions/2753563/attachments/1567554/2471183/PS_Impedance_Status_Injection_Losses_Meeting_301117.pdf)

## Credits

To be completed

