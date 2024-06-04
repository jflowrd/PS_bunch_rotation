# Impedance of KFA79

The KFA79 is a fast extraction kicker for extraction to the SPS.

Note that KFA13 & KFA21 are identical and that KFA13, KFA21, KFA71 and KFA79 all contain identical modules.

The first estimation was done using a Tsutsui model available on the following [PS impedance webpage](http://impedance.web.cern.ch/impedance/)
and on the [IRIS PS repository](https://gitlab.cern.ch/IRIS/PS_IW_model/-/tree/master/Impedances/Longitudinal).

More recent evaluation consisted in modelling the kicker
using CST. The model was built from old drawings available on
DFS (\\cern.ch\dfs\Workspaces\o\Old Drawings\Complexe_PS\PS\)
and CATIA files (reference needed).

The present CST model is stored on
[EDMS](https://edms.cern.ch/ui/#!master/navigator/document?P:100327819:100511972:subDocs).

*Note that the CST model does not include transition pieces like the other kickers. This could
explain the noisy impedance profile.*

<details>
  <summary>KFA79 in the ring</summary>
  <img src="http://cern.ch/psring/psring/pictures/fullsize/ss79.jpg">
</details>

## Subfolders and usage

- Tsutsui

The model was computed using ImpedanceWake2D with the following parameters

    - a=0.0735 half width [m]
    - b=0.0265 half heigth (center-ferrite border) [m]
    - d=0.1 half heigth (center-PEC border including ferrite) [m]
    - L=0.666 Kicker's length [m]

- CST Raw Data

The wakefield is obtained with CST simulations. The model was
built based on 2D drawings and CATIA files.

The output result is the impedance in frequency domain and
can directly be used in particle simulations. The impedance
is noisy, to obtain a smoother response the fit with resonators
can be preferred.

<details>
  <summary>CST model view</summary>
  <img src="cst_raw_data/cst_model_images/KFA79_CST_Model_View1.png">
  <img src="cst_raw_data/cst_model_images/KFA79_CST_Model_View2.png">
  <img src="cst_raw_data/cst_model_images/KFA79_CST_Model_View3.png">
  <img src="cst_raw_data/cst_model_images/KFA79_CST_Model_View4.png">
  <img src="cst_raw_data/cst_model_images/KFA79_CST_Model_View5.png">
  <img src="cst_raw_data/cst_model_images/KFA79_CST_Model_View6.png">
  <img src="cst_raw_data/cst_model_images/KFA79_CST_Model_View7.png">
</details>

- Resonators

The result of the CST simulation is fitted with resonators to get
the more global features of the impedance.

The broadband component of the impedance is fitted with 2 resonators
while resonances at low frequency are fitted with 3 resonators.

The minimum Q factor was bound to 0.5, and the fit with extra resonances
was ran taking into account all resonators for the fit.

## Plots

<details>
  <summary>Comparing Tsutsui model with CST</summary>
  <img src="Tsutsui/comparison_tsu_cst.png">
</details>

<details>
  <summary>Resonator fit</summary>
  <img src="Resonators/fitted_broadband.png">
  <img src="Resonators/fitted_broadband_realimag.png">
  <img src="Resonators/fitted_broadband_extra_resonances.png">
  <img src="Resonators/fitted_broadband_extra_resonances_zoom.png">
</details>

## History

To be completed

## Further work

To be completed

- Determine whether transition pieces are installed (see [KFA45](https://indico.cern.ch/event/830722/contributions/3479602/attachments/1870503/3077852/IWG_PS_KFA45_17_MeasSims_FINAL.pdf))
- Confirm 8C11 material data

## References

To be completed

### Simulations

- [S. Persichelli, The beam coupling impedance model of CERN Proton Synchrotron](https://cds.cern.ch/record/2027523)

- [B. Popovic, Longitudinal limitations with LIU-PS RF upgrades and mitigation strategy](https://indico.cern.ch/event/750790/contributions/3108016/attachments/1719965/2776247/RG_PS_Impedance_Model_Meeting.pdf)
- [B. Popovic, SPS Injection Losses Review](https://indico.cern.ch/event/672967/contributions/2753563/attachments/1567554/2471183/PS_Impedance_Status_Injection_Losses_Meeting_301117.pdf)
- [B. Popovic, LIU-PS Beam Dynamics WG meeting \#6](https://indico.cern.ch/event/678530/contributions/2779037/attachments/1555749/2446435/Update_on_PS_Impedance_Model_091117_Final.pdf)
- [B. Popovic, LIU-PS Beam Dynamics WG meeting \#1](https://indico.cern.ch/event/662292/contributions/2704441/attachments/1516648/2367108/BP_PS_Impedance_Model_Update_310817_Final.pdf)

## Credits

To be completed

- Tsutsui: BE-ABP-HSC Serena Persichelli, Elias Metral, Nicolo Biancacci or Benoit Salvant
- Model: B. Popovic
- 8C11 Material: S. Persichelli
