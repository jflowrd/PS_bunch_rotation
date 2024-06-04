# Impedance of Wire Scanners

There are 4 wire scanners located in straight sections 54, 64, 75, 85.

The impedance model correspond to the new wire scanners, installed during LS2.
The older versions were not modeled.

The impedance was calculated using CST with the eigenmode solver
for two positions for the wire:

- 0 deg (parking)
- 135 deg (wire in beam)

The model is identical to the SPS one (the eigenmode computed for the
SPS impedance model are used).

<details>
  <summary>WS in the ring</summary>
  <img src="http://cern.ch/psring/psring/pictures/fullsize/ss54.jpg">
  <img src="http://cern.ch/psring/psring/pictures/fullsize/ss64.jpg">
  <img src="http://cern.ch/psring/psring/pictures/fullsize/ss75.jpg">
  <img src="http://cern.ch/psring/psring/pictures/fullsize/ss85.jpg">
</details>

## Subfolders and usage

- CST Raw Data

The wakefield and eigenmodes are obtained using CST.

The eigenmodes are saved in a text file. The modes can be loaded
as resonators. The case where the wire is parked (0 deg) is the more
relevant cases for instability simulations.

The eigenmodes can directly be used as resonators.

<details>
  <summary>CST model view</summary>
  <img src="cst_raw_data/WS_CST_Model_View1.png">
</details>


## Plots

<details>
  <summary>Eigenmodes parked (0 deg)</summary>
  <img src="cst_raw_data/PS_SPS_wire_scanner.dat.png">
</details>

<details>
  <summary>Eigenmodes wire in beam (135 deg)</summary>
  <img src="cst_raw_data/PS_SPS_ws_135deg_fork_rot_fRQRoQ_EIG.txt.png">
</details>

## History

To be completed

## Further work

To be completed

- Availability of the model?
- Model identical to the SPS one, but eigenmode files different, to be clarified

## References

To be completed

- [T. Kaltenbacher, LIU BWS Mechanics Production Readiness Review](https://indico.cern.ch/event/618262/contributions/2509953/attachments/1435689/2207628/CVollingerLIUWSReview29-Mar-2017.pdf)
- [T. Kaltenbacher, Impedance Working Group meeting #7](https://indico.cern.ch/event/613751/contributions/2474294/attachments/1414608/2165243/17022017_Kaltenbacher_LIU_WS.pdf)

## Credits

To be completed
