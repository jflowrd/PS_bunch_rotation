# Impedance of Flanges and Step Transitions

The impedance of the vacuum flanges and step transitions from one vacuum chamber
to another.

The flanges are used to preserve the continuity of the vacuum chamber and
are distributed around the ring. Large gaps/steps can lead to high frequency resonances
and inductive impedance at low frequency.

In \[1\] the counting of the flanges was the following (259 flanges in total):

- PS 195, 179 elements, external diameter 195 mm; +2mm gap;
- PS 250, 56 elements, external diameter 250 mm;
- SPS 273, 9 elements, external diameter 332 mm; +15mm gap;
- SPS 159, 8 elements, external diameter 206 mm;
- SPS 219, 7 elements, external diameter 266 mm;

The internal diameter is 148mm corresponding to the large axis of the elliptic vacuum chamber.
The default gap width is 2.5mm on which additional hollow is present for some flanges (in the list above).
The PS/SPS types of flanges can be recognized by their collars (CDD ref 8095.0277.1B).

**The impedance of the gaps is partially included in the model (flange included in the model of the standard downstream assembly).**

**The impedance of steps was evaluated in \[2\] to be ImZ/n=0.96 Ohm with f0=476.82 kHz.**

Around each main magnet units, the flanges are isolated to reduce the beam induced current in ground loops. The insulation
lead to resonances at low frequency, which are mitigated using RF-bypasses. All flanges on the
upstream/downstream ends of the main Magnet Units are isolated and equipped with an RF-bypass (actual counting to be checked).
The measurements are stored on [EDMS](https://edms.cern.ch/ui/#!master/navigator/document?D:100063527:100063527:subDocs).

**The other flanges are not isolated (named "metallic flanges"). These are only partially included in the model.**

**The resonance at >40MHz from the circuit model of the rf bypasses needs to be checked.**

## Subfolders and usage

- Ground loops

The impedance of the ground loops on the isolated flanges is modeled with a circuit model
with the parameters

    - R = 100  # Ohm
    - C = 1000e-12  # Farad
    - L = 10e-6  # Henri

The RF-bypass is modeled as a circuit in series with the following parameters

    - R1 = 1  # Ohm
    - C1 = 0.4e-6  # Farad
    - L1 = 12e-9  # Henri

## Plots

<details>
  <summary>Ground loops and RF-bypasses</summary>
  <img src="ground_loops/circuit_model/abs_impedance.png">
  <img src="ground_loops/circuit_model/abs_impedance_zoom.png">
  <img src="ground_loops/circuit_model/abs_impedance_z_over_n.png">
  <img src="ground_loops/circuit_model/abs_impedance_z_over_n_zoom.png">
</details>

<details>
  <summary>Resonator fit of ground loops and RF-bypasses</summary>
  <img src="ground_loops/Resonators/fitted_impedance_0_bypasses.png">
  <img src="ground_loops/Resonators/fitted_impedance_1_bypasses.png">
  <img src="ground_loops/Resonators/fitted_impedance_2_bypasses.png">
</details>

## History

To be completed

## Further work

To be completed

- Complete a listing of the flanges and identify standalone flanges that wouldn't be included in the model
- Adjust the inductance of the modeled rf-bypass to match the measured resonant frequency
- Amplitude and frequency of the resonance at ~45MHz with rf bypass to clarify

## References

To be completed

- \[1\] [S. Persichelli, The beam coupling impedance model of CERN Proton Synchrotron](https://cds.cern.ch/record/2027523/files/CERN-THESIS-2015-076.pdf)
- \[2\] [M. Migliorati et al., Beam-wall interaction in the CERN Proton Synchrotron for the LHC upgrade](https://journals.aps.org/prab/pdf/10.1103/PhysRevSTAB.16.031001)
- [H. Damerau, Impedance Working Group meeting #8](https://indico.cern.ch/event/616762/contributions/2489490/attachments/1429011/2194449/PSRFBypassesForImpedanceWorkingGroupMeeting.pdf)
- [R. Cappi, RF bypass on the proton synchrotron vacuum chamber flanges](https://cds.cern.ch/record/196669/files/CM-P00059176.pdf)

## Credits

To be completed

