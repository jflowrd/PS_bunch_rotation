# Impedance of Pumping Manifolds with CODD inserted

The PS orbit system (CODD) presently has 40 pick-ups (PUs) dedicated to
trajectory or orbit measurement. These are installed in the pumping manifolds at
the downstream end of 40 of the 100 bending magnets, effectively in the straight
sections with numbers ending in 0, 3, 5 and 7. \[1\]

The resonance near 500 MHz belongs to the pick-up while the resonances above 1 GHz belong to the pumping manifold.

Damping resistors are inserted in the arm of the pumping manifold \[2\]. Their contribution is not included in the model.

## Subfolders and usage

- CST Raw Data

The model is based on the main Magnet Unit pumping manifold,
standard downstream assembly model.

The CODD pick-up is inserted in the model with a simplified geometry to avoid convergence issues.

Damping resistors are not included and may lead to increased Q of the resonances.

The output result is the impedance in frequency domain and
can directly be used in particle simulations. The impedance
is noisy, to obtain a smoother response the fit with resonators
can be preferred.

<details>
  <summary>CST model view</summary>
  <img src="cst_raw_data/cst_model_images/CODD_CST_Model_View1.png">
  <img src="cst_raw_data/cst_model_images/CODD_CST_Model_View2.png">
  <img src="cst_raw_data/cst_model_images/CODD_CST_Model_View3.png">
</details>

- Resonators

The CST output wakefield is fitted with 9 resonators and a contribution of constant ImZ/f.

Smaller resonances were ignored for better overall fit convergence. The rise of impedance suggest that a resonance is present at >1600MHz.

## Plots

<details>
  <summary>Resonator fit</summary>
  <img src="Resonators/fitted_impedance.png">
</details>

## History

To be completed

## Further work

To be completed

- Modelling the remaining pick-ups
- Draw a list of the damping resistors insertion and evaluate their effect

## References

To be completed

- \[1\] [J. Belleman, A Proposal for New Head-Electronics for the PS Orbit Measurement System](http://jeroen.web.cern.ch/jeroen/reports/coddhead.pdf)
- \[2\] [D. Boussard, Observation of microwave longitudinal instabilities in the CPS](https://cds.cern.ch/record/872559/files/cer-002556319.pdf)

## Credits

To be completed

- CST model: B. Popovic
