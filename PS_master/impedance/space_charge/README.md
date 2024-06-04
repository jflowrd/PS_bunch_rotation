# Space charge impedance

Model for longitudinal space charge impedance along the ramp.

<details>
  <summary>Space charge along the ramp for the proton beam</summary>
  <img src="figures/im_z_over_n_along_ramp.png">
</details>

Takes into account the aperture information as gathered in 2013 and obtained from
/eos/project/l/liu/PS/Aperture

The aperture is by default assumed rectangular (see ref. about geometrical factors),
which is a close assumption to the actual elliptical aperture in the PS.

The optics is obtained from the https://gitlab.cern.ch/injector-optics/PS repository (input from 2017)

The load_model function requires the normalized transverse emittance (in m rad),
the particle mass, momentum and momentum spread of the beam (can be passed as arrays to
get the ImZ/n along the ramp). 

## Further work

- Refresh the model with more recent aperture info
- Refresh the model with more recent optics info
- Improve the modeling with different aperture geometry
- Improve the modeling by accounting for particle trajectory with dispersion

## References

### Geometrical factors
- [L. Wang and Y. Li, Analysis of the longitudinal space charge impedance of a round uniform beam inside parallel plates and rectangular chambers](https://journals.aps.org/prab/abstract/10.1103/PhysRevSTAB.18.024201)

