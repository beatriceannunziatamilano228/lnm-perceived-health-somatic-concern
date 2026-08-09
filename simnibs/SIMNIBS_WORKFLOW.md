# SimNIBS electric-field workflow

The electric-field panels were generated with **SimNIBS 4.1.0**, the distributed `m2m_MNI152` example head model, isotropic conductivities, and the `MagVenture_MCF-B65_new.ccd` coil model.

Two final runs were documented:

| Target | Placement | Recorded coil–cortex distance | E-field 99.9th percentile |
|---|---|---:|---:|
| medial frontal / dACC | Fpz placement matrix | 19.26 mm | 0.888 V/m |
| lateral OFC | Fp1 placement matrix | 18.12 mm | 0.963 V/m |

The exact 4 × 4 `matsimnibs` matrices are stored in `parameters.json` and `reproduce_simulations.py`. The recovered Fp1 session explicitly recorded `dI/dt = 1,000,000 A/s` and a session distance of 4 mm. The final Fpz log did not explicitly retain `dI/dt`; the reproduction script uses the same 1,000,000 A/s setting as the contemporaneous Fp1 session. This affects absolute field magnitude but not the saved coil-placement geometry.

The original raw log and session files are not distributed because they contain personal local filesystem paths. The public parameter files retain the scientifically relevant settings without those paths.
