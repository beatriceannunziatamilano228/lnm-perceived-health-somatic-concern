#!/usr/bin/env python3
"""Recreate the two reported SimNIBS 4.1 TMS placement simulations.

Requires a local SimNIBS 4.1 installation and its MNI152 example head model.
The script uses exact saved coil-placement matrices. It does not require the
participant-level lesion data or the normative connectome.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np

FPZ_DACC = np.array([
    [0.0, 0.99926, 0.038343, -0.70611],
    [0.079782, 0.038221, -0.99608, 89.465],
    [-0.99681, 0.0030591, -0.079723, -4.6399],
    [0.0, 0.0, 0.0, 1.0],
])
FP1_OFC = np.array([
    [-1.38777878e-17, 0.897386211, 0.441245949, -31.2153131],
    [0.220737822, 0.430361815, -0.875250548, 84.1104250],
    [-0.975333181, 0.0973996697, -0.198087078, -8.39749825],
    [0.0, 0.0, 0.0, 1.0],
])

def make_session(sim_struct, subject_dir: Path, output_dir: Path, coil_file: Path,
                 name: str, matrix: np.ndarray):
    session = sim_struct.SESSION()
    session.subpath = str(subject_dir)
    session.fnamehead = str(subject_dir / "MNI152.msh")
    session.pathfem = str(output_dir / name)
    session.map_to_vol = True
    session.map_to_MNI = True
    session.map_to_surf = False
    session.map_to_fsavg = False
    session.tissues_in_niftis = 2
    session.fields = "eED"
    tms = session.add_tmslist()
    tms.fnamecoil = str(coil_file)
    tms.anisotropy_type = "scalar"
    position = tms.add_position()
    position.name = name
    position.matsimnibs = matrix
    position.didt = 1_000_000.0
    return session

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--subject-dir", type=Path, required=True,
                        help="Path to the SimNIBS m2m_MNI152 directory.")
    parser.add_argument("--coil-file", type=Path, required=True,
                        help="Path to MagVenture_MCF-B65_new.ccd.")
    parser.add_argument("--output-dir", type=Path, default=Path("simnibs_outputs"))
    args = parser.parse_args()
    try:
        from simnibs import sim_struct, run_simnibs
    except ImportError as exc:
        raise SystemExit("Run this script inside a SimNIBS 4.1 environment.") from exc
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for name, matrix in [("Fpz_dACC", FPZ_DACC), ("Fp1_OFC", FP1_OFC)]:
        session = make_session(sim_struct, args.subject_dir.resolve(),
                               args.output_dir.resolve(), args.coil_file.resolve(),
                               name, matrix)
        run_simnibs(session)

if __name__ == "__main__":
    main()
