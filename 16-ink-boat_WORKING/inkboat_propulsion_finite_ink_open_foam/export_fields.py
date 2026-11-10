"""Export the OpenFOAM Marangoni case to a compact NPZ for the Streamlit app.

Usage:
    python export_fields.py /path/to/run/case  [-o data/marangoni_fields.npz]

Requires only numpy. The case must be reconstructed, written in ASCII, and
contain cell centres at time 0 (postProcess -func writeCellCentres -time 0).

Contents (float16 for fields):
    x, y, z            cell-centre coordinates of the structured tank grid [m]
    t                  written times [s]
    surf_c, surf_ux, surf_uy      top cell layer (free surface), shape (T, nx, ny)
    slice_c, slice_ux, slice_uz   mid-plane y = 0, shape (T, nx, nz)
    vol_c, vx, vy, vz             2x2x2-averaged concentration, shape (T, nx/2, ny/2, nz/2)
    force_t, X_boat, U_boat, a_boat       boat surge motion (boat frame solver)
    F_cl_x, F_pressure_x, F_viscous_x, F_net_x   forces on the boat [N]
    frame_X, frame_U                      boat position/velocity at the written times
Velocities are stored in the BOAT frame; add frame_U to u_x for the lab frame.
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np

GAMMA = 1.0e-3        # -dgamma/dc [N/m], same as constant/fvOptions
INK_DURATION = 3.0    # [s], same as system/controlDict
BOAT = dict(x=(-0.015, 0.015), y=(-0.0075, 0.0075), z=(-0.003, 0.0))
BOAT_MASS = 1.347e-3  # [kg], same as system/controlDict (boatMotion)


def read_internal_field(path: Path):
    txt = path.read_text()
    m = re.search(r"internalField\s+(nonuniform|uniform)", txt)
    if m.group(1) == "uniform":
        val = re.search(r"internalField\s+uniform\s+([^;]+);", txt).group(1)
        return None, val
    start = m.end()
    n_match = re.search(r"\s*(\d+)\s*\(", txt[start:])
    n = int(n_match.group(1))
    body_start = start + n_match.end()
    body_end = txt.index("\n)", body_start)
    body = txt[body_start:body_end]
    if "(" in body:  # vector
        arr = np.array(body.replace("(", " ").replace(")", " ").split(), dtype=float).reshape(n, 3)
    else:
        arr = np.array(body.split(), dtype=float)
    assert len(arr) == n, (path, len(arr), n)
    return arr, None


def axis_index(v, tol=1e-6):
    """Cluster coordinates that differ by < tol (ASCII output has ~6 significant digits)."""
    s = np.sort(np.unique(v))
    starts = np.concatenate([[True], np.diff(s) > tol])
    group = np.cumsum(starts) - 1
    centres = np.bincount(group, weights=s) / np.bincount(group)
    idx = group[np.searchsorted(s, v)]
    return centres, idx


def to_grid(vals, idx, shape):
    g = np.full(shape, np.nan, dtype=float)
    g[idx] = vals
    return g


def coarsen(a):
    """nan-aware 2x2x2 block mean over the last three axes."""
    T, nx, ny, nz = a.shape
    b = a[:, : nx // 2 * 2, : ny // 2 * 2, : nz // 2 * 2]
    b = b.reshape(T, nx // 2, 2, ny // 2, 2, nz // 2, 2)
    with np.errstate(invalid="ignore"):
        return np.nanmean(b, axis=(2, 4, 6))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case", type=Path)
    ap.add_argument("-o", "--out", type=Path, default=Path(__file__).with_name("data") / "marangoni_fields.npz")
    args = ap.parse_args()
    case = args.case

    C, _ = read_internal_field(case / "0" / "C")
    xs, ix = axis_index(C[:, 0])
    ys, iy = axis_index(C[:, 1])
    zs, iz = axis_index(C[:, 2])
    shape = (len(xs), len(ys), len(zs))
    idx = (ix, iy, iz)
    print("grid", shape, "cells", len(C))

    times = sorted(
        (d for d in case.iterdir() if d.is_dir() and re.fullmatch(r"[0-9.eE+-]+", d.name) and float(d.name) > 0),
        key=lambda d: float(d.name),
    )
    t = np.array([0.0] + [float(d.name) for d in times])

    jy = np.argsort(np.abs(ys))[:2]  # two rows straddling y = 0
    c_all, ux_s, uy_s, ux_m, uz_m, c_m, vol_u = [], [], [], [], [], [], []
    blank = np.zeros(len(C))
    for d in [None] + times:
        if d is None:
            c, U = blank, np.zeros((len(C), 3))
        else:
            c, cu = read_internal_field(d / "c")
            if c is None:
                c = np.full(len(C), float(cu))
            U, _ = read_internal_field(d / "U")
        cg = to_grid(c, idx, shape)
        Ug = [to_grid(U[:, k], idx, shape) for k in range(3)]
        c_all.append(cg)
        ux_s.append(Ug[0][:, :, -1])
        uy_s.append(Ug[1][:, :, -1])
        c_m.append(np.nanmean(cg[:, jy, :], axis=1))
        ux_m.append(np.nanmean(Ug[0][:, jy, :], axis=1))
        uz_m.append(np.nanmean(Ug[2][:, jy, :], axis=1))
        vol_u.append(np.stack(Ug))
        print(f"t = {t[len(c_all)-1]:.2f}  max c = {np.nanmax(cg):.3f}  max|U| = {np.nanmax(np.linalg.norm(np.stack(Ug), axis=0)):.4f}")

    c_all = np.array(c_all)
    surf_c = c_all[:, :, :, -1]
    vol_c = coarsen(c_all)
    vol_u = np.array(vol_u)  # (T, 3, nx, ny, nz)
    vu = [coarsen(vol_u[:, k]) for k in range(3)]
    pair_mean = lambda a: a[: len(a) // 2 * 2].reshape(-1, 2).mean(1)
    xs2, ys2, zs2 = pair_mean(xs), pair_mean(ys), pair_mean(zs)

    # ---- boat dynamics, integrated inside the solver (boatMotion function object):
    # t, X, U, a, F_contactLine, F_pressure, F_viscous   (surge, SI units)
    bm = np.loadtxt(case / "postProcessing" / "boatMotion.dat", comments="#")
    force_t, X_boat, U_boat, a_boat, F_cl, F_p, F_v = bm.T
    frame_X = np.interp(t, force_t, X_boat, left=0.0)
    frame_U = np.interp(t, force_t, U_boat, left=0.0)
    print(f"boat: max U = {100*U_boat.max():.2f} cm/s, distance = {100*X_boat[-1]:.2f} cm")

    f16 = lambda a: np.nan_to_num(a, nan=np.nan).astype(np.float16)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.out,
        x=xs, y=ys, z=zs, t=t, x2=xs2, y2=ys2, z2=zs2,
        surf_c=f16(surf_c), surf_ux=f16(np.array(ux_s)), surf_uy=f16(np.array(uy_s)),
        slice_c=f16(np.array(c_m)), slice_ux=f16(np.array(ux_m)), slice_uz=f16(np.array(uz_m)),
        vol_c=f16(vol_c), vol_ux=f16(vu[0]), vol_uy=f16(vu[1]), vol_uz=f16(vu[2]),
        force_t=force_t, X_boat=X_boat, U_boat=U_boat, a_boat=a_boat,
        F_cl_x=F_cl, F_pressure_x=F_p, F_viscous_x=F_v, F_net_x=F_cl + F_p + F_v,
        frame_X=frame_X, frame_U=frame_U, boat_mass=BOAT_MASS,
        boat=np.array([*BOAT["x"], *BOAT["y"], *BOAT["z"]]),
        gamma_coeff=GAMMA, ink_duration=INK_DURATION,
    )
    print("wrote", args.out, f"{args.out.stat().st_size/1e6:.1f} MB")


if __name__ == "__main__":
    main()
