# Ink Boat — Marangoni propulsion (reduced model)

Run `pip install -r requirements.txt` and `streamlit run app.py`.

**Scope:** top-view effective surface tension around a circular hull, integrated to obtain net capillary force and one-dimensional translation. Angular concentration evolves by imposed rear source, angular mixing and relaxation; it is NOT an ink-spreading solver. The three interfacial tensions are represented individually and used to calculate the thermodynamic spreading coefficient `S = gamma_wa - gamma_ia - gamma_iw`. They are NOT sufficient to determine the force without a contact-line geometry and spatially varying tension. Thus `delta_gamma` is a separate *phenomenological* parameter. No claim of a validated three-phase Navier–Stokes solver.

Equations:
- `gamma_eff(theta) = gamma_wa - delta_gamma * c(theta)`
- `F = R integral_0^{2pi} gamma_eff(theta) (cos theta, sin theta) dtheta`
- `m dv/dt = Fx - b v`
- `dc/dt = q(theta)(1-c) - c/tau + Dtheta * d²c/dtheta²`

The line integral assumes a circular contact line and constant geometric projection; contact angles, hull wettability, wave drag, wake hydrodynamics, rotation and recoil are neglected. The computed `Fy` is displayed nowhere in the 1D dynamics. Use only for order-of-magnitude and qualitative tests, pending calibration and validation.
