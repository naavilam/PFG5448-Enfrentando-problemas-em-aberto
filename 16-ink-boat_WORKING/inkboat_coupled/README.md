# Ink Boat — coupled surface transport + Stokes subphase (research prototype)

Run: `streamlit run inkboat_coupled/app.py` from repository root.

## What is actually solved
- A **2D insoluble surface species** Gamma on a flat water–air interface with conservative upwind advection, diffusion, prescribed ink release, and removal.
- **Depth-resolved 3D Stokes half-space response** reduced analytically to a nonlocal Fourier mobility at the interface: `u_s(k) = [I - kk/(2k²)] tau_s(k)/(mu |k|)`, with `tau_s = grad_s gamma` and `gamma=gamma_wa-beta Gamma`.
- Contact-line capillary force sampled around a circular boat; Newtonian translation with **approximate** Stokes drag `6*pi*mu*R*c_d`.

## What is NOT solved
- Full time-dependent Navier–Stokes; water inertia, finite depth, walls, and no-slip boat boundary.
- Separate moving water–ink, ink–air and water–air interfaces (gamma_ia and gamma_iw only enter spreading coefficient diagnostic).
- Full hydrodynamic stress integral on submerged boat; rotation, contact angles and triple line dynamics.
- A chemically validated equation of state for real ballpoint ink. Gamma is an *effective surface-active species*, not necessarily ink mass.

## Numerical notes
- Periodic horizontal domain, deep subphase; zero Fourier mode fixed to zero.
- Positive upwind advection, centered explicit diffusion, CFL-controlled steps.
- Source is localized outside the circular boat. Boat is not excluded from the Stokes flow solve, and contact-line motion is not coupled as a no-slip immersed boundary.
- If `Re` is not much less than one, the Stokes approximation is inconsistent. Use this version for code validation and model exploration, not quantitative prediction.
- Validation tests: zero species implies zero force, surface source changes mass consistently, positive surface tension gradient induces surface velocity.
