MATLAB path-shadow used when generating `test/reference/raycon_reference.json`:
`dmagnetic.m` is the upstream file with the finite-difference step changed from
1e-8 to 1e-5, matching the Julia port. The upstream 1e-8 step puts
~4·eps·|b|/h² ≈ 10–20% roundoff noise on the second derivatives of |B|, which
propagates into the mode-conversion Hessian and shifts η²/τ by ~15% at the
C-Mod reference point (measured 2026-07-05: η² = 0.2435 noisy vs 0.2804
matched-step, with T2 agreeing to 17 digits once the step is matched).
Add this directory to the MATLAB path BEFORE the raycon directory.
