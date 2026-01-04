# Equation verification checklist (Pickett/Simmons tunnel-barrier memristor)

This note is meant to accompany the code implementation and make the mapping
between paper equations and the SPICE/PySpice implementation explicit.

## What the paper provides (and what it points to)
The paper states the device i–v curves are modeled by an **Ohmic resistor in
series with a Simmons tunnel barrier**, and that the *full set of expressions*
for the equivalent circuit are **included in the supplementary information**.

That means: to verify *every* term used in the Simmons/mean-barrier expressions,
you must treat the **supplementary equations** as the authoritative source.

## Mapping used by the implementation

### Equivalent circuit
- Series resistance `Rseries` in series with the tunnel element
- Tunnel current is implemented as a behavioral current source `Bi` using
  `I = sgn(V) * |I(V,w)|`

### State variable
- `x ≡ w` (tunnel gap width, meters), clamped to `[a_off, a_on]`
- Output node `x` exposes `w(t)` directly (V(x) = w)
- Output node `dxdt` exposes `dw/dt` directly (V(dxdt) = dw/dt)

### Simmons / mean barrier height terms
The implementation follows the common Pickett-style structure:
- `lambda_e(w) = lambdaw_e / w`
- `w2(w,v)` with the voltage-dependent denominator
- `dw(w,v) = w2 - w1`
- `Phi1_e(w,v)` mean barrier height expression with a logarithmic term
- `I(w,v) = j0Ae/dw^2 * ( term1 - term2 )`

## Numerical stability fixes (intentional deviations)
ngspice can produce non-physical values from roundoff (e.g., negative log args),
which then explode due to ln/sqrt/division.

The implementation therefore uses:
- `safe_ln(x) = ln(max(x, eps))` (avoid ln(0))
- `safe_sqrt(x) = sqrt(max(x, eps))` (avoid sqrt(negative))
- `safe_div(a,b) = a/max(|b|, eps)` (avoid division by ~0)
- explicit floors for barrier-height terms (`eps_phi`)

These are **numerical guards** and should not change results in the normal range.

## What to verify in your report
1. Cite the paper statements about the equivalent circuit and the use of
   Simmons tunneling (plus the supplement pointer).
2. Cite the supplement for the full expressions (Eq. 20–27 style set).
3. Show that `w(t)` stays within `[a_off,a_on]` and that the model reproduces:
   - pinched hysteresis I–V
   - monotonic potentiation/depression under write–read trains
   - plausible DRM curves (|dw/dt| vs w) with polarity dependence
