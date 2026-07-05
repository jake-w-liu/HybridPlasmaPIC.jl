# RAYCON → HybridPlasmaPIC.jl Porting Notes (working reference)

Source: `/Users/jake/PlasmaWorkspace/raycon` (MATLAB, Jaun–Kaufman–Tracy, v7.0 2006,
plus Steve Richardson's modifications). All 20 `.m` files read end-to-end on
2026-07-05. This file is the durable extraction of every algorithm needed for the
Julia port in `src/raycon/`, including upstream quirks. MATLAB R2026a is available at
`/Applications/MATLAB_R2026a.app/bin/matlab` for reference-data generation
(`matlab -batch`).

References: Tracy, Kaufman & Jaun, Phys. Lett. A 290 (2001) 309; Jaun, Tracy &
Kaufman, PPCF 49 (2007) 43; Tracy, Kaufman & Jaun, PoP 14 (2007) 082102.

## Architecture (verified data flow)

- Driver (`main.m`): TYPE='Con', odeDim=4, state z=(R, Z, kR, kZ) in SI, toroidal
  wavenumber kφ = kant(2) held CONSTANT (upstream convention). Antenna IC via
  `adjust_disp_m([r,z,kant(1),kant(3)], m=0)`. Time span 5e-2 (units 1/freq… it is
  the ODE σ parameter span; rays traced as functions of ray parameter σ, not t).
- Propagation: `ode45('trajectory', tspan, y, odeset(AbsTol=1e-7(Con)/[1e-8 1e-8 1e-6
  1e-6](Trj), RelTol=1e-6, Events=on, InitialStep=1e-7*timespan, Refine=16))`.
  RHS = `disp_eig(y,'Trj')`. After each segment: conversion detection, ray split,
  repeat until σ>timespan or monErrAbort or numConvert≥3.
- Conversion analysis uses `dispertok` (det-based U) at the stopped point with
  velocity/acceleration from 3-point divided differences of the stored trajectory.

## Units and constants (initCnst — keep EXACT for 1:1 parity)

c=2.9979e8, e=1.6022e-19, mp=1.6726e-27, eps0=8.8542e-12 (SI). Frequencies rad/s;
k in 1/m; B Tesla; n 1/m³; T keV (vth² = 2000·T·e/(amass·mp) — the 2000 = keV→eV
×2? actually vth²=2eT/m with T in keV → 2·1000·T·e/m; upstream writes 2000·T·e).

## data('cmod') parameters (primary test case)

EQ='Solovev', MODEL='cld2x2', PROBL='tok', b0=7.9 T, r0=0.67 m, q0=2.0,
iaspr=0.22/0.67, elong=1.6, amass=[1/1836, 2, 3], acharge=[-1, 1, 2],
n0=[10, 5.2, 2.4]·1e19 m⁻³, na=[1, .7, .7], nb=[3,3,3], t0=[3,3,3] keV,
ta=[1,1,1], tb=[1,1,1], sant=0.4, thant=[-.61, .6], freq=80e6 Hz,
kant=[-31.5, 10, 0] (main.m overrides kant=[-31.5, -10, 0], thant=[-.5,.5],
th0=0.001, timespan=5e-2, aThres=0.05, NRAY=1).
Derived: psin = 0.5·b0/q0·elong·(r0·iaspr)², omega = 2π·freq.

## solovev.m (equilibrium, flux label s)

Inputs (ρ, θ, r0, ε=iaspr, E=elong, sflxa). r = ρcosθ + r0, z = −ρsinθ.
fac = psin/(ε·r0²)²; ψ = fac·(r²z²/E² + ¼(r²−r0²)²); s = √(ψ/psin) − sflxa.
- 1st: dψdr = fac·r·(2z²/E² + (r²−r0²)); dψdz = fac·z·(2r²/E²); dψds = 2·s·psin
  (NB uses s AFTER subtracting sflxa only for root finding; for derivatives
  sflxa=0 in all magnetic.m calls); ds/dr = dψdr/dψds, ds/dz = dψdz/dψds.
- 2nd: dψdr2 = fac(2z²/E² + 3r² − r0²); dψdrz = fac·4rz/E²; dψdz2 = fac·2r²/E².
  dsdr2 = dsdp(dψdr2 − dsdp/s·dψdr²), dsdrz = dsdp(dψdrz − dsdp/s·dψdr·dψdz),
  dsdz2 = dsdp(dψdz2 − dsdp/s·dψdz²) with dsdp = 1/dψds.
- 3rd: dψdr3=fac·6r, dψdr2z=fac·4z/E², dψdrz2=fac·4r/E², dψdz3=0;
  dLNsdp=dsdp/s; dsdrp=−dLNsdp·dsdp·dψdr; dLNsdrp=dLNsdp·(dsdrp/dsdp−dsdr/s);
  dsdzp=−dLNsdp·dsdp·dψdz; dLNsdzp=dLNsdp·(dsdzp/dsdp−dsdz/s);
  dsdr3 = dsdr2·(dsdrp/dsdp) + dsdp·(dψdr3 − dLNsdrp·dψdr² − 2·dLNsdp·dψdr·dψdr2)
  dsdr2z= dsdr2·(dsdzp/dsdp) + dsdp·(dψdr2z − dLNsdzp·dψdr² − 2·dLNsdp·dψdr·dψdrz)
  dsdrz2= dsdz2·(dsdrp/dsdp) + dsdp·(dψdrz2 − dLNsdrp·dψdz² − 2·dLNsdp·dψdz·dψdrz)
  dsdz3 = dsdz2·(dsdzp/dsdp) + dsdp·(dψdz3 − dLNsdzp·dψdz² − 2·dLNsdp·dψdz·dψdz2)

θ convention EVERYWHERE: ρ=√((r−r0)²+z²), θ = atan2(r0−r, −z) + π/2
(so r=r0+ρcosθ, z=−ρsinθ).

## mapFlux.m

fzero of solovev(ρ;θ,…, sflxa=s) over ρ∈[0, 1.5·ρa], ρa = elong·iaspr·r0, TolX 1e-4.
Julia: bisection/Brent to 1e-12 (upgrade, documented). Returns ρ, r, z.

## magnetic.m (+dmagnetic.m)

From solovev derivs (sflxa=0): ψ=psin·s², dψds=2s·psin, dψds2=2psin, t0=b0·r0
(toroidal field function, dt0ds=0).
- dtdr=−sinθ/ρ, dtdz=−cosθ/ρ, dtdr2=2sinθcosθ/ρ², dtdrz=(cos²−sin²)/ρ²,
  dtdz2=−2sinθcosθ/ρ².
- jac = r/(dsdz·dtdr − dsdr·dtdz); drds=−jac/r·dtdz; dzds=jac/r·dtdr;
  drdt=jac/r·dsdz; dzdt=−jac/r·dsdr.
- h11=dsdr²+dsdz²; dh11dr=2(dsdr·dsdr2+dsdz·dsdrz); dh11dz=2(dsdr·dsdrz+dsdz·dsdz2);
  dh11ds/dt via chain with drds…; gp2=dψds²·h11; dgp2ds=dψds²(dh11ds+2h11/s);
  dgp2dt=dψds²·dh11dt.
- b = √((t0/r)² + gp2/r²); a2=t0²+gp2; dbds = da2ds/(2b r²) − b·drds/r (da2ds=dgp2ds);
  dbdt analogous.
- Basis vectors (n=∇s direction, b=B̂, p=b̂×n̂ — components in cylindrical (r,φ,z)):
  gs=√h11; ener=dsdr/gs, enef=0, enez=dsdz/gs;
  eber=t0·dsdz/(r·b·gs), ebef=−dψds·gs/(r·b), ebez=−t0·dsdr/(r·b·gs);
  eper=dψds·dsdz/(r·b), epef=t0/(r·b), epez=−dψds·dsdr/(r·b).
- dbdr,dbdz,dbdr2,dbdrz,dbdz2 from dmagnetic: CENTRAL FD of b(R,Z) with
  dr=dz=1e-8 on a 3×3 stencil (scalar). Julia: same FD (parity) — analytic possible
  but keep FD for 1:1, step 1e-8 hardcoded upstream. NOTE 1e-8 with b~8 →
  ~1e-2..1e-4 relative accuracy on 2nd derivs; acceptable upstream choice.
  Then dbds2, dbdt2, dbdst by chaining with drds2 etc.:
  drds2=−(djacds·dtdz − jac·drds·dtdz/r + jac·dtdzs)/r, etc. (see magnetic.m
  lines 96-110; djacdr=jac/r·(1−jac·(dsdrz·dtdr+dsdz·dtdr2−dsdr2·dtdz−dsdr·dtdrz)),
  djacdz=jac/r·(−jac·(dsdz2·dtdr+dsdz·dtdrz−dsdrz·dtdz−dsdr·dtdz2))).
- Basis-vector 1st derivatives wrt (r,z) [Steve R]: with ζ=1/(rb),
  dζdr=−ζ(1/r+dbdr/b), dζdz=−ζ·dbdz/b; mag=1/gs, dmagdr=−½mag³·dh11dr, …
  denerdr = dsdr2/gs + dsdr·dmagdr; … (full list magnetic.m lines 152-186;
  dψds derivative: dpdsr = dsdr·dψds2 with dψds2=2psin).
- Basis-vector 2nd derivatives (lines 199-297): needed only for disp_eig eval2nd
  (dU_mat, i.e. 20-dim symplectic tracing). Port complete.
- Output convention groups by nargout: 18 / 31 / 49 / 76 outputs. Julia: single
  struct with levels (:first, :second) to avoid recomputation.

## Profiles (dispertok/disp_eig identical)

p_k = 1 − s²·na_k; n_k = n0_k·p_k^nb_k; T_k = t0_k·(1 − s²·ta_k)^nb_k
⚠ UPSTREAM QUIRK: temperature exponent uses **nb**, not tb (both files) — port
as-is with a comment (tb is dead config).
dLNn/ds = −2s·na·nb/p; dLNn/ds² (their dLNpds2) = dLNpds/s − dLNpds²/nb.
ωp²_k = (acharge_k·e)²·n_k/(amass_k·mp·eps0); ωc_k = b·acharge_k·e/(amass_k·mp)
(signed!); caoc2 = 1/Σ(ωp²/ωc²).
dLNωc ds/dt = (dbds/b, dbdt/b); 2nd: dLNomcds2=(dbds2/b−(dbds/b)²)…

## Stix elements (cold), ALL models

omc2Mom2 = ωc²−ω² (iomceps=0 in practice); S = 1+Σ Sᵢ, Sᵢ=ωp²/omc2Mom2;
D = Σ Dᵢ, Dᵢ=(ωc/ω)·Sᵢ (signed ωc!); P = 1−Σ Pᵢ, Pᵢ=ωp²/ω².
1st derivatives (log form): dLNSids = dLNnds − 2ωc²/omc2Mom2·dLNomcds;
dLNSidt = −2ωc²/omc2Mom2·dLNomcdt; dLNSidom = 2ω/omc2Mom2;
dLNDids = dLNSids+dLNomcds; dLNDidom = (3ω²−ωc²)/(ω·omc2Mom2);
dLNPids = dLNnds; dPds = −Σ Pᵢ·dLNPids; dPdom = +Σ Pᵢ·2/ω  ⚠ SIGN QUIRK:
dispertok has dLNPidom=2/ω, dPdom=+Σ(Pi·dLNPidom) [WRONG SIGN vs P=1−ΣPi ⇒
dPdom = +Σ2Pi/ω is CORRECT since dPi/dom=−2Pi/ω → dP/dom=+2ΣPi/ω. OK both agree];
disp_eig cld3x3 has dLNPidom=−2/ω AND dPdom=−sum(Pi·dLNPidom) = +2ΣPi/ω — same
result. Fine.
2nd derivatives: omOM=ω²/omc2Mom2, ocOM=ωc²/omc2Mom2;
dLNSids2 = 2·ocOM·(2·omOM·dLNomcds² − dLNomcds2) + dLNnds2; (st, t2 analogous, no
dLNnds2 term for t); dSids2 = Sᵢ(dLNSids2 + dLNSids²); dLNDids2 = dLNomcds2 +
dLNSids2; dDids2 = Dᵢ(dLNDids2 + dLNDids·dLNDids); P: dLNPids2 = dLNnds2,
dPids2 = Pᵢ(dLNPids2+dLNPids²), dPds2 = −Σ (sign per P=1−ΣPᵢ; disp_eig cld3x3 has
dPds2=−sum(dPids2) ✓).

## Wave vector projections

coom=c/ω; kn = kr·ener + kf·enef + kz·enez (kf=kant(2) const); kb, kp analogous;
Nn=coom·kn etc.; Nr=coom·kr, Nf=coom·kf, Nz=coom·kz.

## Tensor (rotating? NO — Stix frame (n,b,p) with n=∇s/|∇s|, b=B̂, p):

D11 = Nb²+Np²−S; D12 = −Nn·Nb − iD; D22 = Nn²+Np²−S; D13 = −Nn·Np;
D23 = −Nb·Np; D33 = Nn²+Nb²−P. DD Hermitian (D21=cD12 etc.).
k-derivatives (in kn,kb,kp): dD11dkn=0, dD11dkb=2Nb·coom, dD11dkp=2Np·coom;
dD12dkn=−Nb·coom, dD12dkb=−Nn·coom, dD12dkp=0; dD22dkn=2Nn·coom, dD22dkb=0,
dD22dkp=2Np·coom; dD13dkn=−Np·coom, dD13dkb=0, dD13dkp=−Nn·coom; dD23dkn=0,
dD23dkb=−Np·coom, dD23dkp=−Nb·coom; dD33dkn=2Nn·coom, dD33dkb=2Nb·coom, dD33dkp=0.
ω-derivatives: dD11dom=−(2/ω)(Nb²+Np²)−dSdom; dD12dom=−(2/ω)(NnNb)−i·dDdom;
dD22dom=−(2/ω)(Nn²+Np²)−dSdom; dD13dom=−(2/ω)(NnNp); dD23dom=−(2/ω)(NbNp);
dD33dom=−(2/ω)(Nn²+Nb²)−dPdom.
s,t-derivatives: dD11ds=−dSds, dD12ds=−i·dDds, dD13ds=0, dD23ds=0, dD33ds=−dPds…

## disp_eig (ODE Hamiltonian) — PRIMARY for tracing

U = eigenvalue of Hermitian DD nearest 0; pol = its eigenvector.
Monitors: mon1=0 (odeDim 4); mon2 = |D11+D22| (2x2) or |sum principal 2×2 minors|
(3x3). Output 'Mon' = [mon1, mon2].
Exact eigenvalue gradient (2x2): with subsum = D11+D22−2U,
dU_vec = [dD11_vec·(D22−U) + dD22_vec·(D11−U) − 2Re(cD12·dD12_vec)]/subsum
where each dDij_vec = [dom, dr, dz, dkr, dkz] and SPATIAL derivatives include
curvature corrections:
dNndr = Nr·denerdr + Nf·denefdr + Nz·denezdr (etc.);
dD11dr = dD11ds·dsdr + dD11dt·dtdr + 2(Nb·dNbdr + Np·dNpdr); (mirror z)
dD12dr = dD12ds·dsdr + dD12dt·dtdr − Nn·dNbdr − Nb·dNndr;
dD22dr = … + 2(Nn·dNndr + Np·dNpdr).
k-space: dD11dkr = dD11dkn·ener + dD11dkb·eber + dD11dkp·eper (etc. — NOTE these
use e{n,b,p}{r} as ∂k{n,b,p}/∂kr).
3x3: U1=det(DD2[23,23]), U2=det(DD2[13,13]), U3=det(DD2[12,12]) with DD2=DD−U·I,
subsum=U1+U2+U3;
dU_vec = [dD11_vec·U1 − (D11−U)·2Re(cD23·dD23_vec) + dD22_vec·U2 −
(D22−U)·2Re(cD13·dD13_vec) + dD33_vec·U3 − (D33−U)·2Re(cD12·dD12_vec)]/subsum
+ 2Re(D12·D23·conj(dD13_vec) + D12·dD23_vec·cD13 + dD12_vec·D23·cD13)/subsum.
2nd derivative dU_mat (4×4 over (r,z,kr,kz)): full expressions in disp_eig lines
280-320 (2x2) and 560-700 (3x3), including dD11_mat etc. built from products of
basis-vector 1st+2nd derivatives; symmetry self-check (3x3) tol 1e-10.
'Trj' output: J=[0 I;−I 0] (2×2 blocks); dz/dσ = J·dU(2:5) — REAL part; if
20-dim state, also dS/dσ = J·dU_mat·S (tangent flow). σ is ray parameter;
physical time direction from sign of dUdom ('Sgn').

## dispertok — CONVERSION machinery (det-based U, cld2x2 path to port exactly)

U = D11·D22 − |D12|²; V = ½(D11+D22); eig2 = V − sign(V)·√(V²−U); T=U.
mon2 (dispertok) = log|V| (Steve's version).
1st derivs: dUdx = dD11dx·D22 + D11·dD22dx − 2Re(D12·conj(dD12dx)) for
x∈{s,t,kn,kb,kp,om}; then to (r,z,kr,kf,kz) via chain (NO curvature corrections in
dispertok! dUdr = dUds·dsdr + dUdt·dtdr — upstream difference vs disp_eig. Port
as-is for conversion parity; conversion quantities are evaluated near the saddle
where this level of approximation is upstream's).
2nd derivs (eval2nd), needed keys: dUds2, dUdst, dUdt2, dUdkn2, dUdknkb, dUdknkp,
dUdkb2, dUdkbkp, dUdkp2, dUdskn…dUdtkp (mixed) — formulas dispertok lines 377-391;
then chain to (r,z,kr,kz): dUdrs=dUds2·dsdr+dUdst·dtdr, …, dUdr2=dUdrs·dsdr+
dUdrt·dtdr (lines 755-787) ⚠ note upstream's dUdr2 CHAIN OMITS 2nd-derivative
terms of s,t wrt r,z (no dsdr2·dUds terms) — this is upstream behavior; PORT
AS-IS (documented divergence from exact math; the disp_eig dU_mat includes them
properly).
Osculating plane: A = ẏ1·ÿ3 + ẏ2·ÿ4 − ẏ3·ÿ1 − ẏ4·ÿ2 (symplectic product of
velocity & acceleration 4-vectors); eq = ẏ/√|A|, ep = ÿ/√|A|;
dDijdq = eq·gDij, dDijdp = ep·gDij with gDij = [dDijdr, dDijdz, dDijdkr, dDijdkz].
H11 = 2dD11dq·dD22dq − 2dD12dq·conj(dD12dq); H22 same with p;
H12 = dD11dq·dD22dp + dD11dp·dD22dq − 2Re(dD12dq·conj(dD12dp));
dUdq = dD11dq·D22 + D11·dD22dq − 2Re(D12·conj(dD12dq)); dUdp analogous.
Conversion estimates: detH=H11·H22−H12²; Hm = inv; 
eta2 = 0.5/√|detH| · |Hm11·dUdq² + 2Hm12·dUdq·dUdp + Hm22·dUdp²| (Tracy eq 26);
qst = Hm11·dUdq + Hm12·dUdp; pst = Hm21·dUdq + Hm22·dUdp;
zinzst = −(qst·eq + pst·ep) (4-vector); saddle guess z* = z + zinzst ('Sdl'),
ITERATED from ray.m ≤30 times, tol=1e-4 (relative ∞-norm), with FIXED eq,ep
(zdot,zddot from the 3 last trajectory points, divided differences:
dm2=(tm2−tm1)(tm2−t0), dm1=(tm1−tm2)(tm1−t0), d0=(t0−tm2)(t0−tm1);
zdot = zm2/dm2·(t0−tm1) + zm1/dm1·(t0−tm2) + z0/d0·(2t0−tm1−tm2);
zddot = 2(zm2/dm2 + zm1/dm1 + z0/d0)).
Guard: if not converged → abort; ⚠ upstream QUIRK: the "unlikely conversion
point" branch is UNREACHABLE as coded (elseif commented out) — the displayed
Estimate/Iterated and abort happen only on non-convergence. Port: converged →
proceed; not converged → abort conversion (skip split, continue ray).
Normal form at z* ('Mch'/'Cnv'/'Trs' eval): J4 = [0 I2; −I2 0] BUT written
component-wise [0 0 1 0; 0 0 0 1; −1 0 0 0; 0 −1 0 0];
T2 = 4×4 Hessian [dUdr2 dUdrz dUdrkr dUdrkz; …] (det-U based, at z*);
[M,val]=eig(J4·T2), sort by real part ascending; opposite=val4/val1,
separate=|val4/val3|; hyperbola OK iff |1+opposite|<0.1 AND separate>4;
vp=M(:,ind4), vm=M(:,ind1);
gdv rows: [gD11·vp, gD12·vp, conj(gD12)·vp, gD22·vp] and same with vm → reshape
2×2; uncoupled polarizations by power iteration (tol 1e-4, ≤100 it, start ones);
gD = [gD11ᵀ gD12ᵀ cgD12ᵀ gD22ᵀ] (4×4: rows=phase-space dirs, cols=ij);
gdalf = gD·vec(conj(pol1·pol1ᴴ)); gdlam = gD·vec(conj(pol2·pol2ᴴ));
braket = gdalfᴴ·J4·gdlam; eta = pol1ᴴ·reshape(DD,2,2)·pol2/√braket — ⚠ NB
MATLAB's column-major reshape makes this the TRANSPOSED tensor, consistent
with the transposed gd matrices the polarizations are iterated on (a
conjugate-polarization convention runs through the whole upstream pipeline;
substituting the un-transposed Hermitian DD changes |η|² by ~15% at the C-Mod
reference saddle and is WRONG — the port matches the original to 7+ digits
once the convention is honored, with arg β defined mod π by eigenvector sign).
eta2=|eta|²; tau = exp(−π·eta2); beta = √(2π·tau)/(eta·cgamma(−i·eta2)).
If hyperbola malformed → tau=0 (upstream sets tau=0 and continues; beta left
stale ⚠ — port: mark conversion invalid, no split).
'Trs' (transmitted): direction zinzst = z*(1:4) − yalf0(1:4); guess −2.2;
fact = fzero(σ ↦ U(yalf0 + σ·zinzst), start 2.0); y_trs = yalf0 + fact·zinzst.
⚠ fzero from scalar start (MATLAB expands bracket automatically) — Julia:
implement scalar root find with automatic bracket expansion around start
(mirror fzero's search: geometric expansion ±, then Brent).
'Cnv' (converted): y_cnv = yalf0 (position/k unchanged — continues on the
eigenvalue-followed branch); τ, β recorded.
Amp extras (odeDim≥9): Salf/Slam matching (lines 987-998), transmitted lnE²
−= 2π·eta2, converted lnE² += log|β|², phase += arg β.

## trajectory.m events (port semantics)

Monitors recorded per accepted step (upstream: inside RHS with step heuristics;
port: after accepted steps — documented improvement). MON row = [σ, mon1, mon2]
from disp_eig('Mon'). Warmup: no events for first 15 recorded rows.
Event logic per check: last 6 rows; rescale times by first; quadratic LSQ fit
A=[1 t t²]; caustic: fitted VALUE of mon1 at current σ crosses 0 (col jcaust=1);
conversion: fitted DERIVATIVE of mon2 (col jconvt=2) crosses 0. Normalizations
adjcnv/adjctc = 1/first-value; event value = adj·mon (sign change = trigger),
terminal. convert_which (ray.m): refit last ≤5 rows of |mon2| (cf=A\|TrD|);
convert iff |moncnv|<1e-8-normalized event fired AND z2>0 (fit curvature: min)
AND sign(z1(first))≠sign(z1(last)). NOTE moncnv is the events-fit derivative —
in port: trigger conversion when the events fit says derivative crossed zero
with positive curvature (minimum of |trace| ↔ closest approach).

## cgamma.m (complex Γ, Stirling + recurrence)

If Im(z)<0: z→conj(z), conjugate result at end. While Re(z)<9: lnΓ −= log z,
z += 1 (loop count floor(9−x)). Stirling: lnΓ += (z−½)log z − z + ½log 2π +
1/(12z) − 1/(360z³) + 1/(1260 z⁵) − 1/(1680 z⁷) + 1/(1188 z⁹) −
691/(360360 z¹¹) + 1/(156 z¹³) − 3617/(122400 z¹⁵). Return exp(lnΓ).
Tests: Γ(1)=1, Γ(½)=√π, |Γ(iy)|²=π/(y sinh πy), |Γ(1+iy)|²=πy/sinh πy,
Γ(z+1)=zΓ(z).

## adjust_disp_m (launch on dispersion surface)

Given (r,z,kr,kz), poloidal mode m: kθ=m/ρ fixed; kρ0 = −√(kr²+kz²−(m/ρ)²);
solve kρ: U(r,z, kρcosθ−(m/ρ)sinθ, −(kρsinθ+(m/ρ)cosθ)) = 0 (disp_eig 'Dsp',
fzero from kρ0). Return adjusted [r,z,kr,kz].

## Maslov/caustic transform (ray.m caustic_list; Amp only)

GG=[dWr2 dWrz; dWrz dWz2]; invert; y(5:7) ← [inv11, inv12, inv22];
idx=π/4·Σsign(eig GG); k→x: f=1/(2π√|det|), sgn=+1; x→k: f=2π/√|det|, sgn=−1;
lnE² += log f²; phase += idx + sgn·(x·k); toggle inKspace; in k-space use short
intervals 0.05·timespan.

## Julia file layout (src/raycon/)

- `raycon_types.jl` — RayconConstants, TokamakEquilibrium (Solovev params+psin),
  PlasmaProfiles (species), RayconProblem (assembled, ω, kφ, MODEL), presets
  (raycon_cmod).
- `raycon_solovev.jl` — solovev_flux (levels 1/3/6/10 outputs), map_flux (Brent),
  flux_surface_mesh.
- `raycon_magnetic.jl` — magnetic_geometry(ρ,θ; level) incl. dmagnetic FD.
- `raycon_dispersion.jl` — stix_elements + derivs, tensor+derivs, disp_eig
  equivalents: raycon_U ('Dsp'), raycon_rhs ('Trj'), monitors, polarization,
  frequencies ('Frq'), dUdom ('Sgn'), msw.
- `raycon_conversion.jl` — det-U + 1st/2nd derivs (dispertok path), osculating
  Hessian, saddle iteration, normal form, η/τ/β, cgamma, transmitted/converted.
- `raycon_trace.jl` — DP45 adaptive with events + monitor recording.
- `raycon_driver.jl` — trace_rays (main-equivalent, conversion splitting),
  antenna launch (adjust + Msw/fzero variant), Maslov transform, deposition
  accumulation (calcFlux equivalents).
Tests: test/test_raycon.jl (FD oracles + MATLAB reference data
test/reference/raycon_reference.json generated by tools/raycon_reference.m).

## Deliberate port decisions (document in code + design doc)

1. No globals: explicit problem/config/result structs. No plotting/GUI/beep.
2. ODE: own Dormand–Prince 4(5) with rtol/atol per upstream; events at accepted
   steps (upstream approximates the same via step heuristics inside RHS).
3. mapFlux root tolerance 1e-12 (upstream 1e-4). fzero → bracket-expand + Brent.
4. cld3x3: full tracing support (disp_eig complete); conversion machinery
   REFUSES cld3x3 (upstream sgn_fix known bug, lines 476-477) with clear error.
5. Temperature-profile nb/tb quirk preserved; flagged in docstring.
6. dispertok's non-curvature spatial chain + 2nd-deriv chain omissions preserved
   in the conversion path (parity), NOT in the tracing path (disp_eig exact).
7. Same SI constants as upstream (not CODATA-2018) for 1:1 comparability.
10. UPSTREAM BUGS FIXED in the port (found by finite-difference oracles,
    RCN-005/006/007; kept out of the parity claims):
    a. `dD12dom` (and `dD13dom`, `dD23dom`) sign: elements are −Nn·Nb etc. and
       N ∝ 1/ω, so ∂(−NnNb)/∂ω = +2NnNb/ω; both dispertok.m and disp_eig.m
       write −2NnNb/ω (pattern-copied from the diagonal). Affects ∂U/∂ω
       ('Sgn' and any dω normalization) by ~10% at cmod parameters.
    b. disp_eig.m swaps the dkz² and dkr·dkz second-derivative expressions of
       the off-diagonal elements (dD12/dD13/dD23) — their "dD12dkz2" holds the
       mixed derivative and vice versa. Affects the 20-dim tangent-map
       (symplectic S-matrix) evolution only.
11. dmagnetic FD step 1e-5 instead of upstream 1e-8: the upstream step puts
   ~4·eps·|b|/h² ≈ 10-20% roundoff noise on dbdr2/dbdrz/dbdz2 (which feeds the
   conversion Hessian via dbds2 → Stix 2nd derivatives); 1e-5 reduces this to
   ~1e-5 relative. Reference comparisons of 2nd-derivative-dependent
   quantities therefore use tolerances that absorb the UPSTREAM noise.
8. 'Amp' machinery (focusing, lnE², phase, deposition, Maslov): ported where
   load-bearing for 'Con'; full Amp evolution not claimed — upstream ships
   with the Amp RHS call commented out.
12. UNIFIED UNITS INTERFACE (raycon_normalized.jl): every entry point has a
    `PlasmaUnits`-first method speaking the package's Ω_ci normalization
    (lengths d_i, k in 1/d_i, ω in Ω_ci, B in B0, n in n0, T in m_i·v_A²);
    conversion happens at the boundary and the SI core underneath is the
    MATLAB-pinned engine, so the rescaling is exact (RCN-013 pins normalized ==
    rescaled-SI to ~1e-12). σ is scale-invariant (U is dimensionless).
    Default constants for `PlasmaUnits`-built problems are the package's
    CODATA values (self-consistency); pass `cnst = RayconConstants()` for
    upstream-constant parity (~1e-5 shifts).
