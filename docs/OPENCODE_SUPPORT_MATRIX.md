# OpenCode support and evidence matrix

Version 2026-11-09 (P00 initial inventory). This matrix makes the required
model combinations explicit. It is a planning and evidence index, not a claim
that a row is accepted. Status vocabulary follows the handoff:

- **source** — canonical implementation exists and has had source review.
- **focused** — an existing focused harness exists; it may not have run on the
  current exact head.
- **exact-head** — minimum/latest Julia, focused/neighbors, Aqua, concrete JET,
  whole-tree JET, formatter and required workflows pass on the recorded SHA.
- **science** — independent oracle/refinement/failure/restart evidence passes.
- **resource** — predeclared cost/accuracy evidence passes.
- **hardware** — real ranks/devices/nodes pass; emulation or one-rank calls do
  not count.

`open` means at least one required gate is missing. A source row never closes
because a file exists or because an unrelated primitive test is green.

## Cross-cutting required variants

Every applicable family row must be exercised, where the model admits it, with:

| Dimension | Required variants | Owner packages |
|---|---|---|
| Species | one, many, empty, unequal charges/masses/weights, mobile/fixed/tracer, added/deleted/grown populations | P02, P29, P31, P35–P39 |
| Geometry | 1D/2D/3D; odd/even grids; periodic seam; thin/empty/unequal partitions; coarse/fine interface; axis; material wall | P05, P08–P11, P19–P24 |
| Shape/kinematics | NGP, CIC, TSC; nonrelativistic, relativistic and mixed species; original-factor round trips | P02, P12–P14, P35 |
| Closure | isothermal, adiabatic/polytropic, evolving scalar, transport (S1/S2/S3b), field-aligned conduction, ideal CGL, evolving tensor/inertia where promised | P17, P33, P35, P37, P40, P41 |
| Collisions/sources | ion–ion self/cross, electron–ion thermal exchange NGP/CIC/TSC, ionization/growth, antenna/Raycon/radiation, external reservoirs | P02–P03, P35–P39 |
| Boundaries | periodic, PEC/conducting, grounded/material, impact/absorbing, PML/conductive layer, open shock reservoirs | P04–P07, P20, P32–P34 |
| Calendar/controller | fixed step, subcycled, S1/S2/generic dividing, persistent arbitrary schedules, coupled CIC/TSC error control, pending events | P29–P31 |
| Backend | serial CPU, real pencil MPI, CUDA, Metal, multi-device; host-only extension load is not a device pass | P08–P15, P34, P39 |
| Precision | Float32, Float64, Float16/BigFloat only where the family documents them | P02, P35, P40 |
| Persistence | raw and archive, schema versions/migrations, fresh-process continuation, corruption/stale aliasing, full phase/RNG/mesh/source state | P05–P07, P09–P11, P13–P15, P29–P31, P35–P41 |
| Failure | invalid/nonfinite/range, aliases, budget exhaustion, late solve/IO failure, collective rejection, rollback before commit | all packages |

## Model-family rows

| ID | Family / canonical owner | Required geometry, shape, closure and species | Boundary / backend / persistence | Packages | Current evidence state |
|---|---|---|---|---|---|
| F01 | Periodic hybrid `HybridStepper` (`src/models/hybrid_pic.jl`, `src/integrators/`) | 1D3V/2D3V/3D3V; NGP/CIC/TSC; isothermal, adiabatic, evolving scalar, transport, field-aligned and CGL closures; arbitrary ion collections | periodic Fourier; serial CPU F32/F64; raw/archive checkpoints | P35, P40, P41 | source + focused; exact-head/science/resource open |
| F02 | Collisional hybrid `CollisionalHybridStepper` | scalar, evolving, transport and CGL electrons; ion–ion cell-local collisions; unequal species/weights | periodic serial; RNG/source ledger/whole-step rollback; collisional schema 1–4 | P02, P35–P37 | source + focused; extreme-occupancy/resource/order open |
| F03 | Electron–ion thermal exchange | NGP reduced and CIC/TSC overlap; scalar/collection consumers | serial CPU; policy 2/schema 3 continuation; zero-flux reduction | P35 | source + focused; native energy/momentum/entropy and mixed-precision gates open |
| F04 | Bounded longitudinal electrostatic `YeeLongitudinalPIC` | 1D3V grid, CIC cloud, event-free finite cloud; fixed/background plates | PEC/closed longitudinal; serial F32/F64; schema-2 force policy | P04–P05 | source; independent cavity/trajectory/work and event scope open |
| F05 | Bounded 1D3V EM `YeeElectromagneticPIC` | both transverse polarizations; Boris/Yee synchronized stages; mobile/fixed collections | PEC/periodic mask; serial F32/F64; schema 2 | P04–P07 | source; general nonlinear, material/event and broader bounded EM open |
| F06 | Impact/capture longitudinal `YeeImpactLongitudinalPIC` | finite-cloud reflection, repeated crossings/capture | closed longitudinal; serial CPU; pending-event restart | P04–P05 | source; D2/3 crossed-cell/face/edge/corner and material laws open |
| F07 | Absorbing longitudinal `YeeAbsorbingLongitudinalPIC` | finite cloud; absorbing layer | layer state/restart | P05–P06 | source; reflection/angle/frequency/thickness and coupled plasma exit open |
| F08 | PML / conductive layer (`yee_pml.jl`, `yee_conductive_layer.jl`) | vacuum conductivity/ADE PML; layer profiles | face/edge/corner; layer memory/restart | P06 | source; CPML generality, corner overlap, late-time and plasma coupling open |
| F09 | Electrostatic PIC `ElectrostaticPIC` | 1D–3D periodic Fourier, CIC/NGP/TSC, signed/multiple species | serial CPU F32/F64; ES checkpoint | P02, P29 | source + focused; broad species/schedule intersections open |
| F10 | EMPIC collections `EMPIC`, `EMPIC1D` | charge-conserving summed currents; relativistic/fixed/mobile; specialized S1/S2 and generic D1–3; persistent schedules | serial CPU; schema 1–3/policy 1–2; pending deadlines | P29–P31 | source + focused; controller adoption, variable-macro and complete failure/restart open |
| F11 | Full-PIC `full_pic.jl` | full electromagnetic particle/field collection paths | serial CPU; checkpoint | P29–P31 | source; parser/JET/full headed acceptance in progress |
| F12 | Implicit hybrid `ImplicitHybridStepper` | periodic D1–3; isothermal and evolving pressure; CIC; arbitrary mobile species | serial F32/F64; schema 1/2 | P16–P17 | source + focused; true residual/phase/order/science/resource open |
| F13 | Implicit spectral EM `ImplicitEMPICStepper` | periodic D1–3; exact-neutral all-mobile CIC collections | serial F32/F64; schema restart | P16, P18 | source + focused; coupled equation/continuity/order/restart open |
| F14 | Implicit Yee EM `ImplicitYeeEMPICStepper` | compatible periodic Yee maps, transverse paths | serial; schema restart | P16, P18 | source; D2/3, wider orbits, boundaries and backend open |
| F15 | Sheared periodic hybrid `ShearedPeriodicHybridPIC` | mapped logical positions/Cartesian velocities; odd grids; coupled RK4 ions/fields | fixed nonorthogonal map; serial; checkpoint | P19 | source + focused; weighted adjoint, physical boundary and resource open |
| F16 | Toroidal annulus `ToroidalAnnulusHybridPIC` | logical/Cartesian ions, metric flux, adiabatic pressure, electrical/thermal reservoirs | fixed contacts; serial; full checkpoint | P19–P20 | source + focused; solid axis/material/equilibria/imports/backends open |
| F17 | Slab gyrokinetic `SlabGyrokineticPIC` | full-f ring-TSC, modified adiabatic electrons, uniform B | serial; weight/equilibrium/RNG checkpoint | P21–P22 | source + focused; collective RK mode/FLR/free-energy/science open |
| F18 | Gyrokinetic utilities/geometry (`gyrokinetics.jl`) | guiding centres, drifts, gyroaverage | serial | P22 | source; toroidal/mirror/bounce, collisions, kinetic electrons/EM open |
| F19 | Static composite ES `StaticCompositeElectrostaticPIC` | 1D hierarchy, triangular cloud, physical walls/Neumann, signed ions | serial F32/F64; schema 1 | P23–P24 | source + focused; independent quadrature/self-force/convergence open |
| F20 | Dynamic composite ES `DynamicCompositeElectrostaticPIC` | completed-boundary 1D hierarchy, regrid/adapt, mesh energy | serial; full epoch/history restart | P23–P24 | source + focused; arbitrary topology/refinement/hybrid-EM open |
| F21 | Adaptive composite ES `AdaptiveCompositeElectrostaticPIC` | Cartesian AMR, grounded/periodic multidimensional, composite Poisson | serial; hierarchy/solve history | P23–P24 | source + focused; material/electrode/Neumann/mixed-axis and subcycling open |
| F22 | Hybrid/EM AMR | compatible face/edge currents, flux/field synchronization, local time | partition/device hierarchy | P24 | not implemented; explicitly unsupported rather than silently routed |
| F23 | Open Hall hybrid `OpenHallHybrid` | D2/3 CIC isothermal, mixed/SBP/SAT operators, upstream/downstream reservoirs | numerical reservoir boundaries; serial; full histories/RNG | P32–P34 | source + focused; physical shock regime, pressure/transport/inertia and rank/device open |
| F24 | Shock/initial-condition families (`shock_cn.jl`, `PerpShock`, `LeroyBoundary`) | 1D/2D perpendicular/parallel shock drivers; CN and explicit paths | fixed/periodic boundaries; serial | P16, P25–P27 | source + validation histories; exact-head and independent datasets open |
| F25 | Hall-MHD/CAM-CL (`hall_mhd.jl`, `camcl.jl`) | periodic fluid hybrid, RK4/CAM-CL | serial; checkpoint | P16, P35 | source + focused; broader closure/source intersections open |
| F26 | Antenna/Raycon/radiation sources | standalone antenna field, geometric-optics rays, reduced radiation reaction | serial; no accounted plasma reservoir | P38–P39 | source components; accounted power/impulse/coupling/combined restart open |
| F27 | Reaction/ionization/growth sources | ionization, neutral reservoirs, newborn phase/current, monotonic IDs | serial/rank/device RNG/source clocks | P36–P39 | primitives/source partially present; complete composition open |
| F28 | Distributed pencil ES/EM/hybrid/Maxwell | periodic D3 pencil, charged multiple species, nonuniform modes | real 2/4/8 ranks; partition-aware checkpoint | P08–P11 | primitives focused + some 2-rank model controls; production families/scaling open |
| F29 | Distributed collisional ES/EM | migrated populations, canonical cell ownership, full per-rank RNG/state | real ranks; restart | P03, P10 | source; actual-rank statistical/full composition open |
| F30 | Pencil checkpoint/repartition | changed layout/rank, IDs/weights/RNG redistribution, multi-node | real ranks/nodes | P11 | partial source; rebalance and scaling open |
| F31 | CUDA primitives/resident hybrid | D1 CIC F64 isothermal lifecycle; device primitives | actual NVIDIA GPU; host-only load is not acceptance | P12–P15 | source/host admission; actual kernel/device/science/resource open |
| F32 | Metal primitives/resident solver | primitive Boris/gather/deposit paths | actual Apple GPU; resident lifecycle required | P12–P15 | primitives source; genuine resident solver and device acceptance open |
| F33 | Multi-GPU/multi-device | domain ownership + migration + device arrays | actual multiple devices/nodes | P15, P39 | not implemented; explicitly open |
| F34 | Validation references (`validation/`) | pinned NHDS/Hybrid-VPIC provenance, production dispersion, oblique Hall/kinetic/thermal/EM, nonlinear shocks/reconnection/AMR/GK | named hardware and raw data | P25–P28 | reference infrastructure and some campaigns source/pass; broad independent and comparative campaigns open |

## Explicit non-claims

- No row above is `complete`; the 23 legacy checked milestones remain the only
  historically checked rows in the 79-row tracker.
- Host-only CUDA extension load, `COMM_SELF`, a one-rank run, a serial replay or
  a CPU conversion does not satisfy any hardware row.
- The broad Cartesian product is not promised. A combination is required only
  when the family guide admits it; unsupported combinations must be rejected
  before state mutation and listed here as unsupported, not as validated.
- Missing promised physics is not an applicability exclusion.
- This matrix was built before P01 exact-head execution; P01 and later packages
  must replace `source`/`focused` entries with dated raw evidence, never with a
  file-presence statement.
