"""
    toroidal_annulus_diagnostics(state, species=state.owners)

Detached physical reservoirs, source/flux/metric observations and actual signed
contact/pressure-controller ledgers. Source review is not boundary stability,
axis convergence, finite-step conservation or resource acceptance.
"""
function toroidal_annulus_diagnostics(
    st::ToroidalAnnulusHybridPIC{T},
    species = st.owners,
) where {T}
    _tah_acquire(st)
    try
        observed = _tah_ready(st, species)
        Q, mesh, f = _exact_rational, st.mesh, st.fields
        reservoirs = _tah_reservoirs(st, species, st.beta, f.pe)
        number = [_exact_float_sum(p.weight) for p in species]
        charge = [Q(p.q) * n for (p, n) in zip(species, number)]
        momentum = ntuple(
            c -> sum(
                (
                    Q(p.m) * Q(p.weight[i]) * Q(p.v[c][i]) for p in species for
                    i in eachindex(p.weight)
                );
                init = zero(_CompositeQ),
            ),
            3,
        )
        axial_momentum = zero(_CompositeQ)
        for p in species, marker in eachindex(p.weight)
            point = _tah_point_geometry(mesh, ntuple(c -> p.x[c][marker], 3))
            su, cu, sv, cv = map(Q, point.trig)
            position = _toroidal_cartesian_point(Q(mesh.R0), Q(point.r), su, cu, sv, cv)
            axial_momentum +=
                Q(p.m) *
                Q(p.weight[marker]) *
                (position[1] * Q(p.v[2][marker]) - position[2] * Q(p.v[1][marker]))
        end
        volume = grid_charge = snapshot_magnetic = zero(_CompositeQ)
        volumes = Array{_CompositeQ}(undef, mesh.n)
        for I in CartesianIndices(f.n)
            V = _tah_volume(mesh, I)
            volumes[I] = V
            volume += V
            grid_charge += V * Q(f.n[I])
            snapshot_magnetic += V * sum(abs2, _tah_row(f.B, I)) / 2
        end
        radial_flux = fill(zero(_CompositeQ), mesh.n[1])
        poloidal_flux = fill(zero(_CompositeQ), mesh.n[2])
        toroidal_flux = fill(zero(_CompositeQ), mesh.n[3])
        for I in CartesianIndices(f.n)
            radial_flux[I[1]] += Q(st.beta[1][I]) * Q(mesh.dx[2]) * Q(mesh.dx[3])
            poloidal_flux[I[2]] += Q(st.beta[2][I]) * Q(mesh.op.sbp.H[I[1]]) * Q(mesh.dx[3])
            toroidal_flux[I[3]] += Q(st.beta[3][I]) * Q(mesh.op.sbp.H[I[1]]) * Q(mesh.dx[2])
        end
        ion_trace, current_trace, electron_trace =
            ntuple(_ -> [zero(_CompositeQ), zero(_CompositeQ)], 3)
        for (face, side) in enumerate((1, mesh.n[1])), index in CartesianIndices(Base.tail(mesh.n))
            I = CartesianIndex((side, Tuple(index)...))
            metric, _ = _tah_metric(mesh, I)
            normal = face == 1 ? -1 : 1
            area = metric[2] * metric[3] * Q(mesh.dx[2]) * Q(mesh.dx[3])
            radial = _tah_basis(mesh, I)[1]
            ion_trace[face] += normal * area * _tah_dot(_tah_row(f.Ji, I), radial)
            current_trace[face] += normal * area * _tah_dot(_tah_row(f.J, I), radial)
            electron_trace[face] += normal * area * Q(f.n[I]) * _tah_dot(_tah_row(f.ue, I), radial)
        end
        fields = NamedTuple{keys(f)}(map(values(f)) do value
            value isa Tuple ? map(copy, value) : copy(value)
        end)
        geometry = (
            n = mesh.n,
            R0 = mesh.R0,
            r_inner = mesh.r_inner,
            r_outer = mesh.r_outer,
            tau = mesh.tau,
            orientation = -1,
            coordinates = map(copy, mesh.op.coordinates),
            r = copy(mesh.r),
            R = copy(mesh.R),
            volume = volumes,
        )
        nominal = Q(mesh.tau)^2 * Q(mesh.R0) * (Q(mesh.r_outer)^2 - Q(mesh.r_inner)^2) / 2
        energy = sum(reservoirs)
        return (
            time = st.time,
            step = st.step,
            last_dt = st.last_dt,
            phase = :synchronous_integer,
            velocity_time_offset = zero(T),
            number = number,
            charge = charge,
            momentum = momentum,
            mechanical_axial_angular_momentum = axial_momentum,
            energy = (
                kinetic = reservoirs[1],
                magnetic = reservoirs[2],
                electron = reservoirs[3],
                total = energy,
            ),
            snapshot_magnetic_energy = snapshot_magnetic,
            snapshot_metric_energy_defect = snapshot_magnetic - reservoirs[2],
            initial_energy = deepcopy(st.initial_energy),
            cumulative_outward_work = deepcopy(st.outward_work),
            cumulative_controller_work = deepcopy(st.controller_work),
            cumulative_stage_defect = deepcopy(st.stage_defect),
            cumulative_rk_defect = deepcopy(st.rk_defect),
            accounting_defect = energy - st.initial_energy + st.outward_work - st.controller_work -
                                st.stage_defect - st.rk_defect,
            physical_volume = volume,
            represented_nominal_volume = nominal,
            quadrature_volume_defect = volume - nominal,
            deposited_charge = grid_charge,
            deposited_charge_defect = grid_charge - sum(charge),
            minimum_density = minimum(f.n),
            geometry = geometry,
            fields = fields,
            beta = map(copy, st.beta),
            cross_section_flux = (
                radial = radial_flux,
                poloidal = poloidal_flux,
                toroidal = toroidal_flux,
            ),
            grid_boundary_trace = (
                ion_charge = Tuple(ion_trace),
                total_current = Tuple(current_trace),
                electron_number = Tuple(electron_trace),
                current_closure_defect = ntuple(
                    i -> ion_trace[i] - current_trace[i] - electron_trace[i],
                    2,
                ),
            ),
            pressure_reservoir = st.pressure_reservoir,
            normal_flux_data = deepcopy(st._normal_flux),
            last_report = deepcopy(st.last_report),
            last_endpoint = deepcopy(st.last_endpoint),
            observed_endpoint = observed,
        )
    finally
        _tah_release(st)
    end
end
