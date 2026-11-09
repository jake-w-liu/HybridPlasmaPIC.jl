# Full-time phases and exact represented receipts remain distinct from the
# canonical rounded physical moment views returned by the common diagnostic.
_empic_collection_velocity_offsets(schedule::_EMPICCoupledSchedule, mobile, offset) =
    fill(zero(_EMPICCoupledQ), length(mobile))
function _empic_collection_schedule_diagnostics(schedule::_EMPICCoupledSchedule, result)
    phase = schedule.phase
    report = phase.report
    measured = phase.field_defect + phase.particle_defect + phase.coupling_defect
    return deepcopy(
        merge(
            result,
            (
                schedule_policy = 3,
                phase_policy = :full_time_coupled,
                magnetic_time_offset = zero(_EMPICCoupledQ),
                source_time = phase.time,
                actual_elapsed = phase.elapsed,
                public_clock_defect = _exact_rational(schedule.shadow.clock[1]) - phase.time,
                public_minus_actual = _exact_rational(schedule.shadow.clock[1]) - phase.elapsed,
                actual_minus_source = phase.elapsed - phase.time,
                nominal_H = phase.interval,
                last_deadline = phase.last,
                next_deadline = phase.next,
                species_dt = schedule.proposal.cap,
                species_adaptive = schedule.proposal.adaptive,
                species_dt_min = schedule.proposal.minimum,
                schedule_limits = schedule.proposal.limits,
                event_budget = schedule.proposal.budget,
                error_control = schedule.control,
                admission_control = phase.admission_control,
                controller_next_h = phase.next_h,
                spectral_signature = schedule.signature,
                initial_reported_energy = phase.initial_energy,
                field_defect = phase.field_defect,
                particle_defect = phase.particle_defect,
                coupling_defect = phase.coupling_defect,
                physical_work_defect = measured,
                cumulative_costs = phase.costs,
                last_error_control_report = report,
            ),
        ),
    )
end
function _empic_coupled_diagnostics(es, species)
    schedule = es.collection.schedule
    _ec_acquire(schedule)
    try
        meter = _ec_meter(schedule.phase.admission_control)
        charges = _ec_ready(es, species, meter)
        result = _empic_collection_moment_diagnostics(es, species, charges, es.collection)
        K, U = _ec_energy(es, species)
        defect =
            K + U - schedule.phase.initial_energy - schedule.phase.field_defect -
            schedule.phase.particle_defect - schedule.phase.coupling_defect
        return merge(
            result,
            (
                reported_kinetic_energy = K,
                represented_field_energy = U,
                reported_total_energy = K + U,
                accounting_defect = defect,
                observation_costs = _ec_costs(meter),
            ),
        )
    finally
        _ec_release(schedule)
    end
end
