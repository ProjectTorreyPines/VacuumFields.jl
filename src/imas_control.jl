function boundary_control_points(dd::IMAS.dd{T};
    x_points_weight::Real=1.0,
    strike_points_weight::Real=1.0,
    active_x_points::AbstractVector{Int}=Int[]) where {T<:Real}

    eqt = dd.equilibrium.time_slice[]

    # Boundary control points
    iso_cps = IsoControlPoints(eqt.boundary.outline.r, eqt.boundary.outline.z)

    strike_weight = strike_points_weight / length(eqt.boundary.strike_point)
    strike_cps =
        IsoControlPoint{T}[IsoControlPoint{T}(strike_point.r, strike_point.z, iso_cps[1].R2, iso_cps[1].Z2, 0.0, strike_weight) for strike_point in eqt.boundary.strike_point]
    append!(iso_cps, strike_cps)

    # Saddle control points
    saddle_weight = x_points_weight / length(eqt.boundary.x_point)
    saddle_cps = SaddleControlPoint{T}[SaddleControlPoint{T}(x_point.r, x_point.z, saddle_weight) for x_point in eqt.boundary.x_point]
    if !isempty(active_x_points)
        xiso_cps = Vector{IsoControlPoint{T}}(undef, length(active_x_points))
        for (k, ax) in enumerate(active_x_points)
            xiso_cps[k] = IsoControlPoint{T}(eqt.boundary.x_point[ax].r, eqt.boundary.x_point[ax].z, iso_cps[1].R2, iso_cps[1].Z2, 0.0, saddle_weight)
        end
        append!(iso_cps, xiso_cps)
    end
    return iso_cps, saddle_cps
end

function magnetic_control_points(dd::IMAS.dd{T};
    reference_flux_loop_index::Int=1,
    flux_loop_weights::AbstractVector{<:Real}=T[],
    magnetic_probe_weights::AbstractVector{<:Real}=T[]) where {T<:Real}

    mags = dd.magnetics
    pavgs = IMAS.getavg(Val(:b_field_pol_probe), dd, dd.global_time, 0.005)
    pdata = pavgs.data

    # Probe control points (Note: type and size ignored)
    min_σ = minimum(probe.err for probe in pdata)
    field_cps = VacuumFields.FieldControlPoint{T}[]
    for (k, probe) in enumerate(mags.b_field_pol_probe)
        !isempty(probe.field.validity) && probe.field.validity < 0 && continue
        weight = isempty(magnetic_probe_weights) ? 1.0 : magnetic_probe_weights[k]
        if !isempty(probe.field.data_σ)
            pdata[k].err < eps() && continue
            weight /= pdata[k].err / min_σ
        end
        #IMAS allows probes to mix toroidal and poloidal fields so that needs to be accounted for
        bpol_field = isempty(probe.toroidal_angle) ? pdata[k].val : pdata[k].val * (1 - cos(probe.poloidal_angle) * sin(probe.toroidal_angle))

        push!(field_cps, VacuumFields.FieldControlPoint{T}(probe.position.r, probe.position.z, probe.poloidal_angle, bpol_field, weight))
    end

    # Flux loop control points (Note: type and size ignored)

    iref = reference_flux_loop_index
    loops = mags.flux_loop
    flavgs = IMAS.getavg(Val(:flux_loop), dd, dd.global_time, 0.005)
    fldata = flavgs.data

    flux_cps = VacuumFields.FluxControlPoint{T}[]
    loop_cps = VacuumFields.IsoControlPoint{T}[]
    if iref > 0
        # Fit a single reference flux value and the difference between the other flux_loops and this value (this matches how data is collected on DIII-D and fit with EFIT)
        N = length(loops)
        @assert iref <= N
        # Convert from total flux uncertainty (IMAS convention) to differential
        relative_σ(k) = sqrt(fldata[k].err^2 - fldata[iref].err^2)
        min_σ = 1 / eps()
        for k in (iref+1):(iref+N-1) # exclude i_ref explicitly
            ck = IMAS.index_circular(N, k)
            if !isempty(loops[ck].flux.data_σ)
                fldata[ck].err < eps() && continue
                min_σ = relative_σ(ck) < min_σ ? relative_σ(ck) : min_σ
            end
        end
        for k in (iref+1):(iref+N-1) # exclude i_ref explicitly
            ck = IMAS.index_circular(N, k)

            !isempty(loops[ck].flux.validity) && loops[ck].flux.validity < 0 && continue
            weight = isempty(flux_loop_weights) ? 1.0 : flux_loop_weights[ck]
            if !isempty(loops[ck].flux.data_σ)
                fldata[ck].err < eps() && continue
                weight /= relative_σ(ck) / min_σ
            end

            push!(
                loop_cps,
                VacuumFields.IsoControlPoint{T}(
                    loops[ck].position[1].r,
                    loops[ck].position[1].z,
                    loops[iref].position[1].r,
                    loops[iref].position[1].z,
                    fldata[ck].val - fldata[iref].val,
                    weight
                )
            )
        end
        weight = isempty(flux_loop_weights) ? 1.0 : flux_loop_weights[iref]
        if !isempty(loops[iref].flux.validity)
            @assert loops[iref].flux.validity >= 0
        end

        push!(flux_cps, VacuumFields.FluxControlPoint{T}(loops[iref].position[1].r, loops[iref].position[1].z, fldata[iref].val, weight))
    else
        # Fit each flux loop separately (this is how data is collected on NSTX and fit with EFIT)
        min_σ = minimum(loop.err for loop in fldata)
        for (k, loop) in enumerate(loops)
            !isempty(loop.flux.validity) && loop.flux.validity < 0 && continue
            weight = isempty(flux_loop_weights) ? 1.0 : flux_loop_weights[k]
            if !isempty(loop.flux.data_σ)
                fldata[k].err < eps() && continue
                weight /= fldata[k].err / min_σ
            end

            push!(flux_cps, VacuumFields.FluxControlPoint{T}(loop.position[1].r, loop.position[1].z, fldata[k].val, weight))
        end
    end
    return flux_cps, loop_cps, field_cps
end