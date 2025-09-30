function boundary_control_points(
    eqt::IMAS.equilibrium__time_slice{T};
    x_points_weight::Real,
    strike_points_weight::Real,
    active_x_points::AbstractVector{Int}=Int[]) where {T<:Real}

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
    return (iso_cps=iso_cps, saddle_cps=saddle_cps)
end

function equilibrium_control_points(
    eqt::IMAS.equilibrium__time_slice{T},
    pc::IMAS.pulse_schedule__position_control{T};
    x_points_weight::Float64,
    strike_points_weight::Float64,
    magnetic_axis_weight::Float64=0.0
) where {T<:Real}

    # boundary
    if ismissing(eqt.global_quantities, :ip) # field nulls
        iso_cps = FluxControlPoints(eqt.boundary.outline.r, eqt.boundary.outline.z, eqt.global_quantities.psi_boundary)
    elseif !isempty(pc.boundary_outline)
        # we favor taking the boundary from the pulse schedule, if available
        iso_cps = IsoControlPoints(IMAS.boundary(pc)...)
    else
        # solutions with plasma
        iso_cps = IsoControlPoints(eqt.boundary.outline.r, eqt.boundary.outline.z)
    end

    # x points
    saddle_cps = SaddleControlPoint{T}[]
    if x_points_weight == 0.0
        # pass
    elseif !isempty(pc.x_point)
        # we favor taking the x-points from the pulse schedule, if available
        saddle_weight = x_points_weight / length(pc.x_point)
        for x_point in pc.x_point
            r = IMAS.get_time_array(x_point.r, :reference, :constant)
            if r == 0.0 || isnan(r)
                continue
            end
            z = IMAS.get_time_array(x_point.z, :reference, :constant)
            push!(saddle_cps, SaddleControlPoint{T}(r, z, saddle_weight))
        end
    else
        saddle_weight = x_points_weight / length(eqt.boundary.x_point)
        for x_point in eqt.boundary.x_point
            push!(saddle_cps, SaddleControlPoint{T}(x_point.r, x_point.z, saddle_weight))
        end
    end

    # strike points
    flux_cps = FluxControlPoint{T}[]
    if strike_points_weight == 0.0
        # pass
    elseif !isempty(pc.strike_point)
        # we favor taking the strike points from the pulse schedule, if available
        flux_cps = FluxControlPoint{T}[]
        strike_weight = strike_points_weight / length(pc.strike_point)
        for strike_point in pc.strike_point
            r = IMAS.get_time_array(strike_point.r, :reference, :constant)
            if r == 0.0 || isnan(r)
                continue
            end
            z = IMAS.get_time_array(strike_point.z, :reference, :constant)
            push!(flux_cps, FluxControlPoint{T}(r, z, eqt.global_quantities.psi_boundary, strike_weight))
        end
    else
        strike_weight = strike_points_weight / length(eqt.boundary.strike_point)
        for strike_point in eqt.boundary.strike_point
            push!(flux_cps, FluxControlPoint{T}(strike_point.r, strike_point.z, eqt.global_quantities.psi_boundary, strike_weight))
        end
    end

    # magnetic axis
    if magnetic_axis_weight > 0.0
        push!(saddle_cps, SaddleControlPoint{T}(eqt.global_quantities.magnetic_axis.r, eqt.global_quantities.magnetic_axis.z, magnetic_axis_weight))
        push!(
            flux_cps,
            FluxControlPoint{T}(
                eqt.global_quantities.magnetic_axis.r,
                eqt.global_quantities.magnetic_axis.z,
                eqt.global_quantities.psi_axis,
                magnetic_axis_weight
            )
        )
    end

    return (iso_cps=iso_cps, flux_cps=flux_cps, saddle_cps=saddle_cps)
end

function magnetic_control_points(
    mags::IMAS.magnetics{T};
    reference_flux_loop_index::Int=1,
    flux_loop_weight::T=1.0,
    flux_loop_weights::AbstractVector{T}=T[],
    magnetic_probe_weight::T=1.0,
    magnetic_probe_weights::AbstractVector{T}=T[]) where {T<:Real}

    field_cps = FieldControlPoint{T}[]
    flux_cps = FluxControlPoint{T}[]
    iso_cps = IsoControlPoint{T}[]

    # Probe control points (Note: type and size ignored)
    if !isempty(mags.b_field_pol_probe)
        min_σ = minimum(IMAS.@ddtime(probe.field.data_σ) for probe in mags.b_field_pol_probe)
        for (k, probe) in enumerate(mags.b_field_pol_probe)
            if (!isempty(probe.field.validity) && probe.field.validity < 0) || !isempty(probe.field.data_σ) && IMAS.@ddtime(probe.field.data_σ) == 0.0
                continue
            end
            weight = (isempty(magnetic_probe_weights) ? 1.0 : magnetic_probe_weights[k]) * magnetic_probe_weight / length(mags.b_field_pol_probe) * min_σ / IMAS.@ddtime(probe.field.data_σ)
            if weight > 0.0
                #IMAS allows probes to mix toroidal and poloidal fields so that needs to be accounted for
                bpol_field = isempty(probe.toroidal_angle) ? IMAS.@ddtime(probe.field.data) : IMAS.@ddtime(probe.field.data) * (1 - cos(probe.poloidal_angle) * sin(probe.toroidal_angle))
                push!(field_cps, FieldControlPoint{T}(probe.position.r, probe.position.z, probe.poloidal_angle, bpol_field, weight))
            end
        end
    end

    # Flux loop control points (Note: type and size ignored)
    if !isempty(mags.flux_loop)
        iref = reference_flux_loop_index
        loops = mags.flux_loop
        @assert 0 <= iref <= length(loops)
        if iref > 0
            @assert isempty(loops[iref].flux.validity) || loops[iref].flux.validity >= 0
            # Fit a single reference flux value and the difference between the other flux_loops and this value (this matches how data is collected on DIII-D and fit with EFIT)

            # Convert from total flux uncertainty (IMAS convention) to differential
            relative_σ(k) = sqrt(IMAS.@ddtime(loops[k].flux.data_σ)^2 + IMAS.@ddtime(loops[iref].flux.data_σ)^2)
            min_σ = minimum(relative_σ(ck) for ck in 1:length(loops) if ck != iref)

            for ck in 1:length(loops)
                if (ck == iref) || (!isempty(loops[ck].flux.validity) && loops[ck].flux.validity < 0) || (!isempty(loops[ck].flux.data_σ) && IMAS.@ddtime(loops[ck].flux.data_σ) == 0.0)
                    continue
                end
                weight = (isempty(flux_loop_weights) ? 1.0 : flux_loop_weights[ck]) * flux_loop_weight / length(loops) * min_σ / relative_σ(ck)

                if weight > 0.0
                    push!(
                        iso_cps,
                        IsoControlPoint{T}(
                            loops[ck].position[1].r,
                            loops[ck].position[1].z,
                            loops[iref].position[1].r,
                            loops[iref].position[1].z,
                            IMAS.@ddtime(loops[ck].flux.data) - IMAS.@ddtime(loops[iref].flux.data),
                            weight
                        )
                    )
                end
            end

            weight = (isempty(flux_loop_weights) ? 1.0 : flux_loop_weights[iref]) * flux_loop_weight / length(loops)
            if weight > 0.0
                push!(flux_cps, FluxControlPoint{T}(loops[iref].position[1].r, loops[iref].position[1].z, IMAS.@ddtime(loops[iref].flux.data), weight))
            end

        else
            # Fit each flux loop separately (this is how data is collected on NSTX and fit with EFIT)
            min_σ = minimum(IMAS.@ddtime(loop.flux.data_σ) for loop in loops)
            for (k, loop) in enumerate(loops)
                if (!isempty(loop.flux.validity) && loop.flux.validity < 0) || (!isempty(loop.flux.data_σ) && IMAS.@ddtime(loop.flux.data_σ) == 0.0)
                    continue
                end
                weight = (isempty(flux_loop_weights) ? 1.0 : flux_loop_weights[k]) * flux_loop_weight / length(loops) * min_σ / IMAS.@ddtime(loop.flux.data_σ)
                if weight > 0.0
                    push!(flux_cps, FluxControlPoint{T}(loop.position[1].r, loop.position[1].z, IMAS.@ddtime(loop.flux.data), weight))
                end
            end
        end
    end

    return (iso_cps=iso_cps, flux_cps=flux_cps, field_cps=field_cps)
end