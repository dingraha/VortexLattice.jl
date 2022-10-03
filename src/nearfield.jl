"""
    near_field_forces!(properties, surfaces, wakes, reference, freestream, Γ;
        dΓdt, additional_velocity, Vh, Vv, symmetric, nwake, surface_id,
        wake_finite_core, wake_shedding_locations, trailing_vortices, xhat)

Calculate local panel forces in the body frame.
"""
@inline function near_field_forces!(props, surfaces, wakes, ref, fs, Γ;
    dΓdt, additional_velocity, Vh, Vv, symmetric, nwake, surface_id,
    wake_finite_core, wake_shedding_locations, trailing_vortices, xhat, prandtl_glauert)

    # number of surfaces
    nsurf = length(surfaces)

    # loop through receiving surfaces
    iΓ = 0 # index for accessing Γ
    for isurf = 1:nsurf

        receiving = surfaces[isurf]
        nr = length(receiving)
        nr1, nr2 = size(receiving)
        cr = CartesianIndices(receiving)

        # loop through receiving panels
        for i = 1:length(receiving)

            # get panel cartesian index
            I = cr[i]

            # --- Calculate forces on the panel bound vortex --- #

            # bound vortex location
            rc = top_center(receiving[i])

            # freestream velocity
            Vi = freestream_velocity(fs)

            # rotational velocity
            Vi += rotational_velocity(rc, fs, ref)

            # additional velocity field
            if !isnothing(additional_velocity)
                Vi += additional_velocity(rc)
            end

            # velocity due to surface motion
            if !isnothing(Vh)
                Vi += Vh[isurf][i]
            end

            # induced velocity from surfaces and wakes
            jΓ = 0 # index for accessing Γ
            for jsurf = 1:nsurf

                # number of panels on sending surface
                sending = surfaces[jsurf]
                Ns = length(sending)

                # see if wake panels are being used
                wake_panels = nwake[jsurf] > 0

                # check if we need to shift shedding locations
                if isnothing(wake_shedding_locations)
                    shedding_locations = nothing
                else
                    shedding_locations = wake_shedding_locations[jsurf]
                end

                # extract circulation values corresonding to the sending surface
                vΓ = view(Γ, jΓ+1:jΓ+Ns)

                # induced velocity from this surface
                if isurf == jsurf
                    # induced velocity on self
                    Vi += induced_velocity(I, surfaces[jsurf], vΓ;
                        finite_core = surface_id[isurf] != surface_id[jsurf],
                        wake_shedding_locations = shedding_locations,
                        symmetric = symmetric[jsurf],
                        trailing_vortices = trailing_vortices[jsurf] && !wake_panels,
                        xhat = xhat)
                else
                    # induced velocity on another surface
                    Vi += induced_velocity(rc, surfaces[jsurf], vΓ;
                        finite_core = surface_id[isurf] != surface_id[jsurf],
                        wake_shedding_locations = shedding_locations,
                        symmetric = symmetric[jsurf],
                        trailing_vortices = trailing_vortices[jsurf] && !wake_panels,
                        xhat = xhat)
                end

                # induced velocity from corresponding wake
                if wake_panels
                    Vi += induced_velocity(rc, wakes[jsurf];
                        finite_core = wake_finite_core[jsurf] || (surface_id[isurf] != surface_id[jsurf]),
                        symmetric = symmetric[jsurf],
                        nc = nwake[jsurf],
                        trailing_vortices = trailing_vortices[jsurf],
                        xhat = xhat)
                end

                jΓ += Ns # increment Γ index for sending panels
            end

            # steady part of Kutta-Joukowski theorem
            Γi = I[1] == 1 ? Γ[iΓ+i] : Γ[iΓ+i] - Γ[iΓ+i-1] # net circulation
            Δs = top_vector(receiving[I]) # bound vortex vector
            tmp = cross(Vi, Δs)
            Fbi = RHO*Γi*tmp

            if !isnothing(dΓdt)
                # unsteady part of Kutta-Joukowski theorem

                #TODO: decide whether to divide by perpindicular velocity like
                # Drela does in ASWING?

                dΓdti = I[1] == 1 ? dΓdt[iΓ+i] : (dΓdt[iΓ+i] + dΓdt[iΓ+i-1])/2
                c = receiving[I].chord
                Fbi += RHO*dΓdti*c*tmp/norm(Vi)

            end

            if prandtl_glauert
                # Prandtl-Glauert compressibility correction:
                mach = norm(Vi)/fs.speedofsound
                PG = 1/sqrt(1 - mach^2)
                Fbi *= PG
            end

            # --- Calculate forces on the left bound vortex --- #

            # bound vortex location
            rc = left_center(receiving[I])

            # freestream velocity
            Veff = freestream_velocity(fs)

            # rotational velocity
            Veff += rotational_velocity(rc, fs, ref)

            # additional velocity field
            if !isnothing(additional_velocity)
                Veff += additional_velocity(rc)
            end

            # velocity due to surface motion
            if !isnothing(Vv)
                Veff += Vv[isurf][I[1], I[2]]
            end

            # NOTE: We don't include induced velocity in the effective velocity
            # for the vertical segments because its influence is likely negligible
            # once we take the cross product with the bound vortex vector. This
            # is also assumed in AVL. This could change in the future.

            # steady part of Kutta-Joukowski theorem
            Γli = Γ[iΓ+i]
            Δs = left_vector(receiving[I])
            Fbli = RHO*Γli*cross(Veff, Δs)

            if prandtl_glauert
                # Prandtl-Glauert compressibility correction:
                mach = norm(Veff)/fs.speedofsound
                PG = 1/sqrt(1 - mach^2)
                Fbli *= PG
            end

            # --- Calculate forces on the right bound vortex --- #

            rc = right_center(receiving[I])

            # freestream velocity
            Veff = freestream_velocity(fs)

            # rotational velocity
            Veff += rotational_velocity(rc, fs, ref)

            # additional velocity field
            if !isnothing(additional_velocity)
                Veff += additional_velocity(rc)
            end

            # velocity due to surface motion
            if !isnothing(Vv)
                Veff += Vv[isurf][I[1], I[2]+1]
            end

            # NOTE: We don't include induced velocity in the effective velocity
            # for the vertical segments because its influence is likely negligible
            # once we take the cross product with the bound vortex vector. This
            # is also assumed in AVL. This could change in the future.

            # steady part of Kutta-Joukowski theorem
            Γri = Γ[iΓ+i]
            Δs = right_vector(receiving[I])
            Fbri = RHO*Γri*cross(Veff, Δs)

            if prandtl_glauert
                # Prandtl-Glauert compressibility correction:
                mach = norm(Veff)/fs.speedofsound
                PG = 1/sqrt(1 - mach^2)
                Fbri *= PG
            end

            # store panel circulation, velocity, and forces
            q = 1/2*RHO*ref.V^2

            props[isurf][i] = PanelProperties(Γ[iΓ+i]/ref.V, Vi/ref.V,
                Fbi/(q*ref.S), Fbli/(q*ref.S), Fbri/(q*ref.S))
        end

        # increment Γ index for receiving panels
        iΓ += nr
    end

    return props
end

"""
    near_field_forces_derivatives!(properties, dproperties, surfaces, reference,
        freestream, Γ, dΓ; dΓdt, additional_velocity, Vh, Vv, symmetric, nwake,
        surface_id, wake_finite_core, wake_shedding_locations, trailing_vortices, xhat)

Version of [`near_field_forces!`](@ref) that also calculates the derivatives of
the local panel forces with respect to the freestream variables.
"""
near_field_forces_derivatives!

@inline function near_field_forces_derivatives!(props, dprops, surfaces, wakes,
    ref, fs, Γ, dΓ; dΓdt, additional_velocity, Vh, Vv, symmetric, nwake,
    surface_id, wake_finite_core, wake_shedding_locations, trailing_vortices, xhat, prandtl_glauert)

    # unpack derivatives
    props_a, props_b, props_p, props_q, props_r = dprops
    Γ_a, Γ_b, Γ_p, Γ_q, Γ_r = dΓ

    # number of surfaces
    nsurf = length(surfaces)

    # loop through receiving surfaces
    iΓ = 0 # index for accessing Γ
    for isurf = 1:nsurf

        receiving = surfaces[isurf]
        nr = length(receiving)
        nr1, nr2 = size(receiving)
        cr = CartesianIndices(receiving)

        # loop through receiving panels
        for i = 1:length(receiving)

            I = cr[i]

            # --- Calculate forces on the panel bound vortex -- #
            rc = top_center(receiving[i])

            # freestream velocity
            Vi, dVi = freestream_velocity_derivatives(fs)
            Vi_a, Vi_b = dVi

            # rotational velocity
            Vrot, dVrot = rotational_velocity_derivatives(rc, fs, ref)
            Vi += Vrot
            Vi_p, Vi_q, Vi_r = dVrot

            # additional velocity field
            if !isnothing(additional_velocity)
                Vi += additional_velocity(rc)
            end

            # velocity due to surface motion
            if !isnothing(Vh)
                Vi += Vh[isurf][i]
            end

            # induced velocity from surfaces and wakes
            jΓ = 0 # index for accessing Γ
            for jsurf = 1:nsurf

                # number of panels on sending surface
                sending = surfaces[jsurf]
                Ns = length(sending)

                # see if wake panels are being used
                wake_panels = nwake[jsurf] > 0

                # check if we need to shift shedding locations
                if isnothing(wake_shedding_locations)
                    shedding_locations = nothing
                else
                    shedding_locations = wake_shedding_locations[jsurf]
                end

                # extract circulation values corresonding to the sending surface
                vΓ = view(Γ, jΓ+1:jΓ+Ns)

                vΓ_a = view(Γ_a, jΓ+1:jΓ+Ns)
                vΓ_b = view(Γ_b, jΓ+1:jΓ+Ns)
                vΓ_p = view(Γ_p, jΓ+1:jΓ+Ns)
                vΓ_q = view(Γ_q, jΓ+1:jΓ+Ns)
                vΓ_r = view(Γ_r, jΓ+1:jΓ+Ns)

                vdΓ = (vΓ_a, vΓ_b, vΓ_p, vΓ_q, vΓ_r)

                # induced velocity from this surface
                if isurf == jsurf
                    # induced velocity on self
                    Vind, dVind = induced_velocity_derivatives(I, surfaces[jsurf], vΓ, vdΓ;
                        finite_core = surface_id[isurf] != surface_id[jsurf],
                        wake_shedding_locations = shedding_locations,
                        symmetric = symmetric[jsurf],
                        trailing_vortices = trailing_vortices[jsurf] && !wake_panels,
                        xhat = xhat)
                else
                    # induced velocity on another surface
                    Vind, dVind = induced_velocity_derivatives(rc, surfaces[jsurf], vΓ, vdΓ;
                        finite_core = surface_id[isurf] != surface_id[jsurf],
                        wake_shedding_locations = shedding_locations,
                        symmetric = symmetric[jsurf],
                        trailing_vortices = trailing_vortices[jsurf] && !wake_panels,
                        xhat = xhat)
                end

                Vind_a, Vind_b, Vind_p, Vind_q, Vind_r = dVind

                Vi += Vind

                Vi_a += Vind_a
                Vi_b += Vind_b
                Vi_p += Vind_p
                Vi_q += Vind_q
                Vi_r += Vind_r

                # induced velocity from corresponding wake
                if wake_panels
                    Vi += induced_velocity(rc, wakes[jsurf];
                        finite_core = wake_finite_core[jsurf] || surface_id[isurf] != surface_id[jsurf],
                        symmetric = symmetric[jsurf],
                        nc = nwake[jsurf],
                        trailing_vortices = trailing_vortices[jsurf],
                        xhat = xhat)
                end

                jΓ += Ns # increment Γ index for sending panels
            end

            # steady part of Kutta-Joukowski theorem
            if I[1] == 1
                Γi = Γ[iΓ+i]

                Γi_a = Γ_a[iΓ+i]
                Γi_b = Γ_b[iΓ+i]
                Γi_p = Γ_p[iΓ+i]
                Γi_q = Γ_q[iΓ+i]
                Γi_r = Γ_r[iΓ+i]
            else
                Γi = Γ[iΓ+i] - Γ[iΓ+i-1]

                Γi_a = Γ_a[iΓ+i] - Γ_a[iΓ+i-1]
                Γi_b = Γ_b[iΓ+i] - Γ_b[iΓ+i-1]
                Γi_p = Γ_p[iΓ+i] - Γ_p[iΓ+i-1]
                Γi_q = Γ_q[iΓ+i] - Γ_q[iΓ+i-1]
                Γi_r = Γ_r[iΓ+i] - Γ_r[iΓ+i-1]
            end

            # bound vortex vector
            Δs = top_vector(receiving[I])

            tmp = cross(Vi, Δs)

            Fbi = RHO*Γi*tmp

            Fbi_a = RHO*(Γi_a*tmp + Γi*cross(Vi_a, Δs))
            Fbi_b = RHO*(Γi_b*tmp + Γi*cross(Vi_b, Δs))
            Fbi_p = RHO*(Γi_p*tmp + Γi*cross(Vi_p, Δs))
            Fbi_q = RHO*(Γi_q*tmp + Γi*cross(Vi_q, Δs))
            Fbi_r = RHO*(Γi_r*tmp + Γi*cross(Vi_r, Δs))

            if !isnothing(dΓdt)
                # unsteady part of Kutta-Joukowski theorem

                #TODO: decide whether to divide by perpindicular velocity like
                # Drela does in ASWING?

                dΓdti = I[1] == 1 ? dΓdt[iΓ+i] : (dΓdt[iΓ+i] + dΓdt[iΓ+i-1])/2
                c = receiving[I].chord
                Fbi += RHO*dΓdti*c*tmp/norm(Vi)

            end

            if prandtl_glauert
                # Prandtl-Glauert compressibility correction:
                mach = norm(Vi)/fs.speedofsound
                PG = 1/sqrt(1 - mach^2)
                Fbi *= PG
            end

            # --- Calculate forces for the left bound vortex --- #

            rc = left_center(receiving[I])

            # freestream velocity
            Vfs, dVfs = freestream_velocity_derivatives(fs)
            Veff = Vfs
            Veff_a, Veff_b = dVfs

            # rotational velocity
            Vrot, dVrot = rotational_velocity_derivatives(rc, fs, ref)
            Veff += Vrot
            Veff_p, Veff_q, Veff_r = dVrot

            # additional velocity field
            if !isnothing(additional_velocity)
                Veff += additional_velocity(rc)
            end

            # velocity due to surface motion
            if !isnothing(Vv)
                Veff += Vv[isurf][I[1], I[2]]
            end

            # NOTE: We don't include induced velocity in the effective velocity
            # for the vertical segments because its influence is likely negligible
            # once we take the cross product with the bound vortex vector. This
            # is also assumed in AVL. This could change in the future.

            # steady part of Kutta-Joukowski theorem
            Γli = Γ[iΓ+i]

            Γli_a = Γ_a[iΓ+i]
            Γli_b = Γ_b[iΓ+i]
            Γli_p = Γ_p[iΓ+i]
            Γli_q = Γ_q[iΓ+i]
            Γli_r = Γ_r[iΓ+i]

            Δs = left_vector(receiving[I])

            tmp = cross(Veff, Δs)

            Fbli = RHO*Γli*tmp

            Fbli_a = RHO*(Γli_a*tmp + Γli*cross(Veff_a, Δs))
            Fbli_b = RHO*(Γli_b*tmp + Γli*cross(Veff_b, Δs))
            Fbli_p = RHO*(Γli_p*tmp + Γli*cross(Veff_p, Δs))
            Fbli_q = RHO*(Γli_q*tmp + Γli*cross(Veff_q, Δs))
            Fbli_r = RHO*(Γli_r*tmp + Γli*cross(Veff_r, Δs))

            if prandtl_glauert
                # Prandtl-Glauert compressibility correction:
                mach = norm(Veff)/fs.speedofsound
                PG = 1/sqrt(1 - mach^2)
                Fbli *= PG
            end

            # --- Calculate forces on the right bound vortex --- #

            rc = right_center(receiving[I])

            # freestream velocity
            Vfs, dVfs = freestream_velocity_derivatives(fs)
            Veff = Vfs
            Veff_a, Veff_b = dVfs

            # rotational velocity
            Vrot, dVrot = rotational_velocity_derivatives(rc, fs, ref)
            Veff += Vrot
            Veff_p, Veff_q, Veff_r = dVrot

            # additional velocity field
            if !isnothing(additional_velocity)
                Veff += additional_velocity(rc)
            end

            # velocity due to surface motion
            if !isnothing(Vv)
                Veff += Vv[isurf][I[1], I[2]+1]
            end

            # NOTE: We don't include induced velocity in the effective velocity
            # for the vertical segments because its influence is likely negligible
            # once we take the cross product with the bound vortex vector. This
            # is also assumed in AVL. This could change in the future.

            # steady part of Kutta-Joukowski theorem
            Γri = Γ[iΓ+i]

            Γri_a = Γ_a[iΓ+i]
            Γri_b = Γ_b[iΓ+i]
            Γri_p = Γ_p[iΓ+i]
            Γri_q = Γ_q[iΓ+i]
            Γri_r = Γ_r[iΓ+i]

            Δs = right_vector(receiving[I])

            tmp = cross(Veff, Δs)

            Fbri = RHO*Γri*tmp

            Fbri_a = RHO*(Γri_a*tmp + Γri*cross(Veff_a, Δs))
            Fbri_b = RHO*(Γri_b*tmp + Γri*cross(Veff_b, Δs))
            Fbri_p = RHO*(Γri_p*tmp + Γri*cross(Veff_p, Δs))
            Fbri_q = RHO*(Γri_q*tmp + Γri*cross(Veff_q, Δs))
            Fbri_r = RHO*(Γri_r*tmp + Γri*cross(Veff_r, Δs))

            if prandtl_glauert
                # Prandtl-Glauert compressibility correction:
                mach = norm(Veff)/fs.speedofsound
                PG = 1/sqrt(1 - mach^2)
                Fbri *= PG
            end

            # store panel circulation, velocity, and forces
            q = 1/2*RHO*ref.V^2

            props[isurf][I] = PanelProperties(Γ[iΓ+i]/ref.V, Vi/ref.V, Fbi/(q*ref.S),
                Fbli/(q*ref.S), Fbri/(q*ref.S))

            props_a[isurf][I] = PanelProperties(Γ_a[iΓ+i]/ref.V, Vi_a/ref.V, Fbi_a/(q*ref.S),
                Fbli_a/(q*ref.S), Fbri_a/(q*ref.S))
            props_b[isurf][I] = PanelProperties(Γ_b[iΓ+i]/ref.V, Vi_b/ref.V, Fbi_b/(q*ref.S),
                Fbli_b/(q*ref.S), Fbri_b/(q*ref.S))
            props_p[isurf][I] = PanelProperties(Γ_p[iΓ+i]/ref.V, Vi_p/ref.V, Fbi_p/(q*ref.S),
                Fbli_p/(q*ref.S), Fbri_p/(q*ref.S))
            props_q[isurf][I] = PanelProperties(Γ_q[iΓ+i]/ref.V, Vi_q/ref.V, Fbi_q/(q*ref.S),
                Fbli_q/(q*ref.S), Fbri_q/(q*ref.S))
            props_r[isurf][I] = PanelProperties(Γ_r[iΓ+i]/ref.V, Vi_r/ref.V, Fbi_r/(q*ref.S),
                Fbli_r/(q*ref.S), Fbri_r/(q*ref.S))
        end

        # increment Γ index for receiving panels
        iΓ += nr
    end

    return props, dprops
end

"""
    lifting_line_forces!(lifting_line_properties, lifting_lines, props, surfaces, ref)

Calculate the force and moment coefficients (per unit span) for each spanwise
segment of a lifting line representation of the geometry. The loading is written
to the `lifting_line_properties` argument and returned.

This function requires that a near-field analysis has been performed on `system`
to obtain panel forces.

# Arguments
 - `lifting_line_properties`: Vector{Vector{LiftingLineProperties}} representing
    the loading along a lifting line for each lifting surface (see [`LiftingLineProperties`](@ref))
 - `lifting_lines`: Vector{Vector{LiftingLineSegment}} representing the
    lifting line geometry of each lifting surface (see [`LiftingLineSegment`](@ref))
 - `props`: Vector of Matrix of surface panel-specific properties (see [`PanelProperties`](@ref))
 - `surfaces`: Vector of matrices of surface panels (see [`SurfacePanel`](@ref))
 - `ref`: Reference parameters (see [`Reference`](@ref))
"""
function lifting_line_forces!(lifting_line_properties, lifting_lines, props, surfaces, ref)
    # number of surfaces
    nsurf = length(surfaces)

    TF = promote_type(eltype(eltype(eltype(lifting_line_properties))),
                      eltype(eltype(eltype(lifting_lines))),
                      eltype(eltype(eltype(props))),
                      eltype(eltype(eltype(surfaces))))

    # iterate through each lifting surface
    for isurf = 1:nsurf
        nc, ns = size(surfaces[isurf])
        # extract current surface panels and panel properties
        panels = surfaces[isurf]
        properties = props[isurf]
        lifting_line = lifting_lines[isurf]
        for j = 1:ns
            # calculate segment length
            rls = lifting_line[j].rl
            rrs = lifting_line[j].rr
            ds = norm(rrs - rls)
            # calculate reference location
            rs = (rls + rrs)/2
            # calculate reference chord
            # cs = (c[isurf][j] + c[isurf][j+1])/2
            cs = (lifting_line[j].chord_l + lifting_line[j].chord_r)/2

            # calculate section force and moment coefficients
            cfj = @SVector zeros(TF, 3)
            cmj = @SVector zeros(TF, 3)
            for i = 1:nc
                # add influence of bound vortex
                rb = top_center(panels[i,j])
                cfb = properties[i,j].cfb
                cfj += cfb
                cmj += cross(rb - rs, cfb)
                # add influence of left vortex leg
                rl = left_center(panels[i,j])
                cfl = properties[i,j].cfl
                cfj += cfl
                cmj += cross(rl - rs, cfl)
                # add influence of right vortex leg
                rr = right_center(panels[i,j])
                cfr = properties[i,j].cfr
                cfj += cfr
                cmj += cross(rr - rs, cfr)
            end
            # update normalization
            cfj *= ref.S/(ds*cs)
            cmj *= ref.S/(ds*cs^2)

            # set the viscous coefficients to zero
            cfvj = zero(cfj)
            cmvj = zero(cmj)
            # Not calculating all the velocities here (we don't need it for the inviscid forces), so set them to zero.
            V = zero(cfj)
            V_fs = zero(cfj)
            V_rot = zero(cfj)
            V_add = zero(cfj)
            V_sm = zero(cfj)
            alpha = zero(eltype(cfj))
            alpha_from_twist = zero(eltype(cfj))
            phi = zero(eltype(cfj))
            V_airfoil = zero(eltype(cfj))
            cl = zero(eltype(cfj))
            cd = zero(eltype(cfj))

            # Get a unit vector pointing from the trailing edge to the leading edge.
            u_chord = (lifting_line[j].u_chord_l + lifting_line[j].u_chord_r)/2
            u_chord /= norm(u_chord)
            lifting_line_properties[isurf][j] = LiftingLineProperties(rs, cs, ds, cfj, cfvj, cmj, cmvj, u_chord, V, V_fs, V_rot, V_add, V_sm, alpha, alpha_from_twist, phi, V_airfoil, cl, cd)
        end
    end
    return lifting_line_properties
end

function lifting_line_viscous_forces!(lifting_line_properties, lifting_lines, drag_polar, clmax, props, surfaces, wakes, ref, fs, Γ; additional_velocity, Vh, symmetric, nwake, surface_id, wake_finite_core, wake_shedding_locations, trailing_vortices, xhat, re_correction, rotation_correction)
    # number of surfaces
    nsurf = length(lifting_line_properties)

    # iterate through each lifting surface
    for isurf = 1:nsurf
        ns = length(lifting_line_properties[isurf])
        # loop through each chordwise set of panels
        for j = 1:ns
            llp = lifting_line_properties[isurf][j]
            # get the invsicid lifting line coefficient for this spanwise station.
            cfj = llp.cf

            # get segment length and reference location
            ds = llp.ds
            # Get the unit vector pointing from the trailing edge to the leading edge.
            u_chord = llp.u_chord
            rs = llp.r
            # rs = llp.r .+ 1.0*llp.c*[-1, 0, 0]
            # rs = llp.r .+ 0.5*llp.c*[-1, 0, 0]
            # rs = llp.r + 1000.0*llp.c*u_chord
            # rs = llp.r .+ 1000.0*llp.c*[-1, 0, 0]
            # @show rs

            # Get a unit vector in the direction of the lifting line's rotation.
            # Ah, actually, I want the opposite of that, since the velocity
            # `lifting_line_velocity` gives me is in the rotor frame, sort of.
            TF = eltype(llp.r)
            radii = sqrt(llp.r[2]^2 + llp.r[3]^2)
            theta_hat = -SVector{3,TF}(zero(TF), -llp.r[3], llp.r[2])/radii

            # Get a unit vector in the direction of the freestream:
            # fs_hat = SVector{3,TF}(one(TF), zero(TF), zero(TF))
            # cos(alpha)^2*cos(beta)^2 + sin(beta)^2 + sin(alpha)^2*cos(beta)^2 = cos(beta)^2*(cos(alpha)^2 + sin(alpha)^2) + sin(beta)^2 = cos(beta)^2 + sin(beta)^2 = 1
            fs_hat = freestream_velocity(fs)/fs.Vinf

            # println("radii should be $radii m")

            # calculate the local velocity for the current chordwise set of panels at the reference location
            V, V_fs, V_rot, V_add, V_sm = lifting_line_velocity(isurf, j, rs, surfaces, wakes, ref, fs, Γ; additional_velocity, Vh, symmetric, nwake, surface_id, wake_finite_core, wake_shedding_locations, trailing_vortices, xhat)
            # @show V - (V_fs + V_sm)

            # normalized direction along the lifting line's length
            rls = lifting_lines[isurf][j].rl
            rrs = lifting_lines[isurf][j].rr
            span_dir = (rrs - rls)/ds
            # @show norm(span_dir)

            # Now I want to know the lift force, which is the force in the
            # direction that's normal to both the local velocity and span_dir
            # I don't like this, because I'm not sure which order I should cross
            # the vectors.
            # lift_dir = cross(span_dir, V/norm(V))
            # V_dir = V/norm(V)
            # @show norm(V_dir)
            # lift_dir_tmp = cross(span_dir, V_dir)
            # lift_dir = lift_dir_tmp/norm(lift_dir_tmp)
            lift_dir_tmp = cross(span_dir, V)
            lift_dir = lift_dir_tmp/norm(lift_dir_tmp)
            # @show span_dir lift_dir
            # @show norm(lift_dir) sqrt(lift_dir[1]^2 + lift_dir[2]^2 + lift_dir[3]^2)
            # Now, cfj is normalized by what? ds and cs. That's good.
            # And by what else?
            # 1/2*RHO*ref.V^2
            # So I want the 1/2 rho, but I think I need it to be nondimenialized
            # by the velocity component normal to the span and lift directions.
            # This isn't the same as norm(V) since that might have a component
            # along the span direction, and I don't want that.
            drag_dir_tmp = cross(lift_dir, span_dir)
            drag_dir = drag_dir_tmp/norm(drag_dir_tmp)
            # @show norm(drag_dir)
            V_airfoil = dot(drag_dir, V)

            # Get the local angle of attack.
            # Use the definition of the dot product to get the angle between the -u_chord vector and the drag direction.
            # alpha = acos(dot(-u_chord, drag_dir))
            # radii = sqrt(rs[2]^2 + rs[3]^2)
            # omega = sqrt(V[2]^2 + V[3]^2)/radii
            # @show V[1] omega/(2*pi)*60
            # Find the angle the chord line makes with the rotor rotation plane.
            # Need to switch the sign on u_chord[1] since positive is downstream.
            twist = atan(-u_chord[1], sqrt(u_chord[2]^2 + u_chord[3]^2))
            # @show dot(V, span_dir)

            # So now I should be able to get the axial and tangential velocities.
            Vx_lame = dot(V, fs_hat)
            Vy_lame = dot(V, theta_hat)
            omega = Vy_lame/radii
            # Rtip = 0.30479999999999996
            # println("radii/Rtip = $(radii/Rtip), Vx_lame = $(Vx_lame) m/s, Vy_lame = $(Vy_lame) m/s, twist = $(twist*180/pi) deg")
            alpha_from_twist = twist - atan(Vx_lame, Vy_lame)

            # So remove the ref.V^2 and replace it with V_airfoil^2. Now I should
            # have the force divided by 1/2*RHO*V_airfoil^2*ds*cs
            cfj_airfoil = cfj * ref.V^2/V_airfoil^2
            # Now I want the lift coefficient, which is the force in the lift
            # direction.
            clj_airfoil = dot(cfj_airfoil, lift_dir)

            # Now I want the drag coefficient.
            cdj_airfoil = drag_polar(clj_airfoil)

            # Limit the lift coefficient in such a way that when added to the inviscid value, it will be at most clmax (or at least -clmax).
            # But this will be possibly modified by the rotation correction.
            if clj_airfoil > clmax
                clj_limited = (clmax - clj_airfoil)
            elseif clj_airfoil < -clmax
                clj_limited = (-clmax - clj_airfoil)
            else
                clj_limited = zero(clj_airfoil)
            end

            # Reynolds number correction.
            # First need a local Reynolds number based on chord.
            Re = V_airfoil*llp.c/fs.viscocity
            # Then do the correction.
            cdj_airfoil = re_correction(cdj_airfoil, Re)

            # Need these things for the rotational stall correction:
            #
            #   * radial location along the blade.
            #   * Rtip
            #   * Vy, the velocity along the blade's motion, without induction.
            #   * Vx, the velocity normal to the blade's motion and span, without induction.
            #   * alpha, the local angle of attack.
            #   * phi, the direction the velocity makes with the rotation plane of the rotor.

            # How do I get things like r and Rtip?
            # `r` is the distance from the rotation axis to the current lifting line segment point `rs`.
            rot_axis = lifting_lines[isurf][j].rot_axis
            # Do I know where the rotation axis is, though?
            # No, but it would be easiy to keep track of.
            # I guess I'll add that too.
            rot_origin = lifting_lines[isurf][j].rot_origin

            # Now this is a line from the rotor origin to the current point (span vector) and the norm of that is the distance I want.
            svec = rs - rot_origin
            r = norm(svec)

            # Now, how do I get Rtip?
            # I guess that would be the same thing as r, but for the tip of blade.
            # But where is the tip of the blade?
            # It would be the position of the last lifting line segment.
            # Oh, but is it the "right" or "left?"
            # And am I sure it's at the end?
            # Well, the tip would be either the inner-most or outer-most
            # segment.
            # The one I want is of course farthest from the origin.
            # And it looks like I want the right point of the last segment,
            # and the left point of the first segment.
            to_Rtip_maybe1 = right(lifting_lines[isurf][end]) - rot_origin
            to_Rtip_maybe2 = left(lifting_lines[isurf][begin]) - rot_origin
            Rtip_maybe1 = norm(to_Rtip_maybe1)
            Rtip_maybe2 = norm(to_Rtip_maybe2)
            if Rtip_maybe1 > Rtip_maybe2
                Rtip = Rtip_maybe1
            else
                Rtip = Rtip_maybe2
            end

            # println("Rtip = $Rtip, should be $(0.5*24*0.0254)")
            Rhub = 0.2*Rtip
            # println("r = $r, should be $((Rtip - Rhub)/ns*(j - 1 + 0.5) + Rhub)")

            # Now I'll work on Vy and Vx.
            # We eventually want to get the velocity normal to the lifting line span and it's rotation direction, without induction.
            # So the first step is to get the velocity without induction:
            V_no_ind_tmp = V_fs + V_rot + V_add + V_sm
            # Do I need to take out any velocity along the span?
            # I guess I will.
            V_no_ind = V_no_ind_tmp - dot(V_no_ind_tmp, span_dir)*span_dir

            # `t_hat` is a tangential unit normal vector.
            # It points opposite the direction of the lifting line's rotation, since the grid motion in VortexLattice.jl is actually relative to the lifting surface.
            # And that's what I want anyways.
            t_hat = normalize(cross(svec, rot_axis))

            # Now we can find the velocity in the tangential direction.
            Vy = dot(V_no_ind, t_hat)

            # Also need Vx, the velocity normal to the lifting line and the t_hat direction.
            # I think that direction will usually be this:
            #
            #     n_hat = normalize(cross(t_hat, svec))
            #
            # But to be sure, I will check if it is in the same direction as the chord line going from leading to trailing edge.
            # The dot product of two vectors is |a||b|cos(θ), where θ is the angle between them.
            # If that angle is greater than 90° or less than -90°, then they're not in the same direction.
            # And that would imply that cos(θ) is < 0.
            # So maybe I can just multiply by the sign of that.
            n_hat_maybe = normalize(cross(t_hat, svec))
            n_hat = n_hat_maybe * sign(dot(n_hat_maybe, -u_chord))

            # So now Vx is just this:
            Vx = dot(V_no_ind, n_hat)
            # println("Re_induction = $(Re), Re_fs = $(sqrt(Vx^2 + Vy^2)*llp.c/fs.viscocity), factor_induction = $(re_correction(1.0, Re)), cdj_airfoil = $(cdj_airfoil)")

            # println("Vx = $Vx, Vx_lame = $(Vx_lame), should be 38.274895154399985")
            # println("Vy = $Vy, Vy_lame = $(Vy_lame), should be $(7100*2*pi/60*r)")
            omega = Vy / r
            # println("omega from Vy and r = $(Vy/r*60/(2*pi)) rev/min")

            # Now need alpha and phi.
            # alpha, the angle of attack, is just the angle the total velocity makes with the chord line.
            # We already know that the total velocity direction: `drag_dir`.
            # So, let's find that.
            alpha = acos(dot(drag_dir, -u_chord))

            # For phi, we want the direction the total velocity makes with the rotation plane.
            # That should be the angle the drag_dir makes with the tangential direction.
            phi = acos(dot(drag_dir, t_hat))

            # alpha = twist - phi
            phi_from_twist = twist - alpha
            # println("phi = $(phi*180/pi)°, phi_from_twist = $(phi_from_twist*180/pi)°")
            # println("alpha = $(alpha*180/pi)°, alpha_from_twist = $(alpha_from_twist*180/pi)°")

            # Now we have everything except the chord.
            chord = llp.c
            # To do the rotation correction properly, we need the total lift.
            clj_total = clj_airfoil + clj_limited
            # println("before rc: clj_total = $(clj_total), cdj_airfoil = $(cdj_airfoil)")
            if ! isnothing(rotation_correction)
                # # println("r/Rtip = $(r/Rtip), clj_total = $clj_total, cdj_airfoil = $cdj_airfoil, Vx = $Vx, Vy = $Vy, chord/r = $(chord/r), alpha = $(alpha*180/pi) deg, phi = $(phi*180/pi) deg")
                # println("r/Rtip = $(r/Rtip), clj_total = $clj_total, cdj_airfoil = $cdj_airfoil, chord/r = $(chord/r), r/Rtip = $(r/Rtip), Vy/Vx*Rtip/r = $(Vy/Vx*Rtip/r), alpha = $(alpha), phi = $(phi)")
                # clj_total, cdj_airfoil = CCBlade.rotation_correction(rotation_correction, clj_total, cdj_airfoil, chord/r, r/Rtip, Vy/Vx*Rtip/r, alpha, phi)

                # Copying this from CCBlade.jl/src/airfoils.jl
                # Du-Selig correction for lift
                tsr = Vy/Vx*Rtip/r
                rR = r/Rtip
                cr = chord/r
                Lambda = tsr / sqrt(1 + tsr^2)
                expon = rotation_correction.d / (Lambda * rR)
                fcl = 1.0/rotation_correction.m*(1.6*cr/0.1267*(rotation_correction.a-cr^expon)/(rotation_correction.b+cr^expon)-1)

                mach2 = (Vx^2 + Vy^2)/fs.speedofsound^2
                if clj_airfoil > clmax
                    cl_linear_minus_cl = clj_airfoil - clmax
                elseif clj_airfoil < -clmax
                    cl_linear_minus_cl = clj_airfoil - (-clmax)
                else
                    # We'll assume that the only difference between the linear lift slope and actual is due to the PG correction (which is what happens in CCBlade.jl when using a linear lift slope).
                    cl_linear_minus_cl = clj_airfoil*(1 - 1/sqrt(1 - mach2))
                end

                # adjustment for max correction
                amax = atan(1/0.12) - 5*pi/180  # account for singularity in Eggers (not pi/2)
                alpha_max_corr = 30*pi/180  # default from CCBlade.jl
                if abs(alpha) >= amax 
                    adj = 0.0
                elseif abs(alpha) > alpha_max_corr
                    adj = ((amax-abs(alpha))/(amax-alpha_max_corr))^2
                else
                    adj = 1.0
                end

                # increment in cl
                deltacl = adj*fcl*(cl_linear_minus_cl)
                clj_total += deltacl
    
                # Eggers correction for drag
                deltacd = deltacl * (sin(phi) - 0.12*cos(phi))/(cos(phi) + 0.12*sin(phi))  # note that we can actually use phi instead of alpha as is done in airfoilprep.py b/c this is done at each iteration
                # cd += deltacd
                cdj_airfoil += deltacd

                clj_limited = clj_total - clj_airfoil
                # println("after rc:  clj_total = $(clj_total), cdj_airfoil = $(cdj_airfoil)")
            end

            # Get the total viscous force.
            cfvj = (cdj_airfoil*drag_dir + clj_limited*lift_dir) * V_airfoil^2/ref.V^2

            # What do I do about the moment coefficient?
            # Nothing, since I'm assuming that the viscous load passes
            # through the center of the lifting line element, and so the
            # moment arm is zero.
            cmvj = zero(cfvj)

            lifting_line_properties[isurf][j] = LiftingLineProperties(llp.r, llp.c, llp.ds, llp.cf, cfvj, llp.cm, cmvj, u_chord, V/ref.V, V_fs/ref.V, V_rot/ref.V, V_add/ref.V, V_sm/ref.V, alpha, alpha_from_twist, phi, V_airfoil/ref.V, clj_total, cdj_airfoil)
        end
    end
end

"""
    body_forces(system; kwargs...)

Return the body force coefficients given the panel properties for `surfaces`

Note that this function assumes that a near-field analysis has already been
performed to obtain the panel forces.

# Arguments
 - `system`: Object of type [`System`](@ref) which holds system properties

# Keyword Arguments
 - `frame`: frame in which to return `CF` and `CM`, options are [`Body()`](@ref) (default),
   [`Stability()`](@ref), and [`Wind()`](@ref)`
"""
function body_forces(system::System{TF}; frame = Body()) where TF

    @assert system.near_field_analysis[] "Near field analysis required"

    # unpack parameters stored in `system`
    surfaces = system.surfaces # surface panels defining each surface
    properties = system.properties # panel properties
    ref = system.reference[] # reference parameters
    fs = system.freestream[] # freestream parameters
    symmetric = system.symmetric # symmetric flag for each surface

    return body_forces(surfaces, properties, ref, fs, symmetric, frame)
end

"""
    body_forces(surfaces, properties, reference, freestream, symmetric; kwargs...)

Return the body force coefficients given the panel properties for `surfaces`

Note that this function assumes that a near-field analysis has already been
performed to obtain the panel forces.

# Arguments:
 - `surfaces`: Collection of surfaces, where each surface is represented by a
    matrix of surface panels (see [`SurfacePanel`](@ref)) of shape (nc, ns)
    where `nc` is the number of chordwise panels and `ns` is the number of
    spanwise panels
 - `properties`: Surface properties for each surface, where surface
    properties for each surface are represented by a matrix of panel properties
    (see [`PanelProperties`](@ref)) of shape (nc, ns) where `nc` is the number
    of chordwise panels and `ns` is the number of spanwise panels
 - `reference`: Reference parameters (see [`Reference`](@ref))
 - `freestream`: Freestream parameters (see [`Freestream`]@ref)
 - `symmetric`: (required) Flag for each surface indicating whether a mirror image
   (across the X-Z plane) was used when calculating induced velocities
 - `frame`: frame in which to return `CF` and `CM`, options are [`Body()`](@ref) (default),
   [`Stability()`](@ref), and [`Wind()`](@ref)
"""
function body_forces(surfaces, properties, ref, fs, symmetric, frame)

    TF = eltype(eltype(eltype(properties)))

    # initialize body force coefficients
    CF = @SVector zeros(TF, 3)
    CM = @SVector zeros(TF, 3)

    # loop through all surfaces
    for isurf = 1:length(surfaces)

        # initialize surface contribution to body force coefficients
        CFi = @SVector zeros(TF, 3)
        CMi = @SVector zeros(TF, 3)

        # loop through all panels on this surface
        for i = 1:length(surfaces[isurf])

            # top bound vortex
            rc = top_center(surfaces[isurf][i])
            Δr = rc - ref.r
            cf = properties[isurf][i].cfb
            CFi += cf
            CMi += cross(Δr, cf)

            # left bound vortex
            rc = left_center(surfaces[isurf][i])
            Δr = rc - ref.r
            cf = properties[isurf][i].cfl
            CFi += cf
            CMi += cross(Δr, cf)

            # right bound vortex
            rc = right_center(surfaces[isurf][i])
            Δr = rc - ref.r
            cf = properties[isurf][i].cfr
            CFi += cf
            CMi += cross(Δr, cf)
        end

        # adjust forces from this surface to account for symmetry
        if symmetric[isurf]
            CFi = SVector(2*CFi[1], 0.0, 2*CFi[3])
            CMi = SVector(0.0, 2*CMi[2], 0.0)
        end

        # add to body forces
        CF += CFi
        CM += CMi

    end

    # add reference length in moment normalization
    reference_length = SVector(ref.b, ref.c, ref.b)
    CM = CM ./ reference_length

    # switch to specified frame
    CF, CM = body_to_frame(CF, CM, ref, fs, frame)

    return CF, CM
end


"""
    body_viscous_forces(system; kwargs...)

Return the viscous body force coefficients for the lifting surfaces defined in `system`

Note that this function assumes that a near-field and viscous lifting line
analyses have already been performed.

# Arguments
 - `system`: Object of type [`System`](@ref) which holds system properties

# Keyword Arguments
 - `frame`: frame in which to return `CF` and `CM`, options are [`Body()`](@ref) (default),
   [`Stability()`](@ref), and [`Wind()`](@ref)`
"""
function body_viscous_forces(system; frame = Body())
    @assert system.viscous_lifting_line[] "Viscous lifting line analysis required"

    # unpack parameters stored in `system`
    lifting_line_properties = system.lifting_line_properties # lifting line properties
    ref = system.reference[] # reference parameters
    fs = system.freestream[] # freestream parameters
    symmetric = system.symmetric # symmetric flag for each surface

    return body_viscous_forces(lifting_line_properties, ref, fs, symmetric, frame)
end

"""
    body_viscous_forces(lifting_line_properties, reference, freestream, symmetric, frame)

Return the viscous body force coefficients acting on the lifting lines

 - `lifting_line_properties`: Vector of vector of [`LiftingLineProperties`](@ref)
 - `reference`: Reference parameters (see [`Reference`](@ref))
 - `freestream`: Freestream parameters (see [`Freestream`]@ref)
 - `symmetric`: (required) Flag for each surface indicating whether a mirror image
   (across the X-Z plane) was used when calculating induced velocities
 - `frame`: frame in which to return `CF` and `CM`, options are [`Body()`](@ref) (default),
   [`Stability()`](@ref), and [`Wind()`](@ref)
"""
function body_viscous_forces(lifting_line_properties, ref, fs, symmetric, frame)
    TF = promote_type(eltype(eltype(eltype(lifting_line_properties))), eltype(ref), eltype(fs))

    # initialize body force coefficients
    CF = @SVector zeros(TF, 3)
    CM = @SVector zeros(TF, 3)

    # loop through all surfaces
    for isurf = 1:length(lifting_line_properties)

        # initialize surface contribution to body force coefficients
        CFi = @SVector zeros(TF, 3)
        CMi = @SVector zeros(TF, 3)

        # loop through all the lifting line elements for this surface
        # ns = size(r[isurf], 2)-1
        ns = length(lifting_line_properties[isurf])
        for j = 1:ns
            # calculate segment length
            # rls = SVector(r[isurf][1,j], r[isurf][2,j], r[isurf][3,j])
            # rrs = SVector(r[isurf][1,j+1], r[isurf][2,j+1], r[isurf][3,j+1])
            # ds = norm(rrs - rls)
            ds = lifting_line_properties[isurf][j].ds
            # calculate reference location
            # rs = (rls + rrs)/2
            rs = lifting_line_properties[isurf][j].r
            # calculate reference chord
            # cs = (c[isurf][j] + c[isurf][j+1])/2
            cs = lifting_line_properties[isurf][j].c
            # add viscous loading for this spanwise station
            Δr = rs - ref.r
            # cf = @view cfv[isurf][:, j]
            cf = lifting_line_properties[isurf][j].cfv
            # adjust non-dimensionalization to match what body_forces uses
            cf *= ds*cs/ref.S  # LiftingLineProperties.cfv is a SVector{3,TF}, so this won't mutate lifting_line_properties
            CFi += cf
            CMi += cross(Δr, cf)
        end

        # adjust forces from this surface to account for symmetry
        if symmetric[isurf]
            CFi = SVector(2*CFi[1], 0.0, 2*CFi[3])
            CMi = SVector(0.0, 2*CMi[2], 0.0)
        end

        # add to body forces
        CF += CFi
        CM += CMi

    end

    # add reference length in moment normalization
    reference_length = SVector(ref.b, ref.c, ref.b)
    CM = CM ./ reference_length

    # switch to specified frame
    CF, CM = body_to_frame(CF, CM, ref, fs, frame)

    return CF, CM
end

"""
    body_forces_derivatives(system)

Return the body force coefficients for the `system` and their derivatives with
respect to the freestream variables

Note that this function assumes that a near-field analysis has already been
performed to obtain the panel forces.

# Arguments:
 - `system`: Object of type `System` which holds system properties
"""
@inline function body_forces_derivatives(system::System)

    # float number type
    TF = eltype(system)

    @assert system.near_field_analysis[] "Near field analysis required"
    @assert system.derivatives[] "Derivative computations required"

    # unpack parameters stored in `system`
    surfaces = system.surfaces # surface panels defining each surface
    ref = system.reference[] # reference parameters
    fs = system.freestream[] # freestream parameters
    symmetric = system.symmetric # symmetric flag for each surface
    properties = system.properties
    props_a, props_b, props_p, props_q, props_r = system.dproperties

    # initialize body force coefficients
    CF = @SVector zeros(TF, 3)
    CM = @SVector zeros(TF, 3)

    CF_a = @SVector zeros(TF, 3)
    CF_b = @SVector zeros(TF, 3)
    CF_p = @SVector zeros(TF, 3)
    CF_q = @SVector zeros(TF, 3)
    CF_r = @SVector zeros(TF, 3)

    CM_a = @SVector zeros(TF, 3)
    CM_b = @SVector zeros(TF, 3)
    CM_p = @SVector zeros(TF, 3)
    CM_q = @SVector zeros(TF, 3)
    CM_r = @SVector zeros(TF, 3)

    # loop through all surfaces
    for isurf = 1:length(surfaces)

        # initialize surface contribution to body force coefficients
        CFi = @SVector zeros(TF, 3)
        CMi = @SVector zeros(TF, 3)

        CFi_a = @SVector zeros(TF, 3)
        CFi_b = @SVector zeros(TF, 3)
        CFi_p = @SVector zeros(TF, 3)
        CFi_q = @SVector zeros(TF, 3)
        CFi_r = @SVector zeros(TF, 3)

        CMi_a = @SVector zeros(TF, 3)
        CMi_b = @SVector zeros(TF, 3)
        CMi_p = @SVector zeros(TF, 3)
        CMi_q = @SVector zeros(TF, 3)
        CMi_r = @SVector zeros(TF, 3)

        for i = 1:length(surfaces[isurf])

            # top bound vortex
            rc = top_center(surfaces[isurf][i])
            Δr = rc - ref.r
            cf = properties[isurf][i].cfb
            CFi += cf
            CMi += cross(Δr, cf)

            cf_a = props_a[isurf][i].cfb
            CFi_a += cf_a
            CMi_a += cross(Δr, cf_a)

            cf_b = props_b[isurf][i].cfb
            CFi_b += cf_b
            CMi_b += cross(Δr, cf_b)

            cf_p = props_p[isurf][i].cfb
            CFi_p += cf_p
            CMi_p += cross(Δr, cf_p)

            cf_q = props_q[isurf][i].cfb
            CFi_q += cf_q
            CMi_q += cross(Δr, cf_q)

            cf_r = props_r[isurf][i].cfb
            CFi_r += cf_r
            CMi_r += cross(Δr, cf_r)

            # left bound vortex
            rc = left_center(surfaces[isurf][i])
            Δr = rc - ref.r
            cf = properties[isurf][i].cfl
            CFi += cf
            CMi += cross(Δr, cf)

            cf_a = props_a[isurf][i].cfl
            CFi_a += cf_a
            CMi_a += cross(Δr, cf_a)

            cf_b = props_b[isurf][i].cfl
            CFi_b += cf_b
            CMi_b += cross(Δr, cf_b)

            cf_p = props_p[isurf][i].cfl
            CFi_p += cf_p
            CMi_p += cross(Δr, cf_p)

            cf_q = props_q[isurf][i].cfl
            CFi_q += cf_q
            CMi_q += cross(Δr, cf_q)

            cf_r = props_r[isurf][i].cfl
            CFi_r += cf_r
            CMi_r += cross(Δr, cf_r)

            # right bound vortex
            rc = right_center(surfaces[isurf][i])
            Δr = rc - ref.r
            cf = properties[isurf][i].cfr
            CFi += cf
            CMi += cross(Δr, cf)

            cf_a = props_a[isurf][i].cfr
            CFi_a += cf_a
            CMi_a += cross(Δr, cf_a)

            cf_b = props_b[isurf][i].cfr
            CFi_b += cf_b
            CMi_b += cross(Δr, cf_b)

            cf_p = props_p[isurf][i].cfr
            CFi_p += cf_p
            CMi_p += cross(Δr, cf_p)

            cf_q = props_q[isurf][i].cfr
            CFi_q += cf_q
            CMi_q += cross(Δr, cf_q)

            cf_r = props_r[isurf][i].cfr
            CFi_r += cf_r
            CMi_r += cross(Δr, cf_r)
        end

        # adjust forces from this surface to account for symmetry
        if symmetric[isurf]
            CFi = SVector(2*CFi[1], 0.0, 2*CFi[3])
            CMi = SVector(0.0, 2*CMi[2], 0.0)

            CFi_a = SVector(2*CFi_a[1], 0.0, 2*CFi_a[3])
            CFi_b = SVector(2*CFi_b[1], 0.0, 2*CFi_b[3])
            CFi_p = SVector(2*CFi_p[1], 0.0, 2*CFi_p[3])
            CFi_q = SVector(2*CFi_q[1], 0.0, 2*CFi_q[3])
            CFi_r = SVector(2*CFi_r[1], 0.0, 2*CFi_r[3])

            CMi_a = SVector(0.0, 2*CMi_a[2], 0.0)
            CMi_b = SVector(0.0, 2*CMi_b[2], 0.0)
            CMi_p = SVector(0.0, 2*CMi_p[2], 0.0)
            CMi_q = SVector(0.0, 2*CMi_q[2], 0.0)
            CMi_r = SVector(0.0, 2*CMi_r[2], 0.0)
        end

        # add to body forces
        CF += CFi
        CM += CMi

        CF_a += CFi_a
        CF_b += CFi_b
        CF_p += CFi_p
        CF_q += CFi_q
        CF_r += CFi_r

        CM_a += CMi_a
        CM_b += CMi_b
        CM_p += CMi_p
        CM_q += CMi_q
        CM_r += CMi_r

    end

    # add reference length in moment normalization
    reference_length = SVector(ref.b, ref.c, ref.b)
    CM = CM ./ reference_length

    CM_a = CM_a ./ reference_length
    CM_b = CM_b ./ reference_length
    CM_p = CM_p ./ reference_length
    CM_q = CM_q ./ reference_length
    CM_r = CM_r ./ reference_length

    # pack up derivatives
    dCF = (CF_a, CF_b, CF_p, CF_q, CF_r)
    dCM = (CM_a, CM_b, CM_p, CM_q, CM_r)

    return CF, CM, dCF, dCM
end

"""
    body_forces_history(system, surface_history, property_history; frame=Body())

Return the body force coefficients `CF`, `CM` at each time step in `property_history`.

# Arguments:
 - `system`: Object of type [`System`](@ref) which holds system properties
 - `surface_history`: Vector of surfaces at each time step, where each surface is
    represented by a matrix of surface panels (see [`SurfacePanel`](@ref)) of shape
    (nc, ns) where `nc` is the number of chordwise panels and `ns` is the number
    of spanwise panels
 - `property_history`: Vector of surface properties for each surface at each
    time step, where surface properties are represented by a matrix of panel
    properties (see [`PanelProperties`](@ref)) of shape (nc, ns) where `nc` is
    the number of chordwise panels and `ns` is the number of spanwise panels

# Keyword Arguments
 - `frame`: frame in which to return `CF` and `CM`, options are [`Body()`](@ref) (default),
   [`Stability()`](@ref), and [`Wind()`](@ref)`
"""
function body_forces_history(system, surface_history::AbstractVector{<:AbstractVector{<:AbstractMatrix}},
    property_history; frame=Body())

    # unpack system parameters
    symmetric = system.symmetric
    ref = system.reference[]
    fs = system.freestream[]

    # float type
    TF = eltype(system)

    # number of time steps
    nt = length(property_history)

    # convert single freestream input to vector
    if isa(fs, Freestream)
        fs = fill(fs, nt)
    end

    # initialize time history coefficients
    CF = Vector{SVector{3, TF}}(undef, nt)
    CM = Vector{SVector{3, TF}}(undef, nt)

    # populate time history coefficients
    for it = 1:nt
        CF[it], CM[it] = body_forces(surface_history[it], property_history[it],
            ref, fs[it], symmetric, frame)
    end

    return CF, CM
end

"""
    body_viscous_forces_history(system, lifting_line_property_history; frame=Body())

Return the body viscous force coefficients `CFv`, `CMv` at each time step in `lifting_line_property_history`.

# Arguments:
 - `system`: Object of type [`System`](@ref) which holds system properties
 - `lifting_line_property_history`: Vector of lifting line properties for each surface at each
    time step, where lifting line properties are represented by a vector of lifting line
    properties (see [`LiftingLineProperties`](@ref)) of shape `ns` where `ns` is
    the number of spanwise panels

# Keyword Arguments
 - `frame`: frame in which to return `CF` and `CM`, options are [`Body()`](@ref) (default),
   [`Stability()`](@ref), and [`Wind()`](@ref)`
"""
function body_viscous_forces_history(system, lifting_line_property_history::AbstractVector{<:AbstractVector{<:AbstractVector{<:LiftingLineProperties}}}; frame=Body())
    @assert system.viscous_lifting_line[] "Viscous lifting line analysis required"

    # unpack system parameters
    symmetric = system.symmetric
    ref = system.reference[]
    fs = system.freestream[]

    # float type
    TF = eltype(system)

    # number of time steps
    nt = length(lifting_line_property_history)

    # convert single freestream input to vector
    if isa(fs, Freestream)
        fs = fill(fs, nt)
    end

    # initialize time history coefficients
    CFv = Vector{SVector{3, TF}}(undef, nt)
    CMv = Vector{SVector{3, TF}}(undef, nt)

    # populate time history coefficients
    for it = 1:nt
        CFv[it], CMv[it] = body_viscous_forces(lifting_line_property_history[it],
            ref, fs[it], symmetric, frame)
    end

    return CFv, CMv
end

"""
    lifting_line_coefficients(system; frame=Body())

Return the force and moment coefficients (per unit span) for each spanwise segment
of a lifting line representation of the geometry.

This function requires that a lifting line analysis has been performed on `system`
to obtain lifting line forces.

# Arguments
 - `system`: Object of type [`System`](@ref) that holds precalculated
    system properties.

# Keyword Arguments
 - `frame`: frame in which to return `cf` and `cm`, possible options are
    [`Body()`](@ref) (default), [`Stability()`](@ref), and [`Wind()`](@ref)`

# Return Arguments:
 - `cf`: Vector with length equal to the number of surfaces, with each element
    being a matrix with size (3, ns) which contains the x, y, and z direction
    force coefficients (per unit span) for each spanwise segment.
 - `cm`: Vector with length equal to the number of surfaces, with each element
    being a matrix with size (3, ns) which contains the x, y, and z direction
    moment coefficients (per unit span) for each spanwise segment.
"""
function lifting_line_coefficients(system; frame=Body())
    # check that a near field analysis has been performed
    @assert system.lifting_line_analysis[] "Lifting line analysis required"

    # extract reference properties
    ref = system.reference[]
    fs = system.freestream[]

    # number of surfaces
    nsurf = length(system.lifting_line_properties)

    TF = eltype(system)
    cf = Vector{Matrix{TF}}(undef, nsurf)
    cm = Vector{Matrix{TF}}(undef, nsurf)
    for isurf = 1:nsurf
        ns = length(system.lifting_line_properties[isurf])
        cf[isurf] = Matrix{TF}(undef, 3, ns)
        cm[isurf] = Matrix{TF}(undef, 3, ns)
        for j = 1:ns
            cfj = system.lifting_line_properties[isurf][j].cf
            cmj = system.lifting_line_properties[isurf][j].cm
            cfj, cmj = body_to_frame(cfj, cmj, ref, fs, frame)
            cf[isurf][:,j] = cfj
            cm[isurf][:,j] = cmj
        end
    end

    return cf, cm
end

function lifting_line_coefficients_history(system, surface_history, property_history, r, c; frame=Body())
    # unpack system parameters
    symmetric = system.symmetric
    ref = system.reference[]
    fs = system.freestream[]
    surfaces = system.surfaces

    # float type
    TF = eltype(system)

    # number of time steps
    nt = length(property_history)

    # convert single freestream input to vector
    if isa(fs, Freestream)
        fs = fill(fs, nt)
    end

    # initialize time history coefficients
    cf = Vector{Vector{Matrix{TF}}}(undef, nt)
    cm = Vector{Vector{Matrix{TF}}}(undef, nt)
    for it in 1:nt
        cf[it], cm[it] = lifting_line_coefficients(surface_history[it], property_history[it], ref, fs[it], r[it], c[it]; frame=frame)
    end

    return cf, cm
end

"""
    lifting_line_viscous_coefficients(system; frame=Body())

Return the viscous force and moment coefficients (per unit span) for each
spanwise segment of a lifting line representation of the geometry.

This function requires that a viscous lifting line analysis has been performed
on `system` to obtain lifting line forces.

# Arguments
 - `system`: Object of type [`System`](@ref) that holds precalculated
    system properties.

# Keyword Arguments
 - `frame`: frame in which to return `cf` and `cm`, possible options are
    [`Body()`](@ref) (default), [`Stability()`](@ref), and [`Wind()`](@ref)`

# Return Arguments:
 - `cfv`: Vector with length equal to the number of surfaces, with each element
    being a matrix with size (3, ns) which contains the x, y, and z direction
    force coefficients (per unit span) for each spanwise segment.
 - `cmv`: Vector with length equal to the number of surfaces, with each element
    being a matrix with size (3, ns) which contains the x, y, and z direction
    moment coefficients (per unit span) for each spanwise segment.
"""
function lifting_line_viscous_coefficients(system; frame=Body())
    # check that a near field analysis has been performed
    @assert system.viscous_lifting_line[] "Viscous lifting line analysis required"

    # extract reference properties
    ref = system.reference[]
    fs = system.freestream[]

    # number of surfaces
    nsurf = length(system.lifting_line_properties)

    TF = eltype(system)
    cf = Vector{Matrix{TF}}(undef, nsurf)
    cm = Vector{Matrix{TF}}(undef, nsurf)
    for isurf = 1:nsurf
        ns = length(system.lifting_line_properties[isurf])
        cf[isurf] = Matrix{TF}(undef, 3, ns)
        cm[isurf] = Matrix{TF}(undef, 3, ns)
        for j = 1:ns
            cfj = system.lifting_line_properties[isurf][j].cfv
            cmj = system.lifting_line_properties[isurf][j].cmv
            cfj, cmj = body_to_frame(cfj, cmj, ref, fs, frame)
            cf[isurf][:,j] = cfj
            cm[isurf][:,j] = cmj
        end
    end

    return cf, cm
end

# """
#     lifting_line_coefficients(system, r, c; frame=Body())

# Return the force and moment coefficients (per unit span) for each spanwise segment
# of a lifting line representation of the geometry.

# This function requires that a near-field analysis has been performed on `system`
# to obtain panel forces.

# # Arguments
#  - `system`: Object of type [`System`](@ref) that holds precalculated
#     system properties.
#  - `r`: Vector with length equal to the number of surfaces, with each element
#     being a matrix with size (3, ns+1) which contains the x, y, and z coordinates
#     of the resulting lifting line coordinates
#  - `c`: Vector with length equal to the number of surfaces, with each element
#     being a vector of length `ns+1` which contains the chord lengths at each
#     lifting line coordinate.

# # Keyword Arguments
#  - `frame`: frame in which to return `cf` and `cm`, possible options are
#     [`Body()`](@ref) (default), [`Stability()`](@ref), and [`Wind()`](@ref)`

# # Return Arguments:
#  - `cf`: Vector with length equal to the number of surfaces, with each element
#     being a matrix with size (3, ns) which contains the x, y, and z direction
#     force coefficients (per unit span) for each spanwise segment.
#  - `cm`: Vector with length equal to the number of surfaces, with each element
#     being a matrix with size (3, ns) which contains the x, y, and z direction
#     moment coefficients (per unit span) for each spanwise segment.
# """
# function lifting_line_coefficients(system, r, c; frame=Body())
#     TF = promote_type(eltype(system), eltype(eltype(r)), eltype(eltype(c)))
#     nsurf = length(system.surfaces)
#     cf = Vector{Matrix{TF}}(undef, nsurf)
#     cm = Vector{Matrix{TF}}(undef, nsurf)
#     for isurf = 1:nsurf
#         ns = size(system.surfaces[isurf], 2)
#         cf[isurf] = Matrix{TF}(undef, 3, ns)
#         cm[isurf] = Matrix{TF}(undef, 3, ns)
#     end
#     return lifting_line_coefficients!(cf, cm, system, r, c; frame)
# end

# """
#     lifting_line_coefficients!(cf, cm, system, r, c; frame=Body())

# In-place version of [`lifting_line_coefficients`](@ref)
# """
# function lifting_line_coefficients!(cf, cm, system, r, c; frame=Body())

#     # number of surfaces
#     nsurf = length(system.surfaces)

#     # check that a near field analysis has been performed
#     @assert system.near_field_analysis[] "Near field analysis required"

#     # extract reference properties
#     ref = system.reference[]
#     fs = system.freestream[]

#     # iterate through each lifting surface
#     for isurf = 1:nsurf
#         nc, ns = size(system.surfaces[isurf])
#         # extract current surface panels and panel properties
#         panels = system.surfaces[isurf]
#         properties = system.properties[isurf]
#         # loop through each chordwise set of panels
#         for j = 1:ns
#             # calculate segment length
#             rls = SVector(r[isurf][1,j], r[isurf][2,j], r[isurf][3,j])
#             rrs = SVector(r[isurf][1,j+1], r[isurf][2,j+1], r[isurf][3,j+1])
#             ds = norm(rrs - rls)
#             # calculate reference location
#             rs = (rls + rrs)/2
#             # calculate reference chord
#             cs = (c[isurf][j] + c[isurf][j+1])/2

#             # calculate section force and moment coefficients
#             cfj = @SVector zeros(eltype(cf[isurf]), 3)
#             cmj = @SVector zeros(eltype(cm[isurf]), 3)
#             for i = 1:nc
#                 # add influence of bound vortex
#                 rb = top_center(panels[i,j])
#                 cfb = properties[i,j].cfb
#                 cfj += cfb
#                 cmj += cross(rb - rs, cfb)
#                 # add influence of left vortex leg
#                 rl = left_center(panels[i,j])
#                 cfl = properties[i,j].cfl
#                 cfj += cfl
#                 cmj += cross(rl - rs, cfl)
#                 # add influence of right vortex leg
#                 rr = right_center(panels[i,j])
#                 cfr = properties[i,j].cfr
#                 cfj += cfr
#                 cmj += cross(rr - rs, cfr)
#             end

#             # update normalization
#             cfj *= ref.S/(ds*cs)
#             cmj *= ref.S/(ds*cs^2)

#             # change coordinate frame
#             cfj, cmj = body_to_frame(cfj, cmj, ref, fs, frame)
#             # save coefficients
#             cf[isurf][:,j] = cfj
#             cm[isurf][:,j] = cmj
#         end
#     end

#     return cf, cm
# end

# """
#     lifting_line_viscous_coefficients!(cf, cm, system, r, drag_polar; frame=Body(), unsteady=true, additional_velocity=nothing)

# Return the viscous force and moment coefficients (per unit span) for each spanwise segment
# of a lifting line representation of the geometry.

# The inviscid force and moment coefficients `cf` and `cm` must be expressed in
# the [`Body()`](@ref) frame, but will be modified and returned in the frame
# specified by the `frame` argument.

# This function requires that a near-field analysis has been performed on `system`
# to obtain panel forces.

# # Arguments
#  - `cf`: Vector with length equal to the number of surfaces, with each element
#     being a matrix with size (3, ns) which contains the x, y, and z direction
#     inviscid force coefficients (per unit span) for each spanwise segment,
#     expressed in the [`Body()`](@ref) frame.
#  - `cm`: Vector with length equal to the number of surfaces, with each element
#     being a matrix with size (3, ns) which contains the x, y, and z direction
#     inviscid moment coefficients (per unit span) for each spanwise segment,
#     expressed in the [`Body()`](@ref) frame.
#  - `system`: Object of type [`System`](@ref) that holds precalculated
#     system properties.
#  - `r`: Vector with length equal to the number of surfaces, with each element
#     being a matrix with size (3, ns+1) which contains the x, y, and z coordinates
#     of the resulting lifting line coordinates
#  - `drag_polar`: Function defining the drag coefficient as a function of lift
#     coefficient for each spanwise section

# # Keyword Arguments
#  - `frame`: frame in which to return `cf`, `cm`, `cfv`, `cfm`, possible options are
#     [`Body()`](@ref) (default), [`Stability()`](@ref), and [`Wind()`](@ref)`
#  - `unsteady`: Flag indicating if `system` contains the result of an unsteady calculation, and thus
#     should incorperate possible grid motion and the wake shedding locations in the
#     lifting line velocity
#  - `additional_velocity`: Function defining additional velocity field

# # Return Arguments:
#  - `cf`: Vector with length equal to the number of surfaces, with each element
#     being a matrix with size (3, ns) which contains the x, y, and z direction
#     inviscid force coefficients (per unit span) for each spanwise segment,
#     expressed in the reference frame indicated by the `frame` argument.
#  - `cm`: Vector with length equal to the number of surfaces, with each element
#     being a matrix with size (3, ns) which contains the x, y, and z direction
#     inviscid moment coefficients (per unit span) for each spanwise segment,
#     expressed in the reference frame indicated by the `frame` argument.
#  - `cfv`: Vector with length equal to the number of surfaces, with each element
#     being a matrix with size (3, ns) which contains the x, y, and z direction
#     viscous force coefficients (per unit span) for each spanwise segment,
#     expressed in the reference frame indicated by the `frame` argument.
#  - `cmv`: Vector with length equal to the number of surfaces, with each element
#     being a matrix with size (3, ns) which contains the x, y, and z direction
#     viscous moment coefficients (per unit span) for each spanwise segment,
#     expressed in the reference frame indicated by the `frame` argument.
# """
# function lifting_line_viscous_coefficients!(cf, cm, system, r, drag_polar; frame=Body(), unsteady=true, additional_velocity=nothing)
#     TF = promote_type(eltype(eltype(cf)), eltype(eltype(cm)), eltype(system), eltype(eltype(r)))
#     nsurf = length(system.surfaces)
#     cfv = Vector{Matrix{TF}}(undef, nsurf)
#     cmv = Vector{Matrix{TF}}(undef, nsurf)
#     for isurf = 1:nsurf
#         ns = size(system.surfaces[isurf], 2)
#         cfv[isurf] = Matrix{TF}(undef, 3, ns)
#         cmv[isurf] = Matrix{TF}(undef, 3, ns)
#     end
#     return lifting_line_viscous_coefficients!(cfv, cmv, cf, cm, system, r, drag_polar; frame, unsteady, additional_velocity)
# end


# """
#     lifting_line_viscous_coefficients!(cfv, cmv, cf, cm, system, r, drag_polar; frame=Body(), unsteady=true, additional_velocity=nothing)

# In-place version of [`lifting_line_viscous_coefficients`](@ref)
# """
# function lifting_line_viscous_coefficients!(cfv, cmv, cf, cm, system, r, drag_polar; frame=Body(), unsteady=true, additional_velocity=nothing)
#     # number of surfaces
#     nsurf = length(system.surfaces)

#     # check that a near field analysis has been performed
#     @assert system.near_field_analysis[] "Near field analysis required"

#     # extract reference properties
#     ref = system.reference[]
#     fs = system.freestream[]

#     # iterate through each lifting surface
#     for isurf = 1:nsurf
#         nc, ns = size(system.surfaces[isurf])
#         # extract current surface panels and panel properties
#         panels = system.surfaces[isurf]
#         properties = system.properties[isurf]
#         # loop through each chordwise set of panels
#         for j = 1:ns

#             # Get the invsicid lifting line coefficient for this spanwise
#             # station.
#             cfj = SVector{3,eltype(cf[isurf])}(cf[isurf][:, j])
#             cmj = SVector{3,eltype(cm[isurf])}(cm[isurf][:, j])

#             # calculate segment length
#             rls = SVector(r[isurf][1,j], r[isurf][2,j], r[isurf][3,j])
#             rrs = SVector(r[isurf][1,j+1], r[isurf][2,j+1], r[isurf][3,j+1])
#             ds = norm(rrs - rls)
#             # calculate reference location
#             rs = (rls + rrs)/2

#             # calculate the local velocity for the current chordwise set of panels
#             V = lifting_line_velocity(isurf, j, rs, system, unsteady, additional_velocity)

#             # normalized direction along the lifting line's length
#             span_dir = (rrs - rls)/ds

#             # Now I want to know the lift force, which is the force in the
#             # direction that's normal to both the local velocity and span_dir
#             # I don't like this, because I'm not sure which order I should cross
#             # the vectors.
#             lift_dir = cross(span_dir, V./norm(V))
#             # Now, cfj is normalized by what? ds and cs. That's good.
#             # And by what else?
#             # 1/2*RHO*ref.V^2
#             # So I want the 1/2 rho, but I think I need it to be nondimenialized
#             # by the velocity component normal to the span and lift directions.
#             # This isn't the same as norm(V) since that might have a component
#             # along the span direction, and I don't want that.
#             drag_dir = cross(lift_dir, span_dir)
#             Vairfoil = dot(drag_dir, V)
#             # So remove the ref.V^2 and replace it with Vairfoil^2. Now I should
#             # have the force divided by 1/2*RHO*Vairfoil^2*ds*cs
#             cfj_airfoil = cfj * ref.V^2/Vairfoil^2
#             # Now I want the lift coefficient, which is the force in the lift
#             # direction.
#             clj = dot(cfj_airfoil, lift_dir)

#             # Now I want the drag coefficient.
#             cdj = drag_polar(clj)

#             # Update the dimensionalization to match cfj.
#             cdj *= Vairfoil^2/ref.V^2

#             # The force needs to be in the drag direction.
#             cfvj = cdj*drag_dir
#             # What do I do about the moment coefficient?
#             # Nothing, since I'm assuming that the viscous load passes
#             # through the center of the lifting line element, and so the
#             # moment arm is zero.
#             cmvj = zero(cfvj)
            
#             # change coordinate frame
#             cfj, cmj = body_to_frame(cfj, cmj, ref, fs, frame)
#             cfvj, cmvj = body_to_frame(cfvj, cmvj, ref, fs, frame)

#             # save coefficients
#             cf[isurf][:,j] = cfj
#             cm[isurf][:,j] = cmj
#             cfv[isurf][:,j] = cfvj
#             cmv[isurf][:,j] = cmvj
#         end
#     end

#     return cfv, cmv, cf, cm
# end

# """
#     lifting_line_viscous_coefficients!(cfv, cmv, cf, cm, surfaces, properties, reference, freestream, r, drag_polar; frame=Body(), unsteady=true, additional_velocity=nothing)

# In-place version of [`lifting_line_viscous_coefficients`](@ref)
# """
# function lifting_line_viscous_coefficients!(cfv, cmv, cf, cm, surfaces, properties, ref, fs, Vh, r, drag_polar; frame=Body(), unsteady=true, additional_velocity=nothing)
#     # number of surfaces
#     nsurf = length(surfaces)

#     # iterate through each lifting surface
#     for isurf = 1:nsurf
#         nc, ns = size(surfaces[isurf])
#         # extract current surface panels and panel properties
#         panels = surfaces[isurf]
#         properties = properties[isurf]
#         # loop through each chordwise set of panels
#         for j = 1:ns

#             # Get the invsicid lifting line coefficient for this spanwise
#             # station.
#             cfj = SVector{3,eltype(cf[isurf])}(cf[isurf][:, j])
#             cmj = SVector{3,eltype(cm[isurf])}(cm[isurf][:, j])

#             # calculate segment length
#             rls = SVector(r[isurf][1,j], r[isurf][2,j], r[isurf][3,j])
#             rrs = SVector(r[isurf][1,j+1], r[isurf][2,j+1], r[isurf][3,j+1])
#             ds = norm(rrs - rls)
#             # calculate reference location
#             rs = (rls + rrs)/2

#             # calculate the local velocity for the current chordwise set of panels
#             V = lifting_line_velocity(isurf, j, rs, system, unsteady, additional_velocity)

#             # normalized direction along the lifting line's length
#             span_dir = (rrs - rls)/ds

#             # Now I want to know the lift force, which is the force in the
#             # direction that's normal to both the local velocity and span_dir
#             # I don't like this, because I'm not sure which order I should cross
#             # the vectors.
#             lift_dir = cross(span_dir, V./norm(V))
#             # Now, cfj is normalized by what? ds and cs. That's good.
#             # And by what else?
#             # 1/2*RHO*ref.V^2
#             # So I want the 1/2 rho, but I think I need it to be nondimenialized
#             # by the velocity component normal to the span and lift directions.
#             # This isn't the same as norm(V) since that might have a component
#             # along the span direction, and I don't want that.
#             drag_dir = cross(lift_dir, span_dir)
#             Vairfoil = dot(drag_dir, V)
#             # So remove the ref.V^2 and replace it with Vairfoil^2. Now I should
#             # have the force divided by 1/2*RHO*Vairfoil^2*ds*cs
#             cfj_airfoil = cfj * ref.V^2/Vairfoil^2
#             # Now I want the lift coefficient, which is the force in the lift
#             # direction.
#             clj = dot(cfj_airfoil, lift_dir)

#             # Now I want the drag coefficient.
#             cdj = drag_polar(clj)

#             # Update the dimensionalization to match cfj.
#             cdj *= Vairfoil^2/ref.V^2

#             # The force needs to be in the drag direction.
#             cfvj = cdj*drag_dir
#             # What do I do about the moment coefficient?
#             # Nothing, since I'm assuming that the viscous load passes
#             # through the center of the lifting line element, and so the
#             # moment arm is zero.
#             cmvj = zero(cfvj)
            
#             # change coordinate frame
#             cfj, cmj = body_to_frame(cfj, cmj, ref, fs, frame)
#             cfvj, cmvj = body_to_frame(cfvj, cmvj, ref, fs, frame)

#             # save coefficients
#             cf[isurf][:,j] = cfj
#             cm[isurf][:,j] = cmj
#             cfv[isurf][:,j] = cfvj
#             cmv[isurf][:,j] = cmvj
#         end
#     end

#     return cfv, cmv, cf, cm
# end

# function lifting_line_viscous_coefficients_history(cf, cm, system, surface_history, property_history, r, drag_polar; frame=Body(), unsteady=true, additional_velocity=nothing)
#     # unpack system parameters
#     symmetric = system.symmetric
#     ref = system.reference[]
#     fs = system.freestream[]
#     surfaces = system.surfaces

#     # float type
#     TF = eltype(system)

#     # number of time steps
#     nt = length(property_history)

#     # convert single freestream input to vector
#     if isa(fs, Freestream)
#         fs = fill(fs, nt)
#     end

#     # initialize time history coefficients
#     cfv = Vector{Vector{Matrix{TF}}}(undef, nt)
#     cmv = Vector{Vector{Matrix{TF}}}(undef, nt)
#     for it in 1:nt
#         cfv[it], cmv[it] = lifting_line_viscous_coefficients!(cf[it], cm[it], surface_history[it], property_history[it], ref, fs[it], r[it], c[it]; frame=frame)
#     end

#     return cf, cm
# end


"""
    lifting_line_velocity(isurf, is, rc, system, unsteady=true, additional_velocity=nothing)

Return the velocity for a spanwise segment of a lifting line representation
of the geometry associated with surface `isurf` at spanwise segment `is`.

# Arguments
 - `isurf`: Index of the lifting surface this spanwise segment is associated with.
 - `is`: Spanwise index on surface `isurf` this segment is associated with.
 - `rc`: Length-3 vector containing the x, y, and z coordinates of the center of the
    lifting line segment.
 - `system`: Object of type [`System`](@ref) that holds precalculated system properties.
 - `unsteady`: Flag indicating if `system` contains the result of an unsteady
    calculation, and thus should incorperate possible grid motion and the wake
    shedding locations in the lifting line velocity
 - `additional_velocity`: Function defining additional velocity field

# Return Arguments:
 - `V`: Length-3 SVector containing the x, y, and z components of the velocity for the spanwise segment.
"""
@inline function lifting_line_velocity(isurf, is, rc, surfaces, wakes, ref, fs, Γ; additional_velocity, Vh, symmetric, nwake, surface_id, wake_finite_core, wake_shedding_locations, trailing_vortices, xhat)
    # number of surfaces
    nsurf = length(surfaces)

    # size of the current surface
    nc, ns = size(surfaces[isurf])

    # freestream velocity
    V_fs = freestream_velocity(fs)
    V = zero(V_fs)
    V += V_fs

    # rotational velocity
    # V += rotational_velocity(rc, fs, ref)
    V_rot = rotational_velocity(rc, fs, ref)
    V += V_rot

    # additional velocity field
    # if !isnothing(additional_velocity)
    #     V += additional_velocity(rc)
    # end
    if !isnothing(additional_velocity)
        V_add = additional_velocity(rc)
    else
        V_add = zero(V)
    end
    V += V_add

    # velocity due to surface motion by averaging the surface motion at the horizontal bound vorticies for this chordwise set of panels.
    if !isnothing(Vh)
        # V_sm = sum(Vh[isurf][:, is])/nc
        # @show size(Vh[isurf]) nc
        V_sm = sum(Vh[isurf][:, is])/length(Vh[isurf][:, is])
        # for ic in 1:nc
        #     p = top_center(surfaces[isurf][ic, is])
        #     println("ic = $ic, is = $is, r = $(sqrt(p[2]^2 + p[3]^2))")
        # end
        # p = bottom_center(surfaces[isurf][nc, is])
        # println("ic = $(nc+1), is = $is, r = $(sqrt(p[2]^2 + p[3]^2))")
    else
        V_sm = zero(V)
    end
    V += V_sm

    # induced velocity from surfaces and wakes
    jΓ = 0 # index for accessing Γ
    for jsurf = 1:nsurf

        # number of panels on sending surface
        sending = surfaces[jsurf]
        Ns = length(sending)

        # see if wake panels are being used
        wake_panels = nwake[jsurf] > 0

        # check if we need to shift shedding locations
        if isnothing(wake_shedding_locations)
            shedding_locations = nothing
        else
            shedding_locations = wake_shedding_locations[jsurf]
        end

        # extract circulation values corresonding to the sending surface
        vΓ = view(Γ, jΓ+1:jΓ+Ns)

        # OK, so what to do about this?
        # `I` is the CartesianIndex for the recieving panel in
        # `near_field_forces!`.
        # And how does that affect the calculation? Ah, it just changes which
        # vortex elements should be fixed.
        # Don't want that, but I do need to make sure that the finite core model
        # is always used, since the induced velocity calculation would blow up
        # if that wasn't the case and the lifting line location was on one of
        # the vortex elements.
        V += induced_velocity(rc, surfaces[jsurf], vΓ;
            finite_core = true,
            wake_shedding_locations = shedding_locations,
            symmetric = symmetric[jsurf],
            trailing_vortices = trailing_vortices[jsurf] && !wake_panels,
            xhat = xhat)

        # induced velocity from corresponding wake
        if wake_panels
            V += induced_velocity(rc, wakes[jsurf];
                finite_core = wake_finite_core[jsurf] || (surface_id[isurf] != surface_id[jsurf]),
                symmetric = symmetric[jsurf],
                nc = nwake[jsurf],
                trailing_vortices = trailing_vortices[jsurf],
                xhat = xhat)
        end

        jΓ += Ns # increment Γ index for sending panels
    end

    return V, V_fs, V_rot, V_add, V_sm
end

"""
    body_to_frame(CF, CM, reference, freestream, frame)

Transform the coefficients `CF` and `CM` from the body frame to the frame
specified in `frame`
"""
body_to_frame

@inline body_to_frame(CF, CM, ref, fs, ::Body) = CF, CM

@inline function body_to_frame(CF, CM, ref, fs, ::Stability)
    R = body_to_stability(fs)
    return R*CF, R*CM
end

@inline function body_to_frame(CF, CM, ref, fs, ::Wind)
    # remove reference lengths
    reflen = SVector(ref.b, ref.c, ref.b)
    CM = CM .* reflen

    # rotate
    R = body_to_wind(fs)
    CF = R*CF
    CM = R*CM

    # reapply reference lengths
    CM = CM ./ reflen

    return CF, CM
end
