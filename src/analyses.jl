"""
    steady_analysis(grids, reference, freestream; kwargs...)
    steady_analysis(surfaces, lifting_lines, reference, freestream; kwargs...)

Perform a steady vortex lattice method analysis.  Return an object of type
[`System`](@ref) containing the system state.

# Arguments
 - `grids`: Vector of grids of shape (3, nc+1, ns+1) which represent lifting surfaces
 - `surfaces`:
   - Vector of matrices of shape (nc, ns) containing surface panels (see
   [`SurfacePanel`](@ref))
   where `nc` is the number of chordwise panels and `ns` is the number of
   spanwise panels
 - `lifting_lines`: Vector of length-`ns` vectors of
    [`LiftingLineSegment`](@ref), containing a lifting line representation of each
    surface
 - `reference`: Reference parameters (see [`Reference`](@ref))
 - `freestream`: Freestream parameters (see [`Freestream`](@ref))

# Keyword Arguments
 - `symmetric`: Flag for each surface indicating whether a mirror image across
    the X-Z plane should be used when calculating induced velocities. Defaults to
    `false` for each surface
 - `wakes`: Matrix of wake panels (see [`WakePanel`](@ref)) for each surface.  Each
    matrix has shape (nw, ns) where `nw` is the number of chordwise wake panels
    and `ns` is the number of spanwise panels for each surface, defaults to no
    wake panels for each surface
 - `nwake`: Number of chordwise wake panels to use from each wake in `wakes`,
    defaults to all wake panels for each surface
 - `surface_id`: Surface ID for each surface.  The finite core model is disabled
    when calculating the influence of surfaces/wakes that share the same ID.
    By default, all surfaces are assigned their own IDs
 - `wake_finite_core`: Flag for each wake indicating whether the finite core
    model should be enabled when calculating a wake's influence on itself and
    surfaces/wakes with the same surface ID.  Defaults to `true` for each surface.
 - `trailing_vortices`: Flags to enable/disable trailing vortices for each surface,
    defaults to `true` for each surface
 - `xhat`: Direction in which to shed trailing vortices, defaults to `[1, 0, 0]`
 - `additional_velocity`: Function which defines additional velocity as a
    function of location.
 - `fcore`: function which sets the finite core size for each surface based on
    the chord length and/or the panel width. Defaults to `(c, Δs) -> 1e-3`.
    Only used for grid inputs.
 - `calculate_influence_matrix`: Flag indicating whether the aerodynamic influence
    coefficient matrix has already been calculated.  Re-using the same AIC matrix
    will reduce calculation times when the underlying geometry has not changed.
    Defaults to `true`.  Note that this argument is only valid for the pre-allocated
    version of this function.
 - `near_field_analysis`: Flag indicating whether a near field analysis should be
    performed to obtain panel velocities, circulation, and forces. Defaults to `true`.
 - `lifting_line_analysis`: Flag indicating whether force coefficients should be
    calculated on a lifting line representation of the geometry. Defaults to
    `false`.
 - `viscous_lifting_line`: Flag indicating whether a viscous loading model
    should be used when calculating the forces acting on the lifting line. Defaults to `false`.
 - `drag_polar`: Function defining the drag coefficient as a function of lift
    coefficient for each spanwise section. Used by the viscous lifting line model.
 - `clmax`: Maximum lift coefficient for each spanwise section. Used by the viscous lifting line model.
 - `derivatives`: Flag indicating whether the derivatives with respect
    to the freestream variables should be calculated. Defaults to `true`. 
 - `prandtl_glauert`: Flag indicating whether the Prandtl-Glauert compressibility correction should be used. Defaults to `false`.
    to the freestream variables should be calculated. Defaults to `true`.
 - `xc`: Normalized chordwise location of the lifting line from the leading edge.
    Defaults to `0.25`, indicating the lifting line will be constructed along the quarter chord
 - `re_correction`: function of the form `cd_corrected = re_correction(cd, Re)`, where `cd` is a local drag coefficient and `Re` a local Reynolds number based on chord.
"""
steady_analysis

# grid input
function steady_analysis(grids::AbstractVector{<:AbstractArray{<:Any, 3}}, reference, freestream;
        fcore = (c, Δs) -> 1e-3, xc = 0.25, kwargs...)

    # pre-allocate system storage
    system = System(grids)

    # generate surface panels
    # grid_to_surface_panels returns an array of the panel corners *and* the
    # Matrix of SurfacePanels, but we just want the latter, hence the [2].
    surfaces = [grid_to_surface_panels(grid; fcore)[2] for grid in grids]

    # create lifting lines from the grids
    lifting_lines = lifting_line_geometry(grids, xc)

    return steady_analysis!(system, surfaces, lifting_lines, reference, freestream; kwargs...,
        calculate_influence_matrix = true)
end

# surface panel and lifting line input
function steady_analysis(surfaces::AbstractVector{<:AbstractMatrix{<:SurfacePanel}}, lifting_lines, reference, freestream; kwargs...)

    # pre-allocate system storage
    system = System(surfaces)

    return steady_analysis!(system, surfaces, lifting_lines, reference, freestream; kwargs...,
        calculate_influence_matrix = true)
end

"""
    steady_analysis!(system, surfaces, reference, freestream; kwargs...)

Pre-allocated version of `steady_analysis`.
"""
function steady_analysis!(system, surfaces::AbstractVector{<:AbstractMatrix{<:SurfacePanel}}, lifting_lines, ref, fs;
    symmetric = fill(false, length(surfaces)),
    wakes = [Matrix{WakePanel{Float64}}(undef, 0, size(surfaces[i],2)) for i = 1:length(surfaces)],
    nwake = size.(wakes, 1),
    surface_id = 1:length(surfaces),
    wake_finite_core = fill(true, length(surfaces)),
    trailing_vortices = fill(true, length(surfaces)),
    xhat = SVector(1, 0, 0),
    additional_velocity = nothing,
    calculate_influence_matrix = true,
    near_field_analysis = true,
    lifting_line_analysis = false,
    viscous_lifting_line = false,
    drag_polar = nothing,
    drag_alpha = nothing,
    prandtl_glauert = false,
    derivatives = true,
    re_correction = (cd, Re) -> cd,
    cp_offset = 0.0)

    # This probably isn't the "right" way to throw an error.
    if lifting_line_analysis && !near_field_analysis
        @error "lifting line analysis requires nearfield analysis"
    end
    if viscous_lifting_line && !lifting_line_analysis
        @error "viscous lifting line analysis requires lifting_line_analysis"
    end
    if viscous_lifting_line && (drag_polar === nothing)
        @error "viscous lifting line requires a drag polar"
    end
    if viscous_lifting_line && (drag_alpha === nothing)
        @error "viscous lifting line requires a drag alpha polar"
    end

    # number of surfaces
    nsurf = length(surfaces)

    # if only one value is provided, use it for all surfaces
    symmetric = isa(symmetric, Number) ? fill(symmetric, nsurf) : symmetric
    nwake = isa(nwake, Number) ? fill(nwake, nsurf) : nwake
    surface_id = isa(surface_id, Number) ? fill(surface_id, nsurf) : surface_id
    wake_finite_core = isa(wake_finite_core, Number) ? fill(wake_finite_core, nsurf) : wake_finite_core
    trailing_vortices = isa(trailing_vortices, Number) ? fill(trailing_vortices, nsurf) : trailing_vortices

    # store provided surface panels and lifting lines
    for isurf = 1:nsurf
        system.surfaces[isurf] .= surfaces[isurf]
        system.lifting_lines[isurf] .= lifting_lines[isurf]
    end

    # update other parameters stored in `system`
    system.reference[] = ref
    system.freestream[] = fs
    system.symmetric .= symmetric
    system.wakes .= wakes
    system.nwake .= nwake
    system.surface_id .= surface_id
    system.wake_finite_core .= wake_finite_core
    system.trailing_vortices .= trailing_vortices
    system.xhat[] = xhat
    system.prandtl_glauert[] = prandtl_glauert

    # unpack variables stored in `system`
    surfaces = system.surfaces
    lifting_lines = system.lifting_lines
    AIC = system.AIC
    w = system.w
    Γ = system.Γ
    dw = system.dw
    dΓ = system.dΓ
    properties = system.properties
    dproperties = system.dproperties
    lifting_line_properties = system.lifting_line_properties

    # see if wake panels are being used
    wake_panels = nwake .> 0

    # calculate/re-calculate AIC matrix (if necessary)
    if calculate_influence_matrix
        influence_coefficients!(AIC, surfaces;
            symmetric = symmetric,
            surface_id = surface_id,
            trailing_vortices = trailing_vortices .& .!wake_panels,
            xhat = xhat)
    end

    # calculate RHS
    if derivatives
        normal_velocity_derivatives!(w, dw, surfaces, wakes, ref, fs;
            additional_velocity = additional_velocity,
            Vcp = nothing, # no velocity at control points due to surface motion
            symmetric = symmetric,
            surface_id = surface_id,
            nwake = nwake,
            wake_finite_core = wake_finite_core,
            trailing_vortices = trailing_vortices,
            xhat = xhat)
    else
        normal_velocity!(w, surfaces, wakes, ref, fs;
            additional_velocity = additional_velocity,
            Vcp = nothing, # no velocity at control points due to surface motion
            symmetric = symmetric,
            surface_id = surface_id,
            nwake = nwake,
            wake_finite_core = wake_finite_core,
            trailing_vortices = trailing_vortices,
            xhat = xhat)
    end

    # solve for the circulation distribution
    if derivatives
        circulation_derivatives!(Γ, dΓ, AIC, w, dw)
    else
        circulation!(Γ, AIC, w)
    end

    if near_field_analysis
        # perform a near field analysis to obtain panel properties
        if derivatives
            near_field_forces_derivatives!(properties, dproperties, surfaces, wakes,
                ref, fs, Γ, dΓ;
                dΓdt = nothing, # no unsteady forces
                additional_velocity = additional_velocity,
                Vh = nothing, # no velocity due to surface motion
                Vv = nothing, # no velocity due to surface motion
                symmetric = symmetric,
                nwake = nwake,
                surface_id = surface_id,
                wake_finite_core = wake_finite_core,
                wake_shedding_locations = nothing, # shedding location at trailing edge
                trailing_vortices = trailing_vortices,
                xhat = xhat,
                prandtl_glauert = prandtl_glauert)
        else
            near_field_forces!(properties, surfaces, wakes, ref, fs, Γ;
                dΓdt = nothing, # no unsteady forces
                additional_velocity = additional_velocity,
                Vh = nothing, # no velocity due to surface motion
                Vv = nothing, # no velocity due to surface motion
                symmetric = symmetric,
                nwake = nwake,
                surface_id = surface_id,
                wake_finite_core = wake_finite_core,
                wake_shedding_locations = nothing, # shedding location at trailing edge
                trailing_vortices = trailing_vortices,
                xhat = xhat,
                prandtl_glauert = prandtl_glauert)
        end
        if lifting_line_analysis
            # Do the lifting line stuff.
            lifting_line_forces!(lifting_line_properties, lifting_lines, properties, surfaces, ref)
            if viscous_lifting_line
                lifting_line_viscous_forces!(lifting_line_properties, lifting_lines, drag_polar, drag_alpha, clmax, properties, surfaces, wakes,
                    ref, fs, Γ;
                    additional_velocity = additional_velocity,
                    Vh = nothing,
                    symmetric = symmetric,
                    nwake = nwake,
                    surface_id = surface_id,
                    wake_finite_core = wake_finite_core,
                    wake_shedding_locations = wake_shedding_locations,
                    trailing_vorticies = trailing_vorticies,
                    xhat = xhat,
                    re_correction = re_correction,
                    rotation_correction = nothing,
                    cp_offset = cp_offset)
            end
        end
    end

    # save flags indicating whether certain analyses have been performed
    system.near_field_analysis[] = near_field_analysis
    system.lifting_line_analysis[] = lifting_line_analysis
    system.viscous_lifting_line[] = viscous_lifting_line
    system.derivatives[] = derivatives

    # return the modified system
    return system
end

"""
    unsteady_analysis(grids, reference, freestream, dt; kwargs...)
    unsteady_analysis(surfaces, lifting_lines, reference, freestream, dt; kwargs...)

Perform a unsteady vortex lattice method analysis.  Return an object of type
[`System`](@ref) containing the final system state, a matrix of surface panels
(see [`SurfacePanel`](@ref) for each surface at each time step, a matrix of surface
panel properties (see [`PanelProperties`](@ref)) for each surface at each time step,
a matrix of wake panels (see [`WakePanel`](@ref)) for each surface at each time
step, a vector of lifting line geometries (see [`LiftingLineSegment`](@ref)) for
each surface at each time step, and a vector of lifting line properties (see
[`LiftingLineProperties`](@ref)) for each surface at each time step.

# Arguments
 - `grids`: Vector of grids of shape (3, nc+1, ns+1) which represent lifting
    surfaces. Alternatively, a vector containing surface shapes/positions
    at each time step (including at `t=0`) may be provided to model
    moving/deforming lifting surfaces.
 - `surfaces`: Matrices of surface panels (see [`SurfacePanel`](@ref)) of shape
    `(nc, ns)` where `nc` is the number of chordwise panels and `ns` is the number of
    spanwise panels. Alternatively, a vector containing surface shapes/positions
    at each time step (including at `t=0`) may be provided to model
    moving/deforming lifting surfaces.
 - `lifting_lines`: Vector of length-`ns` vectors of lifting line segments (see
    [`LiftingLineSegment`](@ref)), containing a lifting line representation of each
    surface
 - `reference`: Reference parameters (see [`Reference`](@ref))
 - `freestream`: Freestream parameters for each time step (see [`Freestream`](@ref))
 - `dt`: Time step vector

# Keyword Arguments
 - `symmetric`: Flags indicating whether a mirror image (across the X-Z plane) should
    be used when calculating induced velocities, defaults to `false` for each surface
 - `initial_wakes`: Vector of initial wakes corresponding to each surface, represented
    by matrices of wake panels (see [`WakePanel`](@ref)) of shape (nw, ns) where
    `nw` is the number of chordwise wake panels and `ns` is the number of
    spanwise panels. Defaults to no wake panels for each surface
 - `initial_circulation`: Vector containing the initial circulation of all surface
    panels in the system.  Defaults to `zeros(N)` where `N` is the total number
    of surface panels in `surfaces`.
 - `nwake`: Maximum number of wake panels in the chordwise direction for each
    surface.  Defaults to `length(dx)` for all surfaces.
 - `surface_id`: Surface ID for each surface.  The finite core model is disabled
    when calculating the influence of surfaces/wakes that share the same ID.
 - `wake_finite_core`: Flag for each wake indicating whether the finite core
    model should be enabled when calculating the wake's influence on itself and
    surfaces/wakes with the same surface ID.  Defaults to `true` for each surface.
 - `save`: Time indices at which to save the time history, defaults to `1:length(dx)`
 - `calculate_influence_matrix`: Flag indicating whether the aerodynamic influence
    coefficient matrix needs to be calculated.  Re-using the same AIC matrix
    will (slightly) reduce calculation times when the underlying geometry has not changed.
    Defaults to `true`.  Note that this argument only affects the pre-allocated
    version of this function.
 - `near_field_analysis`: Flag indicating whether a near field analysis should be
    performed to obtain panel velocities, circulation, and forces for the final
    time step. Defaults to `true`.
 - `lifting_line_analysis`: Flag indicating whether force coefficients should be
    calculated on a lifting line representation of the geometry. Defaults to
    `false`.
 - `viscous_lifting_line`: Flag indicating whether a viscous loading model
    should be used when calculating the forces acting on the lifting line. Defaults to `false`.
 - `drag_polar`: Function defining the drag coefficient as a function of lift
    coefficient for each spanwise section. Used by the viscous lifting line model.
 - `clmax`: Maximum lift coefficient for each spanwise section. Used by the viscous lifting line model.
 - `derivatives`: Flag indicating whether the derivatives with respect to the
    freestream variables should be calculated for the final time step. Defaults
    to `true`.
 - `prandtl_glauert': Flag indicating whether the Prandtl-Glauert
    compressibility correction should be used. Defaults to `false`.
 - `unsteady_kj': Flag indicating whether the unsteady part of 
    the Kutta-Joukowski theorem should be used for the nearfield loading calculation.
    Defaults to `true`.
 - `xc`: Normalized chordwise location of the lifting line from the leading edge.
    Defaults to `0.25`, indicating the lifting line will be constructed along the quarter chord
 - `re_correction`: function of the form `cd_corrected = re_correction(cd, Re)`, where `cd` is a local drag coefficient and `Re` a local Reynolds number based on chord.
 - `rotation_correction`: `CCBlade.RotationCorrection` object, or `nothing`.
 - `oseen`: `length(system.surfaces)`-Vector of oseen coefficients. Set to zero to disable the finite core length growth model.
 - `a1`: `length(system.surfaces)`-Vector of non-dimensional scaling parameter for circulation's contribution to finite core length growth.
 - `bq_s`: `length(system.surfaces)`-Vector of circulation exponential decay factor. Set to zero for no decay.
 - `cp_offset`: factor displacing the control point for each spanwise station used in the viscous loading calculation.
"""
unsteady_analysis

# same geometry at each time step, grid input
function unsteady_analysis(grids::AbstractVector{<:AbstractArray{T, 3}}, ref, fs, dt;
    nwake = fill(length(dt), length(grids)), fcore = (c, Δs) -> 1e-3, xc=0.25, kwargs...) where {T}

    # pre-allocate system storage
    system = System(grids; nw = nwake)

    # generate surface panels
    # grid_to_surface_panels returns an array of the panel corners *and* the
    # Matrix of SurfacePanels, but we just want the latter, hence the [2].
    surfaces = [grid_to_surface_panels(grid; fcore)[2] for grid in grids]

    # create lifting lines from the grid
    lifting_lines = lifting_line_geometry(grids, xc)

    return unsteady_analysis!(system, surfaces, lifting_lines, ref, fs, dt;
        kwargs..., nwake, calculate_influence_matrix = true)
end

# same geometry at each time step, surface and lifting line input
function unsteady_analysis(surfaces::AbstractVector{<:AbstractMatrix{<:SurfacePanel}}, lifting_lines, ref, fs, dt;
    nwake = fill(length(dt), length(surfaces)), kwargs...)

    # pre-allocate system storage
    system = System(surfaces; nw = nwake)

    return unsteady_analysis!(system, surfaces, lifting_lines, ref, fs, dt;
        kwargs..., nwake, calculate_influence_matrix = true)
end

# # different grids/surfaces at each time step
# function unsteady_analysis(surfaces::AbstractVector{<:AbstractVector{<:AbstractArray}},
#     ref, fs, dt; nwake = fill(length(dt), length(surfaces[1])), kwargs...)

#     # pre-allocate system storage
#     system = System(surfaces[1]; nw = nwake)

#     return unsteady_analysis!(system, surfaces, ref, fs, dt;
#         kwargs..., nwake, calculate_influence_matrix = true)
# end

# different geometry at each time step, grid input
function unsteady_analysis(grids::AbstractVector{<:AbstractVector{<:AbstractArray{T, 3}}},
    ref, fs, dt; nwake = fill(length(dt), length(grids[1])), fcore = (c, Δs) -> 1e-3, xc = 0.25, kwargs...) where {T}

    # pre-allocate system storage
    system = System(grids[1]; nw = nwake)

    # generate surface panels
    # grid_to_surface_panels returns an array of the panel corners *and* the
    # Matrix of SurfacePanels, but we just want the latter, hence the [2].
    surfaces = [[grid_to_surface_panels(grid; fcore)[2] for grid in grids_t] for grids_t in grids]

    # create lifting lines from the grid
    lifting_lines = [lifting_line_geometry(grids_t, xc) for grids_t in grids]

    return unsteady_analysis!(system, surfaces, lifting_lines, ref, fs, dt;
        kwargs..., nwake, calculate_influence_matrix = true)
end

# different geometry at each time step, surface and lifting line input
function unsteady_analysis(surfaces::AbstractVector{<:AbstractVector{<:AbstractMatrix{<:SurfacePanel}}}, lifting_lines,
    ref, fs, dt; nwake = fill(length(dt), length(surfaces[1])), kwargs...) where {T}

    # pre-allocate system storage
    system = System(surfaces[1]; nw = nwake)

    return unsteady_analysis!(system, surfaces, lifting_lines, ref, fs, dt;
        kwargs..., nwake, calculate_influence_matrix = true)
end

"""
    unsteady_analysis!(system, surfaces, reference, freestream, dt; kwargs...)

Pre-allocated version of `unsteady_analysis`.
"""
function unsteady_analysis!(system, surfaces::Union{AbstractVector{<:AbstractMatrix{<:SurfacePanel}}, AbstractVector{<:AbstractVector{<:AbstractMatrix{<:SurfacePanel}}}}, lifting_lines, ref, fs, dt;
    symmetric = fill(false, length(system.surfaces)),
    initial_circulation = zero(system.Γ),
    initial_wakes = [Matrix{WakePanel{Float64}}(undef, 0, size(surfaces[i], 2)) for i = 1:length(system.surfaces)],
    nwake = fill(length(dt), length(system.surfaces)),
    surface_id = 1:length(system.surfaces),
    wake_finite_core = fill(true, length(system.surfaces)),
    additional_velocity = nothing,
    fcore = (c, Δs) -> 1e-3,
    save = 1:length(dt),
    calculate_influence_matrix = true,
    near_field_analysis = true,
    lifting_line_analysis = false,
    viscous_lifting_line = false,
    drag_polar = nothing,
    drag_alpha = nothing,
    clmax = Inf,
    derivatives = true,
    prandtl_glauert = false,
    unsteady_kj=true,
    re_correction = (cd, Re) -> cd,
    rotation_correction = nothing,
    oseen = fill(zero(eltype(system)), length(system.surfaces)),
    a1 = fill(zero(eltype(system)), length(system.surfaces)),
    bq_s = fill(zero(eltype(system)), length(system.surfaces)),
    cp_offset = zero(eltype(system)))

    # This probably isn't the "right" way to throw an error.
    if lifting_line_analysis && !near_field_analysis
        @error "lifting line analysis requires nearfield analysis"
    end
    if viscous_lifting_line && !lifting_line_analysis
        @error "viscous lifting line analysis requires lifting_line_analysis"
    end
    if (drag_polar === nothing) && viscous_lifting_line
        @error "viscous lifting line requires a drag polar"
    end
    if (drag_alpha === nothing) && viscous_lifting_line
        @error "viscous lifting line requires a drag alpha polar"
    end

    # float number type
    TF = eltype(system)

    # number of surfaces
    nsurf = length(system.surfaces)

    # surface motion?
    surface_motion = eltype(surfaces) <: AbstractVector

    # --- Input Pre-Processing --- #

    # convert scalar inputs to appropriately sized vectors
    symmetric = isa(symmetric, Number) ? fill(symmetric, nsurf) : symmetric
    nwake = isa(nwake, Number) ? fill(nwake, nsurf) : nwake
    surface_id = isa(surface_id, Number) ? fill(surface_id, nsurf) : surface_id
    wake_finite_core = isa(wake_finite_core, Number) ? fill(wake_finite_core, nsurf) : wake_finite_core
    fs = isa(fs, Freestream) ? fill(fs, length(dt)) : fs

    # extract initial surface panels from input
    if surface_motion
        # surface moves, store initial surface panel locations
        initial_surfaces = surfaces[1]
        initial_lifting_lines = lifting_lines[1]
    else
        # surface doesn't move
        initial_surfaces = surfaces
        initial_lifting_lines = lifting_lines
    end

    # --- Update System Parameters --- #

    system.reference[] = ref
    system.symmetric .= symmetric
    system.surface_id .= surface_id
    system.wake_finite_core .= wake_finite_core
    system.trailing_vortices .= false
    system.prandtl_glauert[] = prandtl_glauert

    # check if existing wake panel storage is sufficient, replace if necessary
    for isurf = 1:nsurf
        if size(system.wakes[isurf], 1) < nwake[isurf]
            # get surface/wake dimensions
            nc, ns = size(surfaces[isurf])
            nw = nwake[isurf]

            # update wake panel storage
            system.wakes[isurf] = Matrix{WakePanel{TF}}(undef, nwake[isurf], size(surfaces[isurf], 2))
        end
    end

    # find intersecting surfaces
    repeated_points = repeated_trailing_edge_points(initial_surfaces)

    # --- Set Initial Simulation Variables --- #

    # store initial surface panels in `system`
    # if eltype(initial_surfaces) <: AbstractArray{<:Any, 3}
    #     # initial surfaces are input as a grid, convert to surface panels
    #     for isurf = 1:nsurf
    #         update_surface_panels!(system.surfaces[isurf], initial_surfaces[isurf]; fcore)
    #     end
    # else
    #     # initial surfaces are input as matrices of surface panels
    #     for isurf = 1:nsurf
    #         system.surfaces[isurf] .= initial_surfaces[isurf]
    #     end
    # end
    # initial surfaces are *always* input as matrices of surface panels now
    for isurf = 1:nsurf
        system.surfaces[isurf] .= initial_surfaces[isurf]
        system.lifting_lines[isurf] .= initial_lifting_lines[isurf]
    end

    # store initial wake panels in `system`
    for isurf = 1:nsurf
        for I in CartesianIndices(initial_wakes[isurf])
            system.wakes[isurf][I] = initial_wakes[isurf][I]
        end
    end

    # store initial freestream parameters in `system`
    system.freestream[] = fs[1]

    # store initial circulation parameters in `system`
    system.Γ .= initial_circulation

    # set the initial number of wake panels for each surface
    iwake = [min(size(initial_wakes[isurf], 1), nwake[isurf]) for isurf = 1:nsurf]

    # --- Begin Simulation --- #

    # initialize solution history for each time step
    surface_history = Vector{Vector{Matrix{SurfacePanel{TF}}}}(undef, length(save))
    property_history = Vector{Vector{Matrix{PanelProperties{TF}}}}(undef, length(save))
    wake_history = Vector{Vector{Matrix{WakePanel{TF}}}}(undef, length(save))
    lifting_line_history = Vector{Vector{Vector{LiftingLineSegment{TF}}}}(undef, length(save))
    lifting_line_property_history = Vector{Vector{Vector{LiftingLineProperties{TF}}}}(undef, length(save))
    isave = 1

    # # loop through all time steps
    for it = 1 : length(dt)

        first_step = it == 1
        last_step = it == length(dt)

        near_field_analysis_it = it in save || (last_step && near_field_analysis)
        lifting_line_analysis_it = lifting_line_analysis && near_field_analysis_it
        viscous_lifting_line_it = viscous_lifting_line && lifting_line_analysis_it
        if surface_motion
            propagate_system!(system, surfaces[1+it], lifting_lines[1+it], fs[it], dt[it];
                additional_velocity,
                repeated_points,
                nwake = iwake,
                eta = 0.25,
                calculate_influence_matrix = first_step && calculate_influence_matrix,
                near_field_analysis = near_field_analysis_it,
                lifting_line_analysis = lifting_line_analysis_it,
                viscous_lifting_line = viscous_lifting_line_it,
                derivatives = last_step && derivatives,
                prandtl_glauert = prandtl_glauert,
                unsteady_kj = unsteady_kj,
                drag_polar = drag_polar,
                drag_alpha = drag_alpha,
                clmax = clmax,
                re_correction = re_correction,
                rotation_correction = rotation_correction,
                oseen = oseen,
                a1 = a1,
                bq_s = bq_s,
                cp_offset = cp_offset)
        else
            propagate_system!(system, fs[it], dt[it];
                additional_velocity,
                repeated_points,
                nwake = iwake,
                eta = 0.25,
                calculate_influence_matrix = first_step && calculate_influence_matrix,
                near_field_analysis = near_field_analysis_it,
                lifting_line_analysis = lifting_line_analysis_it,
                viscous_lifting_line = viscous_lifting_line_it,
                derivatives = last_step && derivatives,
                prandtl_glauert = prandtl_glauert,
                unsteady_kj = unsteady_kj,
                drag_polar = drag_polar,
                drag_alpha = drag_alpha,
                clmax = clmax,
                re_correction = re_correction,
                rotation_correction = rotation_correction,
                oseen = oseen,
                a1 = a1,
                bq_s = bq_s,
                cp_offset = cp_offset)
        end

        # increment wake panel counter for each surface
        for isurf = 1:nsurf
            if iwake[isurf] < nwake[isurf]
                iwake[isurf] += 1
            end
        end

        # save the surface shape, properties, and resulting shed wake
        if it in save
            surface_history[isave] = [copy(system.surfaces[isurf]) for isurf = 1:nsurf]
            property_history[isave] = [copy(system.properties[isurf]) for isurf = 1:nsurf]
            wake_history[isave] = [system.wakes[isurf][1:iwake[isurf], :] for isurf = 1:nsurf]
            lifting_line_history[isave] = [copy(system.lifting_lines[isurf]) for isurf = 1:nsurf]
            lifting_line_property_history[isave] = [copy(system.lifting_line_properties[isurf]) for isurf = 1:nsurf]
            isave += 1
        end

    end

    # return the modified system and time history
    return system, surface_history, property_history, wake_history, lifting_line_history, lifting_line_property_history
end

"""
    propagate_system!(system, [surfaces, ] freestream, dt; kwargs...)

Propagate the state variables in `system` forward one time step using the
unsteady vortex lattice method system of equations.

# Arguments
 - `system`: Object of type `system` which contains the current system state
 - `surfaces`: Surface locations at the end of this time step. If omitted,
   surfaces are assumed to be stationary.
 - `freestream`: Freestream parameters corresponding to this time step.
 - `dt`: Time increment

# Keyword Arguments
 - `additional_velocity`: Function which defines additional velocity as a
    function of location.
 - `repeated_points`: Dictionary of the form `Dict((isurf, i) => [(jsurf1, j1),
    (jsurf2, j2)...]` which defines repeated trailing edge points.  Trailing edge
    point `i` on surface `isurf` is repeated on surface `jsurf1` at point `j1`,
    `jsurf2` at point `j2`, and so forth. See [`repeated_trailing_edge_points`](@ref)
 - `nwake`: Number of wake panels in the chordwise direction for each surface.
 - `eta`: Time step fraction used to define separation between trailing
    edge and wake shedding location.  Typical values range from 0.2-0.3.
 - `calculate_influence_matrix`: Flag indicating whether the aerodynamic influence
    coefficient matrix needs to be calculated. If argument `surfaces` is provided
    the influence matrix will always be recalculated.
 - `near_field_analysis`: Flag indicating whether a near field analysis should be
    performed to obtain panel velocities, circulation, and forces.
 - `lifting_line_analysis`: Flag indicating whether force coefficients should be
    calculated on a lifting line representation of the geometry. Defaults to
    `false`.
 - `viscous_lifting_line`: Flag indicating whether a viscous loading model
    should be used when calculating the forces acting on the lifting line. Defaults to `false`.
 - `drag_polar`: Function defining the drag coefficient as a function of lift
    coefficient for each spanwise section. Used by the viscous lifting line model.
 - `derivatives`: Flag indicating whether the derivatives with respect to the
    freestream variables should be calculated.
 - `prandtl_glauert`: Flag indicating whether the Prandtl-Glauert
    compressibility correction should be used.
 - `unsteady_kj': Flag indicating whether the unsteady part of 
    the Kutta-Joukowski theorem should be used for the nearfield loading calculation.
 - `re_correction`: function of the form `cd_corrected = re_correction(cd, Re)`, where `cd` is a local drag coefficient and `Re` a local Reynolds number based on chord.
 - `rotation_correction`: `CCBlade.RotationCorrection` object, or `nothing`.
 - `oseen`: `length(wake)`-Vector of oseen coefficients. Set to zero to disable the finite core length growth model.
 - `a1`: `length(wake)`-Vector of non-dimensional scaling parameter for circulation's contribution to finite core length growth.
 - `bq_s`: `length(wake)`-Vector of circulation exponential decay factor. Set to zero for no decay.
 - `cp_offset`: factor displacing the control point for each spanwise station.
"""
propagate_system!

# stationary surfaces
propagate_system!(system, fs, dt; kwargs...) = propagate_system!(system,
    nothing, nothing, fs, dt; kwargs...)

# moving/deforming surfaces
function propagate_system!(system, surfaces, lifting_lines, fs, dt;
    additional_velocity,
    repeated_points,
    nwake,
    eta,
    calculate_influence_matrix,
    near_field_analysis,
    lifting_line_analysis,
    viscous_lifting_line = false,
    drag_polar = nothing,
    drag_alpha = nothing,
    clmax = Inf,
    derivatives,
    prandtl_glauert,
    unsteady_kj,
    re_correction = (cd, Re) -> cd,
    rotation_correction = nothing,
    oseen = fill(zero(eltype(system)), length(surfaces)),
    a1 = fill(zero(eltype(system)), length(surfaces)),
    bq_s = fill(zero(eltype(system)), length(surfaces)),
    cp_offset = zero(eltype(system)))

    # This probably isn't the "right" way to throw an error.
    if lifting_line_analysis && !near_field_analysis
        @error "lifting line analysis requires nearfield analysis"
    end
    if viscous_lifting_line && !lifting_line_analysis
        @error "viscous lifting line analysis requires lifting_line_analysis"
    end
    if viscous_lifting_line && (drag_polar === nothing)
        @error "viscous lifting line requires a drag polar"
    end
    if viscous_lifting_line && (drag_alpha === nothing)
        @error "viscous lifting line requires a drag alpha polar"
    end

    # NOTE: Each step models the transition from `t = t[it]` to `t = [it+1]`
    # (e.g. the first step models from `t = 0` to `t = dt`).  Properties are
    # assumed to be constant during each time step and resulting transient
    # forces/moments are provided at `t[it+1]` (e.g. the properties returned
    # for the first step correspond to `t = dt`).

    nsurf = length(system.surfaces)

    # unpack constant system parameters
    ref = system.reference[]
    symmetric = system.symmetric
    surface_id = system.surface_id
    wake_finite_core = system.wake_finite_core
    trailing_vortices = system.trailing_vortices
    xhat = system.xhat[]

    # unpack system storage (including state variables)
    previous_surfaces = system.previous_surfaces
    previous_lifting_lines = system.previous_lifting_lines
    current_surfaces = system.surfaces
    current_lifting_lines = system.lifting_lines
    properties = system.properties
    lifting_line_properties = system.lifting_line_properties
    dproperties = system.dproperties
    wakes = system.wakes
    wake_velocities = system.V
    wake_shedding_locations = system.wake_shedding_locations
    AIC = system.AIC
    w = system.w
    dw = system.dw
    Γ = system.Γ
    dΓ = system.dΓ
    dΓdt = system.dΓdt
    Vcp = system.Vcp
    Vh = system.Vh
    Vv = system.Vv
    Vte = system.Vte
    prandtl_glauert = system.prandtl_glauert[]

    # check if the surfaces are moving
    surface_motion = !isnothing(surfaces)

    # move geometry and calculate velocities for this time step
    if surface_motion

        # # check if grids are used to represent the new surfaces
        # grid_input = eltype(surfaces) <: AbstractArray{<:Any, 3}

        for isurf = 1:nsurf

            # save current surfaces as previous surfaces
            previous_surfaces[isurf] .= current_surfaces[isurf]
            previous_lifting_lines[isurf] .= current_lifting_lines[isurf]

            # set new surface shape...
            # if grid_input
            #     # ...based on grid inputs
            #     update_surface_panels!(current_surfaces[isurf], surfaces[isurf]; fcore)
            # else
            #     # ...based on surface panels
            #     current_surfaces[isurf] .= surfaces[isurf]
            # end
            # only allowing surface panel input now
            current_surfaces[isurf] .= surfaces[isurf]
            current_lifting_lines[isurf] .= lifting_lines[isurf]

            # calculate surface velocities
            get_surface_velocities!(Vcp[isurf], Vh[isurf], Vv[isurf], Vte[isurf],
                current_surfaces[isurf], previous_surfaces[isurf], dt)
        end
    else
        # zero out surface motion
        for isurf = 1:nsurf
            system.Vcp[isurf] .= Ref(SVector(0.0, 0.0, 0.0))
            system.Vh[isurf] .= Ref(SVector(0.0, 0.0, 0.0))
            system.Vv[isurf] .= Ref(SVector(0.0, 0.0, 0.0))
            system.Vte[isurf] .= Ref(SVector(0.0, 0.0, 0.0))
        end
    end

    # update stored freestream parameters for this time step
    system.freestream[] = fs

    # update number of wake panels for each surface for this time step
    system.nwake .= nwake

    # update the wake shedding location for this time step
    update_wake_shedding_locations!(wakes, wake_shedding_locations,
        current_surfaces, ref, fs, dt, additional_velocity, Vte,
        nwake, eta)

    # calculate/re-calculate AIC matrix (if necessary)
    if surface_motion || calculate_influence_matrix
        influence_coefficients!(AIC, current_surfaces;
            symmetric = symmetric,
            wake_shedding_locations = wake_shedding_locations,
            surface_id = surface_id,
            trailing_vortices = trailing_vortices,
            xhat = xhat)
    end

    # update the AIC matrix to use the new wake shedding locations
    update_trailing_edge_coefficients!(AIC, current_surfaces;
        symmetric = symmetric,
        wake_shedding_locations = wake_shedding_locations,
        trailing_vortices = trailing_vortices)

    # calculate RHS
    if derivatives
        normal_velocity_derivatives!(w, dw, current_surfaces, wakes,
            ref, fs; additional_velocity, Vcp, symmetric, nwake,
            surface_id, wake_finite_core, trailing_vortices, xhat)
    else
        normal_velocity!(w, current_surfaces, wakes, ref, fs;
            additional_velocity, Vcp, symmetric, nwake, surface_id,
            wake_finite_core, trailing_vortices, xhat)
    end

    if unsteady_kj
        # save (negative) previous circulation in dΓdt
        dΓdt .= .-Γ
    end

    # solve for the new circulation
    if derivatives
        circulation_derivatives!(Γ, dΓ, AIC, w, dw)
    else
        circulation!(Γ, AIC, w)
    end

    if unsteady_kj
        # solve for dΓdt using finite difference `dΓdt = (Γ - Γp)/dt`
        dΓdt .+= Γ # add newly computed circulation
        dΓdt ./= dt # divide by corresponding time step
    else
        fill!(dΓdt, 0)
    end

    # compute transient forces on each panel (if necessary)
    if near_field_analysis
        if derivatives
            near_field_forces_derivatives!(properties, dproperties,
                current_surfaces, wakes, ref, fs, Γ, dΓ; dΓdt,
                additional_velocity, Vh, Vv, symmetric, nwake,
                surface_id, wake_finite_core, wake_shedding_locations,
                trailing_vortices, xhat, prandtl_glauert)
        else
            near_field_forces!(properties, current_surfaces, wakes,
                ref, fs, Γ; dΓdt, additional_velocity, Vh, Vv,
                symmetric, nwake, surface_id, wake_finite_core,
                wake_shedding_locations, trailing_vortices, xhat, prandtl_glauert)
        end
        if lifting_line_analysis
            # Do the lifting line stuff.
            lifting_line_forces!(lifting_line_properties, current_lifting_lines, properties, current_surfaces, ref)
            if viscous_lifting_line
                lifting_line_viscous_forces!(lifting_line_properties,
                    current_lifting_lines, drag_polar, drag_alpha, clmax, properties,
                    current_surfaces, wakes, ref, fs, Γ; additional_velocity,
                    Vh, symmetric, nwake, surface_id, wake_finite_core,
                    wake_shedding_locations, trailing_vortices, xhat, re_correction, rotation_correction, cp_offset)
            end
        end

        # save flag indicating that various analyses has been performed
        system.near_field_analysis[] = near_field_analysis
        system.lifting_line_analysis[] = lifting_line_analysis
        system.viscous_lifting_line[] = viscous_lifting_line
    end

    # save flag indicating that derivatives wrt freestream variables have been obtained
    system.derivatives[] = derivatives

    # update wake velocities
    get_wake_velocities!(wake_velocities, current_surfaces,
        wakes, ref, fs, Γ, additional_velocity, Vte, symmetric,
        repeated_points, nwake, surface_id, wake_finite_core,
        wake_shedding_locations, trailing_vortices, xhat)

    # shed additional wake panel (and translate existing wake panels)
    shed_wake!(wakes, wake_shedding_locations, wake_velocities,
        dt, current_surfaces, Γ, nwake, fs.viscosity, oseen, a1, bq_s)

    return system
end
