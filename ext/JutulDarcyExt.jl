module JutulDarcyExt

using ConfigurationsJutulDarcy
using JutulDarcy
using JutulDarcy.Jutul

function Jutul.CartesianMesh(options::MeshOptions)
    CartesianMesh(options.n, options.n .* options.d)
end

function create_field(mesh, options::FieldOptions)
    if options.type == "constant"
        if options.pad_boundary
            field = fill(options.value, mesh.dims)
            if mesh.dims[2] == 1
                field = dropdims(field; dims=2)
            end
            field[1, :] .= options.pad_value
            field[end, :] .= options.pad_value
            field[:, end] .= options.pad_value
            return field
        else
            return options.value
        end
    end
    error("Unknown field type: '$(options.type)'")
end

function JutulDarcy.reservoir_domain(mesh, options::JutulOptions)
    porosity = create_field(mesh, options.porosity)
    temperature = create_field(mesh, options.temperature)
    permeability = create_field(mesh, options.permeability)
    rock_density = create_field(mesh, options.rock_density)
    rock_heat_capacity = create_field(mesh, options.rock_heat_capacity)
    rock_thermal_conductivity = create_field(mesh, options.rock_thermal_conductivity)
    fluid_thermal_conductivity = create_field(mesh, options.fluid_thermal_conductivity)
    component_heat_capacity = create_field(mesh, options.component_heat_capacity)
    domain = reservoir_domain(mesh; permeability, porosity, temperature,
        rock_density,
        rock_heat_capacity,
        rock_thermal_conductivity,
        fluid_thermal_conductivity,
        component_heat_capacity,
    )
end

function JutulDarcy.setup_well(D::DataDomain, options::WellOptions; kwarg...)
    start_loc = [options.loc..., options.search_zrange[1]]
    stop_loc = [options.loc..., options.search_zrange[2]]

    start_dist, start_idx = findmin(x -> sum((x .- start_loc).^2), eachcol(D[:cell_centroids])) 
    stop_dist, stop_idx = findmin(x -> sum((x .- stop_loc).^2), eachcol(D[:cell_centroids])) 
    if start_idx != stop_idx
        error("I didn't figure out how to implement this yet.")
    end
    reservoir_cells = [start_idx]
    return setup_well(D, reservoir_cells; name=Symbol(options.name), kwarg...)
end

function JutulDarcy.InjectorControl(options::WellRateOptions)
    rate_kg_s = options.rate_mtons_year * 1e9 / (365.2425 * 24 * 60 * 60)
    rate_m3_s = rate_kg_s / options.fluid_density
    rate_target = TotalRateTarget(rate_m3_s)
    I_ctrl = InjectorControl(rate_target, [0.0, 1.0], density = options.fluid_density)
end

function setup_control(options)
    if options.type == "injector"
        return InjectorControl(options)
    end
    error("Not implemented for type: $type")
end

function JutulDarcy.setup_reservoir_forces(model, options::Vector{<:TimeDependentOptions})
    nsteps = sum(x->x.steps, options)
    dt = fill(0.0, nsteps)
    forces = Vector{Any}(undef, nsteps)

    start_step = 1
    for opt in options
        stop_step = start_step + opt.steps - 1
        control = Dict(Symbol(o.name) => setup_control(o) for o in opt.controls)
        forces_opt = setup_reservoir_forces(model; control)
        dt[start_step:stop_step] .= (365.2425 * 24 * 60 * 60) * opt.years / opt.steps
        fill!(view(forces, start_step:stop_step), forces_opt)
        start_step = stop_step + 1
    end
    return dt, forces
end

end # module