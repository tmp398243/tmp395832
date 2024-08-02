using ConfigurationsJutulDarcy
using JutulDarcy
using JutulDarcy.Jutul

Darcy, bar, kg, meter, day = si_units(:darcy, :bar, :kilogram, :meter, :day)
options = JutulOptions(;
    mesh=MeshOptions(; n=(10, 1, 5), d=(1e2, 1e0, 1e1)),
    porosity=FieldOptions(; value=0.3),
    permeability=FieldOptions(; value=1.0Darcy),
    temperature=FieldOptions(; value=convert_to_si(30.0, :Celsius)),
    rock_density=FieldOptions(; value=30.0),
    rock_heat_capacity=FieldOptions(; value=900.0),
    rock_thermal_conductivity=FieldOptions(; value=3.0),
    fluid_thermal_conductivity=FieldOptions(; value=0.6),
    component_heat_capacity=FieldOptions(; value=4184.0),
    injection=WellOptions(;
        loc=[645.0, 0.5], search_zrange=[43.5, 43.5], length=37.5, name=:Injector
    ),
    time=[
        TimeDependentOptions(;
            years=25.0,
            steps=25,
            controls=[
                WellRateOptions(;
                    type="injector",
                    name=:Injector,
                    fluid_density=9e2,
                    rate_mtons_year=8.14e-4,
                ),
            ],
        ),
        TimeDependentOptions(; years=25.0, steps=25, controls=[]),
    ],
)

# ## Set up a 2D aquifer model
# We set up a Cartesian mesh that is then transformed into an unstructured mesh.
# We can then modify the coordinates to create a domain with a undulating top
# surface. CO2 will flow along the top surface and the topography of the top
# surface has a large impact on where the CO2 migrates.
mesh = UnstructuredMesh(CartesianMesh(options.mesh))

points = mesh.node_points
for (i, pt) in enumerate(points)
    x, y, z = pt
    x_u = 2 * π * x / 1000.0
    w = 0.2
    dz =
        0.05 * x +
        w * (
            30 * sin(2.0 * x_u) +
            20 * sin(5.0 * x_u) +
            10 * sin(10.0 * x_u) +
            5 * sin(25.0 * x_u)
        )
    points[i] = pt + [0, 0, dz]
end
# ## Set up simulation model
# We set up a domain and a single injector. We pass the special :co2brine
# argument in place of the system to the reservoir model setup routine. This
# will automatically set up a compositional two-component CO2-H2O model with the
# appropriate functions for density, viscosity and miscibility.
#
# Note that this model can be run with a thermal mode by setting
domain = reservoir_domain(mesh, options)
Injector = setup_well(domain, options.injection)
model, parameters = setup_reservoir_model(domain, :co2brine; wells=Injector);
# ## Find the boundary and set increased volume
# We find the left and right boundary of the model and increase the volume of
# those cells. This mimicks a constant pressure boundary condition.
boundary = Int[]
for cell in 1:number_of_cells(mesh)
    I, J, K = cell_ijk(mesh, cell)
    if I == 1 || I == options.mesh.n[1]
        push!(boundary, cell)
    end
end
parameters[:Reservoir][:FluidVolume][boundary] *= 1000;
## Plot the model
using GLMakie
plot_reservoir(model)
# ## Set up schedule
# We set up 25 years of injection and 25 years of migration where the well is
# shut. The density of the injector is set to 900 kg/m^3, which is roughly
# the density of CO2 at in-situ conditions.
dt, forces = setup_reservoir_forces(model, options.time)

# ## Set up initial state
state0 = setup_reservoir_state(model; Pressure=200bar, OverallMoleFractions=[1.0, 0.0])
# ## Simulate the schedule
# We set a maximum internal time-step of 30 days to ensure smooth convergence
# and reduce numerical diffusion.
wd, states, t = simulate_reservoir(
    state0, model, dt; parameters=parameters, forces=forces, max_timestep=30day
)
# ## Plot the density of brine
## The density of brine depends on the CO2 concentration and gives a good
## visualization of where the mass of CO2 exists.
function plot_co2!(fig, ix, x, title="")
    ax = Axis3(
        fig[ix, 1];
        zreversed=true,
        azimuth=-0.51π,
        elevation=0.05,
        aspect=(1.0, 1.0, 0.3),
        title=title,
    )
    plt = plot_cell_data!(ax, mesh, x; colormap=:seaborn_icefire_gradient)
    return Colorbar(fig[ix, 2], plt)
end
fig = Figure(; size=(900, 1200))
nstep = options.time[1].steps
nstep_shut = options.time[2].steps
for (i, step) in enumerate([1, 5, nstep, nstep + nstep_shut])
    plot_co2!(
        fig,
        i,
        states[step][:PhaseMassDensities][1, :],
        "Brine density report step $step/$(nstep+nstep_shut)",
    )
end
fig
## Plot result in interactive viewer
plot_reservoir(model, states)
