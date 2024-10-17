
# From https://github.com/sintefmath/JutulDarcy.jl/blob/v0.2.35/examples/co2_sloped.jl
using GLMakie
# using HYPRE
using ConfigurationsJutulDarcy
using JutulDarcy
using JutulDarcy.Jutul

Darcy, bar, kg, meter, day, yr = si_units(:darcy, :bar, :kilogram, :meter, :day, :year)
injection_well_trajectory = [
    645.0 0.5 75;    # First point
    660.0 0.5 85;    # Second point
    710.0 0.5 100.0  # Third point
]

options = JutulOptions(;
    mesh=MeshOptions(; n=(100, 1, 50), d=(1e1, 1e0, 1e0)),
    system=CO2BrineOptions(co2_physics=:immiscible, thermal=false),
    porosity=FieldOptions(; value=-1.0),
    permeability=FieldOptions(; value=-1.0),
    temperature=FieldOptions(; value=convert_to_si(30.0, :Celsius)),
    rock_density=FieldOptions(; value=30.0),
    rock_heat_capacity=FieldOptions(; value=900.0),
    rock_thermal_conductivity=FieldOptions(; value=3.0),
    fluid_thermal_conductivity=FieldOptions(; value=0.6),
    component_heat_capacity=FieldOptions(; value=4184.0),
    injection=WellOptions(;
        trajectory=injection_well_trajectory, name=:Injector
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
                    rate_mtons_year=2.05e-5,
                ),
            ],
        ),
        TimeDependentOptions(; years=475.0, steps=475, controls=[]),
    ],
)

# ## Set up a 2D aquifer model
# We set up a Cartesian mesh that is then transformed into an unstructured mesh.
# We can then modify the coordinates to create a domain with a undulating top
# surface. CO2 will flow along the top surface and the topography of the top
# surface has a large impact on where the CO2 migrates.
cart_mesh = CartesianMesh(options.mesh)
mesh = UnstructuredMesh(cart_mesh)

points = mesh.node_points
for (i, pt) in enumerate(points)
    x, y, z = pt
    x_u = 2 * π * x / 1000.0
    w = 0.2
    dz =
        0.05 * x +
        0.05 * abs(x - 500.0) +
        w * (
            30 * sin(2.0 * x_u) +
            20 * sin(5.0 * x_u)
        )
    points[i] = pt + [0, 0, dz]
end

# ## Define permeability and porosity
# We loop over all cells and define three layered regions by the K index of each
# cell. We can then set a corresponding diagonal permeability tensor (3 values)
# and porosity (scalar) to introduce variation between the layers.

nc = number_of_cells(mesh)
perm = zeros(3, nc)
poro = fill(0.3, nc)
region = zeros(Int, nc)
for cell in 1:nc
    I, J, K = cell_ijk(mesh, cell)
    if K < 0.3*options.mesh.n[3]
        reg = 1
        permxy = 0.3*Darcy
        phi = 0.2
    elseif K < 0.7*options.mesh.n[3]
        reg = 2
        permxy = 1.2*Darcy
        phi = 0.35
    else
        reg = 3
        permxy = 0.1*Darcy
        phi = 0.1
    end
    permz = 0.5*permxy
    perm[1, cell] = perm[2, cell] = permxy
    perm[3, cell] = permz
    poro[cell] = phi
    region[cell] = reg
end

plot_cell_data(mesh, poro)
# ## Set up simulation model
# We set up a domain and a single injector. We pass the special :co2brine
# argument in place of the system to the reservoir model setup routine. This
# will automatically set up a compositional two-component CO2-H2O model with the
# appropriate functions for density, viscosity and miscibility.
#
# Note that this model by default is isothermal, but we still need to specify a
# temperature when setting up the model. This is because the properties of CO2
# strongly depend on temperature, even when thermal transport is not solved.
domain = reservoir_domain(mesh, options; permeability = perm, porosity = poro)


# ## Plot cells intersected by the deviated injector well
# We place a single injector well. This well was unfortunately not drilled
# completely straight, so we cannot directly use `add_vertical_well` based on
# logical indices. We instead define a matrix with three columns x, y, z that
# lie on the well trajectory and use utilities from `Jutul` to find the cells
# intersected by the trajectory.
import JutulDarcy.Jutul: plot_mesh_edges

Injector = setup_well(domain, options.injection)

wc = Injector.perforations.reservoir
fig, ax, plt = plot_mesh_edges(mesh, z_is_depth = true)
plot_mesh!(ax, mesh; cells = wc, transparency = true, alpha = 0.4)
# View from the side
ax.azimuth[] = 1.5*π
ax.elevation[] = 0.0
lines!(ax, options.injection.trajectory', color = :red)
display(fig)
fig
# Make model
model = setup_reservoir_model(domain, options.system; wells=Injector, extra_out=false)

# ## Customize model by adding relative permeability with hysteresis
# We define three relative permeability functions: kro(so) for the brine/liquid
# phase and krg(g) for both drainage and imbibition. Here we limit the
# hysteresis to only the non-wetting gas phase, but either combination of
# wetting or non-wetting hysteresis is supported.
#
# Note that we import a few utilities from JutulDarcy that are not exported by
# default since hysteresis falls under advanced functionality.
import JutulDarcy: table_to_relperm, add_relperm_parameters!, brooks_corey_relperm
so = range(0, 1, 10)
krog_t = so.^2
krog = PhaseRelativePermeability(so, krog_t, label = :og)

# Higher resolution for second table
sg = range(0, 1, 50)

# Evaluate Brooks-Corey to generate tables
tab_krg_drain = brooks_corey_relperm.(sg, n = 2, residual = 0.1)
tab_krg_imb = brooks_corey_relperm.(sg, n = 3, residual = 0.25)

krg_drain  = PhaseRelativePermeability(sg, tab_krg_drain, label = :g)
krg_imb  = PhaseRelativePermeability(sg, tab_krg_imb, label = :g)

fig, ax, plt = lines(sg, tab_krg_drain, label = "krg drainage")
lines!(ax, sg, tab_krg_imb, label = "krg imbibition")
lines!(ax, 1 .- so, krog_t, label = "kro")
axislegend()
display(fig)
fig
## Define a relative permeability variable
# JutulDarcy uses type instances to define how different variables inside the
# simulation are evaluated. The `ReservoirRelativePermeabilities` type has
# support for up to three phases with w, ow, og and g relative permeabilities
# specified as a function of their respective phases. It also supports
# saturation regions.
#
# Note: If regions are used, all drainage curves come first followed by equal
# number of imbibition curves. Since we only have a single (implicit) saturation
# region, the krg input should have two entries: One for drainage, and one for
# imbibition.
#
# We also call `add_relperm_parameters` to the model. This makes sure that when
# hysteresis is enabled, we track maximum saturation for hysteresis in each
# reservoir cell.
import JutulDarcy: KilloughHysteresis, ReservoirRelativePermeabilities
krg = (krg_drain, krg_imb) 
H_g = KilloughHysteresis() # Other options: CarlsonHysteresis, JargonHysteresis
relperm = ReservoirRelativePermeabilities(g = krg, og = krog, hysteresis_g = H_g)
replace_variables!(model, RelativePermeabilities = relperm)
add_relperm_parameters!(model);
# ## Define approximate hydrostatic pressure and set up initial state
# The initial pressure of the water-filled domain is assumed to be at
# hydrostatic equilibrium. If we use an immiscible model, we must provide the
# initial saturations. If we are using a compositional model, we should instead
# provide the overall mole fractions. Note that since both are fractions, and
# the CO2 model has correspondence between phase ordering and component ordering
# (i.e. solves for liquid and vapor, and H2O and CO2), we can use the same input
# value.
nc = number_of_cells(mesh)
p0 = zeros(nc)
depth = domain[:cell_centroids][3, :]
g = Jutul.gravity_constant
@. p0 = 200bar + depth*g*1000.0
state0 = setup_reservoir_state(model, options.system; Pressure=p0)
parameters = setup_parameters(model);

# ## Find the boundary and apply a constant pressureboundary condition
# We find cells on the left and right boundary of the model and set a constant
# pressure boundary condition to represent a bounding aquifer that retains the
# initial pressure far away from injection.

boundary = Int[]
for cell in 1:number_of_cells(mesh)
    I, J, K = cell_ijk(mesh, cell)
    if I == 1 || I == options.mesh.n[1]
        push!(boundary, cell)
    end
end
bc = flow_boundary_condition(boundary, domain, p0[boundary], fractional_flow = [1.0, 0.0])
println("Boundary condition added to $(length(bc)) cells.")
## Plot the model
fig = plot_reservoir(model)
display(fig)
fig
# ## Set up schedule
# We set up 25 years of injection and 475 years of migration where the well is
# shut. The density of the injector is set to 900 kg/m^3, which is roughly the
# density of CO2 at these high-pressure in-situ conditions.
dt, forces = setup_reservoir_forces(model, options.time; bc)
nstep = options.time[1].steps
nstep_shut = options.time[2].steps
println("$nstep report steps with injection, $nstep_shut report steps with migration.")
# ## Add some more outputs for plotting
rmodel = reservoir_model(model)
push!(rmodel.output_variables, :RelativePermeabilities)
push!(rmodel.output_variables, :PhaseViscosities)
# ## Simulate the schedule
# We set a maximum internal time-step of 90 days to ensure smooth convergence
# and reduce numerical diffusion.
wd, states, t = simulate_reservoir(
    state0, model, dt; parameters=parameters, forces=forces, max_timestep=90day
)
# ## Plot the CO2 mole fraction
# We plot the overall CO2 mole fraction. We scale the color range to account for
# the fact that the mole fraction in cells made up of only the aqueous phase is
# much smaller than that of cells with only the gaseous phase, where there is
# almost just CO2.
function plot_co2!(fig, ix, x, title="")
    ax = Axis3(
        fig[ix, 1];
        zreversed=true,
        azimuth=-0.51π,
        elevation=0.05,
        aspect=(1.0, 1.0, 0.3),
        title=title,
    )
    plt = plot_cell_data!(ax, mesh, x, colormap = :seaborn_icefire_gradient, colorrange = (0.0, 0.1))
    return Colorbar(fig[ix, 2], plt)
end
fig = Figure(; size=(900, 1200))
for (i, step) in enumerate([1, 5, nstep, nstep + nstep_shut])
    if options.system.co2_physics == :immiscible
        plot_co2!(fig, i, states[step][:Saturations][2, :], "CO2 plume saturation at report step $step/$(nstep+nstep_shut)")
    else
        plot_co2!(fig, i, states[step][:OverallMoleFractions][2, :], "CO2 mole fraction at report step $step/$(nstep+nstep_shut)")
    end
end
display(fig)
fig
# ## Plot all relative permeabilities for all time-steps
# We can plot all relative permeability evaluations. This both verifies that the
# hysteresis model is active, but also gives an indication to how many cells are
# exhibiting imbibition during the simulation.
kro_val = Float64[]
krg_val = Float64[]
sg_val = Float64[]
for state in states
    kr_state = state[:RelativePermeabilities]
    s_state = state[:Saturations]
    for c in 1:nc
        push!(kro_val, kr_state[1, c])
        push!(krg_val, kr_state[2, c])
        push!(sg_val, s_state[2, c])
    end
end

fig = Figure()
ax = Axis(fig[1, 1], title = "Relative permeability during simulation")
fig, ax, plt = scatter(sg_val, kro_val, label = "kro", alpha = 0.3)
scatter!(ax, sg_val, krg_val, label = "krg", alpha = 0.3)
axislegend()
display(fig)
fig
# ## Plot result in interactive viewer
# If you have interactive plotting available, you can explore the results
# yourself.
plot_reservoir(model, states)
## Calculate and display inventory of CO2
# We can classify and plot the status of the CO2 in the reservoir. We use a
# fairly standard classification where CO2 is divided into:
#
# - dissolved CO2 (dissolution trapping)
# - residual CO2 (immobile due to zero relative permeability, residual trapping)
# - mobile CO2 (mobile but still inside domain)
# - outside domain (left the simulation model and migrated outside model)
#
# We also note that some of the mobile CO2 could be considered to be
# structurally trapped, but this is not classified in our inventory.

inventory = co2_inventory(model, wd, states, t)
fig = JutulDarcy.plot_co2_inventory(t, inventory)
display(fig)
fig
# ## Pick a region to investigate the CO2
# We can also specify a region to the CO2 inventory. This will introduce
# additional categories to distinguish between outside and inside the region of
# interest.
cells = findall(region .== 2)
inventory = co2_inventory(model, wd, states, t, cells = cells)
fig = JutulDarcy.plot_co2_inventory(t, inventory)
display(fig)
fig

# ## Define a region of interest using geometry
# Another alternative to determine a region of interest is to use geometry. We
# pick all cells within an ellipsoid a bit away from the injection point.
is_inside = fill(false, nc)
centers = domain[:cell_centroids]
for cell in 1:nc
    x, y, z = centers[:, cell]
    is_inside[cell] = sqrt((x - 720.0)^2 + 20*(z-70.0)^2) < 75
end
plot_cell_data(mesh, is_inside)
# ## Plot inventory in ellipsoid
# Note that a small mobile dip can be seen when free CO2 passes through this region.
inventory = co2_inventory(model, wd, states, t, cells = findall(is_inside))
fig = JutulDarcy.plot_co2_inventory(t, inventory)
display(fig)
fig
# ## Plot the average pressure in the ellipsoid region
# Now that we know what cells are within the region of interest, we can easily
# apply a function over all time-steps to figure out what the average pressure
# value was.
using Statistics
p_avg = map(
    state -> mean(state[:Pressure][is_inside])./bar,
    states
)
fig = lines(t./yr, p_avg,
    axis = (
        title = "Average pressure in region",
        xlabel = "Years", ylabel = "Pressure (bar)"
    )
)
display(fig)
fig
# ## Make a composite plot to correlate CO2 mass in region with spatial distribution
# We create a pair of plots that combine both 2D and 3D plots to simultaneously
# show the ellipsoid, the mass of CO2 in that region for a specific step, and
# the time series of the CO2 in the same region.

stepno = clamp(100, nstep, length(states)) 
co2_mass_in_region = map(
    state -> sum(state[:TotalMasses][2, is_inside])/1e3,
    states
)
fig = Figure(size = (1200, 600))
ax1 = Axis(fig[1, 1],
    title = "Mass of CO2 in region",
    xlabel = "Years",
    ylabel = "Tonnes CO2"
)
lines!(ax1, t./yr, co2_mass_in_region)
scatter!(ax1, t[stepno]./yr, co2_mass_in_region[stepno], markersize = 12, color = :red)
ax2 = Axis3(fig[1, 2], zreversed = true)
plot_cell_data!(ax2, mesh, states[stepno][:TotalMasses][2, :])
plot_mesh!(ax2, mesh, cells = findall(is_inside), alpha = 0.5)
ax2.azimuth[] = 1.5*π
ax2.elevation[] = 0.0
display(fig)
fig
