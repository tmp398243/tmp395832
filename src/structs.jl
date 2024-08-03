using Configurations

export JutulOptions, MeshOptions, SystemOptions
export FieldOptions, FluidOptions
export WellOptions, WellRateOptions, WellPressureOptions
export TimeDependentOptions

@option struct JutulOptions
    mesh = MeshOptions()

    "number of time steps stored in one file"
    nt = 25

    system = SystemOptions(; label=:co2brine)

    "time interval between 2 adjacent time steps (in days)"
    dt = 73.0485

    "number of files, each of which has nt timesteps."
    nbatches = 1

    "ratio of vertical permeability over horizontal permeability."
    kv_over_kh = 0.36

    sat0_radius_cells = 4
    sat0_range = [0.2, 0.8]

    fluid1 = FluidOptions(;
        name="H₂O", viscosity=1e-3, density=1.053e3, compressibility=3.6563071e-10
    )

    fluid2 = FluidOptions(;
        name="CO₂", viscosity=1e-4, density=7.766e2, compressibility=8e-9
    )

    "m/s^2"
    g = 9.81

    "Pascals"
    reference_pressure = 1.5e7

    porosity = FieldOptions(; value=0.3)
    permeability = FieldOptions(; value=1e-12)
    temperature = FieldOptions(; value=3e2)
    rock_density = FieldOptions(; value=30.0)
    rock_heat_capacity = FieldOptions(; value=900.0)
    rock_thermal_conductivity = FieldOptions(; value=3.0)
    fluid_thermal_conductivity = FieldOptions(; value=0.6)
    component_heat_capacity = FieldOptions(; value=4184.0)

    injection = WellOptions(;
        loc=[1875.0, 50.0], search_zrange=[1693.75, 1812.5], length=37.5, name=:Injector
    )

    production = WellOptions(;
        active=false,
        loc=[2875.0, 50.0],
        search_zrange=[1693.75, 1756.25],
        length=37.5,
        name=:Producer,
    )

    time = [
        TimeDependentOptions(;
            years=1.0,
            steps=10,
            controls=[
                WellRateOptions(;
                    type="injector", name=:Injector, fluid_density=5e2, rate_mtons_year=1e-3
                ),
            ],
        ),
        TimeDependentOptions(; years=1.0, steps=10, controls=[]),
    ]
end

@option struct MeshOptions
    "Grid dimensions"
    n = (30, 1, 20)

    "Grid cell size"
    d = (12.5, 100.0, 6.25)
end

@option struct SystemOptions
    label::Symbol
end

@option struct FluidOptions
    "Identifier for viewing options"
    name::Any

    "Pascal seconds (decapoise) Reference: https://github.com/lidongzh/FwiFlow.jl"
    viscosity::Any

    "kg/m^3"
    density::Any

    "1 / Pascals (should be reciprocal of bulk modulus)"
    compressibility::Any
end

@option struct FieldOptions
    type = "constant"
    value::Any
    pad_boundary = false
    pad_value = 0.0
end

@option struct WellOptions
    active = true
    loc::Any
    search_zrange::Any

    "meters"
    length::Any

    name::Symbol
end

@option struct WellRateOptions
    type
    name
    "kg/m^3"
    fluid_density::Any
    rate_mtons_year::Any
end

@option struct WellPressureOptions
    "Pa"
    bottom_hole_pressure_target::Any
end

@option struct TimeDependentOptions
    years::Any
    steps::Any
    controls::Any
end
