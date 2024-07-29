module JutulDarcyExt

using JutulModelConfigurations
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
    domain = reservoir_domain(mesh; permeability, porosity, temperature)
end

end # module