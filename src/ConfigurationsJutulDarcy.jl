module ConfigurationsJutulDarcy

include("structs.jl")

using PackageExtensionCompat
function __init__()
    @require_extensions
end

end # module ConfigurationsJutulDarcy
