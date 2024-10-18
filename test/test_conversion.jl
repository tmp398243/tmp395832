
using Test
using Configurations

@testset "Conversion to and from dictionary" begin
    opt_orig = JutulOptions()
    T = typeof(opt_orig)
    opt_conv = from_dict(T, to_dict(opt_orig))
    @test typeof(opt_orig) == typeof(opt_conv)
    for f in fieldnames(T)
        @test getfield(opt_orig, f) == getfield(opt_conv, f)
    end
end
