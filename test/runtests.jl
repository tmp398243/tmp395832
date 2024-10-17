using Pkg: Pkg
using ConfigurationsJutulDarcy
using Test
using TestReports
using Aqua
using Documenter

ts = @testset ReportingTestSet "" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(ConfigurationsJutulDarcy; ambiguities=false)
        Aqua.test_ambiguities(ConfigurationsJutulDarcy)
    end

    # Set metadata for doctests.
    DocMeta.setdocmeta!(
        ConfigurationsJutulDarcy,
        :DocTestSetup,
        :(using ConfigurationsJutulDarcy, Test);
        recursive=true,
    )

    # Run doctests.
    doctest(ConfigurationsJutulDarcy; manual=true)
    if ConfigurationsJutulDarcy.HAS_NATIVE_EXTENSIONS
        using JutulDarcy
        doctest(
            ConfigurationsJutulDarcy.get_extension(
                ConfigurationsJutulDarcy, :JutulDarcyExt
            );
            manual=true,
        )
    end

    # Run examples.
    examples_dir = joinpath(@__DIR__, "..", "examples")
    for example in readdir(examples_dir)
        example_path = joinpath(examples_dir, example)
        @show example_path
        orig_project = Base.active_project()
        @testset "Example: $(example)" begin
            if isdir(example_path)
                Pkg.activate(example_path)
                Pkg.develop(; path=joinpath(@__DIR__, ".."))
                Pkg.instantiate()
            end
            script_path = joinpath(example_path, "main.jl")
            try
                include(script_path)
                println("Included script_path")
            finally
                Pkg.activate(orig_project)
            end
        end
    end
end

outputfilename = joinpath(@__DIR__, "..", "report.xml")
open(outputfilename, "w") do fh
    print(fh, report(ts))
end
exit(any_problems(ts))
