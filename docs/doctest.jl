
using Documenter: DocMeta, doctest
using ConfigurationsJutulDarcy
DocMeta.setdocmeta!(
    ConfigurationsJutulDarcy,
    :DocTestSetup,
    :(using ConfigurationsJutulDarcy, Test);
    recursive=true,
)
doctest(ConfigurationsJutulDarcy)
