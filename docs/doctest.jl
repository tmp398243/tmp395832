
using Documenter: DocMeta, doctest
using JutulModelConfigurations
DocMeta.setdocmeta!(
    JutulModelConfigurations,
    :DocTestSetup,
    :(using JutulModelConfigurations, Test);
    recursive = true,
)
doctest(JutulModelConfigurations)
