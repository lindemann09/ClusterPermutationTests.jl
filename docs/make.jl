push!(LOAD_PATH,"../src/")

using Documenter
using ClusterPermutationTests

makedocs(
    sitename = "ClusterPermutationTests.jl",
    doctest = false,
    authors = "Oliver Lindemann",
    pages = [
        "index.md",
    ],
)

#deploydocs(; repo = "github.com/lindemann09/ClusterPermutationTests.jl", push_preview = true)