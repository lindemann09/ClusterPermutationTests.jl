using ClusterPermutationTests
using Aqua
using Test

Aqua.test_all(ClusterPermutationTests; ambiguities=false, deps_compat=false)

@testset "ClusterPermutationTests.jl" begin
    # Write your tests here.
end
