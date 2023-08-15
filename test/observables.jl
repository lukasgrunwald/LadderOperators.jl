#=
Verify that full-partial trace and reduced partial trace give the same result
=#
using LadderOperators
using Test

# ————————————————————————————————— Partial-Trace tests ———————————————————————————————— #
@testset "Partial Trace" begin
    for dim in [6, 28]
        x = rand(ComplexF64, dim, dim)
        x_full = map_full_fock(x)

        pr_red = partial_trace(x)
        pr_full = partial_trace(x_full)

        @test pr_red == pr_full
    end
end
