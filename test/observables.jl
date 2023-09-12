#=
Verify that full-partial trace and reduced partial trace give the same result
=#
using LadderOperators
using LinearAlgebra
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

# Check that eentropy is the same for tracing out dimer1 or dimer2
@testset "Eentropy" begin
    for dim in [6, 16, 28, 256]
        # Create random, pure density matrix
        x = rand(ComplexF64, dim)
        rho = x * x'
        rho *= 1 / tr(rho)

        rho1 = partial_trace(rho; trace_out = :dimer1)
        rho2 = partial_trace(rho; trace_out = :dimer2)

        # The entanglement entropy depends on the eigenvalues of the
        # density matrices. If this tests true, the eentropies will also be equivalent
        @test eigvals(rho1) ≈ eigvals(rho2)
    end
end
