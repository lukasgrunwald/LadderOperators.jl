#=
Commutator tests for the local ladder-operatos. If the commutators are correct, the representation
of the operators must be correct!
=#
using LadderOperators
using LinearAlgebra
using Test

# Convenience functions for checking the functions
commutator(A, B) = A*B - B * A
anti_commutator(A, B) = A*B + B * A

# ——————————————————————————————— Test spinless Fermions ——————————————————————————————— #
Cdag, C, Z, N = local_spinless_fermions()

@testset "local_spinless" begin
    @test N == [0 0; 0 1]
    @test anti_commutator(Cdag, Cdag) == 0I
    @test anti_commutator(Cdag, Cdag) == 0I
    @test anti_commutator(C, Cdag) == I
    @test Z^2 == I
end

# ———————————————————————————————————— Test Fermions ——————————————————————————————————— #
Cdag, C, Z, N = local_fermions()

@testset "local_fermion" begin
    for i in 1:2
        for j in 1:2
            @test anti_commutator(Cdag[i], Cdag[j]) == 0I
            @test anti_commutator(C[i], C[j]) == 0I

            @test anti_commutator(Cdag[i], C[j]) == (i == j ? I : 0I)
        end
    end
    @test Z^2 == I
end

# —————————————————————————————————— Test multiorbital ————————————————————————————————— #
Nb = 2
Cdag, C, Z, N = local_fermions_multiorbital(Nb)

@testset "fermion multiorbital" begin
    for b1 in 1:Nb, b2 in 1:Nb
        for i in 1:2, j in 1:2
            @test anti_commutator(Cdag[b1, i], Cdag[b2, j]) == 0I
            @test anti_commutator(C[b1, i], C[b2, j]) == 0I
            @test anti_commutator(Cdag[b1, i], C[b2, j]) == ((i == j && b1 == b2) ? I : 0I)
        end
    end
    @test Z^2 == I # Jordan-Wigner factor only gives phase
end

# —————————————————————————————————— Bosonic Operators ————————————————————————————————— #
@testset "bosonic" begin
    dim = 10
    idl = (1:dim-1, 1:dim-1)

    a_dagger, a = ladder_operators(:bosons, 1, dim)
    n = Diagonal(0:(dim-1))

    @test a_dagger * a ≈ n
    # Need to remove the last entry from a, a_dagger commutator, since this
    # can never be fullfilled in a finite hilbert space.
    @test (a * a_dagger - a_dagger * a)[idl...] ≈ I(dim-1)
    @test (n * a_dagger - a_dagger * n) ≈ a_dagger
    @test (n * a - a * n) ≈ -a
end
