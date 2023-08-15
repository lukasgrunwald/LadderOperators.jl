#=
Test embedding of operators into the lattice by verifying the commutators.
=#
using LadderOperators
using Test

# —————————————————————————————————— Spinless Fermions ————————————————————————————————— #
Nsites = 3
Cdag, C, Z, N = local_spinless_fermions()

Cdag = embed_fock(Nsites, Cdag, Z)
C = embed_fock(Nsites, C, Z)
N = embed_fock(Nsites, N) # directly translate to fock space

@testset "fock_local_spin" begin
    for i in 1:Nsites, j in 1:Nsites
        @test anti_commutator(Cdag[i], Cdag[j]) == 0I
        @test anti_commutator(Cdag[i], Cdag[j]) == 0I
        @test anti_commutator(C[i], Cdag[j]) == (i == j ? I : 0I)
        @test N == Cdag .* C
    end
end

# —————————————————————————————————— Spinful Fermions —————————————————————————————————— #
Nsites = 2
Cdag, C, Z, N = local_fermions()

Cdag = embed_fock(Nsites, Cdag, Z)
C = embed_fock(Nsites, C, Z)
N = embed_fock(Nsites, N) # directly translate to fock space

@testset "fock_local_fermions" begin
for x1 in 1:Nsites, x2 in 1:Nsites
    for i in 1:2, j in 1:2
        @test anti_commutator(Cdag[x1, i], Cdag[x2, j]) == 0I
        @test anti_commutator(Cdag[x1, i], Cdag[x2, j]) == 0I
        @test anti_commutator(C[x1, i], Cdag[x2, j]) == ((i == j && x1 == x2) ? I : 0I)
    end
end
@test N == Cdag .* C
end

# ———————————————————————————————— Multiorbital Fermions ——————————————————————————————— #
Nb = 2
Nsites = 2
Cdag, C, Z, N = local_fermions_multiorbital(Nb)

Cdag = embed_fock(Nsites, Cdag, Z)
C = embed_fock(Nsites, C, Z)
N = embed_fock(Nsites, N) # directly translate to fock space

@testset "fock_fermion_multiorbital" begin
for x1 in 1:Nsites, x2 in 1:Nsites
    for b1 in 1:Nb, b2 in 1:Nb
        for i in 1:2, j in 1:2
            @test anti_commutator(Cdag[x1, b1, i], Cdag[x2, b2, j]) == 0I
            @test anti_commutator(C[x1, b1, i], C[x2, b2, j]) == 0I
            @test anti_commutator(Cdag[x1, b1, i], C[x2, b2, j]) == ((i == j && b1 == b2 && x1 == x2) ? I : 0I)
        end
    end
end
@test N == Cdag .* C
end

# ——————————————————————————— Move to and from U(1) subspace ——————————————————————————— #
Base.:*(x::Integer, word::String) =  x == 0 ? "" : word
Base.:+(st1::String, st2::String) = st1 * st2
Base.zero(st::String) = ""

@testset "U(1)→Full Fock space" begin
    for x in [:fermi, :fermi_multi]
        Cdag, C = ladder_operators(x, 2)
        N = Cdag .* C
        idx = idx_u1(sum(N))

        d = sqrt(size(Cdag[1], 2)) |> Int
        trans_mat = LadderOperators.trans_matrix_fock(idx, d^2)
        temp = trans_mat * (LadderOperators.show_basis(N) |> Base.Fix2(map_u1, idx))
        # display(temp)
        @test LadderOperators.show_basis(N)[idx] == temp[idx]
    end
end
