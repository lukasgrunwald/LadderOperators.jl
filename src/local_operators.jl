#=
Site local ladder operators, represented on single site using Wigner-Strings
=#

# ———————————————————————————————— site local operators ———————————————————————————————— #
"""
    local_bosons(ncut)

Matrix representation of the site local ladder operators a and a_dagger
in the harmonic oscillator basis |n⟩ with cutoff `ncut`.

`ncut:` Dimension of the matrix rep
`return:` a_dagger, a, n as `ncut⋅ncut` matrices
"""
function local_bosons(ncut)
    # Container for ladder operator
    a = zeros(ncut, ncut)

    for i in range(1, ncut-1)
        # <i|a|j> = sqrt(j) δ(i, j - 1)
        a[i, i + 1] = sqrt(i) # + 1, but labeling starts at 1 and not at 0!
    end

    a_dagger = a' |> collect # Else mapping to Fock space does not work so well
    n = Diagonal(0:(ncut-1))
    return a_dagger, a, n
end

"""
    local_spinless_fermions()

Representation of local spinless fermion operators (on one lattice site) using Jordan-Wigner
string Z
`- returns:` C†, C, Zjw, N
"""
function local_spinless_fermions()
    # Ladder and Wigner string for spinless fermions
    C = Int8[0 1; 0 0] # c
    Cdag = Int8[0 0; 1 0] # c†
    N = Cdag * C
    Zjw = Int8[1 0; 0 -1] # Jordan-Wigner string operator

    return Cdag, C, Zjw, N
end

"""
    local_fermions()

Representation of local spinful fermion operators (on one lattice site) using Jordan-Wigner
string Z. Ordered, such that ↑ always comes before ↓.
`- returns:` C†, C, Zjw, N. Ladder operators are vectors in the spin indices, with 1 => ↑, 2 => ↓
"""
function local_fermions()
    Cdag, C, Zjw, _ = local_spinless_fermions()

    # Spinful Fermions
    Cdagup = kron(Cdag, I(2))
    Cup = kron(C, I(2))

    Cdagdn = kron(Zjw, Cdag)
    Cdn = kron(Zjw, C)

    # Number operators
    Nup = Cdagup * Cup
    Ndn = Cdagdn * Cdn

    Zspin = kron(Zjw, Zjw)

    return [Cdagup, Cdagdn], [Cup, Cdn], Zspin, [Nup, Ndn]
end

"""
    local_fermions_multiorbital(Nb = 2)

Representation of local Nb-band spinful fermion operators (on one lattice site) using
Jordan-Wigner string Z. Ordered, such that ↑ always comes before ↓ and band i comes before band i + 1.
`- returns:` C†, C, Zjw, N. Index-structure of operators is (#band, 2 = {↑, ↓}).
"""
function local_fermions_multiorbital(Nb = 2)
    Cdag, C, Z, _ = local_fermions()

    # Multi-Orbital Spinful-Fermions
    Z_band = kron(fill(Z, Nb)...)
    C_band, Cdag_band = [Matrix{Matrix{Int8}}(undef, Nb, 2) for _ in 1:2]
    N = Matrix{Matrix{Int8}}(undef, Nb, 2) # Number operators

    # Implementation of A in Hilbert space as
    # Z ⊗ ... ⊗ Z ⊗ A ⊗ 1 ... ⊗ 1
    for i in 1:Nb
        for j in 1:2
            Cdag_band[i, j] =  kron(fill(Z, i - 1)..., Cdag[j], fill(I(4), Nb - i)...)
            C_band[i, j] =  kron(fill(Z, i - 1)..., C[j], fill(I(4), Nb - i)...)
            N[i, j] = Cdag_band[i, j] * C_band[i, j]
        end
    end

    return Cdag_band, C_band, Z_band, N
end
