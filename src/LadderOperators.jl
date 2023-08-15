module LadderOperators

# Matrix rep of Ladder operators
export local_spinless_fermions, local_fermions, local_fermions_multiorbital, local_bosons
export ladder_operators

# Basis transformations
export embed_fock
export idx_u1, map_u1, map_full_fock # Mapping between half-filled and Full hilbert space

# Observables
export partial_trace, projection_operators

using LinearAlgebra

include("local_operators.jl")
include("basis_embedding.jl")
include("obersvables.jl")
include("show_basis.jl")

"""
    ladder_operators(type, Nsites, args...)

Matrix representation of Ladder operator of full antisymmetric Fock-space. Options are
{:spinless, :fermi, :fermi_multi, :boson}. The return operator will have shape (Nsites, Ninternal...)
"""
function ladder_operators(type, Nsites, args...)
    # Load the appropriate opeartors
    if type === :spinless
        Cdag, C, Z, _ = local_spinless_fermions()
    elseif type === :fermi
        Cdag, C, Z, _ = local_fermions()
    elseif type === :fermi_multi
        Cdag, C, Z, _ = local_fermions_multiorbital(args...)
    elseif type === :bosons
        Cdag, C, _ = local_bosons(args...)
        ncut = size(Cdag, 1)
        Z = I(ncut) # No Wigner string needed for bosonic operators
    end

    # Embed into full Fock-space.
    Cdag = embed_fock(Nsites, Cdag, Z)
    C = embed_fock(Nsites, C, Z)

    return Cdag, C
end

end # module
