module LadderOperators

# Matrix rep of Ladder operators
export local_spinless_fermions, local_fermions, local_fermions_multiorbital
export local_bosons, local_spin
export ladder_operators, spin_operators

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

Generate matrix represantion of ladder operators for lattice models in full (anti-) symmetric
Fock space based on the specified type and system size.

# Arguments
- `type{:spinless, :fermi, :fermi_multi, :bosons}`: The type of ladder operators to generate.
- `Nsites::Int`: The number of sites in the system.
- `args...`: Additional arguments required for specific operator types.

# Returns
- A tuple containing two matrices: `Cdag` and `C`. These matrices represent the creation and
annihilation operators for the specified ladder operators and have shape (Nsites, Ninternal...)
"""
function ladder_operators(type, Nsites, args...)
    # Load the appropriate local opeartors
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
    else
        error("$type not implemented!")
    end

    # Embed into full Fock-space.
    Cdag = embed_fock(Nsites, Cdag, Z)
    C = embed_fock(Nsites, C, Z)

    return Cdag, C
end
"""
    spin_operators(type, Nsites, args...)

Generate multi-site spin operators based on the specified type and system size.

# Arguments
- `type{:heisenberg, :fermi}`: The type of spin representation to generate.
- `Nsites::Int`: The number of sites or particles in the system.
- `args...`: Additional arguments, if needed.

# Returns
- Array containing multi-site spin operators. The dimensions are (Nsites, S_{i=x,y,z}).
"""
function spin_operators(type, Nsites, args...)
    # Generate site-local operators
    s⃗ = local_spin(type) # Spin-1/2 operator
    dim = size(s⃗[1], 1)
    Z = I(dim) # No Wigner strings for spin-operators

    S⃗ = embed_fock(Nsites, s⃗, Z) # Embed into multi-particle fock space

    return S⃗
end

end # module
