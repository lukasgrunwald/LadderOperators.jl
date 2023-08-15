#=
Small implementation for visualizing the generated dimer basis
#! Functions not exported, since do not expect to use this often
=#

# Illustration of the Fock space basis
function __fermion_basis(N, i)
    occp = ""

    if N[1][i, i] == 1
        occp *= "↑"
    end
    if N[2][i, i] == 1
        occp *= "↓"
    end

    # Site not occupied
    if N[1][i, i] == 0 == N[2][i, i]
        occp *= "0"
    end

    return occp
end

"""
    show_basis(N::AbstractArray{T, #}) where T

Return pretty printed illustration of the many body dimer basis for (Multiband) Fermions.
"""
function show_basis(N::AbstractArray{T, 1}) where T
    d = size(N[1], 1)
    states = String[]

    for i in 1:d
        # Occupations
        occp1 = if N[1][i, i] == 1 "1" else "0" end
        occp2 = if N[2][i, i] == 1 "1" else "0"end

        push!(states, "|$occp1⟩⊗|$occp2⟩")
    end

    return states
end

function show_basis(N::AbstractArray{T, 2}) where T
    d = size(N[1], 1)
    states = String[]

    for i in 1:d
        # Occupations
        occp1 = __fermion_basis(N[1, :], i)
        occp2 = __fermion_basis(N[2, :], i)

        push!(states, "|$occp1⟩⊗|$occp2⟩")
    end

    return states
end

"""
    show_basis(N::AbstractArray{T, .})

Visual representation of the chosen basis with ↑, ↓, etc. Useful as consistency check!
"""
function show_basis(N::AbstractArray{T, 3}) where T
    d = size(N[1], 1)
    states = String[]

    for i in 1:d
        # Occupations
        occp1_1 = __fermion_basis(N[1, 1, :], i) # 1st band
        occp1_2 = __fermion_basis(N[1, 2, :], i) # 2nd band

        # occp2 = __fermion_basis(N[2, :], i)
        occp2_1 = __fermion_basis(N[2, 1, :], i) # 1st band
        occp2_2 = __fermion_basis(N[2, 2, :], i) # 2nd band

        push!(states, "|$(occp1_1)₁,$(occp1_2)₂⟩⊗|$(occp2_1)₁,$(occp2_2)₂⟩")
    end

    return states
end
