#=
1.) Embedding of site local operators onto Full Fock space.
2.) Restriction to U(1)-sector and mapping back to full product Hilbert-Space
=#
# ———————————————————————————————— Fock-Space embedding ———————————————————————————————— #
"""
    embed_fock(N, x::AbstractArray{T}, Z = I(size(x[1], 2)); idx = 1:N) where T<:AbstractMatrix
    embed_fock(N, x::AbstractMatrix{T}, Z = I(size(x, 1)); idx = 1:N) where T<:Number

Elevate the local operator(s) `x` to operators on the full antisymmetric Fock-space
H = H_1 ⊗ H_2 ⊗ ... ⊗ H_N. `Z` is the Jordan-Wigner string matrix, which only needs to be
supplied for fermions.

# Arguments
`N:` Number of lattice sites
`x:` Operator(s) to embed into Fock space
`idx(=1:N):` Site indices for which to evaluate operator
`returns:` x with shape (N, size(x)), if x has internal structure, else simply as vector (N, )
"""
function embed_fock(N, x::AbstractArray{T}, Z = I(size(x[1], 1)); idx = 1:N) where T<:AbstractMatrix
    N == 1 && return x

    d = size(x[1], 1) # Local state dimension
    x_fock = Array{eltype(x), 1 + length(size(x))}(undef, length(idx), size(x)...)

    for i in idx
        for j in CartesianIndices(x)
            x_fock[i, j] = kron(fill(Z, i - 1)..., x[j], fill(I(d), N - i)...)
        end
    end

    return x_fock
end

# Embedding without aditional structure
function embed_fock(N, x::AbstractMatrix{T}, Z = I(size(x, 1)); idx = 1:N) where T<:Number
    N == 1 && return x

    d = size(x, 1)
    x_fock = Vector{typeof(x)}(undef, N)

    for i in idx
        x_fock[i] = kron(fill(Z, i - 1)..., x, fill(I(d), N - i)...)
    end

    return x_fock
end


# —————————————————————————————— Restrict to half-filling —————————————————————————————— #
# Indices that map to half-filled system for fermions and 2-band fermions on dimer
const idx_fermi = [4, 6, 7, 10, 11, 13]
const idx_multi = [4, 6, 7, 10, 11, 13, 18, 19, 21, 25, 34, 35, 37, 41, 49, 66, 67, 69,
                   73, 81, 97, 130, 131, 133, 137, 145, 161, 193]


"""
    idx_u1(Ntot; Nfill = 2)

Find indices of U(1) sector with given filling fraction Nfill. `Ntot` is the total number operator
of the given basis.
"""
idx_u1(Ntot; Nfill = 2) = findall(x -> x == Nfill, diag(Ntot))

# Map full matrix to chosen U(1) sector
map_u1(op::AbstractVector, idx) = @view op[idx]
map_u1(op::AbstractMatrix, idx) = @view op[idx, idx]

"""
    map_u1(op)

Map operator `op` in Full-Hilbert space into half-filled U(1) subspace.
"""
function map_u1(op)
    if size(op, 1) == 4^2 # Fermions
        idx = idx_fermi
    elseif size(op, 1) == 16^2 # 2-band Fermions
        idx = idx_multi
    else
        error("No automatic mapping implemented! Call function with idx instead.")
    end
    return map_u1(op, idx)
end

# Map U(1) sector back to full Fock space
"""
    trans_matrix_fock(idx::AbstractVector, dim::Integer)

Transformation matrix from half-filled U(1) space back to full Hilbert space.
"""
function trans_matrix_fock(idx::AbstractVector, dim::Integer)
    # idx are the inidces that were extracted
    mat = zeros(Int8, dim, length(idx))

    @inbounds for i in eachindex(idx)
        mat[idx[i], i] = 1
    end

    return mat
end

map_full_fock(x::AbstractVector, M::AbstractMatrix) = M * x
map_full_fock(X::AbstractMatrix, M::AbstractMatrix) = M * X * M'

"""
    map_full_fock(x)

Map U(1) representation of operator back to full Hilbert space.
"""
function map_full_fock(x)
    if size(x, 1) == 6
        M = trans_matrix_fock(idx_fermi, 4^2)
    elseif size(x, 1) == 28
        M = trans_matrix_fock(idx_multi, (4^2)^2)
    else
        error("No automatic mapping implemented! Call function with M instead.")
    end
    return map_full_fock(x, M)
end
