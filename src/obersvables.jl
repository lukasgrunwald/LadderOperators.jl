#=
1.) Evaluation of partial trace in reduced Hilbert space by matrix multiplications
=#

# ———————————————————————————————————— Partial trace ——————————————————————————————————— #
function projection_operators(d, idx = 1:d^2; trace_out = :dimer1)
    # Containers for evaluation
    projectors = Matrix{Int8}[]
    vec = zeros(Int8, d) # local basis vector

    for i in range(1, d)
        # Reset vector entries
        vec .= 0
        vec[i] = 1.0

        if trace_out == :dimer1
            push!(projectors, kron(vec, I(d)))
        elseif trace_out == :dimer2
            push!(projectors, kron(I(d), vec))
        end
    end

    return [@view x[idx, :] for x in projectors]
end

"""
    partial_trace(O; trace_out = :dimer1)

Evaluate partial trace of O, which is represented as matrix in half-filled basis. Atm only
implemented for the dimer system.
"""
function partial_trace(O; trace_out = :dimer1)
    # Find local dimension
    dim = size(O, 1)
    if dim == 6 || dim == 16 # :fermi
        d = 4
        idx = dim == 6 ? idx_fermi : 1:d^2
    elseif dim == 28 || dim == 256 # :fermi_multi
        d = 4 * 4 # two-bands
        idx = dim == 28 ? idx_multi : 1:d^2
    else
        @error("Automatic method not implemented")
    end

    O_traced = zeros(d, d)
    pjs = projection_operators(d, idx; trace_out)

    for i in range(1, d)
        O_traced += pjs[i]' * O * pjs[i]
    end

    return O_traced
end
