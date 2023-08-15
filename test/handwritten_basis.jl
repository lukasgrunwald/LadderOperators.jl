#=
Map the basis constucted here onto the basis in my handwritten notes as a consistency check.
=#
using LadderOperators
import LadderOperators: show_basis

# Extend some functions such that algebra with matrices works
Base.:*(x::Integer, word::String) =  x == 0 ? "" : word
Base.:+(st1::String, st2::String) = st1 * st2
Base.zero(st::String) = ""

Cdag, C = ladder_operators(:fermi, 2)
N = Cdag .* C
idx = idx_u1(sum(N))

# Basis transformation to handwritten basis
my_basis = [
    0 0 0 0 0 1
    1 0 0 0 0 0
    0 0 0 1 0 0
    0 0 1 0 0 0
    0 0 0 0 1 0
    0 1 0 0 0 0
]

map_mybasis(x::AbstractVector) = my_basis * x
map_mybasis(X::AbstractMatrix) = my_basis * X * my_basis'

# Verify that this is really my basis
map_mybasis(show_basis(N)[idx]) |> display

print("1↑-2↑=")
map_mybasis(Cdag[1, 1] * C[2, 1] |> (Base.Fix2(map_u1, idx))) |> display

print("1↓-2↑=")
map_mybasis(Cdag[1, 2] * C[2, 1] |> (Base.Fix2(map_u1, idx))) |> display

print("1↑-2↓ =")
map_mybasis(Cdag[1, 1] * C[2, 2] |> (Base.Fix2(map_u1, idx))) |> display

print("1↓-2↓ =")
map_mybasis(Cdag[1, 2] * C[2, 2] |> (Base.Fix2(map_u1, idx))) |> display
