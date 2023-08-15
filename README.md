# LadderOperators.jl

Generate matrix representation of fermionic and bosonic ladder operators on the full Fock-space with $N$ lattice sites. Provide additional routines to reduce to U(1)-blocks with fixed particle number for Fermions. Code developed to study Hubbard-Dimer in multimode cavity, hence the focus is not on the performance of the implementation, but rather on ease of use.

The matrix representation is generated using the Jordan-Wigner strings in an MPS like approach. Ladder operators act as
$$
c_i^\dagger |n_1, n_2, \dots \rangle = (-1)^{\sum_{j=1}^{i-1} n_j} |n_1, \dots, (n_i + 1)\mod ~2, \dots \rangle,
$$
so that their matrix representation in a site resolved Slater-Basis reads
$$
c_i^\dagger = (-1)^{n_1} \otimes (-1)^{n_2} \otimes \dots (-1)^{n_{i-1}} \otimes c_i^\dagger \otimes 1 \otimes \dots \otimes 1.
$$
By implementing the site-local representation of $c_i^\dagger$, one can evaluate the representation on the full Fock space by calculating the associated Kronecker products. For spinless-fermions, with site local basis $|0\rangle, |1\rangle$ the site local representation reads
$$
c^\dagger = 
\begin{bmatrix}
0 & 0\\1 & 0
\end{bmatrix}
,
\qquad
c = 
\begin{bmatrix}
0 & 1\\0 & 0
\end{bmatrix}
$$
The representation of additional internal structure, like spins, is performed by treating these degrees of freedom like lattice sites as well.

Implemented are spinless fermions, spinful fermions, multi_band spinful fermions as well as bosons. The corresponding ladder operators can be obtained using the `ladder_operators` function.