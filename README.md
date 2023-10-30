# LadderOperators.jl

[![Tests](https://github.com/lukasgrunwald/LadderOperators.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/lukasgrunwald/LadderOperators.jl/actions/workflows/CI.yml)

Generate matrix representation of fermionic and bosonic ladder operators on the full Fock-space with $N$ lattice sites. Provide additional routines to reduce to U(1)-blocks with fixed particle number for Fermions. Code developed to study Hubbard-Dimer in multimode cavity, hence the focus is not on the performance of the implementation, but rather on ease of use.

The matrix representation is generated using the Jordan-Wigner strings in an MPS like approach. Ladder operators act as
<!-- $$
c_i^\dagger |n_1, n_2, \dots \rangle = (-1)^{\sum_{j=1}^{i-1} n_j} |n_1, \dots, (n_i + 1)\mod ~2, \dots \rangle,
$$ -->
![equation](https://latex.codecogs.com/svg.image?%5Cbg%7Bwhite%7D%7B%5Ccolor%7BWhite%7Dc_i%5E%5Cdagger%7Cn_1,n_2,%5Cdots%5Crangle=(-1)%5E%7B%5Csum_%7Bj=1%7D%5E%7Bi-1%7Dn_j%7D%7Cn_1,%5Cdots,(n_i&plus;1)%5Cmod~2,%5Cdots%5Crangle,%7D)

so that their matrix representation in a site resolved Slater-Basis reads
<!-- $$
c_i^\dagger = (-1)^{n_1} \otimes (-1)^{n_2} \otimes \dots (-1)^{n_{i-1}} \otimes c_i^\dagger \otimes 1 \otimes \dots \otimes 1.
$$ -->

![equation](https://latex.codecogs.com/svg.image?%5Cbg%7Bwhite%7D%7B%5Ccolor%7BWhite%7Dc_i%5E%5Cdagger=(-1)%5E%7Bn_1%7D%5Cotimes(-1)%5E%7Bn_2%7D%5Cotimes%5Cdots(-1)%5E%7Bn_%7Bi-1%7D%7D%5Cotimes%20c_i%5E%5Cdagger%5Cotimes%201%5Cotimes%5Cdots%5Cotimes%201.%7D)

By implementing the site-local representation of $c_i^\dagger$, one can evaluate the representation on the full Fock space by calculating the associated Kronecker products. For spinless-fermions, with site local basis $|0\rangle, |1\rangle$ the site local representation reads
<!-- $$
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
$$ -->
![equation](https://latex.codecogs.com/svg.image?%5Cbg%7Bwhite%7D%7B%5Ccolor%7Bwhite%7Dc%5E%5Cdagger=%5Cbegin%7Bpmatrix%7D0&0%5C%5C1&0%5Cend%7Bpmatrix%7D,%5Cqquad%20c=%5Cbegin%7Bpmatrix%7D0&1%5C%5C0&0%5Cend%7Bpmatrix%7D%7D)

The representation of additional internal structure, like spins, is performed by treating these degrees of freedom like lattice sites as well.

Implemented are spinless fermions, spinful fermions, multi_band spinful fermions as well as bosons. The corresponding ladder operators can be obtained using the `ladder_operators` function.