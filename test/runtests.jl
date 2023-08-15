using Test

@testset "Commutators" verbose = true begin include("local_operators.jl") end
@testset "Embedding" verbose = true begin include("embedding.jl") end
@testset "Observables" verbose = true begin include("observables.jl") end
