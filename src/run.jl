A=rand(5,4)
Σ=A'A
μ=rand(4)
#l=zeros(4)
#u=.6*ones(4)

include("cla.jl")
result = calculate_turningpoints(μ, Σ)
println("result: $result")
