# calculating the turning points given by Markowitz' critical line algorithm
# see http://www.vwl.unibe.ch/papers/dp/dp0701.pdf p. 22 for a description

#using Debug
using Logging
import Base.max

#Logging.configure(level=WARNING)

# extend max to allow for "nothing" as a parameter
max(x::Nothing, y::Nothing) = nothing
max(x, y::Nothing) = x
max(x::Nothing, y) = y

# extend getindex to get x[i] == nothing for i==nothing
getindex(x::Array, i::Nothing) = nothing

argmax(x) = indmax(x)

argmax(x, condition) =
    if any(condition)
        indmax(x-Inf*!convert(Array{Bool},condition))
    else
      nothing
    end

function calculate_turningpoints(μ, Σ)
    j = argmax(μ)
    𝔽 = [j]
    w = zeros(size(μ))
    w[j] = 1
    𝔹 = setdiff(1:length(μ), 𝔽)
    W = {w}
    λcurrent = Inf
    t = 1
    λ = zeros(size(μ))
    while true
        i_inside = nothing
        i_outside = nothing
        λ[:] = -Inf
        # Case a) Free asset moves to its bound
        if length(𝔽) > 1
            for i in 𝔽
                one_F = ones(size(𝔽))
                one_B = ones(size(𝔹))
                i_relative = findfirst(𝔽, i)
                Σ𝔽I = Σ[𝔽,𝔽]^-1
                Ci = first(-(one_F'*Σ𝔽I*one_F)*((Σ𝔽I*μ[𝔽])[i_relative]) +
                           (one_F'*Σ𝔽I*μ[𝔽])*((Σ𝔽I*one_F)[i_relative]))
                λ[i] = (Σ𝔽I*one_F)[i_relative]/Ci
            end
            i_inside = argmax(λ, [λ[i]<λcurrent && (i in 𝔽)
                                  for i=1:length(λ)])
        end
        # Case b) Asset on its bound becomes free
        if length(𝔽) < length(μ)
            for i in 𝔹
                𝔽i = union(𝔽, [i])
                𝔹i = setdiff(1:length(μ), 𝔽i)
                i_relative = findfirst(𝔽i, i)
                one_Fi = ones(size(𝔽i))
                one_Bi = ones(size(𝔹i))
                Σ𝔽iI = Σ[𝔽i,𝔽i]^-1
                temp1 = (one_Fi'*Σ𝔽iI*one_Fi)*(Σ𝔽iI*μ[𝔽i])[i_relative]
                temp2 = (one_Fi'*Σ𝔽iI*μ[𝔽i])*(Σ𝔽iI*one_Fi)[i_relative]
                Ci = first(-temp1 +
                           temp2)
                λ[i] = (Σ𝔽iI*one_Fi)[i_relative]/Ci
            end
            i_outside = argmax(λ, [λ[i]<λcurrent && (i in 𝔹)
                                   for i=1:length(λ)])
        end
        # Find turning point by comparing cases
        if i_inside == nothing && i_outside == nothing
            break
        end
        t = t+1
        λcurrent = max(λ[i_inside], λ[i_outside])
        if λ[i_inside] == max(λ[i_inside], λ[i_outside])
            𝔽 = setdiff(𝔽, [i_inside])
        else
            𝔽 = union(𝔽, [i_outside])
        end
        one_F = ones(size(𝔽))
        𝔹 = setdiff(1:length(μ), 𝔽)
        one_B = ones(size(𝔹))
        Σ𝔽I = Σ[𝔽,𝔽]^-1
        γ = first(-λcurrent*(one_F'*Σ𝔽I*μ[𝔽])./(one_F'*Σ𝔽I*one_F) +
                  1./(one_F'*Σ𝔽I*one_F))
        w = zeros(size(μ))
        w[𝔽] = λcurrent*Σ𝔽I*μ[𝔽] + γ*Σ𝔽I*one_F
        push!(W,w)

    end
    return W
end
