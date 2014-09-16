#using Debug
import Base.max

# extend max to allow for "nothing" as a parameter
max(x::Nothing, y::Nothing) = nothing
max(x, y::Nothing) = x
max(x::Nothing, y) = y


argmax(x, condition) =
    if any(condition)
        indmax(x-Inf*!convert(Array{Bool},condition))
    else
      nothing
    end

function calculate_turningpoints_general(μ, Σ, l, u)
    𝔽, w0 = starting_solution(μ, l, u)
    𝔹 = setdiff(1:length(μ), 𝔽)
    W = {w0}
    λcurrent = Inf
    t = 0
    λcurrent_list = {}
    𝔽list = {}
    while true
        # Case a) Free asset moves to its bound
        i_inside, λ_i_inside, b = asset_moves_to_bound(μ, Σ, l, u,
                                                       𝔽, λcurrent, W[end])
        # Case b) Asset on its bound becomes free
        i_outside, λ_i_outside = asset_becomes_free(μ, Σ,
                                                    𝔽, λcurrent, W[end])
    
        # Find turning points by comparing cases
        if i_inside ≠ nothing || i_outside ≠ nothing
            t = t+1
            push!(W,ones(W[end]))
            W[end][𝔹] = W[end-1][𝔹]
            λcurrent = max(λ_i_inside, λ_i_outside)
            if λ_i_inside == max(λ_i_inside, λ_i_outside)
                𝔽 = setdiff(𝔽, [i_inside])
                W[end][i_inside] = b
            else
                𝔽 = union(𝔽, [i_outside])
            end
            one_F = ones(size(𝔽))
            𝔹 = setdiff(1:length(μ), 𝔽)
            one_B = ones(size(𝔹))
            Σ𝔽I = Σ[𝔽,𝔽]^-1
            γ = -λcurrent*(one_F'*Σ𝔽I*μ[𝔽])./(one_F'*Σ𝔽I*one_F) +
                (1-one_B'*W[end][𝔹]+one_F'*Σ𝔽I*Σ[𝔽,𝔹]*W[end][𝔹])./
                  (one_F'*Σ𝔽I*one_F)
            W[end][𝔽] = -Σ𝔽I*Σ[𝔽,𝔹]*W[end][𝔹] + γ.*Σ𝔽I*one_F +
                        λcurrent*Σ𝔽I*μ[𝔽]
        end
        if i_inside == nothing && i_outside == nothing
            break
        end
        push!(λcurrent_list, λcurrent)
        push!(𝔽list, 𝔽)
    end
    return W[2:end], λcurrent_list, 𝔽list
end

function starting_solution(μ, l, u)
    local i_free
    w = copy(l)
    i = indmax(μ)
    while sum(w)<1
        i_free = i
        w[i] = min(u[i], l[i]+1-sum(w))
        i = argmax(μ, μ.<μ[i])
    end
    𝔽 = [i_free]
    return 𝔽, w
end

function asset_moves_to_bound(μ, Σ, l, u, 𝔽, λcurrent, w)
    # A sole asset cannot move to bound
    if length(𝔽) == 1
        return nothing, nothing, nothing
    end
    λ = -Inf*ones(size(μ))
    b = zeros(size(μ))
    𝔹 = setdiff(1:length(μ), 𝔽)
    for i in 𝔽
        one_F = ones(size(𝔽))
        one_B = ones(size(𝔹))
        i_relative = findfirst(𝔽, i)
        Σ𝔽I = Σ[𝔽,𝔽]^-1
        Ci = first(-(one_F'*Σ𝔽I*one_F)*((Σ𝔽I*μ[𝔽])[i_relative]) +
                   (one_F'*Σ𝔽I*μ[𝔽])*((Σ𝔽I*one_F)[i_relative]))
        if Ci == 0
            @bp
            nothing
        end
        if Ci ≥ 0
            b[i] = u[i]
        else
            b[i] = l[i]
        end
        temp1 = (1-one_B'*w[𝔹]+one_F'*Σ𝔽I*Σ[𝔽,𝔹]*w[𝔹])*((Σ𝔽I*one_F)[i_relative])
        temp2 = (one_F'*Σ𝔽I*one_F)*(b[i]+((Σ𝔽I*Σ[𝔽,𝔹]*w[𝔹])[i_relative]))
        λ[i] = Ci^-1 * first(temp1
                             -temp2)
    end
#    @bp
    correction = if λcurrent == Inf
                    0.0
                 else
                    10eps(λcurrent)
                 end    
    i_inside = argmax(λ, [(λ[i]<λcurrent-correction) && (i in 𝔽) for i=1:length(λ)])
    if i_inside == nothing
        return nothing, nothing, nothing
    end
    return i_inside, λ[i_inside], b[i_inside]
end

function asset_becomes_free(μ, Σ, 𝔽, λcurrent, w)
    # Skip procedure if all assets are free
    if length(𝔽) == length(μ)
        return nothing, nothing
    end
    λ = -Inf*ones(size(μ))
    b = zeros(size(μ))
    𝔹 = setdiff(1:length(μ), 𝔽)
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
        bi = w[i]
        λ[i] = Ci^-1*first((1-one_Bi'*w[𝔹i]+one_Fi'*Σ𝔽iI*Σ[𝔽i,𝔹i]*w[𝔹i])
                                                  *(Σ𝔽iI*one_Fi)[i_relative]
                           -(one_Fi'*Σ𝔽iI*one_Fi)*(bi+(Σ𝔽iI*Σ[𝔽i,𝔹i]*w[𝔹i])[i_relative]))
    end
    correction = if λcurrent == Inf
                    0.0
                 else
                    10eps(λcurrent)
                 end
    i_outside = argmax(λ, [(λ[i]<λcurrent-correction) && (i in 𝔹) for i=1:length(λ)])
    if i_outside == nothing
        return nothing, nothing
    end
    return i_outside, λ[i_outside]
end
