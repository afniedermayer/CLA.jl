using Debug
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

@debug function calculate_turningpoints_general(μ, Σ, l, u)
    #    @bp
    𝔽, w0 = starting_solution(μ, l, u)
    𝔹 = setdiff(1:length(μ), 𝔽)
    W = {w0}
    λcurrent = Inf
    t = 0
    #    WB = {}
    while true
        # Case a) Free asset moves to its bound
        i_inside, λ_i_inside, b = asset_moves_to_bound(μ, Σ, l, u,
                                                       𝔽, λcurrent, W[end])
        # Case b) Asset on its bound becomes free
        i_outside, λ_i_outside = asset_becomes_free(μ, Σ,
                                                    𝔽, λcurrent, W[end])
        println("W: $W")
        println("i_inside: $i_inside, λ_i_inside: $λ_i_inside, b; $b")
        println("i_outside: $i_outside, λ_i_outside: $λ_i_outside, b; $b")
        println("𝔽: $𝔽, 𝔹: $𝔹")
        println("W: $(W[end])")
#        @bp
    
        # Find turning points by comparing cases
        if i_inside ≠ nothing || i_outside ≠ nothing
            t = t+1
#            𝔹 = setdiff(1:length(μ), 𝔽)
#            one_B = ones(size(𝔹))
#            W[end][𝔹] = W[end-1][𝔹]
            push!(W,ones(W[end]))
            W[end][𝔹] = W[end-1][𝔹]
            λcurrent = max(λ_i_inside, λ_i_outside)
            if λ_i_inside == max(λ_i_inside, λ_i_outside)
                deleteat!(𝔽, findin(𝔽, i_inside))
                W[end][i_inside] = b
            else
                𝔽 = union(𝔽, [i_outside])
            end
            one_F = ones(size(𝔽))
            𝔹 = setdiff(1:length(μ), 𝔽)
            one_B = ones(size(𝔹))
            println("𝔽 is ", 𝔽)
            println("Σ[𝔽,𝔽] is ", Σ[𝔽,𝔽])
            println("W is ", W)
            Σ𝔽I = Σ[𝔽,𝔽]^-1
            γ = -λcurrent*(one_F'*Σ𝔽I*μ[𝔽])./(one_F'*Σ𝔽I*one_F) +
                (1-one_B'*W[end][𝔹]+one_F'*Σ𝔽I*Σ[𝔽,𝔹]*W[end][𝔹])./
                  (one_F'*Σ𝔽I*one_F)
            W[end][𝔽] = -Σ𝔽I*Σ[𝔽,𝔹]*W[end][𝔹] + γ.*Σ𝔽I*one_F +
                        λcurrent*Σ𝔽I*μ[𝔽]
            println("w: $(W[end])")
            println("λcurrent: $λcurrent")
            @bp
        end
        if i_inside == nothing && i_outside == nothing
            break
        end
    end
    return W[2:end]
end

function starting_solution(μ, l, u)
    local i_free
    w = copy(l)
    i = indmax(μ)
    while sum(w)<1
        i_free = i
        w[i] = min(u[i], l[i]+1-sum(w))
        i = argmax(μ, μ.<μ[i]) # indmax(μ-Inf*(μ.≥μ[i]))
    end
    𝔽 = [i_free]
    return 𝔽, w
end

@debug function asset_moves_to_bound(μ, Σ, l, u, 𝔽, λcurrent, w)
    # A sole asset cannot move to bound
    @bp
    if length(𝔽) == 1
        return nothing, nothing, nothing
    end
    λ = zeros(size(μ))
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
    i_inside = argmax(λ, [λ[i]<λcurrent && (i in 𝔽) for i=1:length(λ)])
    println("i_inside: $i_inside, condition: $([λ[i]<λcurrent && (i in 𝔽) for i=1:length(λ)])")
    if i_inside == nothing
        @bp
        return nothing, nothing, nothing
    end
    return i_inside, λ[i_inside], b[i_inside]
end

@debug function asset_becomes_free(μ, Σ, 𝔽, λcurrent, w)
    # Skip procedure if all assets are free
    if length(𝔽) == length(μ)
        return nothing, nothing
    end
    λ = zeros(size(μ))
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
        try
        λ[i] = Ci^-1*first((1-one_Bi'*w[𝔹i]+one_Fi'*Σ𝔽iI*Σ[𝔽i,𝔹i]*w[𝔹i])
                                                  *(Σ𝔽iI*one_Fi)[i_relative]
                           -(one_Fi'*Σ𝔽iI*one_Fi)*(bi+(Σ𝔽iI*Σ[𝔽i,𝔹i]*w[𝔹i])[i_relative]))
        catch err
            @bp
            rethrow(err)
        end
    end
    i_outside = argmax(λ, [λ[i]<λcurrent && i in 𝔹 for i=1:length(λ)])
    println("i_outside: $i_outside, condition: $([λ[i]<λcurrent && i in 𝔹 for i=1:length(λ)])")
    if i_outside == nothing
        @bp
        return nothing, nothing
    end
    return i_outside, λ[i_outside]
end
