

@Base.kwdef struct LFSVD
    α = [0, 0.5, 0.5, 0]
    truncated_dims = 8
end


# different typical alphas
# degree only: [0, 0.5, 0.5, 0]
# connectance only: [0 0 0 1]
# hybrid [0 1/3 1/3 1/3]


data = load_data()
M = ColoradoBumblebees.adjacency(metaweb(data))

using LinearAlgebra

function linear_filter(I, α)
    n,m = size(I)

    A = zeros(size(I))

    d = ((m*n)^-1)*sum(I)
    for i in 1:n, j in 1:m
        b = (n^-1)*sum(I[:,j])
        c = (m^-1)*sum(I[i,:])    
        A[i,j] = dot(α, [A[i,j],b,c,d])
    end
    A
end

function lfsvd(I, α, truncated_dims=5)
    A = linear_filter(I, α)
    U,Λ,V = svd(A)
    Λt = Λ
    Λt[truncated_dims:end] .= 0

    U*diagm(Λt)*V'
end


# Simulated annealing on α


α = [0 0.5 0.5 0]


function energy(T₀, λ, t)
    return T₀ * (λ^t)
end

function move_probability(δᵢ, δ₀, ε)
    return exp(-(δᵢ - δ₀) / ε)
end

best_prauc = zeros(Float64, 10_000)
best_prauc[1] = 0.
T₀, λ = 1., 0.9999

alphas = zeros(Float64, 10_000, 4)
alphas[begin,:] .= rand(Dirichlet(ones(4)))
using ProgressMeter

progbar = Progress(10_000)
@showprogress for i in axes(best_prauc, 1)[2:end]
    #α = #rand(Dirichlet(ones(4)))
    α = alphas[i-1,:]  .+ (0.01 .+ 0.02rand(4))
    α = α ./ sum(α)

    imp = lfsvd(M, α)
    δᵢ = computemeasures(M, imp)[:prauc]
    en = energy(T₀, λ, i)
    if δᵢ < best_prauc[i-1]
        best_prauc[i] = δᵢ
    else
        if rand() <= move_probability(δᵢ, best_prauc[i - 1], en)
            best_prauc[i] = δᵢ
            alphas[i,:] .= α
        else
            alphas[i,:] .= alphas[i-1,:]
            best_prauc[i] = best_prauc[i - 1]
        end
    end

    #if i % 100 == 0
    #    ProgressMeter.next!(progbar; showvalues = [(Symbol("Best PR-AUC"), best_prauc[i]), (Symbol("Temp"), en)])
    #end
end


using CairoMakie

scatter(1:length(best_prauc), best_prauc, color=(:dodgerblue, 0.2))



# Steps:

# 1. Take an observed matrix Y to create A 
# A_ij = [Y_ij, (n⁻¹∑ₖ Yₖₗ) , (m⁻¹ ∑ₗ Yᵢₗ) ,  ((mn)⁻¹ ∑Y )] ⋅ α
# ∀ (i,j) ∈ (1:n, 1:m)

# 2. take K = A 
# then we do SVD to get
# A = U ⋅ Σ ⋅ Vᵀ

# Then set all eigenvalues in Σ on the diagonal beyond rank `truncated_dims`
# to 0, creating Σₖ 


# You can recreate K' as 
# K' = U ⋅ Σₖ ⋅ Vᵀ
# then Kᵢⱼ is set to mean of Kₛₜ and Kₜₛ ∀ (s,t)

# but the point is embedding, so take L&R subspaces

# alternative: run at end on matrix predicted by model on species level features
# SVD output prediction mat, and reconstruct with Σₖ

# perhaps hit with FNR correction prior to SVD recontstruction


