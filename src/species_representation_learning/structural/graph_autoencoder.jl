Base.@kwdef struct GraphAutoencoder{S} <: Structural
    num_input_features = 0. 
    embedding_dim = 0.
end

outdim(gae::GraphAutoencoder{S}) where S = gae.embedding_dim
outdim(gae::GraphAutoencoder{S}, ::Union{Type{Bee},Type{Plant}}) where S = outdim(gae)
