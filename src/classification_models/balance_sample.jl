function balance_sample(y::Vector{Bool}, train_idx::Vector{I}, batch_size::I, true_pct::F) where {I<:Integer,F<:AbstractFloat}
    itrue = findall(x -> x == true, y[train_idx])
    ifalse = findall(x -> x == false, y[train_idx])

    ntrue, nfalse = Int32(floor(batch_size * true_pct)), Int32(floor(batch_size * (1 - true_pct)))
    return [shuffle(itrue)[1:ntrue]..., shuffle(ifalse)[1:nfalse]...]
end
