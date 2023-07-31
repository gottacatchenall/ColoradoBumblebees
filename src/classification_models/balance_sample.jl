function balance_sample(y, I, batch_size=64, true_pct=0.5)
    itrue = findall(x -> x == true, y)
    ifalse = findall(x -> x == false, y)

    filter!(x -> x ∈ I, itrue)
    filter!(x -> x ∈ I, ifalse)

    ntrue, nfalse = Int32(floor(batch_size * true_pct)),
    Int32(floor(batch_size * (1 - true_pct)))
    return [shuffle(itrue)[1:ntrue]..., shuffle(ifalse)[1:nfalse]...]
end
