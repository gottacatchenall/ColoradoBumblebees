function cooccurence(data::BeeData)
    cooc = data.cooccurence
    b, p = bees(data), plants(data)
    C = zeros(Bool,length(b),length(p))

    for i in eachindex(b)
        for j in eachindex(p)
            row = findfirst(r-> r.bee == b[i].name && r.plant == p[j].name, eachrow(cooc))
            C[i,j] = cooc[row,:cooccurence]
        end
    end
    C
end