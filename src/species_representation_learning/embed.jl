function embed(data::BeeData, model)
    convert(Dict{Species, Vector}, _embed(data, model))
end