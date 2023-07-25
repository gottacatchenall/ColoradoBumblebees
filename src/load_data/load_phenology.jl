function load_phenology(data)
    species = vcat(bees(data)..., plants(data)...)
    firstdoy, lastdoy = extrema([dayofyear(i.time) for i in interactions(data)])

    dict = Dict()
    for sp in species
        abundances = zeros(Int32, lastdoy - firstdoy + 1)
        ints = interactions(data, sp)
        for i in ints
            abundances[dayofyear(i.time) - firstdoy + 1] += 1
        end
        merge!(dict, Dict(sp => abundances))
    end
    return dict
end
