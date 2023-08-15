using DrWatson
@quickactivate :ColoradoBumblebees

using CairoMakie

fit_dir = joinpath(artifactdir(), "classification_fits")

fits = ColoradoBumblebees.load.([joinpath(fit_dir, x) for x in readdir(fit_dir)])

# we want to turn as set of fits into a dataframe with params and prauc as
# columns 
function fits_to_dataframe(fits)
    colnames = vcat([collect(keys(struct2dict(r.embed_model))) for r in fits[1].fits[1].representation]...)
    df = DataFrame([c=>[] for c in colnames])
    df.prauc = []
    df.rocauc = []
    df.embed_model = []
    df.classifier = []

    for f in fits
        for fit in f.fits
            embed_models = fit.representation

            for e in embed_models
                dict = struct2dict(e.embed_model)
                for (k,v) in dict
                    push!(df[!, k], v)
                end
                push!(df.classifier, fit.model)
                push!(df.embed_model, e.embed_model)
                push!(df.prauc, fit.fit_stats["prauc"])
                push!(df.rocauc, fit.fit_stats["rocauc"])

            end
        end
    end
    df
end

df = fits_to_dataframe(fits)

df.embed_dims = [x[end] for x in df.encoder_dims]
df.hidden_layers = [length(x)-2 for x in df.rnn_dims]
df.hidden_dim = [x[2] for x in df.rnn_dims]

is_variational(::RecurrentAutoencoder{T}) where T = T == Variational

df.variational =[is_variational(x) for x in df.embed_model]

xs = Float32[]
ys = Float32[]
for (i, emb_dim) in enumerate(unique(df.embed_dims))
    for (j, hidden_layers) in enumerate(unique(df.hidden_layers))
        for (k, unit) in enumerate(unique(df.unit))
            for (l, hidden_dim) in enumerate(unique(df.hidden_dim))
                for (m, is_variational) in enumerate(unique(df.variational))

                    this_df = filter(x ->
                        x.embed_dims == emb_dim && 
                        x.hidden_layers== hidden_layers && 
                        x.unit == unit && 
                        x.hidden_dim == hidden_dim &&
                        x.variational == is_variational, 
                        df
                    )
                    push!(xs, fill(i+j+k+l+m, nrow(this_df))...)
                    push!(ys, this_df.rocauc...)
        
                end
            end 
        end
    end
end

rainclouds(xs, ys)

unique(df.hidden_layers)




df

f = Figure()
ax = Axis(f[1,1])

pr = praucs.(fits)

for i in eachindex(pr)
    scatter!(ax, [i for _ in eachindex(pr[i])], pr[i], color=(:blue,0.1))
end

f


representation(fits[end])[1].embed_model