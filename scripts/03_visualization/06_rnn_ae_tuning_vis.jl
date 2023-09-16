# List of figures:

# Supplement:
# ----------------------------------------------------------------------
# - S009_Standard_RNN_tuning.png   
# - S010_Variational_RNN_tuning.png   
# - S011_Standard_LSTM_tuning.png
# - S012_Variational_LSTM_tuning.png
# - S013_Standard_GRU_tuning.png
# - S014_Variational_GRU_tuning.png


using DrWatson
@quickactivate :ColoradoBumblebees

using CairoMakie
using ColorSchemes

CairoMakie.activate!(; px_per_unit=3)
fontsize_theme = Theme(; fontsize=32)
set_theme!(fontsize_theme)

fit_dir = joinpath(artifactdir(), "classification_fits")
fit_paths = filter(x->contains(x, "RecurrentAutoencoder"), readdir(fit_dir))

fits = ColoradoBumblebees.load.([joinpath(fit_dir, x) for x in fit_paths])

reps = representation.(fits)

variationalness(::RecurrentAutoencoder{R}) where {R} = R

function filter_representations(unit, varness)
    fits[findall(r -> varness == variationalness(r[1].embed_model) && r[1].embed_model.unit == unit, reps)]
end

num_rnn_hidden_layers(fit) = length(representation(fit)[1].embed_model.rnn_dims) - 2
hidden_dim_channels(fit) = representation(fit)[1].embed_model.rnn_dims[2]
embedding_dims(fit) = representation(fit)[1].embed_model.encoder_dims[end]


rnn_standard = filter_representations(RNN, Standard)
rnn_variational = filter_representations(RNN, Variational)

lstm_standard = filter_representations(LSTM, Standard)
lstm_variational = filter_representations(LSTM, Variational)

gru_standard = filter_representations(GRU, Standard)
gru_variational = filter_representations(GRU, Variational)

ra_sets = [rnn_standard, rnn_variational, lstm_standard, lstm_variational, gru_standard, gru_variational]
ra_names = ["Standard RNN", "Variational RNN", "Standard LSTM", "Variational LSTM", "Standard GRU", "Variational GRU"]
ra_filenames = ["S009_Standard_RNN_tuning", "S010_Variational_RNN_tuning", "S011_Standard_LSTM_tuning", "S012_Variational_LSTM_tuning", "S013_Standard_GRU_tuning", "S014_Variational_GRU_tuning"]

n_reps = 128

_splat(x, n_reps) = vcat([[i for _ in 1:n_reps] for i in x]...)

function build_plot(ra_set, model_name)

    f = Figure(resolution=(3000,1500))

    pltsettings = (;clouds=nothing,
    side=:right,
    boxplot_width=0.15,
    markersize= 10,
    jitter_width=0.1,
    boxplot_nudge = 0.4)

    markeralpha = 0.4
    num_hidden_layers = _splat(num_rnn_hidden_layers.(ra_set), n_reps)
    hidden_channels = _splat(hidden_dim_channels.(ra_set), n_reps)
    embed_dims = _splat(embedding_dims.(ra_set), n_reps)

    unique_hidden_channels = unique(hidden_channels)[sortperm(unique(hidden_channels))]
    unique_num_hidden_layers = unique(num_hidden_layers)[sortperm(unique(num_hidden_layers))]
    unique_embed_dims = unique(embed_dims)[sortperm(unique(embed_dims))]

    prs = vcat(praucs.(ra_set)...)
    rocs = vcat(rocaucs.(ra_set)...)

    # Goal, each panel is a Embed Dims.
    # Each color is Channel size
    # Each is grouped together based on # Hidden layers
    function make_set(gridpanel, data; ylims=(0.5,1), ylabel="PRAUC", xlabel="", title=false)
        for (i,e) in enumerate(unique_embed_dims)
            ax = Axis(
                gridpanel[1,i],
                subtitle= title ? "Embed dims: $e" : "",
                xticks = ([1.5, 3.675, 5.75], string.(unique_num_hidden_layers)),
                xlabel = xlabel, 
                ylabel = i == 1 ? ylabel : "",
                yticks=0.5:0.025:1
            )
            ylims!(ax, ylims...)
            
            category_labels = Float32[]
            data_matrix = Float32[]

            x_position = 1

            marker_alpha = 0.5
            cols = [(x, marker_alpha) for x in [colorant"#3ebd8c", colorant"#3e8cbd", colorant"#9052a1"]]

            col_vec = []
            for (j,hl) in enumerate(unique_num_hidden_layers)
                for (k, hc) in enumerate(unique_hidden_channels)
                    I = findall(isequal(e), embed_dims)
                    J = findall(isequal(hl), num_hidden_layers[I])
                    K = findall(isequal(hc), hidden_channels[I][J])

                    if ylabel=="PRAUC"
                    @info "Model: $(model_name)"
                    @info "Hidden Layers: $hl"
                    @info "Hidden Channel: $hc"
                    @info "Embed Dims: $e"
                    @info "PRAUC: $(mean(data[I][J][K]))"
                    @info "\n\n"
                    end 
                    category_labels = vcat(category_labels, [x_position for _ in 1:n_reps])
                    x_position += 0.5
                    data_matrix = vcat(data_matrix, data[I][J][K])
                    col_vec = vcat(col_vec, [cols[k] for _ in hidden_channels[I][J][K]])
                end

                x_position += 0.75
            end 
            rainclouds!(ax, category_labels, data_matrix; color=col_vec, pltsettings...)
        end
    end 

    g = GridLayout(f[1,1])
    
    top = g[1,1]
    


    ax = Axis(top, 
        title=model_name, 
        titlealign=:left, 
        titlesize=48, 
        leftspinevisible=false,
        topspinevisible=false,
        rightspinevisible=false,
        bottomspinevisible=false,)
    hidedecorations!(ax)
    prpanel = g[2,1]
    rocpanel = g[3,1]

    els = [PolyElement(color=c) for c in [colorant"#3ebd8c", colorant"#3e8cbd", colorant"#9052a1"]]
    Legend(g[2:3,2], width=220, els, string.(unique_hidden_channels), "Channels per\nhidden node")

    rowsize!(g, 1, Relative(0.02))

    make_set(prpanel, prs; ylims=(0.6,0.85), ylabel="PRAUC", title=true)
    make_set(rocpanel, rocs; ylims=(0.8,0.95), xlabel="Number of hidden layers", ylabel="ROCAUC")
    f
end

f = build_plot(ra_sets[1], ra_names[1])




for i in eachindex(ra_sets)
    f = build_plot(ra_sets[i], ra_names[i])

    save(plotsdir("$(ra_filenames[i]).png"), f)
    save(plotsdir("$(ra_filenames[i]).svg"), f)
end 

