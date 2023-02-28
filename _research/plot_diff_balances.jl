using DrWatson
@quickactivate "ColoradoBumblebees"


using CSV
using DataFrames
using CairoMakie
using Statistics
df = CSV.read(datadir("tree_across_balances.csv"), DataFrame)

using ColorSchemes

batch_sizes = unique(df.batch_size)
ens_sizes = unique(df.ensemble_size)
balances = unique(df.balance)

f = Figure(resolution=(1500, 900))

axisargs = (titlealign=:left, xlabel="Proportion True")

cmap = [get(ColorSchemes.magma, i) for i in range(0,0.8,length(ens_sizes))]

for (i,bs) in enumerate(batch_sizes)
    thisdf = filter(x->x.batch_size == bs, df)
    ax = Axis(f[1,i];title="ROC-AUC",subtitle="Batch size: $bs",yticks=0.5:0.05:0.9, axisargs...)
    ylims!(ax, (0.5,0.9))

    for (j,es) in enumerate(ens_sizes)
        this_ens_df = filter(x->x.ensemble_size==es, thisdf)
        mns = Float32[]
        for b in balances
            thisbalance = filter(x->x.balance==b, this_ens_df)
            push!(mns, mean(thisbalance.rocauc))
        end 
        @info mns, balances
        scatterlines!(ax, balances, mns, color=cmap[j])
    end
end

bslines = []

for (i,bs) in enumerate(batch_sizes)
    thisdf = filter(x->x.batch_size == bs, df)
    ax = Axis(f[2,i]; title="PR-AUC",subtitle="Batch size: $bs", yticks=0.3:0.05:0.8, axisargs...)
    ylims!(ax, (0.3,0.8))

    for (j,es) in enumerate(ens_sizes)
        this_ens_df = filter(x->x.ensemble_size==es, thisdf)
        mns = Float32[]
        for b in balances
            thisbalance = filter(x->x.balance==b, this_ens_df)
            push!(mns, mean(thisbalance.prauc))
        end 
        @info mns, balances
        sl = scatterlines!(ax, balances, mns, color=cmap[j])
        i == 1 && push!(bslines, sl)
    end
end

Legend(f[:, end+1],  bslines, ["$es" for es in ens_sizes], "Ensemble Size")

f

save(plotsdir("tree_balancing.png"), f)