function ColoradoBumblebees.save(overlap::ProjectedOverlap{T,S}; cluster=false) where {T,S}
    outdir = path(overlap; cluster=cluster)
    mkpath(outdir) 

    metadata = Dict(
        :scenario => _ssp_string(S),
        :timespan => string(T),
        :git_commit => gitdescribe(),
        :datetime => string(now()),
    )
    _write_json(joinpath(outdir, "metadata.json"), metadata)

    CSV.write(joinpath(outdir, "cooccurrence.csv"), overlap.cooccurrence_dataframe)
    SpeciesDistributionToolkit.save(joinpath(outdir, "interaction_richness.tif"), overlap.interaction_richness)
    SpeciesDistributionToolkit.save(joinpath(outdir, "interaction_uncertainty.tif"), overlap.interaction_uncertainty)
    SpeciesDistributionToolkit.save(joinpath(outdir, "sdm_uncertainty.tif"), overlap.sdm_uncertainty)
end


_ssp_string(s) = Dict(
    Baseline => "baseline",
    SSP1_26 => "ssp126",
    SSP2_45 => "ssp245",
    SSP3_70 => "ssp370"
)[s]

_ssp_obj(s) = Dict(
    "baseline" => Baseline,
    "ssp126" => SSP1_26,
    "ssp245" => SSP2_45,
    "ssp370" => SSP3_70
)[s]

function _timespan_and_scenario(metadata)
    t, sc = metadata["timespan"], metadata["scenario"]

    yrs = Year.(split(t, "_"))
    Timespan{yrs[1],yrs[2]}, _ssp_obj(sc)
end

function _load_projected_overlap(po_path)
    int_richness = SimpleSDMPredictor(joinpath(po_path, "interaction_richness.tif"))
    int_uncert = SimpleSDMPredictor(joinpath(po_path, "interaction_uncertainty.tif"))
    sdm_uncert = SimpleSDMPredictor(joinpath(po_path, "sdm_uncertainty.tif"))
    cooc_df = CSV.read(joinpath(po_path, "cooccurrence.csv"), DataFrame)

    metadata = JSON.parsefile(    joinpath(po_path, "metadata.json"))

    T, S = _timespan_and_scenario(metadata)
    ProjectedOverlap{T,S}(int_richness, int_uncert, sdm_uncert, cooc_df)    
end