struct ProjectedOverlap{T,S}
    interaction_richness
    interaction_uncertainty
    sdm_uncertainty
    cooccurrence_dataframe
end

function ColoradoBumblebees.path(po::ProjectedOverlap{T,S}; cluster=false) where {T,S}
    lead = cluster ? "/scratch/mcatchen/artifacts/" : artifactdir() 
    joinpath(lead, "projected_overlap", "Scenario_$(string(S))_Timespan_$(string(T))")
end