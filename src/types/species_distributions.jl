abstract type Scenario end 

struct Baseline <: Scenario end 
struct SSP1_26 <: Scenario end 
struct SSP2_45 <: Scenario end 
struct SSP3_70 <: Scenario end 

struct Timespan{Y,Z} end
Base.string(::Type{Timespan{Y,Z}}) where {Y,Z} = string(Y.value)*"_"*string(Z.value)

const TIMESPANS = vcat(Timespan{Year(2000), Year(2015)}, [Timespan{Year(2000+i),Year(2000+i+9)} for i in 20:10:90]...)

struct SpeciesDistribution{L<:SimpleSDMLayer,S<:Species,SC<:Scenario,D<:Dict,T<:Timespan}
    species::S
    probability::L 
    uncertainty::L
    fit_stats::D
    timespan::Type{T}
    scenario::Type{SC}
end



baseline() = Timespan{Year(2000), Year(2015)}


dirname(::Type{SSP1_26}) = "ssp126"
dirname(::Type{SSP2_45}) = "ssp245"
dirname(::Type{SSP3_70}) = "ssp370"


