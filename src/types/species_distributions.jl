abstract type Scenario end 

struct Baseline <: Scenario end 
struct SSP1_26 <: Scenario end 
struct SSP2_45 <: Scenario end 
struct SSP3_70 <: Scenario end 
struct Timespan{Y,Z} end
const TIMESPANS = vcat(Timespan{Year(2000), Year(2015)}, [Timespan{Year(2000+i),Year(2000+i+9)} for i in 20:10:90]...)

struct SpeciesDistribution{L<:SimpleSDMLayer,S<:Species,SC<:Scenario}
    species::Species
    probability::L 
    uncertainty::L
    timespan::Timespan
    scenario::SC
end

baseline() = Timespan{Year(2000), Year(2015)}


dirname(::Type{SSP1_26}) = "ssp126"
dirname(::Type{SSP2_45}) = "ssp245"
dirname(::Type{SSP3_70}) = "ssp370"


