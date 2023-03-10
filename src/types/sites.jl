abstract type Site end
struct PikesPeak <: Site end
struct Gothic <: Site end
struct ElkMeadows <: Site end
struct Metaweb <: Site end

sitename(::Type{PikesPeak}) = "Pikes Peak"
sitename(::Type{Gothic}) = "Gothic"
sitename(::Type{ElkMeadows}) = "Elk Meadows"
