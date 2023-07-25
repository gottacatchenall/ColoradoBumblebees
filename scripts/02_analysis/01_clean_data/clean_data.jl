using DrWatson
@quickactivate "ColoradoBumblebees"

include(srcdir("ColoradoBumblebees.jl"))
using Main.ColoradoBumblebees


clean_interactions(PikesPeak)
clean_interactions(Gothic)
clean_interactions(ElkMeadows)


clean_environmental_covariate_data()