using DrWatson
@quickactivate "ColoradoBumblebees"

# Here you may include files from the source directory
include(srcdir("ColoradoBumblebees.jl"))
using Main.ColoradoBumblebees


# Clean CSVs
clean_interactions(PikesPeak)
clean_interactions(ElkMeadows)
clean_interactions(Gothic)

create_environmental_covariate_data()

create_cooccurence_data()
