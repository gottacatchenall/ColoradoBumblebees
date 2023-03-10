using DrWatson
@quickactivate "ColoradoBumblebees"

# Here you may include files from the source directory
include(srcdir("ColoradoBumblebees.jl"))
using Main.ColoradoBumblebees

data = load_data()

make_sdms(data)
