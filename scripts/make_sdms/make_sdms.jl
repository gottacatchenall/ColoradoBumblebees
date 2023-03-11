using DrWatson
@quickactivate "ColoradoBumblebees"

include(srcdir("ColoradoBumblebees.jl"))
using Main.ColoradoBumblebees

data = load_data()

make_sdms(data, GaussianBRT(); cluster=true)
