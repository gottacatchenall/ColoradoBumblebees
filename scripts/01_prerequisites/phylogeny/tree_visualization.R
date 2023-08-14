install.packages("BiocManager")
install.packages("ggplot2")

BiocManager::install("ggtree")

require(ggtree)
require(treeio)
require(ape)
require(ggplot2)

planttree <- ape::read.nexus(file = "./scripts/01_prerequisites/phylogeny/plants/plants.tree")
beetree <- ape::read.nexus(file = "./scripts/01_prerequisites/phylogeny/bees/bees.tree")


planttree = rename_taxa(planttree, plantdf, Species.Codename, Plants)
p_plant = ggtree(
    planttree,  
    layout="circular",
    branch.length = "none"
) + geom_tiplab(offset=0.3) + geom_text(aes(label=node)) 
p_plant

plantdf = read.csv("data/public/phylogeny/raw_sequences/plant_sequences.csv")
beedf = read.csv("data/public/phylogeny/raw_sequences/bee_sequences.csv")
beetree = rename_taxa(beetree, beedf, Species.Codename, Bee)

p_bee =  
p_plant = ggtree(
    planttree,  
    layout="circular",
    branch.length = "none"
) + geom_tiplab(offset=0.3) + 
    geom_hilight(node=286, fill="#49b883", alpha=.6) +  # lamiids
        geom_hilight(node=164, fill="#2cb8b3", alpha=.4) +    # superas
    geom_hilight(node=165, fill="steelblue", alpha=.5) +    # asterids
    geom_hilight(node=316, fill="#f5992a", alpha=.6) + # monocots 
    geom_hilight(node=282, fill="#d18a98", alpha=.6) + # eudicots 
    geom_hilight(node=240, fill="#d68eed", alpha=.4)  +  # superrosids
    geom_hilight(node=241, fill="#93709e", alpha=.5)  +  # rosid 
    theme(plot.margin = unit(c(1,1,1,1), "in"))
p_plant

ggsave("plants.png", width = 4000, height = 4000,dpi=300, units = "px")




ggtree(
    beetree,  
    branch.length='none'
) + geom_tiplab(align=TRUE, linesize = .7, size=7) + 
    ggplot2::xlim(0, 18)  +
geom_hilight(node=32, fill="#f5992a", alpha=.6) +  # pyro
geom_hilight(node=38, fill="#bd1f5e", alpha=.6) + # Cullumanobombus
geom_hilight(node=27, fill="steelblue", alpha=.6) +  # Thoracobombus
geom_hilight(node=39, fill="#93709e", alpha=.6) # Psithyrus



ggsave("bees.png", width = 3500, height = 4000,dpi=300, units = "px")

