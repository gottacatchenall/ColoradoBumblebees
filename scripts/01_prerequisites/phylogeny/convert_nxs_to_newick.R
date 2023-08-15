library(ape)


bee_tree = read.nexus("./scripts/01_prerequisites/phylogeny/bees/bees.tree")
plant_tree = read.nexus("./scripts/01_prerequisites/phylogeny/plants/plants.tree")

plot(plant_tree)

write.tree(bee_tree, file="bee_tree.newick")
write.tree(plant_tree, file="plant_tree.newick")
