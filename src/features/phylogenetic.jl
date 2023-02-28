using Phylo, Random

planttree = open(parsenewick, datadir("public", "phylogeny", "plant_tree.tre"))

rand!(BrownianTrait(planttree, "Trait"), planttree) 


d = DataFrame(
    nodename=getnodename.(planttree, traversal(planttree, preorder)), 
    trait=getnodedata.(planttree, traversal(planttree, preorder), "Trait"))