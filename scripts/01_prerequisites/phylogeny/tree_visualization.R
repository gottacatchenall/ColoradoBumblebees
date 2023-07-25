install.packages("BiocManager")
install.packages("ggplot2")

BiocManager::install("ggtree")

require(ggtree)
require(treeio)
require(ape)
require(ggplot2)

planttree <- ape::read.tree(file = "plants/plant_tree.newick")
beetree <- read.tree(file = "bees/bee_tree.newick")


planttree <- ape::read.nexus(file = "plant_experiments/plants.newick")
planttree = rename_taxa(planttree, plantdf, Species.Codename, Next.clade)
p_plant = ggtree(
    planttree,  
    layout="circular",
    branch.length = "none"
) + geom_tiplab(offset=0.3) + geom_text(aes(label=node)) 
p_plant




plantdf = read.csv("data/raw_sequences/plant_sequences.csv")
beedf = read.csv("data/raw_sequences/bee_sequences.csv")


planttree = rename_taxa(planttree, plantdf, Species.Codename, Plants)
beetree = rename_taxa(beetree, beedf, Species.Codename, Bee)

p_bee =  ggtree(
    beetree,  
    branch.length='none'
) + geom_tiplab(align=TRUE, linesize = .7, size=7) + 
    ggplot2::xlim(0, 15) 

p_bee


p_plant = ggtree(
    planttree,  
    layout="circular",
    branch.length = "none"
) + geom_tiplab(offset=0.3) + geom_text(aes(label=node)) 
p_plant


ggtree::rotate(p_plant, 285) %>%  
    ggtree::rotate(311)%>%
    flip(263, 236) +
    geom_hilight(node=244, fill="#bd1f5e", alpha=.6) +  # lamiids
    geom_hilight(node=263, fill="#bd1f5e", alpha=.6) +  # lamiids
    geom_hilight(node=237, fill="steelblue", alpha=.6) +    # asterids
    geom_hilight(node=224, fill="steelblue", alpha=.6) +    # asterids
    geom_hilight(node=310, fill="steelblue", alpha=.6) +    # asterids
    geom_hilight(node=313, fill="steelblue", alpha=.6) +    # asterids
    geom_hilight(node=177, fill="steelblue", alpha=.6) +
    geom_hilight(node=312, fill="#7fc87f", alpha=.6) +
    geom_hilight(node=172, fill="orange", alpha=.6) + # monocots 
    geom_hilight(node=167, fill="#dddb61", alpha=.6) + # eudicots 
    geom_hilight(node=8, fill="#dddb61", alpha=.6) + # eudicots 
    geom_hilight(node=9, fill="#dddb61", alpha=.6) + # eudicots 
    geom_hilight(node=168, fill="#db6e9a", alpha=.6) + # rosid 
    geom_hilight(node=279, fill="#db6e9a", alpha=.6) +   # rosid 
    geom_hilight(node=276, fill="#db6e9a", alpha=.6) +   # rosid 
    geom_hilight(node=174, fill="#bd1f5e", alpha=.6) +  # lamiids 
    geom_hilight(node=272, fill="#bd1f5e", alpha=.6) +  # lamiids 
    geom_hilight(node=274, fill="#a876b8", alpha=.6) +  # superrosids
    geom_hilight(node=306, fill="#a876b8", alpha=.6)   # superrosids

p_plant






























p_bee =  ggtree::rotate(p_bee, 23) 

p_bee = p_bee + geom_hilight(node=30, fill="#e5b380", alpha=.6) +  # pyro
geom_hilight(node=25, fill="steelblue", alpha=.6) +  # Thoracobombus
geom_hilight(node=27, fill="#a876b8", alpha=.6) +# Psithyrus
geom_hilight(node=29, fill="#bd1f5e", alpha=.6) # Cullumanobombus

 # + geom_text(aes(label=node))
p_bee

ggsave("bees.png", width = 2500, height = 4000,dpi=300, units = "px")




p = ggtree(
    planttree,  
    layout="circular",
    branch.length = "none"
) + geom_tiplab(offset=0.3)# + geom_text(aes(label=node))

flip(p, 254, 193) %>% 
flip(247,287) %>%
flip(221,190) %>% 
flip(284, 298) %>%
flip(253, 220) %>%
flip(247, 287) %>%
flip(239, 287) %>%
flip(229, 247) %>%
flip(193, 247) %>%
flip(221, 190) %>% 
flip(221, 287) %>%
flip(221, 239) %>%
flip(221, 229) %>% 
flip(239, 287) %>% 
flip(239, 190) %>%
flip(294, 224) %>%
flip(284, 216) %>%
flip(284, 224) %>%
flip(294, 239) %>% 
flip(284, 239) %>%
flip(294, 190) %>%
flip(284, 190) %>%
flip(287, 294) %>%
flip(284, 287) %>%
flip(253, 220) %>%
flip(253, 298) %>%
flip(216, 224) %>%
flip(253, 224) %>%
flip(216, 239) %>%
flip(253, 239) %>%
flip(253, 216) %>%
flip(224, 239) %>%
flip(224, 239) %>%
flip(224, 239) %>%
flip(224, 216) %>%
flip(224, 253) %>%
flip(224,190) %>%
flip(224, 287) %>% 
flip(224, 284) %>%
flip(224, 294) %>%
flip(224, 229) %>%
flip(224, 221) %>%
flip(224, 193) %>%
flip(224, 247) %>%

flip(298, 239) %>%
flip(298, 216) %>%
flip(298, 253) %>%
flip(298, 190) %>%
flip(298, 287) %>%
flip(298, 284) %>%
flip(298, 294) %>%
flip(298, 229) %>%
flip(298, 221) %>%
flip(298, 193) %>%
flip(298, 247) %>%

flip(220, 239) %>%
flip(220, 216) %>%
flip(220, 253) %>%
flip(220, 190) %>%
flip(220, 287) %>%
flip(220,284) %>%
flip(220, 294) %>%
flip(220, 229) %>%
flip(220, 221) %>%
flip(220, 193) %>%
flip(220, 247) %>%
flip(224, 254) %>%
flip(298, 254)  %>%
flip(224, 298) %>%
flip(229, 294) %>%
flip(229, 284) %>%
ggtree::rotate(224) %>%
ggtree::rotate(216) %>%
ggtree::rotate(247) +
geom_tree() + 
theme_tree() + 
geom_rootedge()+
    geom_hilight(node=193, fill="steelblue", alpha=.6) +  # asterid
    geom_hilight(node=252, fill="steelblue", alpha=.6) +  # asterid
    geom_hilight(node=221, fill="steelblue", alpha=.6) +  # asterid
    geom_hilight(node=230, fill="steelblue", alpha=.6) +  #asterids
    geom_hilight(node=294, fill="steelblue", alpha=.6) +  #asterids
    geom_hilight(node=284, fill="steelblue", alpha=.6) +  #asterids
    geom_hilight(node=250, fill="#a876b8", alpha=.6) +    # superasterids
    geom_hilight(node=248, fill="#a876b8", alpha=.6) +    # superasterids
    geom_hilight(node=254, fill="#db6e9a", alpha=.6) +   # rosids
    geom_hilight(node=227, fill="#db6e9a", alpha=.6) +   # rosids
    geom_hilight(node=220, fill="#db6e9a", alpha=.6) +   # rosids
    geom_hilight(node=298, fill="#bd1f5e", alpha=.6) +   # superrosids
    geom_hilight(node=225, fill="#bd1f5e", alpha=.6) +   # superrosids
    geom_hilight(node=287, fill="#7fc87f", alpha=.6) +   # lamiids 
    geom_hilight(node=190, fill="#7fc87f", alpha=.6) +   # lamiids 
    geom_hilight(node=187, fill="#7fc87f", alpha=.6) +   # lamiids 
    geom_hilight(node=301, fill="#7fc87f", alpha=.6) +   # lamiids 
    geom_hilight(node=217, fill="#7fc87f", alpha=.6) +   # lamiids 
    geom_hilight(node=253, fill="#7fc87f", alpha=.6) +   # lamiids 
    geom_hilight(node=239, fill="#e5b380", alpha=.6) +  # eudicots
    geom_hilight(node=218, fill="#dddb61", alpha=.6) + #+  # monocots
    theme(plot.margin = unit(c(4,4,4,4), "cm")) 



    geom_tippoint(node=6, fill="#dddb61", alpha=.6) + # delphinium  -> Eudicots
    geom_hilight(node=7, fill="#dddb61", alpha=.6) + # packera  -> Asterids
    geom_hilight(node=8, fill="#dddb61", alpha=.6) + # symphyotrichum -> Asterids
    geom_hilight(node=93, fill="#dddb61", alpha=.6) + # pedicularis -> Lamiids
    geom_hilight(node=94, fill="#dddb61", alpha=.6) + # monarda -> Lamiids
    geom_hilight(node=95, fill="#dddb61", alpha=.6) + # nepeta -> Lamiids
    geom_hilight(node=156, fill="#dddb61", alpha=.6) + # verbascum  -> Lamiids
    geom_hilight(node=157, fill="#dddb61", alpha=.6) + # agastache  -> Lamiids
    geom_hilight(node=158, fill="#dddb61", alpha=.6) + # prunella -> Lamiids
    geom_hilight(node=178, fill="#dddb61", alpha=.6) + # mentzelia  -> Asterids
 

colors <- c(
    "Asterids" = "steelblue", 
    "Rosids" = "#db6e9a", 
    "Superrosids" = "#a44169",
    "Lamiids" = "#7fc87f",
    "Superasterids" = "#a876b8",
    "Eudicots" = "#e5b380",
    "Monocots" = "#dddb61"
    )


plt

ggsave("test.png", width = 4000, height = 4000,dpi=300, units = "px")

 #+geom_hilight(node=254, fill="darkgreen", alpha=.6)

# superasterid
# 247

# Asterid internal nodes
# 193, 230, 221, 294, 284

# Rosid internal nodes
# 254, 224, 220


