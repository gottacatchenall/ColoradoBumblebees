#import "template.typ": *
#import "@preview/wordometer:0.1.4": word-count, total-words
#import "authors.typ": *
//#show: doc => style(doc)

#let title = "Mapping the seasonal dynamics of the bumble bee pollination network of the Southern Rocky Mountains"
#let abstract = "TK"


#showtitle(title)
#showauthors(authors, institutions, corresponding_email)
#showabstract(abstract)
#show: word-count.with(exclude: (figure, bibliography))


#set heading(numbering: "1.")
#show heading: it => {
  it.body
  v(0em)
}
#show cite: set text(fill: rgb("#528fd1"))

\
*Word Count*: #total-words
\

_Ecology Letters: 5000 words, 6 figures._

= Introduction

The vast majority of flowering plant species are pollinated by animals @Ollerton2011HowMan. Pollination interactions are essential to plant reproduction, and necessary for the persistence of Earth’s biodiversity @Bascompte2007PlaMut. Pollination is also important for humans because pollinator dependent crops contribute to 35% of global crop production by volume and bumblebees commonly pollinate over 25 crops worldwide. Recent models indicate 3%–5% of fruit, vegetable, and nut production is lost due to inadequate pollination, leading to an estimated 427,000 excess deaths annually from lost healthy food consumption and associated diseases @Klein2006ImpPol @Smith2022PolDef. But, the human driven forces of climate and land-use change are altering plant-pollinator interactions, and therefore threatening the stability and function of these ecosystems @Kearns1998EndMut As a result, understanding and predicting how plant-pollinator networks are changing due to anthropogenic pressures is a major question in pollinator research [@Mayer2011PolEco] and conservation biology [@Kearns1998EndMut], and ultimately necessary for climate change adaptation and maintaining food security [@Bailes2015HowCana; @Klein2006ImpPol].

A large amount of research into the effects of climate change on pollinator interactions has focused on phenology [@Ogilvie2017IntBee], as shifts in flowering and hatching times caused by warm temperatures move earlier and earlier in the season. Yet this is not the only factor that will govern changes in plant-pollinator systems in the future—--in particular, the spatial distribution of where particular plant and pollinator species can feasibly live is changing, from small (habitat fragmentation) to large (range shift) spatial scales because of climatic niche shifts. Mountainous ecosystems face a unique threat in this respect. Elevational gradients in niche conditions are being pressured “upward” as average temperatures increase [@Hoiss2015IntEff; @Richman2020AsyRan]. We need to be able to reliably predict species interactions to understand how global change will affect these interactions, the networks they form, and the services they provide people. Why is predicting interactions challenging?—for two reasons: (1) interaction data is filled with false-negatives, where two species interact, even if we do not have a record of it (Jordano 2016) and (2) as climate change causes shifts in species’ ranges and phenologies, these networks are becoming “rewired”, where some interactions are lost (as they cease to coexist in the future) and others are gained among species as they start to co-exist for the first time.

Species interactions intrinsically vary in space and time [@Poisot2015SpeWhy]—–-in addition to directional phenological shifts, there is intrinsic year-to-year variation in the structure of empirical pollination networks [@Alarcon2008YeaVar]---so we need a spatiotemporally explicit approach to interaction prediction [@Strydom2021RoaPre]. Going further, forecasting the rewiring of networks as their species shift under different climate projections (and understanding how this rewiring affects ecosystem functioning, persistence, and the services ecosystems provide to people) is crucial to conservation strategies designed for uncertain climate futures.

Here, we use interaction data from 19 bumble bee (Bombus spp.) and 156 angiosperm species from the Southern Rocky Mountains of Colorado to project how this network of 175 species will become “rewired” as each species’ environmental niche shifts over the coming century. Using novel applications of tools from machine-learning, we develop methods for species representation learning. We use this to (1) predict biologically feasible interactions among the species pool and (2) produce spatially explicit predictions of the bumblebee/wildflower pollination network across the southern Rocky Mountains and (3) project how much the ranges of each bee and plant involved in each interaction will overlap in the future under different climate scenarios. We hope these projections of rewiring can guide future research and monitoring to understand how the interactions expected to have the most extreme change (positive or negative) will impact ecosystem functioning and services.

= Methods

== Data 

=== Interaction Data

The primary data source for this study is interaction data collected across three different sites in the subalpine Rocky Mountains of Colorado across a total of 19 Bombus species and
184 flowering plants. The three sample sites are henceforth referred to as Rocky Mountain Biological Lab (RMBL), Pikes Peak, and University of Colorado Mountain Research Station (MRS) (visible in fig. 2.1, upper left). Observations of bumble bee-plant interactions at the RMBL take place within a long-term bumble bee monitoring project (details in  #note("Michael", "add proper cite here") Ogilvie & CaraDonna 2022). In brief,at six permanent study sites near the RMBL, bumble bee interactions with flowering plants are monitored for one hour at weekly intervals for the entire growing season (May–September). At these same weekly intervals, floral abundances of all flowering plant species visited by bumble bees is surveyed within a series of 15 permanent 20 x 0.5 m transects at each site. Observations of bees, their interactions with flowers, and counts of flowers take place within the three dominant habitat types where bumble bee activity is observed: wet meadows, dry meadows, and aspen forest. We used seven years (2015-2021) of observational data on Bombus visits to flowering plants at these study sites to aggregate interaction networks. The Pikes Peak data consists of three seasons (2019-2021) of interaction data collected along an elevational gradient ranging from the near the summit of Pikes Peak (4285 m) to down near the transition of the Great Plains (1872 m), sampled roughly twice a week (described in @Barthell2025BumbleBee). MRS data consists of seven seasons (2015-2021) of data at a six plot at 2900 meters elevation, sampled weekly over the entire flowering season (described in @Resasco2021PlaPol).

=== Occurrence Data

Occurrence data for each species was curated from GBIF within the region of study. The dataset of all occurrence records can be found here. Each occurrence record is associated with a geospatial coordinate and a timestamp. The GBIF dataset consists of XX occurrences from either human observation or preserved specimens within the bounding box seen FIG XX (BBOX HERE).


=== Climate Data

CHELSA provides 19 bioclimatic variables dating back to 1970 at 1km resolution across Earth’s surfaces [@Karger2017CliHig]. CHELSA also provides these same variables pro-
jected into the future until the end of this century. These projections of bioclimatic variables, and the climate projections they are built on, rely on the framework of Shared-Socioeconomic-Pathways [@ONeill2014NewSce] to describe different scenarios in how humanity responds to climate change. These scenarios of climate response vary on two axes: mitigation and adaptation—how humanity will mitigate climate change (i.e. by reducing and eventually reaching net-zero greenhouse gas emissions) and how we will adapt to climate change (i.e. building infrastructure for a warmer world). We consider SSP1-2.6, SSP2-4.5, SSP3-7.0, representing low, moderate, and extreme warming, respectively. The GFDL-4 Earth System Model  [@Held2019StrPer; GFDL-ESM4] was used as the backend Earth System Model for projections, as it is unique in providing these climate projections from CMIP6 (Coupled Model Intercomparison Project). We used the Python library chelsa-cmip6 @Karger2023ChePyt to construct projected bioclimatic variables for each decade from the 2020s until the 2090s for each SSP. Here we consider the period from 2000-2015 to be the climatic baseline (years between 2016 and 2023 require a specified SSP). These projections are then averaged for each of the 19 bioclimatic variables over each decade from the 2020s to the 2090s.

== Fitting Species Phenologies

Gaussian Mixture Model (GMM). One to three components. All visible in Appendix. 

== Mapping Plant-Pollinator Networks

Species distributions were fit using occurrence data from GBIF and SpeciesDistributionToolkit.jl #note("Michael", "ADD SDT Cite.") and EvoTrees.jl #note("michael", "cite"). The 19 CHELSA bioclimatic variables were used as environmental predictors @Karger2017CliHig. [CITE]. For each species, pseudoabsences are generated using ”background thickening” (#note("michael", "cite"): Vollering et al. 2019), with no pseudoabsences allowed within a XXkm buffer of any GBIF occurrence record for that particular species. To infer the spatial distribution of each species, we used EvoTrees.jl to fit Boosted Regression Trees (BRT) with Gaussian loss metric. The use of this loss metric means each node in the classification tree is estimated as a Gaussian distribution via maximum-likelihood estimation, and therefore after fitting we have explicit variance around each split in each tree, and hence a total uncertainty value associated with each pixel can be constructed by aggregating the variance at each in the classification tree. SDMs were fit using an TODO-fold crossvalidation. The mean ROC-AUC across all SDM fits is XX TK.

= Results

== Species Richness over Space

== Interaction Richness 



= Discussion 

Predicting how interaction networks will change as species’ ranges shift is critical to the management of biodiversity, ecosystem processes, of ecological intactness, and the maintenance of ecosystem services. Our work highlights the potential of synthesizing “small” data from different sources to aid the prediction of interactions [@Todman2023SmaDat]. Further, by placing these interaction networks in a biogeographic context, we can project how these networks will change in the future, and crucially also spatially quantify our uncertainty in these projections. This provides guidance for where future monitoring efforts should target their effort to provide the most information possible on how these networks are changing, so we can detect and attribute the causes of this change [@Gonzalez2023FraDet], and make better decisions about how to manage and mitigate the consequences of this change in these systems [@Chapman2023BriAda]. This also leads to new avenues of research—can we predict the consequences of both extinction and formation of news interactions within species pollination networks, and if so, can we mitigate them? To do this, we should work toward synthesizing predictive models of network structure on biogeographic spatial scales and decadal temporal scales (as done here) with models of dynamics on mutualist networks that reflect processes on smaller scales [@Valdovinos2019MutNet].

---

*Acknowledgements*:  MDC is funded by an IVADO Postdoctoral Fellowship. TP is funded by an NSERC Discovery grant, a Discovery Acceleration Supplement grant, and a Wellcome Trust grant (223764/Z/21/Z).

#showbibliography("refs.bib") 
