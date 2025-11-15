#import "template.typ": *
#import "@preview/wordometer:0.1.4": word-count, total-words
#import "authors.typ": *
//#show: doc => style(doc)

#let title = "Projecting the spatial rewiring of the bumble bee pollination network of the Southern Rocky Mountains"
#let abstract = "TK"


#showtitle(title)
#showauthors(authors, institutions, corresponding_email)
#showabstract(abstract)
#show: word-count.with(exclude: (figure, bibliography))
#set heading(numbering: "1.")
#show heading: it => {
  [#it]
  v(0em)
}
#show cite: set text(fill: rgb("#528fd1"))

\
*Word Count*: #total-words
\

_Ecology Letters: 5000 words, 6 figures._

= Introduction

The vast majority of flowering plant species are pollinated by animals @Ollerton2011HowMan. Pollination interactions are essential to plant reproduction, and necessary for the persistence of Earth’s biodiversity @Bascompte2007PlaMut. Pollinatinator dependent crops contribute to 35% of global crop production by volume and bumble bees commonly pollinate over 25 crops worldwide. Recent models indicate 3%–5% of fruit, vegetable, and nut production is lost due to inadequate pollination, leading to an estimated 427,000 excess deaths annually from lost healthy food consumption and associated diseases @Klein2006ImpPol @Smith2022PolDef. But, the human driven forces of climate and land-use change are altering plant-pollinator interactions, and therefore threatening the stability and function of these ecosystems @Kearns1998EndMut As a result, understanding and predicting how plant-pollinator networks are changing due to anthropogenic pressures is a major question in pollinator research @Mayer2011PolEco and conservation biology @Kearns1998EndMut, and ultimately necessary for climate change adaptation and maintaining food security @Bailes2015HowCana @Klein2006ImpPol.

A large amount of research into the effects of climate change on pollinator interactions has focused on phenology @Ogilvie2017IntBee, as shifts in flowering and hatching times caused by warm temperatures move earlier and earlier in the season. Yet this is not the only factor that will govern changes in plant-pollinator systems in the future—--in particular, the spatial distribution of where particular plant and pollinator species can feasibly live is changing, from small (habitat fragmentation) to large (range shift) spatial scales because of climatic niche shifts. Mountainous ecosystems face a unique threat in this respect. Elevational gradients in niche conditions are being pressured toward higher eleveations as average temperatures increase @Hoiss2015IntEff @Richman2020AsyRan. As these species’ ranges shift, these networks are becoming “rewired”, where some interactions are lost (as the participating species cease to coexist in the future) and others are gained among species as they start to co-exist for the first time.

Species interactions intrinsically vary in space and time @Poisot2015SpeWhy —–- in addition to systematic phenological shifts, there is intrinsic year-to-year variation in the structure of empirical pollination networks @Alarcon2008YeaVar ---so we need a spatiotemporally explicit approach to understand interaction networks and how they are changing @Strydom2021RoaPre. Going further, forecasting the rewiring of networks as their species shift under different climate projections (and understanding how this rewiring affects ecosystem functioning, persistence, and the services ecosystems provide to people) is crucial to conservation strategies designed for uncertain climate futures. Still, spatially explicit predictions of network structure remain rare because thorough spatial coverage of network composition and structure is very labor intensive. However, spatial network predictions can be made by combining a regional metaweb (all possible interactions across the regional species pool) and using species distribution models to predict host composition across space @Dansereau2024SpatiallyExplicit. 

Here, we follow this approach using interaction data from *XX* bumble bee (_Bombus spp._) and *XXX* angiosperm species from the Southern Rocky Mountains of Colorado to project how this network of 175 species will become “rewired” as each species’ environmental niche shifts over the coming century. We use this to (1) produce spatiotemporally explicit predictions of the bumblebee/wildflower pollination network across the southern Rocky Mountains and (2) project how much the ranges of each bee and plant involved in each interaction will overlap in the future under different climate scenarios. We hope these projections of rewiring can guide future research and monitoring to understand how the interactions expected to have the most extreme change will impact ecosystem functioning and services.

= Methods 

*CONCEPT FIG HERE*

== Data 

=== Interaction Data

The primary data source for this study is interaction data collected across three different sites in the subalpine Rocky Mountains of Colorado across a total of XX Bombus species and XXX flowering plants. The three sample sites are henceforth referred to as Rocky Mountain Biological Lab (RMBL), Pikes Peak, and University of Colorado Mountain Research Station (MRS) (visible in fig. 2.1, upper left). Observations of bumble bee-plant interactions at the RMBL take place within a long-term bumble bee monitoring project  @Ogilvie2022ShiftingImportance. In brief, at six permanent study sites near the RMBL, bumble bee interactions with flowering plants are monitored for one hour at weekly intervals for the entire growing season (May–September). At these same weekly intervals, floral abundances of all flowering plant species visited by bumble bees is surveyed within a series of 15 permanent 20 x 0.5 m transects at each site. Observations of bees, their interactions with flowers, and counts of flowers take place within the three dominant habitat types where bumble bee activity is observed: wet meadows, dry meadows, and aspen forest. We used seven years (2015-2021) of observational data on Bombus visits to flowering plants at these study sites to aggregate interaction networks. The Pikes Peak data consists of three seasons (2019-2021) of interaction data collected along an elevational gradient ranging from the near the summit of Pikes Peak (4285 m) to down near the transition of the Great Plains (1872 m), sampled roughly twice a week (described in #cite(<Barthell2025BumbleBee>, form: "prose")). MRS data consists of seven seasons (2015-2021) of data at a six plot at 2900 meters elevation, sampled weekly over the entire flowering season (described in #cite(<Resasco2021PlaPol>, form: "prose")).

=== Occurrence Data

Occurrence data for each species was curated from GBIF within the region of study. The dataset of all occurrence records can be found *here*[DOI]. Each occurrence record is associated with a geospatial coordinate and a timestamp. The GBIF dataset consists of XX occurrences from human observation within the bounding box seen FIG XX (BBOX HERE).

=== Climate Data

CHELSA provides 19 bioclimatic variables dating back to 1970 at 1km resolution across Earth’s surfaces @Karger2017CliHig. CHELSA also provides these same variables projected into the future until the end of this century. These projections of bioclimatic variables, and the climate projections they are built on, rely on the framework of Shared-Socioeconomic-Pathways @ONeill2014NewSce to describe different scenarios in how humanity responds to climate change. These scenarios of climate response vary on two axes: mitigation and adaptation—how humanity will mitigate climate change (i.e. by reducing and eventually reaching net-zero greenhouse gas emissions) and how we will adapt to climate change (i.e. building infrastructure for a warmer world). 
We consider three climate scenarios: SSP1-2.6, SSP2-4.5, SSP3-7.0, representing low, moderate, and extreme warming, respectively. Projections for the global climate for each of these situations are derived from five different Earth System Models (ESMs) --- GFDL-ESM4 @Held2019StrPer, IPSL_CM6A_LR @Boucher2020PresentationEvaluation, MPI_ESM1_2_HR @Muller2018HigherresolutionVersion, MRI_ESM2_0 @Yukimoto2019MeteorologicalResearch, and UKESM1_0_LL @Good2019MohcUkesm10ll --- which are then translated into bioclimatic variable predictions by CHELSA.

== Species Distribution Models 

Species distributions were fit using occurrence data from GBIF and SpeciesDistributionToolkit.jl @Poisot2025JuliaToolkit and EvoTrees.jl @Desgagne-Bouchard2025EvovestEvotreesjl. The 19 CHELSA bioclimatic variables were used as environmental predictors @Karger2017CliHig. For each species, pseudoabsences are generated using "background thickening” @Vollering2019BunchingBackground --- where points are selected in proportion to their minimum distance to an observed prsences --- with no pseudoabsences allowed within a XXkm buffer of any GBIF occurrence record for that particular species. To infer the spatial distribution of each species, we used EvoTrees.jl to fit Boosted Regression Trees (BRT) with Gaussian loss metric. The use of this loss metric means each node in the classification tree is estimated as a Gaussian distribution via maximum-likelihood estimation, and therefore after fitting we have explicit variance around each split in each tree, and hence a total uncertainty value associated with each pixel can be constructed by aggregating the variance at each in the classification tree. 

SDMs were fit using an TODO-fold crossvalidation. Thresholding was done by selecting the value that maximizes Matthew's Correlation Coefficient (MCC). The mean ROC-AUC across all SDM fits is XX TK, the mean value of MCC is XX, and the distribution of both MCC and ROC-AUC for each species are visible in Appendix A4. 

Individual distribution models for each species are available in Appendix A1, and visualization of projected shifts for each species are in Appendix A2. 

== Phenology Estimation

To estimate species phenologies, we take all GBIF occurrence records with a valid timestamp. Then, the total number of observations on each ordinal day of the year is used to fit a Bayesian Gaussian Mixture Model (GMM) to estimate the phenology for each species. 


We denote the total count of observations on a given ordinal day of the year as $y(t)$. This is modeled as 
$
  y(t) ~ cal(N)(sum_(k=1)^K w_k exp(-(mu_k-t)^2/(2sigma_k^2)), epsilon)
$

where the sum over $k$ indicates each Gaussian (each called a _component_), and $cal(N)$ indicates a Normal distribution. Each of the $k$ components has a mean $mu_k$, a standard deviation $sigma_k$, and a weight $w_k$ associated with it, and $epsilon$ accounts for the overall noise causing points to deviate from the mixture's prediction. 

Both the number of observations and the ordinal day of year are scaled to $[0,1]$ to make sampling more efficient. 
We fit the GMM using the No U-Turn Sampling (NUTS; #cite(<Hoffman2014NouturnSampler>, form: "prose")) method for Hamiltonian Monte Carlo (HMC) in Turing.jl @Fjelde2025TuringjlGeneralpurpose, a framework for Bayesian inference in Julia, using the following priors: 

#note("Michael", [*TODO* refit with different $mu_k$ prior])
$
  mu_k ~ "Uniform"(0,1) \
  sigma_k ~ "Truncated"(cal(N)(1,1), [0.05, 1.5]) \ 
  w_k ~ "Truncated"(cal(N)(0,1), [0,1]) \
  epsilon ~ "Truncated"(cal(N)(0, 0.25), [0, infinity))
$

Note that this requires specifying a number of components $K$ --- because $K$ is discrete, it cannot be sampled directly with HMC. Instead, we fit models for K = {1,2,3} and the select the best model using the widely-applicable Akaike Information Criterion (wAIC, #cite(<Watanabe2010AsymptoticEquivalence>, form: "prose")). Fit phenologies with raw data are all visible in Appendix A3 for each species. 


== Rewiring Quantification

In order to quantify the amount of rewiring in the network over the rest of this century, for each interaction in the metaweb we take the baseline distribution for each bee and plant species and compute the amount of area in which their ranges overlap. Then, for each future time-period and climate change scenario, we compute the overlap in the same way, and compute the ratio between the overlap in that time-period and climate scenario and the baseline overlap.

= Results

The predicted species richness and aggregated uncertainty is shown in @richness-uncertainty. 

#figure(
  image("../plots/richness_vs_uncertainty.png"),
  caption: []
) <richness-uncertainty>


== Intra-annual Dynamics

#figure(
  image("../plots/within_season.png"),
  caption: []
) <within-season>



== Range Shift Projections

The projections of gained and lost range for each individual species are contained in Appendix XX. In @range-shift, we see the median elevation shift per time period for each of the climate scnarios.  

#figure(
  image("../plots/range_shift.png"),
  caption: []
) <range-shift>

== Network Rewiring Projections

#figure(
  image("../plots/interaction_loss.png"),
  caption: []
) <overlap-change>

#figure(
  image("../plots/overlap_change.png"),
  caption: []
) <overlap-change>


= Discussion 

By combining in-situ interaction data with community science occurrence records and geospatial projections of bioclimatic variables, we are able to predict the ranges and phenologies of _Bombus_ species and the plants they pollinate in the present, and project the future range shifts. We consistent upward elevational range shifts, and dismantleing of interaction networks over the coming century, with nearly complete collapse by 2100 occurring in the most extreme climate scenario.  

By projecting the distribution of each species, we have an estimate of the particular (bee, plant) interactions that will undergo the most change. This enables a prioritization of interactions likely to affect pollination and community outcomes in the future . This finding is of significant applied importance for the management of agricultural systems, and introduces many questions for future research. In the future, we will use our best projections of the future of a network’s (metaweb) structure and functioning arising from shifting bioclimatic and biogeographical conditions to prioritize conservation action. The ecological consequences of the expected rewiring could be tested in experimental systems in order to best understand the consequences of this rewiring on mutualist networks in the field. We do not envision this work as providing definitive forecasts of ecological networks under change. Instead we believe these results provide a priority-list for what species and particular interactions are expected to undergo the most change, so we can better focus network monitoring to anticipate change and guide conservation action.

Predicting how interaction networks will change as species’ ranges shift is critical to the management of biodiversity, ecosystem processes, of ecological intactness, and the maintenance of ecosystem services. Our work highlights the potential of synthesizing “small” data from different sources to aid the prediction of interactions @Todman2023SmaDat. Further, by placing these interaction networks in a biogeographic context, we can project how these networks will change in the future, and crucially also spatially quantify our uncertainty in these projections. This provides guidance for where future monitoring efforts should target their effort to provide the most information possible on how these networks are changing, so we can detect and attribute the causes of this change @Gonzalez2023FraDet, and make better decisions about how to manage and mitigate the consequences of this change in these systems @Chapman2023BriAda. This also leads to new avenues of research—can we predict the consequences of both extinction and formation of news interactions within species pollination networks, and if so, can we mitigate them? To do this, we should work toward synthesizing predictive models of network structure on biogeographic spatial scales and decadal temporal scales (as done here) with models of dynamics on mutualist networks that reflect processes on smaller scales @Valdovinos2019MutNet.


---

*Acknowledgements*:  MDC is funded by an IVADO Postdoctoral Fellowship. TP is funded by an NSERC Discovery grant, a Discovery Acceleration Supplement grant, and a Wellcome Trust grant (223764/Z/21/Z).

#showbibliography("refs.bib") 
