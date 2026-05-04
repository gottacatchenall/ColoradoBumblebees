
#import "../template.typ": *
#show cite: set text(fill: rgb("#528fd1"))
#set figure(numbering: "S1")


#set table.hline(stroke: .6pt)
#show figure.caption: it => {
  context {
    set align(left)
    set par(leading: 0.4em, hanging-indent: 0pt, justify: false)
    (
      text(10pt, it.supplement + " " + it.counter.display() + "\n", weight: "semibold", font: "Libertinus Sans")
        + text(9pt, it.body, font: "Libertinus Sans", luma(20%))
    )
  }
}
#show figure.where(
  kind: table,
): set figure.caption(position: top)


#show heading.where(
  level: 1,
): it => block(width: 100%)[
  #v(1em)
  #set text(18pt, weight: 600, font: "Libertinus Sans", )
  // #smallcaps(it.body)
  #it.body
  #v(1em)
]

#show heading.where(
  level: 2,
): it => block(width: 100%)[
  #v(1em)
  #set text(14pt, weight: "semibold", fill: black.lighten(15%), font: "Libertinus Sans")
  #it.body
  #v(1.0em)
]


= Projecting the vertical disassembly of the bumble bee pollination network of the Southern Rocky Mountains --- Supplemental Material



#figure(
  image("../../plots/sdm_fit.png"),
  caption: ["Matthew’s Correlation Coefficient (MCC) values for the SDM fits for both bumble bees (left) and plants (right). The median MCC for bumble bees is 0.935, and 0.957 for plants."]
) <sdm-fit>


= Quantifying Intra-annual Variation Across Space


To quantify within-year variation in species abundance, we estimate species phenologies by taking all occurrence records from the iNaturalist data with a valid timestamp. Then, the total number of observations on each ordinal day of the year is used to fit a Bayesian Gaussian Mixture Model (GMM) to estimate the phenology for each species.We denote the total count of observations on a given ordinal day of the year as $y(t)$. This is modeled as

$
  y(t) ~ cal(N)(sum_(k=1)^K w_k exp(-(mu_k-t)^2/(2sigma_k^2)), epsilon)
$

where the sum over $k$ indicates each Gaussian (each called a _component_), and $cal(N)$ indicates a Normal distribution. Each of the $k$ components has a mean $mu_k$, a standard deviation $sigma_k$, and a weight $w_k$ associated with it, and $epsilon$ accounts for the overall noise causing points to deviate from the mixture's prediction. 


Both the number of observations and the ordinal day of year are scaled to $[0,1]$ to make sampling more efficient. 
We fit the GMM using the No U-Turn Sampling (NUTS; #cite(<Hoffman2014NouturnSampler>, form: "prose")) method for Hamiltonian Monte Carlo (HMC) in Turing.jl @Fjelde2025TuringjlGeneralpurpose, a framework for Bayesian inference in Julia, using the following priors: 

$
  mu_k ~ "Uniform"(0,1) \
  sigma_k ~ "Truncated"(cal(N)(1,1), [0.05, 1.5]) \ 
  w_k ~ "Truncated"(cal(N)(0,1), [0,1]) \
  epsilon ~ "Truncated"(cal(N)(0, 0.25), [0, infinity))
$


Note that this requires specifying a number of components $K$ --- because $K$ is discrete, it cannot be sampled directly with HMC. Instead, we fit models for K = {1,2,3} and the select the best model using the widely-applicable Akaike Information Criterion (wAIC, #cite(<Watanabe2010AsymptoticEquivalence>, form: "prose")). Fit phenologies with raw data are all visible in Appendix A3 for each species. 


#figure(
  image("../../plots/within_season.png"),
  caption: [Left: the weighted interaction richness for each month. Each color represents a 20% quantile of the number of interactions. White represents regions with 0 interactions. Top Right: the summed species richness for both bees (blue) and plants (green) across the season. Uncertainty bands represent 5-95% credible intervals sampled from HMC chains. Bottom Right: The average elevation of both bee (blue) and plant (green) species  weighted across months, computed in the same way as weighted richness across space. ]
) <within-season>


= Supplemental Gain and Loss Visualizations

#figure(
  image("../../plots/interaction_gain.png"),
  caption: [The number of gained interactions across space, meaning a pair of species known to interact now both occur at a location in space that they did not co-occur at under the baseline range predictions. Each shade of green represents a 20% quantile of lost interactions in the most extreme scenario, with the cutoffs between quantiles labeled on the color bar. The bar plot in the bottom right of each map is the proportion of area that falls into each quantile. The background (regions in white) have an expected gain of 0 interactions. Bottom Right: the mean number of gained interactions at each elevation in each of the climate scenarios (colors).]
) <gained-interactions>


#figure(
  image("../../plots/novel_cooccurrence.png"),
  caption: [Novel cooccurrence across space. Each heatmap shows quantiles for the number of novel co-occurrences among species that have only ever overlapped across less than 5% of the range (this value is arbitrary, but used because for species that overlap in a very small amount of area, it is unlikely we will have a detected an interaction between them even if it is feasible).]
) <novel_cooccurrence>

= Supplemental Winners and Losers Visualization 

#figure(
  image("../../plots/winners_and_losers.png"),
  caption: [Scatter plots for the amount of the median overlap with interaction partners relative to baseline (x-axis) and the relative range size (y-axis) for every species (points) across each climate scenario and time-period (panels). Point size is proportional to the total number of interactions that species has in the metaweb. Background colors in each region correspond to the status of each species as a relative winner/loser, or absolute winner/loser.]
) <winners-and-losers>


= Baseline Species Distributions 

Here are the baseline species range predictions for each species. Yellow points indicate occurrence records, and red points represent pseudoabsences.


#figure(
  image("../../plots/baseline_sdms/sdm1.png"),
  caption: [Part 1 of the baseline distribution models for each species]
) <sdm1>

#figure(
  image("../../plots/baseline_sdms/sdm2.png"),
  caption: [Part 2 of the baseline distribution models for each species]
) <sdm2>

#figure(
  image("../../plots/baseline_sdms/sdm3.png"),
  caption: [Part 3 of the baseline distribution models for each species]
) <sdm3>
#figure(
  image("../../plots/baseline_sdms/sdm3.png"),
  caption: [Part 4 of the baseline distribution models for each species]
) <sdm4>

= Range Shifts

*Note that this depicts _FIRST_ shift*. In some cases, a given regions is e.g. gained and then lost across the century. Interpret with this caveat.


#let plotfig(i) = {
  figure(
    image("../../plots/range_shifts/" + str(i) + ".png"),
    caption: [Species Range Shifts --- Part #i. Dark Grey background indicates the area is never part of the range.],
  ) 
}

#for i in range(1, 9) {
  plotfig(i)
}


#showbibliography("./refs.bib") 

