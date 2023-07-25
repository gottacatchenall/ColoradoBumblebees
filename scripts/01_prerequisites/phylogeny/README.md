# Building Phylogenies

1. First use `beatui` to configure the settings used for phylogenetic inference.
   Do this by loading in the aligned sequence files in `data/aligned_sequences`
   and make the following configuration changes inside `beatui`: 
    - On the _Sites_ tab, change the Substitution Model to GTR, the Site Heteorgeneity Model to Gamma + Invariant Sites. 
    - On the _Trees_ tab, change the Tree Prior to a Yule Process.
    - Export the XML file with _Generate Beast File_

2. Run BEAST on the two XML configuration files for each species group.

3. Run TreeAnnotater on the `.trees` output file for the `BEAST` run to create a
   MCC tree.

4. Run the R script to convert the NEXUS output files from TreeAnotater into
   Newick format with the codenames replaced with the actual species names.
    - Make sure to use 1,000,000 trees and a burnin, e.g. via `treeannotate
      -burnin 1000000 bees.trees bee_tree.nxs` 
   
Then this can be used by both the embedding scripts and visualization scripts.