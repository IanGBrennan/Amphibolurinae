# Scripts

Basic explanation of which scripts do what to reanalyze data from *Generalists link peaks in the shifting adaptive landscape of Australia's dragon lizards*. Scripts are organized in to X directories based on their utility. 

## Tree Visualization

**CFs_ggtree.R**: plots the concordance factors as pie charts on the genus-level tree of amphibolurines using *ggtree*. 

**MCMCTree_PriorPosterior.R**: Functions for processing MCMCTree mcmc files to combine and plot data from prior, effective prior, and posterior information. Generates plot in Fig.S1. 

**plotTree_arc.R**: plot the arc style tree in Fig.1



## Morphological Analysis

**00_Data_Preparation.R**: General data handling and processing to prepare morphological and phylogenetic data for analysis. Outputs include:

+ Data/Amphibolurinae_Morphology_AllSamples.csv
+ Data/Amphibolurinae_Morphology_spMEANS.csv
+ Data/Amphibolurinae_AllLSR.csv
+ Trees/Amphibolurinae_forR.tre
+ Data/Amphibolurinae_Data.RData

**01_Process_BayesTraits_Output_MultiChain.R**: Script to process fabric results. A version of this script is also available for a single chain *01_Process_BayesTraits_Output_SingleChain.R*. Output is the ancestral trait estimates:

+ Data/Ancestral_Trait_Estimates.RData

**02_fabric.R**: Process and visualize outputs of the fabric model applied to a number of traits. We can visualize directional trends (b) and evolvability (v). Output:

+ Data/Amphibolurinae_fabric_AllChainsSummary.RData

**03_Regimes.R**: Identify morphological regimes using *PhylogeneticEM*. We will extract these state codings for species and apply them to *RandomForests* next. 

**04_randomForests.R**: Use the morphological data and *PhyloEM* regime codings to generate a random forest decision tree to identify morphological regimes. Then apply this to the ancestral trait values to determine the appropriate regime of ancestral taxa.

**05_distance.to.BM.R**: Compare the amount of morphological change occurring on each branch of the tree to a null model (Brownian Motion). This helps visualize periods in time where morphological change was greater than expected. See Fig.3

**06_CentroidDistances.R**: Estimate the morphological distance between each tip/node and the centroid of each genus. 

**07_funspace.R**: Using morphological data to visualize functional spaces using *funspace*. Plots a functional space like in the inset PCA of Fig.1. This script also includes functions to visualize the changes occurring to individual traits along specified branches as seen in Fig.4.

## Morphological Visualization

**CorrPlot.R**: plot the correlations among morphological traits. See Fig.S3

**compare2centroid.R**: visualize the morphological distance between each species and the centroid of all regimes and the amphibolurine MRCA. See Fig.S12

**heatmap.R**: plot the multivariate euclidean morphological distance between all pairs of species. See Fig.S11

**trait_variance.R**: estimate and plot the variance of each individual morphological trait.
 



## Analysis Helper Scripts

**fabric_ControlFile.R**: Generate control files for BayesTraits. Includes the ability to estimate reasonable priors on the root value (alpha) and evolutionary rate (sigma) using *geiger*. 

**fabric_Functions.R**: A number of helper functions for processing BayesTraits *fabric* model outputs.

+ *fabricControl()* will generate a control file for a given trait, including making tags in the control file to track ancestral trait values.
+ *fabricAncestors()* will summarize ancestral trait values for all nodes as estimated by the *fabric* model and tracked through the control file. Might need to adjust the *log.in* line to skip the necessary number of lines in your mcmc file.
+ *inspectBT.log()* will summarize ESS values and optionally return a warning if values are below a minimum threshold.
+ *fabricMD5()* will create a reference key to match the MD5SUM codes from BayesTraits to node numbers in your input tree.
+ *fabricProcess()* helps to process the output of *fabric* runs to extract significant b/v values
+ *fabricProcessMerge()* helps to process the output of Post-Processed *fabric* files to extract significant b/v values. 
+ *fabricPlot()* allows the user to plot shift (b/v) values quickly. 
+ *fabricPlotVB()* allows the user to plot shift (b/v) values in a prettier way.
+ *fabricSuper()* will plot a summary of all shifts across multiple traits and shows shifts in b/v on separate phylogenies.
+ *fabricTraitgram()* generate a traitgram (a la *phytools*) for a given trait to track its evolution through time along the phylogeny. Requires ancestral trait data.
+ *fabricPathgram()* generate a biplot of the path from root to a specified tip and annotate directional trends (arrows) and evolvability changes (colors)
+ *plot.treePath()* plot a phylogeny with the path to a given tip or clade highlighted
+ *getEdges()* this function pulls out all descendant edges from a provided edge or node number

**BayesTraits_ModelComparison.R**: Summarize the fit from a number of BayesTraits models for comparison.

**TimeTree_Additions.R**: script to add three *Hypsilurus* species to our timetree for analysis.

**trait.at.time.R**: functions to extract trait characteristics (variance, range, volume) through time. 

**Visualize_Locus_Sampling.R**: visualize the molecular sampling for each sample, split by locus type (AHE, UCE, gene). 