# ehdn-cardio
Tandem repeat expansions in cardiomyopathy

This directory contains the code for downstream analysis of rare tandem repet expansions in patients with cardiomyopathy.

Prior to the steps described here, I use ExpansionHunterDenovo on WGS data, ExpansionHunterDenovo scripts to combine the output, and tandem-repeat-expansions-in-ASD to define rare expansion outliers. 

The scripts generate the results and plots as follows:

* motifDistribution.R - calculates and plots distribution of GC-composition and motif length for Figures 1A-D.

* expansionAnalysisFunctions.R - contains functions used for burden analysis and epigenetic marks

* correlationBins.R - calculates correlation between the number of rare TREs and different genomic features in 1 kb segments; generates Figure 1E. Output of the script is provided as Supplementary Table S6.
* sizePercentile.R - calculates distribution of repeat size in proband-parent pairs; generates Figure 1F. Repeat size calculated in this script and used to generate the plot is available as Supplementary Table S7.
* tssDistance.R - calculates distance from rare TREs to the nearest TSSs; generates Figure 3A. 
* burdenAnalysis.R - calculates burden in genic regions; generates Figure 2A. Output of the script is provided as Supplementary Table S8.
* constraintScores.R - calculates constraint scores of genes with rare TREs; generates Figure 2B.
* promoterEpigenetic.R - calculates burden in regions overlapping promoter and epigenetic marks; generates Figure 3B. Output of the script is provided as Supplementary Table S9.
* testFactors.R - calculates the number of rare TREs with respect to clinical variables; generates Figure 3C. Output of the script is provided as Supplementary Table S10.
