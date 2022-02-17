# nH_selection
Daya Bay nH Analysis

This is the collection of codes used in the nH Analysis (inspired by the nGd analysis) at the Daya Bay Reactor Neutrino Experiment.

Codes included:
-findIBDs_2000.C: selection and plotting steps of the IBD selection
-findSingles_new.C: selection of single events, pairing into a synthetic accidentals background, and plotting
-finalize2000.C: adding up the plots from the runs done by findIBDs_2000.C
-finalize_singles2000.C: properly scaling the synthetic accidentals background for each run and adding the plots from findSingles_new.C
