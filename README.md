# pipe-ABC

A group repository containing all 4 elements of our project at Imperial College London, completed in 2021.

## Characterizing inter-species kinetics during neuronal differentiation

Authors: Nicholas McQuibban [NM], Abhranil Maiti [AM], Joshua Mayne [JM], Maisey Davidson [MD]
Supervisor: Dr. Ruben Perez-Carrasco

### Abstract

#### Background and Aims: 

The gene regulatory network comprising of Sonic Hedgehog (Shh), Gli, Nkx2.2, Olig2, Irx3, and Pax6controls motor neuron (MN) differentiation during embryonic development and occurs approximately 2.5 times faster in mice than in humans. A previous study attributed this to differences in protein degradation and defined a mathematical model to describe the network. The aim of this study was to re-analyse raw RNA-seq reads from this study; characterise a mathematical relationship between rate of MN differentiation and development stage (time-factor); build and infer the parameters of a mathematical model describing the gene regulatory network.

#### Methods: 

Raw RNA-seq reads for 81 human and mouse samples at various time-points (n=3) were processed using FastQC, STAR, RSEM/StringTie, and DESeq2; Variance-Stabilised-Normalisation (VST) was compared to Transcripts per Million (TPM) using multiple t-tests with Bonferroni correction. Batch effect was assessed using principal component analysis. Regression analysis was carried out on time versus development stage, and the equation derivatives were divided by each other to retrieve the time-factor. Hill Equations were used as the new model and parameterised using Approximate Bayesian Computation Sequential Monte Carlo (ABC SMC). Clustering of similar gene- patterns was carried out using Dirichlet Process Gaussian Process Clustering.

#### Key Results: 

No difference was found between VST-normalisation and TPMs for most data-points (α=0.05). Time-factors ranged from 2.68-7.57 depending on development stage. The ABC SMC algorithm was 7-fold faster than the simple ABC, while parallel multiprocessing further sped it up 8- fold. The proposed Hill function model of mRNA and protein dynamics best reproduced experimental data and the protein production and degradation rates obtained for this model were smaller in humans than mice by 2.5-fold. Ubiquilin 1 was found as a candidate for further investigation into protein degradation differences between human and mouse MN differentiation.

#### Conclusion: 

Following our experimental work, it can be ascertained that the differences in protein stability between species, such as humans and mouse, may depend on both protein production and degradation. Further, it is suggested that protein degradation pathways between these two species could differ significantly, due to findings highlighting inverse levels of expression in orthologous degradation genes between the species. Overall, more work is required to determine the extent of the aforementioned variation in protein degradation pathways.

# Repo Contents

This repository contains all scripts which were executed to run each element of this project. In addition, certain datasets and results relevant to the scripts and project can also be found in their respective directories.

All branches except the ABC section uploaded their scripts in Jupyter Notebook format, using R, Python3 and Bash Kernels to display said code, allowing for easy access and understanding of the various tools at work.

## Logical Order of Contents

1. RNASeqPipeline
2. RNASeqAnalysis
3. ABC
4. Clustering

