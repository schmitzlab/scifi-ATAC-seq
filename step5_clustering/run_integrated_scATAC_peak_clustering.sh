#!/bin/bash
rds=./scifi-merged.peak.raw.soc.rds
meta=./maize_scifi_BayesDoublet.singlet.meta.txt
feature_rate=1
res=0.5
PC=25
name=NAM_8libs


# run
Rscript SVD_analysis_scATAC.integrated.R $rds $meta $feature_rate $res $PC ${name}_${feature_rate}feature_${res}res_${PC}PCs
