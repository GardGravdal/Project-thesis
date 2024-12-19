This repository contains all code used in my Project thesis (20/12/2024).

This file will provide a short description of the methods to derive the results from the project thesis, as well as various vizualisations.
Note that the raw datafiles from house sparrow (genetic-, morphological- and non-genetic data) is not provided, and will have to be saved 
locally to preform the analysis. Files are too large to be stored on GitHub.

Methods of deriving results:

Step XGBoost predictions:
1. dataloader.R -> Pre-processing of the phenotypes to compute the processed phenotype to be used in XGBoost.
2. CreateCV.ipynb -> Create the CV folds to be used for each trait (different for 70k and 180k). NB! Use the same folds in Bayesian animal model.
3. xgBoostOputna.ipynb -> Create XGBoost models, its predictions and correlation measure, which are saved locally.
4. get_boxplot.R -> Plot the correlation measures together with the correlations from GEBV.
5. Note that this procedure is the same for 180k and 70k datasets.

Step Bayesian animal model predictions:
1. Use the CV folds generated in Step XGBoost predictions in this step.
2. AnimalModel_INLA.R -> Generate predictions using INLA on LMM animal model and compute correlations. Store correlations locally.
3. get_boxplot.R -> Plot the correlation measures together with the correlations from XGBoost.

Step SHAP values:
1. This step assumes that XGBoost models from 70k set are generated for each (70k) fold following Step XGBoost predictions.
2. get_shapley.ipynb -> Generate SHAP values, create cumulative importance plot and save it, save mean abs SHAP values locally.
3. Chromosome2Manhattan.R -> Generate Manhattan-style plot for mean abs SHAP used in report. Also links SHAP vals to bp and chromosome.

Step GWAS (GEMMA):
1. qc.R -> Generate qc.raw files from 70k SNP data. The qc.raw files are used by GEMMA package later.
2. h_dataPrep.R -> At the end of this script, a GWAS_$phenotype.txt file is generated to be used later.
3. qc_gard -> Processes BIM,BAM,BED,FAM files and generates file inputGEMMA to be used by GEMMA program.
4. GEMMA-0.98.5 -> Run the "system()" lines in qc_gard in e.g. Command prompt. Make sure GEMMA-0.98.5 is stored locally.
5. GEMMA.R -> Uses the resulting p-values from output of GEMMA-0.98.5 to generate the Manhattan plot.

Step GWAS-SHAP correlation plot:
1. Assumes Step SHAP values and Step GWAS (GEMMA) already done.
2. CorrelationShapGWAS.R -> Generate the correlation plot for each trait.

Step Various visualizations:
1. The various scripts for illustrative visualizations (e.g. for the p-value) is found in folder Illustration_viz.

Note that a lot of the scripts need manual changing of dependencies, for example for which phenotype is analyzed or which dataset/CV folds used.
