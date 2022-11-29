# Lp_rot1129
This is a repository for conducting experiments in Rotation to Sparse Loadings using Lp Losses and Related Inference Problems http://arxiv.org/abs/2206.02263.

*cl1121_Lp_fun.R* contains all the auxiliary functions, such as IRGP algorithm, proximal gradient descent algorithm, etc. It is the source for all experiments in folders Study I, study II and Big 5 Personality Test. 

## Study 1
Table 1, 2 and Figure 4 are produced mainly by *cl1121A_experiment.15S.R*, cl1121A_experiment.30S.R*, where the 15S, 30S indicating the setting and the two files only differ in parameters, not in experiments. 

The initial results are stored in Rdata files and processed by *res.outputfinal1107.Rmd*.

To draw the U-shape curve in Figure 3, we first generate Lasso path solutions by *ushape1011upB.15S.R*, *ushape1011upB.30S.R* and save the results in Rdata form. Then, we use *draw.R* to draw the figure.

## Study 2
Table 3 is produced by *counter1121.R*. The initial results are stored in Rdata files and processed by *res.outputfinal1114.Rmd*

## Big-Five Personality Test
Table 4-7 in the main test and K1-K3 in the supplement are produced by *big5_r_code_1115sgn.Rmd*, which is a rmarkdown file.

## Supplement.H
Here are the files that compare the performance of Lp rotation to other rotation criteria. Table H.1 and Table H.2 are produced by *geomin_1128_chart.R*. Table H.3 is produced by *res.outputfinal1107.Rmd* using the Rdata results from *cl1123A_experiment.21S.R*. The minimizers of the L1 and L0.5 rotation is checked by grid search by *geomin_1128_gridsearch.R*

## Rdata
This folder contains all raw results in the Rdata form which will be further processed by *res.outputfinal1107.Rmd* and *res.outputfinal1114.Rmd*. It also includes the big 5 dataset.

## Tables and Figures
This folder contains all Table and Figures in the main text.
