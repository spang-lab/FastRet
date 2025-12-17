# Old Introduction

⚠️ **Deprecated:** This article was the original introduction shipped
with version 0.99.6 of FastRet. It is deprecated and should no longer be
used. It is retained here for reference purposes only.

Welcome to this introduction vignette for FastRet. This tutorial will
show you how to use FastRet and what to pay attention to when using it.

``` r
library(FastRet)
```

Before starting our analysis lets save some example data just as you
might have it laying around. The resulting excel file does need to have
three columns, RT, SMILES and NAME.

``` r
data(RP)
print(tibble::as_tibble(RP))
#> # A tibble: 442 × 3
#>       RT SMILES                                                    NAME         
#>    <dbl> <chr>                                                     <chr>        
#>  1  3.34 CCC(C(=O)O)O                                              2-HYDROXYBUT…
#>  2  3.35 COC1=C(C=CC(=C1)CCN)O                                     3-METHOXYTYR…
#>  3  2.11 C1=NC2=C(N1)C(=NC=N2)N                                    ADENINE      
#>  4  2.1  C1=NC2=C(C(=N1)N)N=CN2C3C(C(C(O3)COP(=O)(O)O)O)O          ADENOSINE 5'…
#>  5  3.13 C1C2C(C(C(O2)N3C=NC4=C3N=CN=C4N)O)OP(=O)(O1)O             ADENOSINE 3'…
#>  6  2.07 C1=NC2=C(C(=N1)N)N=CN2C3C(C(C(O3)COP(=O)(O)O)OP(=O)(O)O)O ADENOSINE 3'…
#>  7  1.22 [O-]C(=O)C[C@H](O)C[N+](C)(C)C                            L-CARNITINE  
#>  8  1.25 C[N+](C)(CCOP([O-])(OCC(O)CO)=O)C                         SN-GLYCERO-3…
#>  9  1.12 CN1C=C(N=C1)CC(C(=O)O)N                                   N(PAI)-METHY…
#> 10  1.21 C(CC(=O)O)CN                                              4-AMINOBUTAN…
#> # ℹ 432 more rows
openxlsx::write.xlsx(as.data.frame(RP), file="~/RP_data.xlsx", rowNames=FALSE)
```

Starting the interface is done with one line of code. From now on there
is no need to use R anymore.

``` r
FastRet()
```

All objects are additionally explained in more detail if you click on
the question mark.

![FastRet GUI main interface showing the train new model mode with file
upload section and model training options](Old-Intro/01.png)

FastRet GUI main interface showing the train new model mode with file
upload section and model training options

## Train new Model

Staying in the mode “Train new Model” we can upload our previously saved
excel file and start training the regression model. This might take some
time depending on the size of the training set. This is due to the
process of calculating chemical descriptors with rcdk. On the top right
you can see a loading symbol if the process is still running.

Once the training is complete, a boxplot with performance measures as
well as a scatterplot of predicted~actual retention time is shown. This
can be used to get a first impression of the performance that can be
expected.

![Boxplot showing cross-validation performance measures with R-squared
and error metrics for the trained retention time
model](Old-Intro/02.png)![Scatterplot of predicted versus actual
retention times showing model fit quality with diagonal reference
line](Old-Intro/03.png)

If you want to utilize the trained model on new data, you have to first
save the model and then reload it again in the “Utilize Model” category.
The model you download is trained on the whole data set without a
validation set. All performance measures of the boxplot are evaluated
using cross validation.

## Utilize Model to predict on new data

In this mode your previously saved model can be used to predict on new
data. This data can be either provided with the input of a single SMILES
or a whole .xlsx file which.

![FastRet GUI in utilize model mode showing options for uploading models
and entering SMILES for prediction](Old-Intro/04.png)

FastRet GUI in utilize model mode showing options for uploading models
and entering SMILES for prediction

### Adjusting measurements using a linear model

If you have measured some metabolites on your new experiment setup that
were also measured on the original setup, you can use this method to
adjust the models prediction with a linear model. Simply check “Use
measured metabolites to adjust Prediction” and upload an Excel file with
the corresponding SMILES/retention time data. After you click on
“Analyze Linear Model” a scatter plot will be shown, displaying the old
vs the new retention times of your measured metabolites. This plot helps
you decide which components of the linear model to select.

![Scatterplot showing correlation between old and new retention times
for model adjustment, with linear regression line](Old-Intro/05.png)

Scatterplot showing correlation between old and new retention times for
model adjustment, with linear regression line

## Selective Measuring

This mode calculates, for a given data set, the best k molecules to be
measured for a retention time prediction on a new experiment setup. It
uses a combination of Ridge Regression and k-means to determine the best
representatives of your dataset. Representatives as well as their
corresponding clusters can be downloaded afterwards as an excel file.
This step should be used once you have a predictive model and/or data
set and want to use it for new column/gradient/temperature… combination.

By providing an excel file with SMILES, NAME and RT and a custom k you
can simply start the process and export another excel file which will
show your clusters as well as your metabolites.
