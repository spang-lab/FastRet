# Extract predictor names from an 'frm' object

Extracts the predictor names from an 'frm' object.

## Usage

``` r
get_predictors(frm, base = TRUE, adjust = FALSE)
```

## Arguments

- frm:

  An object of class 'frm' from which to extract the predictor names.

- base:

  Logical indicating whether to include base model predictors.

- adjust:

  Logical indicating whether to include adjustment model predictors.

## Value

A character vector with the predictor names.

## Examples

``` r
frm <- read_rp_lasso_model_rds()
get_predictors(frm)
#>   [1] "Fsp3"             "nSmallRings"      "nAromRings"      
#>   [4] "nRingBlocks"      "nAromBlocks"      "nRings5"         
#>   [7] "nRings6"          "tpsaEfficiency"   "Zagreb"          
#>  [10] "XLogP"            "WPATH"            "WPOL"            
#>  [13] "WTPT.1"           "WTPT.2"           "WTPT.3"          
#>  [16] "WTPT.4"           "WTPT.5"           "MW"              
#>  [19] "VAdjMat"          "TopoPSA"          "LipinskiFailures"
#>  [22] "nRotB"            "topoShape"        "PetitjeanNumber" 
#>  [25] "MDEC.12"          "MDEC.13"          "MDEC.22"         
#>  [28] "MDEC.23"          "MDEC.33"          "MDEO.11"         
#>  [31] "MDEO.12"          "MDEN.13"          "MDEN.22"         
#>  [34] "MDEN.23"          "MLogP"            "nAtomP"          
#>  [37] "nAtomLC"          "khs.sCH3"         "khs.ssCH2"       
#>  [40] "khs.dsCH"         "khs.aaCH"         "khs.sssCH"       
#>  [43] "khs.dssC"         "khs.aasC"         "khs.aaaC"        
#>  [46] "khs.sNH2"         "khs.ssNH"         "khs.aaNH"        
#>  [49] "khs.aaN"          "khs.aasN"         "khs.sOH"         
#>  [52] "khs.dO"           "khs.ssO"          "khs.dsssP"       
#>  [55] "Kier1"            "Kier2"            "HybRatio"        
#>  [58] "nHBDon"           "nHBAcc"           "fragC"           
#>  [61] "FMF"              "ECCEN"            "SP.0"            
#>  [64] "SP.1"             "SP.2"             "SP.3"            
#>  [67] "SP.4"             "SP.5"             "SP.6"            
#>  [70] "SP.7"             "VP.1"             "VP.2"            
#>  [73] "VP.3"             "VP.4"             "VP.5"            
#>  [76] "VP.6"             "VP.7"             "SPC.4"           
#>  [79] "SPC.5"            "SPC.6"            "VPC.4"           
#>  [82] "VPC.5"            "VPC.6"            "SC.3"            
#>  [85] "SC.4"             "SC.5"             "VC.3"            
#>  [88] "VC.4"             "VC.5"             "SCH.5"           
#>  [91] "SCH.6"            "SCH.7"            "VCH.5"           
#>  [94] "VCH.6"            "VCH.7"            "C1SP2"           
#>  [97] "C2SP2"            "C3SP2"            "C1SP3"           
#> [100] "C2SP3"            "bpol"             "nB"              
#> [103] "nBase"            "ATSp1"            "ATSp2"           
#> [106] "ATSp3"            "ATSp4"            "ATSp5"           
#> [109] "ATSm1"            "ATSm2"            "ATSm3"           
#> [112] "ATSm4"            "ATSm5"            "nAtom"           
#> [115] "nAromBond"        "naAromAtom"       "apol"            
#> [118] "ALogP"            "ALogp2"           "AMR"             
#> [121] "nAcid"            "nA"               "nG"              
```
