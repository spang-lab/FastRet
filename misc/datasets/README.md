# Files sorted by Category

<!-- no toc -->
## Full Training Datasets

- HILIC-Retip: availble through R function `FastRet::read_retip_hilic_data()`
- HILIC: [Excel](20210702_RT_Prediction_Hilic_Library_FF.xlsx)
- RP-AXMM: Exel missing
- RP: [Excel](RP.xlsx)

## Smaller Datasets for Model Adjustment

- Combined: Sheet 'R8_RT_Medoids' of [Excel](20211022_R8_dif_conditions_Medoids_validSet.xlsx)
- RP-Steep: [Excel](R8P_Medoids_input_Steeper gradient.xlsx)
- RP-Flat: [Excel](R8P_Medoids_input_Flatter gradient.xlsx)
- RP-T25: [Excel](R8P_Medoids_input_Lower Temperature 25.xlsx)
- RP-FR25: [Excel](R8P_Medoids_input_Lower flow rate 025.xlsx)
- RP-T25-FR25  [Excel](R8P_Medoids_input_Lower Temperature 25 + Lower Flow rate 025.xlsx)
- RP-T25-Fr25 Steep: [Excel](R8P_Medoids_input_Steeper gradient+lower Temperature 25+Lower Flow rate 025.xlsx)

## Validation Data Sets

- Combined: Sheet 'R8_RT_Validation_set' of [Excel](20211022_R8_dif_conditions_Medoids_validSet.xlsx)

## Validation Data Set Predictions

- RP-Steep-Val: [Excel](ValidSet_Steeper gradient.xlsx)
- RP-Flat-Val: [Excel](ValidSet_Flatter gradient.xlsx)
- RP-T25-Val: [Excel](ValidSet_Lower Temperature 25.xlsx)
- RP-FR25-Val: [Excel](ValidSet_Lower flow rate 025.xlsx)
- RP-T25-FR25-Val: [Excel](ValidSet_Lower Temperature 25 + Lower Flow rate 025.xlsx)
- RP-T25-Fr25-Steep-Val: [Excel](ValidSet_Steeper gradient+lower Temperature 25+Lower Flow rate 025.xlsx)

## Problems

- Number of metabolites given for HILIC is different between paper (364) and excel (392).
- Excel file for RP-AXMM column is not available.
- Number of metabolites given for RP is different between paper (401) and excel (442).

# Files sorted alphabetically

## 20210702_RT_Prediction_Hilic_Library_FF.xlsx

- Description: **In-House-Library Metabolites** measured on in-house HILIC column
- Dimensions: 393 x 5 (row one contains colnames)
- Columns:
    - Code: To be described by Fadi. Example values: e.g. '7B3' or '7E3'. (TODO: VERIFY)
    - NAME: Name of measured metabolite, e.g. 'MENAQUINONE' or 'PHYLLOQUINONE'
    - InChlKey: Not used
    - SMILES: SMILES representation of metabolite, e.g. 'C1CCCCC1' or 'C=1CC1'
    - RT: Measured retention times in minutes, e.g. '0.91' or '10.2'
- Keep: yes
- Symbol: -

## 20211022_R8_dif_conditions_Medoids_validSet.xlsx

SHEET1: R8_RT_Medoids

- Description: Retention times for **Medoid Metabolites** measured on the original RP column as well as all 6 adjusted reversed phase (ARP) columns
- Dimensions: 25 x 16 (row one contains colnames)
- Columns:
    - PLATE, NROW, NCOL, CODE, M+H, M-H: To be described by Fadi.
    - CNAME, SMILES, FORMULA, HMDB: Identifier for measured Molecule.
    - RT, Steeper, Flatter, T25, FR025, T25_FR025, T25_FR025_Steeper: Retention time measured on Reversed Phase column (see `20210702_RT_Prediction_Hilic_Library_FF.xlsx`) and reversed phase columns with adjusted chromatographic conditions.
- Keep: yes
- Symbol: -

SHEET2: R8_RT_Validation_set

- Description: Retention times for **Validation Set Metabolites** measured on the original RP column as well as all 6 adjusted reversed phase (ARP) columns
- Dimensions: 25 x 16 (row one contains colnames)
- Columns:
    - SMILES, Name: Identifiers for measured Molecule.
    - RT Normal, RT Steep, RT Flatter, T25, FR025, T25_FR025, T25_FR025_Steeper: Retention time measured on Reversed Phase column (see `20210702_RT_Prediction_Hilic_Library_FF.xlsx`) and reversed phase columns with adjusted chromatographic conditions.
- Keep: yes
- Symbol: -

## model-2022-03-10_xxx

- Files:
    - model-2022-03-10_Flatter gradient
    - model-2022-03-10_Lower flow rate 025
    - model-2022-03-10_Lower Temperature 25
    - model-2022-03-10_Lower Temperature 25 + Lower Flow rate 025
    - model-2022-03-10_Steeper gradient
    - model-2022-03-10_Steeper gradient+lower Temperature 25+Lower Flow rate 025
- Description: Models for adjusted for different chromatographic conditions. Adjustment done in 2022-03-10 using an early version of FastRet. Should not be used anymore.
- Keep: no
- Symbol: -

## predictor_set_2024-05-21.xlsx

- Description: 131 Chemical Descriptors obtained using rCDK for 439 SMILES from **In-House-Library**. Chemical Descriptors were obtained using an early version of FastRet (which used an early version of rCDK). Should not be used anymore.
- Dimensions: 440 x 132 (row one contains colnames, col one contains SMILES)
- Columns:
    - (1) Smiles
    - (2) RT: Measured Retention Time
    - (3-132): Chemical Descriptor Names, e.g. 'Fsp3', 'nSmallRings' or 'nAromRings'
- Keep: no
- Symbol: -

## R8P_Medoids_input_xxx

- Files:
    - R8P_Medoids_input_Flatter gradient.xlsx
    - R8P_Medoids_input_Lower flow rate 025.xlsx
    - R8P_Medoids_input_Lower Temperature 25 + Lower Flow rate 025.xlsx
    - R8P_Medoids_input_Lower Temperature 25.xlsx
    - R8P_Medoids_input_Steeper gradient.xlsx
    - R8P_Medoids_input_Steeper gradient+lower Temperature 25+Lower Flow rate 025.xlsx
- Description: TODO

## R8P_Validation_Set_input.xlsx

- Description: Names and Smiles for **Validation Set Metabolites**. This file can be uploaded to the FastRet Web Client to get predictions for the Validation Set. It was used in combination with the models from [model-2022-03-10_xxx](#model-2022-03-10_xxx) to produce the predictions in [ValidSet_xxx](#validset_xxx)
- Dimensions: 26 x 2 (row one contains colnames)
- Columns:
    - NAME, SMILES: Identifiers for measured Molecule.
- Keep: yes
- Symbol: -

## RP_adj.xlsx

- Description: Subset of 25 metabolites of `RP.xlsx` with artifically modified retention times. Used for testing of R functions.
- Dimensions: 26 x 3 (row one contains colnames)
- Columns: RT, NAME, SMILES (see `RP.xlsx`)
- Keep: no
- Symbol: -

## RP_broken.xlsx

- Description: Empty excel file. Used for testing of R functions.
- Dimensions: 0 x 0
- Columns: -
- Keep: no
- Symbol: -

## RP.xlsx

- Description: 442 metabolites measured on in-house HILIC column
- Dimensions: 443 x 3 (row one contains colnames)
- Columns:
    - RT: Measured retention times in minutes, e.g. '0.91' or '10.2'
    - NAME: Name of measured metabolite, e.g. 'MENAQUINONE' or 'PHYLLOQUINONE'
    - SMILES: SMILES representation of metabolite, e.g. 'C1CCCCC1' or 'C=1CC1'
- Keep: yes
- Symbol: `RP`

## ValidSet_xxx

- Files:
    - ValidSet_Flatter gradient.xlsx
    - ValidSet_Lower flow rate 025.xlsx
    - ValidSet_Lower Temperature 25 + Lower Flow rate 025.xlsx
    - ValidSet_Lower Temperature 25.xlsx
    - ValidSet_Steeper gradient.xlsx
    - ValidSet_Steeper gradient+lower Temperature 25+Lower Flow rate 025.xlsx
- Description: TODO
