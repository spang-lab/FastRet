# Raw Files sorted alphabetically

## RP.xlsx

- Link: [RP.xlsx]
- Description: 442 metabolites measured on in-house Reverse-Phase (RP) column
- Dimensions: 443 x 3 (row one contains colnames)
- Columns:
    - RT: Measured retention times in minutes, e.g. '0.91' or '10.2'
    - NAME: Name of measured metabolite, e.g. 'MENAQUINONE' or 'PHYLLOQUINONE'
    - SMILES: SMILES representation of metabolite, e.g. 'C1CCCCC1' or 'C=1CC1'
- Keep: yes
- Symbol: `RP`

## RP-AXMM_FastRet_Input.xlsx

- Link: [RP-AXMM_FastRet_Input.xlsx]
- Description: 438 metabolites measured on in-house Reverse-Phase Anion-Exchange Mixed-Mode (RP-AXMM) Column.
- Dimensions: 439 x 3 (row one contains colnames)
- Columns: RT, NAME, SMILES (like in [RP.xlsx])
- Keep: yes
- Symbol: RP_AXMM

## 20210702_RT_Prediction_Hilic_Library_FF.xlsx

- Link: [20210702_RT_Prediction_Hilic_Library_FF.xlsx]
- Description: 392 Metabolites measured on in-house HILIC column
- Dimensions: 393 x 5 (row one contains colnames)
- Columns:
    - Code: Specifies position on "vendor plate" (i.e. the plate that was sent by the vendor where molecules were ordered). Example values: e.g. '7B3' or '7E3'.
    - NAME: Name of measured metabolite, e.g. 'MENAQUINONE' or 'PHYLLOQUINONE'
    - InChlKey: Not used, i.e. all entries are NA.
    - SMILES: SMILES representation of metabolite, e.g. 'C1CCCCC1' or 'C=1CC1'
    - RT: Measured retention times in minutes, e.g. '0.91' or '10.2'
- Keep: yes
- Symbol: `HILIC`

## 20211022_R8_dif_conditions_Medoids_validSet.xlsx > R8_RT_Medoids

- Link: [20211022_R8_dif_conditions_Medoids_validSet.xlsx]
- Description: Retention times for "Medoid Metabolites" measured on the original RP column under six different chromatographic conditions
- Dimensions: 25 x 16 (row one contains colnames)
- Columns:
    - PLATE, NROW, NCOL, CODE, M+H, M-H: Specifies position on "vendor plate".
    - CNAME, SMILES, FORMULA, HMDB: Identifier for measured Molecule.
    - RT, Steeper, Flatter, T25, FR025, T25_FR025, T25_FR025_Steeper: Retention time of medoid subset measured on Reversed Phase column (see [RP.xlsx]) under varying chromatographic conditions.
- Keep: yes
- Symbols: `RP_Mod`

## 20211022_R8_dif_conditions_Medoids_validSet.xlsx > R8_RT_Validation set

- Link: [20211022_R8_dif_conditions_Medoids_validSet.xlsx]
- Description: Retention times for "Validation Set Metabolites" measured on the original RP column under the original chromatographic conditions as well as six different modified chromatographic conditions
- Dimensions: 25 x 16 (row one contains colnames)
- Columns:
    - SMILES, Name: Identifiers for measured Molecule.
    - RT Normal, RT Steep, RT Flatter, T25, FR025, T25_FR025, T25_FR025_Steeper: Retention time of validation set metabolites measured on Reversed Phase column (see [RP.xlsx] under varying chromatographic conditions.
- Symbol: `RP_Val`

## 20250623_IROA_Library.xlsx

- Link: [20250623_IROA_Library.xlsx]
- Description: List of Metabolites in "In House Library" of Institute of functional Genomics. Created manually.
- Dimensions: 25 x 16 (row one contains colnames)
- Columns:
    - PLATE: Plate identifier of "vendor plate" where metabolite is located.
    - NROW: Row number on the plate.
    - NCOL: Column number on the plate.
    - CNAME: Compound name.
    - SMILES: SMILES representation of the metabolite structure.
    - FORMULA: Chemical formula of the metabolite.
    - MW (as put into the well): Molecular weight as used for well preparation.
    - µmol in well: Amount (micromoles) of compound in the well.
    - µM in 200 uL: Concentration (micromolar) in 200 microliters.
    - Stock dil uM: Stock dilution concentration in micromolar.
    - Formula (as seen in MS): Chemical formula as detected in mass spectrometry.
    - Neutral mass: Neutral (uncharged) mass of the compound.
    - M+H: Mass-to-charge ratio for the protonated molecule ([M+H]+).
    - M-H: Mass-to-charge ratio for the deprotonated molecule ([M-H]−).
    - CHARGE: Charge state of the molecule.
    - PARENT_CID: PubChem Compound ID of the parent compound.
    - CAS: CAS Registry Number.
    - CHEBI: ChEBI (Chemical Entities of Biological Interest) identifier.
    - HMDB: Human Metabolome Database identifier.
    - PC_CID: PubChem Compound ID.
    - PC_SID: PubChem Substance ID.
    - METLIN_ID: METLIN metabolite database identifier.

# Archived Files sorted alphabetically

The following files are no longer used by the package and will be
removed in the future. They are kept for a short while as a
reference.

## model-2022-03-10_xxx

- Links:
    - [archived/model-2022-03-10_Flatter gradient.xlsx]
    - [archived/model-2022-03-10_Lower flow rate 025.xlsx]
    - [archived/model-2022-03-10_Lower Temperature 25.xlsx]
    - [archived/model-2022-03-10_Lower Temperature 25 + Lower Flow rate 025.xlsx]
    - [archived/model-2022-03-10_Steeper gradient.xlsx]
    - [archived/model-2022-03-10_Steeper gradient+lower Temperature 25+Lower Flow rate 025.xlsx]
- Description: Models for adjusted for different chromatographic conditions. Adjustment done in 2022-03-10 using an early version of FastRet. Should not be used anymore.

## predictor_set_2024-05-21.xlsx

- Description: 131 Chemical Descriptors obtained using rCDK for 439 SMILES from **In-House-Library**. Chemical Descriptors were obtained using an early version of FastRet (which used an early version of rCDK). Should not be used anymore.
- Dimensions: 440 x 132 (row one contains colnames, col one contains SMILES)
- Columns:
    - (1) Smiles
    - (2) RT: Measured Retention Time
    - (3-132): Chemical Descriptor Names, e.g. 'Fsp3', 'nSmallRings' or 'nAromRings'

## R8P_Medoids_input_xxx

- Files:
    - [archived/R8P_Medoids_input_Flatter gradient.xlsx]
    - [archived/R8P_Medoids_input_Lower flow rate 025.xlsx]
    - [archived/R8P_Medoids_input_Lower Temperature 25 + Lower Flow rate 025.xlsx]
    - [archived/R8P_Medoids_input_Lower Temperature 25.xlsx]
    - [archived/R8P_Medoids_input_Steeper gradient.xlsx]
    - [archived/R8P_Medoids_input_Steeper gradient+lower Temperature 25+Lower Flow rate 025.xlsx]
- Description: Same as sheet 'R8_RT_Medoids' of [20211022_R8_dif_conditions_Medoids_validSet.xlsx] but split into seperate files so they can be used as input for the FastRet GUI.

## R8P_Validation_Set_input.xlsx

- Link: [archived/R8P_Validation_Set_input.xlsx]
- Description: Names and Smiles for **Validation Set Metabolites**. This file can be uploaded to the FastRet Web Client to get predictions for the Validation Set. It was used in combination with the models from [model-2022-03-10_xxx]
- Dimensions: 26 x 2 (row one contains colnames)
- Columns:
    - NAME, SMILES: Identifiers for measured Molecule.

## RP_adj.xlsx

- Link: [archived/RP_adj.xlsx]
- Description: Subset of 25 metabolites of `RP.xlsx` with artifically modified retention times. Used for testing of R functions. Has been copied to `inst/extdata` and `inst/extdata/RP_adj.xlsx` is the only version that should be used in the future.
- Dimensions: 26 x 3 (row one contains colnames)
- Columns: RT, NAME, SMILES (see `RP.xlsx`)

## RP_broken.xlsx

- Link: [archived/RP_broken.xlsx]
- Description: Empty excel file. Was used for testing of R functions.
- Dimensions: 0 x 0
- Columns: -

## ValidSet_xxx

- Links:
    - [archived/ValidSet_Flatter gradient.xlsx]
    - [archived/ValidSet_Lower flow rate 025.xlsx]
    - [archived/ValidSet_Lower Temperature 25 + Lower Flow rate 025.xlsx]
    - [archived/ValidSet_Lower Temperature 25.xlsx]
    - [archived/ValidSet_Steeper gradient.xlsx]
    - [archived/ValidSet_Steeper gradient+lower Temperature 25+Lower Flow rate 025.xlsx]
- Columns:
    - NA: contains the row-numbers from 1 to 25
    - NAME, SMILES: Identifiers for measured Molecule.
    - pred_RT: Predicted Retention Time in minutes, e.g. '0.91' or '10.2'
- Description: same as Sheet 'R8_RT_Validation set' of [20211022_R8_dif_conditions_Medoids_validSet.xlsx] but without the measured Retention times and split into seperate files so they can be used as input for the FastRet GUI. They predicted retention times observed in 2022 are included as well.

# Files sorted by Category

## Full Training Datasets (~ 400 Metabolites)

- HILIC-Retip: availble through R function `FastRet::read_retip_hilic_data()`
- HILIC: [20210702_RT_Prediction_Hilic_Library_FF.xlsx]
- RP-AXMM: [RP-AXMM_FastRet_Input.xlsx]
- RP: [RP.xlsx]

## Model Adjustment Datasets (Subset of 25 Metabolites)

- Combined: sheet 'R8_RT_Medoids' of [20211022_R8_dif_conditions_Medoids_validSet.xlsx]
- Archived:
    - RP-Steep: [archived/R8P_Medoids_input_Steeper gradient.xlsx]
    - RP-Flat: [archived/R8P_Medoids_input_Flatter gradient.xlsx]
    - RP-T25: [archived/R8P_Medoids_input_Lower Temperature 25.xlsx]
    - RP-FR25: [archived/R8P_Medoids_input_Lower flow rate 025.xlsx]
    - RP-T25-FR25  [archived/R8P_Medoids_input_Lower Temperature 25 + Lower Flow rate 025.xlsx]
    - RP-T25-Fr25-Steep: [archived/R8P_Medoids_input_Steeper gradient+lower Temperature 25+Lower Flow rate 025.xlsx]

## Validation Data Sets (New set of 22 Metabolites)

- Combined: sheet 'R8_RT_Validation' of [20211022_R8_dif_conditions_Medoids_validSet.xlsx]
- Archived:
    - RP-Steep-Val: [archived/ValidSet_Steeper gradient.xlsx]
    - RP-Flat-Val: [archived/ValidSet_Flatter gradient.xlsx]
    - RP-T25-Val: [archived/ValidSet_Lower Temperature 25.xlsx]
    - RP-FR25-Val: [archived/ValidSet_Lower flow rate 025.xlsx]
    - RP-T25-FR25-Val: [archived/ValidSet_Lower Temperature 25 + Lower Flow rate 025.xlsx]
    - RP-T25-Fr25-Steep-Val: [archived/ValidSet_Steeper gradient+lower Temperature 25+Lower Flow rate 025.xlsx]

# How Measurements are done

Note: The following is a general description of the experimental
workflow for measuring metabolite retention times using LC-MS.
Details may vary depending on laboratory protocols and
instrumentation.

1. Compound Acquisition: Metabolite standards are ordered from a
   vendor. They are typically delivered in plates (e.g., 96-well
   or 384-well format), with each well containing a single
   metabolite. Each metabolite can be identified by a unique ID,
   e.g. its HMDBID, CASID, CHEBID, PBCID or INCHIKEY.

2. Preparation of Input Plate: Individual metabolites are
   dissolved in an appropriate solvent to a defined concentration.
   Aliquots of these solutions are transferred to a new plate (the
   "input plate"), which is compatible with the LC-MS autosampler.
   Each well of the input plate contains either a single
   metabolite or a defined mixture, depending on the experimental
   design.

3. LC-MS Analysis: The input plate is loaded into the LC-MS
   system. The autosampler injects samples from each well into the
   liquid chromatography (LC) column, where metabolites are
   separated based on their physicochemical properties. As each
   metabolite elutes from the column, it enters the mass
   spectrometer (MS), which records a mass spectrum for each
   elution event.

4. Data Processing: For each injection, the LC-MS system generates
   chromatograms and mass spectra. For each compound, the average
   retention time is written in a xlsx sheet. If multiple
   compounds are present in a well, they are distinguished based
   on their unique m/z values.

# Problems

- Number of metabolites given for HILIC is different between paper (364) and excel (392).
- Excel file for RP-AXMM column is not available.
- Number of metabolites given for RP is different between paper (401) and excel (442).

<!-- Active Dataset Links -->

[RP.xlsx]: raw/RP.xlsx
[RP-AXMM_FastRet_Input.xlsx]: raw/RP-AXMM_FastRet_Input.xlsx
[20210702_RT_Prediction_Hilic_Library_FF.xlsx]: raw/20210702_RT_Prediction_Hilic_Library_FF.xlsx
[20211022_R8_dif_conditions_Medoids_validSet.xlsx]: raw/20211022_R8_dif_conditions_Medoids_validSet.xlsx
[20250623_IROA_Library.xlsx]: raw/20250623_IROA_Library.xlsx

<!-- Heading Links -->

[model-2022-03-10_xxx]: #model-2022-03-10_xxx

<!-- Archived Dataset Links -->

[archived/ValidSet_Steeper gradient+lower Temperature 25+Lower Flow rate 025.xlsx]: archived/ValidSet_Steeper%20gradient+lower%20Temperature%2025+Lower%20Flow%20rate%20025.xlsx
[archived/ValidSet_Steeper gradient.xlsx]: archived/ValidSet_Steeper%20gradient.xlsx
[archived/ValidSet_Lower Temperature 25.xlsx]: archived/ValidSet_Lower%20Temperature%2025.xlsx
[archived/ValidSet_Lower Temperature 25 + Lower Flow rate 025.xlsx]: archived/ValidSet_Lower%20Temperature%2025%20+%20Lower%20Flow%20rate%20025.xlsx
[archived/ValidSet_Lower flow rate 025.xlsx]: archived/ValidSet_Lower%20flow%20rate%20025.xlsx
[archived/ValidSet_Flatter gradient.xlsx]: archived/ValidSet_Flatter%20gradient.xlsx
[archived/RP_broken.xlsx]: archived/RP_broken.xlsx
[archived/RP_adj.xlsx]: archived/RP_adj.xlsx
[archived/R8P_Validation_Set_input.xlsx]: archived/8P_Validation_Set_input.xlsx
[archived/R8P_Medoids_input_Steeper gradient+lower Temperature 25+Lower Flow rate 025.xlsx]: archived/R8P_Medoids_input_Steeper%20gradient+lower%20Temperature%2025+Lower%20Flow%20rate%20025.xlsx
[archived/R8P_Medoids_input_Steeper gradient.xlsx]: archived/R8P_Medoids_input_Steeper%20gradient.xlsx
[archived/R8P_Medoids_input_Lower Temperature 25.xlsx]: archived/R8P_Medoids_input_Lower%20Temperature%2025.xlsx
[archived/R8P_Medoids_input_Lower Temperature 25 + Lower Flow rate 025.xlsx]: archived/R8P_Medoids_input_Lower%20Temperature%2025%20+%20Lower%20Flow%20rate%20025.xlsx
[archived/R8P_Medoids_input_Lower flow rate 025.xlsx]: archived/R8P_Medoids_input_Lower%20flow%20rate%20025.xlsx
[archived/R8P_Medoids_input_Flatter gradient.xlsx]: archived/R8P_Medoids_input_Flatter%20gradient.xlsx
[archived/model-2022-03-10_Steeper gradient+lower Temperature 25+Lower Flow rate 025.xlsx]: archived/model-2022-03-10_Steeper%20gradient+lower%20Temperature%2025+Lower%20Flow%20rate%20025
[archived/model-2022-03-10_Steeper gradient.xlsx]: archived/model-2022-03-10_Steeper%20gradient
[archived/model-2022-03-10_Lower Temperature 25.xlsx]: archived/model-2022-03-10_Lower%20Temperature%2025
[archived/model-2022-03-10_Lower Temperature 25 + Lower Flow rate 025.xlsx]: archived/model-2022-03-10_Lower%20Temperature%2025%20+%20Lower%20Flow%20rate%20025
[archived/model-2022-03-10_Lower flow rate 025.xlsx]: archived/model-2022-03-10_Lower%20flow%20rate%20025
[archived/model-2022-03-10_Flatter gradient.xlsx]: archived/model-2022-03-10_Flatter%20gradient
