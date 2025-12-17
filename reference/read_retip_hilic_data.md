# Read the HILIC dataset from the Retip package

Reads the `Retip::HILIC` dataset (CC BY 4.0) from the Retip package or,
if Retip is not installed, downloads the dataset directly from the
[Retip GitHub repository](https://github.com/oloBion/Retip). Before
returning the dataset, SMILES strings are canonicalized and the original
`tibble` object is converted to a base R `data.frame`.

## Usage

``` r
read_retip_hilic_data(verbose = 1)
```

## Source

<https://github.com/oloBion/Retip/raw/master/data/HILIC.RData>

## Arguments

- verbose:

  Verbosity. 1 for messages, 0 to suppress them.

## Value

A data frame with 970 rows and the following columns:

- `NAME`: Molecule name

- `INCHIKEY`: InChIKey

- `SMILES`: Canonical SMILES string

- `RT`: Retention time in Minutes

## Details

Attribution as required by CC BY 4.0:  

- Original dataset by: Paolo Bonini, Tobias Kind, Hiroshi Tsugawa,
  Dinesh Kumar Barupal, and Oliver Fiehn as part of the Retip project.  

- Source repository: <https://github.com/oloBion/Retip>  

- Original file:
  <https://github.com/oloBion/Retip/raw/master/data/HILIC.RData>  

- License: CC BY 4.0 (<https://creativecommons.org/licenses/by/4.0/>)  

- Modifications in FastRet:

  - converted tibble to data.frame

  - canonicalized SMILES using
    [`as_canonical()`](https://spang-lab.github.io/FastRet/reference/as_canonical.md)

  - renamed column 'INCHKEY' to 'INCHIKEY'

## References

Retip: Retention Time Prediction for Compound Annotation in Untargeted
Metabolomics  
Paolo Bonini, Tobias Kind, Hiroshi Tsugawa, Dinesh Kumar Barupal, and
Oliver Fiehn  
Analytical Chemistry 2020 92 (11), 7515-7522 DOI:
10.1021/acs.analchem.9b05765

## Examples

``` r
df <- read_retip_hilic_data(verbose = 0)
```
