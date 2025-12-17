# Selective Measuring

The function
[`adjust_frm()`](https://spang-lab.github.io/FastRet/reference/adjust_frm.md)
is used to modify existing FastRet models based on changes in
chromatographic conditions. It requires a set of molecules with measured
retention times on both the original and new column. This function
selects a sensible subset of molecules from the original dataset for
re-measurement. The selection process includes:

1.  Generating chemical descriptors from the SMILES strings of the
    molecules. These are the features used by
    [`train_frm()`](https://spang-lab.github.io/FastRet/reference/train_frm.md)
    and
    [`adjust_frm()`](https://spang-lab.github.io/FastRet/reference/adjust_frm.md).

2.  Standardizing chemical descriptors to have zero mean and unit
    variance.

3.  Training a Ridge Regression model with the standardized chemical
    descriptors as features and the retention times as the target
    variable.

4.  Scaling the chemical descriptors by coefficients of the Ridge
    Regression model.

5.  Clustering the entire dataset, which includes the scaled chemical
    descriptors and the retention times.

6.  Returning the clustering results, which include the cluster
    assignments, the medoid indicators, and the raw data.

## Usage

``` r
selective_measuring(
  raw_data,
  k_cluster = 25,
  verbose = 1,
  seed = NULL,
  rt_coef = "max_ridge_coef"
)
```

## Arguments

- raw_data:

  The raw data to be processed. Must be a dataframe with columns NAME,
  RT and SMILES.

- k_cluster:

  The number of clusters for PAM clustering.

- verbose:

  The level of verbosity.

- seed:

  An optional random seed for reproducibility, set at the beginning of
  the function.

- rt_coef:

  Which coefficient to use for scaling RT before clustering. Options
  are:

  - max_ridge_coef: scale with the maximum absolute coefficient obtained
    in ridge regression. I.e., RT will have approximately the same
    weight as the most important chemical descriptor.

  - 1: do not scale RT any further, i.e., use standardized RT. The
    effect of leaving RT unscaled is kind of unpredictable, as the ridge
    coefficients depend on the dataset. If the maximum absolute
    coefficient is much smaller than 1, RT will dominate the clustering.
    If it is much larger than 1, RT will have little influence on the
    clustering.

  - 0: exclude RT from the clustering.

## Value

A list containing the following elements:

- `clustering`: A data frame with columns RT, SMILES, NAME, CLUSTER and
  IS_MEDOID.

- `clobj`: The clustering object. The object returned by the clustering
  function. Depends on the `method` parameter.

- `coefs`: The coefficients from the Ridge Regression model.

- `model`: The Ridge Regression model.

- `df`: The preprocessed data.

- `dfz`: The standardized features.

- `dfzb`: The features scaled by the coefficients (betas) of the Ridge
  Regression model.

## Examples

``` r
x <- selective_measuring(RP[1:50, ], k = 5, verbose = 0)
# For the sake of a short runtime, only the first 50 rows of the RP dataset
# were used in this example. In practice, you should always use the entire
# dataset to find the optimal subset for re-measurement.
```
