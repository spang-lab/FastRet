# Chemical Descriptors Names

This object contains the names of various chemical descriptors.

## Usage

``` r
CDNames
```

## Format

An object of class `character` of length 45.

## Details

One descriptor can be associated with multiple features, e.g. the BCUT
descriptor corresponds to the following features: BCUTw.1l, BCUTw.1h,
BCUTc.1l, BCUTc.1h, BCUTp.1l, BCUTp.1h. Some descriptors produce
warnings for certain molecules., e.g. "The AtomType null could not be
found" or "Molecule must have 3D coordinates" and return NA in such
cases. Descriptors that produce only NAs in our test datasets will be
excluded. To see which descriptors produce only NAs, run
`analyzeCDNames`. The "LongestAliphaticChain" descriptors sometimes even
produces `Error: segfault from C stack overflow` error, e.g. for SMILES
`c1ccccc1C(Cl)(Cl)Cl` (== `rcdk::bpdata$SMILES[200]`) when using
`OpenJDK Runtime Environment (build 11.0.23+9-post-Ubuntu-1ubuntu122.04.1)`.
Therefore, this descriptor is also excluded.

## See also

[`analyzeCDNames()`](https://spang-lab.github.io/FastRet/reference/analyzeCDNames.md),
[CDFeatures](https://spang-lab.github.io/FastRet/reference/CDFeatures.md)

## Examples

``` r
str(CDNames)
#>  chr [1:45] "org.openscience.cdk.qsar.descriptors.molecular.FractionalCSP3Descriptor" ...
```
