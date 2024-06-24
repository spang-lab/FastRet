## R CMD check results

0 errors | 0 warnings | 1 note

## Changes since initial submission of FastRet (v1.1.2) at Fri 6/21/2024 8:28 AM

*   Fixed issues mentioned by CRAN reviewer:
    *   __Comment 1__: *Please do not modify the .GlobalEnv. This is not allowed by the CRAN policies. -> R/patch.R*
    *   __Solution__: Moved `patch.R` from the `R` folder to `misc/scripts`, which is excluded from the package build using `.Rbuildignore`. The file is conditionally sourced by the private function `start_gui_in_devmode()` if available, allowing its use during development without including it in the package.
    *   __Comment 2__: *Please add \value to .Rd files regarding exported methods and explain the functions' results in the documentation. Please write about the structure of the output (class) and also what the output means. (If a function does not return a value, document that too, e.g., \value{No return value, called for side effects} or similar) -> Missing Rd-tags: adjust_frm.Rd: \value, analyzeCDNames.Rd: \value, getCDs.Rd: \value, getCDsFor1Molecule.Rd: \value, read_rpadj_xlsx.Rd: \value*
    *   __Solution__: Added `\value` tags to the mentioned `.Rd` files describing the functions' return values.
    *   __Comment 3__: *If there are references describing the methods in your package, please add these in the description field of your DESCRIPTION file in the form authors (year) <doi:...> authors (year, ISBN:...) or if those are not available: <https:...> with no space after 'doi:', 'https:' and angle brackets for auto-linking. (If you want to add a title as well please put it in quotes: "Title")*
    *   __Solution__: Added *Bonini et al. (2020) <doi:10.1021/acs.analchem.9b05765>* as reference to the description part of the DESCRIPTION file, listing  it as *Related work*. This reference is used in the documentation for `read_retip_hilic_data()` and `ram_cache`. No additional references are used in the package documentation.
*   Added Fadi Fadil as a contributor. Fadi measured the example datasets shipped with FastRet.
*   Added ORCID IDs for contributors as described in [CRAN's checklist for submissions](https://cran.r-project.org/web/packages/submission_checklist.html).
