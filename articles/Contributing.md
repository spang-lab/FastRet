# Contributing

To contribute to this package, you should follow the below steps:

1.  Create a issue at
    [github.com/spang-lab/FastRet/issues](https://github.com/spang-lab/FastRet/issues)
    describing the problem or feature you want to work on.
2.  Wait until the issue is approved by a package maintainer.
3.  Create a fork of the repository at
    [github.com/spang-lab/FastRet](https://github.com/spang-lab/FastRet)
4.  Make your edits as described in section [Making
    Edits](#making-edits)
5.  Create a pull request at
    [github.com/spang-lab/FastRet/pulls](https://github.com/spang-lab/FastRet/pulls)

## Requirements

FastRet computes chemical descriptors with the `rcdk` package, which
needs a Java Development Kit (JDK) version 11 or higher. See the
[Installation](https://spang-lab.github.io/FastRet/articles/Installation.html)
article for how to set it up.

The GUI end-to-end tests in `tests/testthat/test-gui-e2e.R` drive the
Shiny app through a headless Chrome browser via the suggested packages
`chromote` and `shinytest2`. They are skipped on CRAN (i.e. unless
`NOT_CRAN=true`) and whenever no Chrome or Chromium binary is found. To
run them locally without root permissions, install Chrome for Testing
with the script
[misc/scripts/install-chrome-for-testing.sh](https://github.com/spang-lab/FastRet/blob/main/misc/scripts/install-chrome-for-testing.sh),
point `chromote` at it via the environment variable `CHROMOTE_CHROME`
and set `NOT_CRAN=true`:

``` bash
bash misc/scripts/install-chrome-for-testing.sh
export CHROMOTE_CHROME="$HOME/.local/share/chrome-for-testing/chrome-linux64/chrome"
export NOT_CRAN=true
Rscript -e 'devtools::test(filter = "gui-e2e")'
```

## Making Edits

Things you can update, are:

1.  Function code in folder
    [R](https://github.com/spang-lab/FastRet/tree/main/R)
2.  Function documentation in folder
    [R](https://github.com/spang-lab/FastRet/tree/main/R)
3.  Package documentation in folder
    [vignettes](https://github.com/spang-lab/FastRet/tree/main/vignettes)
4.  Test cases in folder
    [tests](https://github.com/spang-lab/FastRet/tree/main/tests)
5.  Dependencies in file
    [DESCRIPTION](https://github.com/spang-lab/FastRet/blob/main/DESCRIPTION)
6.  Authors in file
    [DESCRIPTION](https://github.com/spang-lab/FastRet/blob/main/DESCRIPTION)

Whenever you update any of those things, you should run the below
commands to check that everything is still working as expected:

``` r

devtools::document() # Build files in man folder
devtools::spell_check() # Check spelling (add false positives to inst/WORDLIST)
urlchecker::url_check() # Check URLs
devtools::test() # Execute tests from tests folder
toscutil::check_pkg_docs() # Check function documentation for missing tags
devtools::check() # Check package formalities
devtools::install() # Install as required by next command
pkgdown::build_site() # Build website in docs folder
```

Every pull request to `main` must increase the version number in
[DESCRIPTION](https://github.com/spang-lab/FastRet/blob/main/DESCRIPTION);
this is enforced by a CI check. Add a matching entry to
[NEWS.md](https://github.com/spang-lab/FastRet/blob/main/NEWS.md)
describing your changes.

After doing these steps, you can push your changes to GitHub.

## Releasing to CRAN

Whenever a package maintainer wants to release a new version of the
package to CRAN, they should:

1.  Check whether the [release
    requirements](https://r-pkgs.org/release.html#sec-release-initial)
    are fulfilled
2.  Make sure the version in `DESCRIPTION` has been bumped and `NEWS.md`
    has an entry for it
3.  Use the following commands to do a final check of the package and
    release it to CRAN

``` r

# Check spelling and URLs. False positive findings of spell check should be
# added to inst/WORDLIST.
devtools::spell_check()
urlchecker::url_check()

# Slower, but more realistic tests than devtools::check()
rcmdcheck::rcmdcheck(
  args = c("--no-manual", "--as-cran"),
  build_args = ("--no-manual"),
  error_on = ("warning"),
  check_dir = "../FastRet-RCMDcheck"
)
devtools::check(
  remote = TRUE,
  manual = TRUE,
  run_dont_test = TRUE
)

# Check reverse dependencies. For details see:
# https://r-pkgs.org/release.html#sec-release-revdep-checks
revdepcheck::revdep_check(num_workers = 8)

# Send your package to CRAN's builder services. You should receive an e-mail
# within about 30 minutes with a link to the check results. Checking with
# check_win_devel is required by CRAN policy and will (also) be done as part
# of CRAN's incoming checks.
devtools::check_win_oldrelease()
devtools::check_win_release()
devtools::check_win_devel()
devtools::check_mac_release()

# Update cran-comments.md with the results of the above checks.

# Use the following command to submit the package to CRAN or submit via the web
# interface available at https://cran.r-project.org/submit.html.
devtools::submit_cran()
```

After CRAN has accepted the submission:

1.  Tag the released commit as `vX.Y.Z`
    (e.g. `git tag v1.5.2 && git push origin v1.5.2`)
2.  Create a GitHub release for the tag, using the corresponding
    `NEWS.md` entry as release notes
