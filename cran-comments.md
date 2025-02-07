## R CMD check results

0 errors | 0 warnings | 0 notes

## Changes since removal of package from CRAN

This package has already been accepted by CRAN at 25.06.2024, but was removed later on due to the following note in the "donttest" log:

> 1 NOTE: Found the following files/directories: ‘~/.cache/R/FastRet’

This note has been addressed by adding an `.onUnload()` function, that ensures that the cache directory is removed if it doesn't contain any files that should persist between R sessions.

Apart from that, the following improvements have been made to the package:

* Added an article about installation details incl. a troubleshooting section
* Improved function docs
* Improved examples by removing `donttest` blocks
* Improved examples & tests by using smaller example datasets to reduce runtime
