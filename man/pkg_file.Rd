% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\name{pkg_file}
\alias{pkg_file}
\title{Get package file}
\usage{
pkg_file(path, mustWork = FALSE)
}
\arguments{
\item{path}{The path to the file within the package.}

\item{mustWork}{If TRUE, an error is thrown if the file does not exist.}
}
\value{
The path to the file.
}
\description{
Returns the path to a file within the FastRet package.
}
\examples{
path <- pkg_file("extdata/RP.xlsx")
}
\keyword{internal}
