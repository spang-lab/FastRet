% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getcds.R
\docType{data}
\name{CDNames}
\alias{CDNames}
\title{Chemical Descriptors Names}
\format{
An object of class \code{character} of length 45.
}
\usage{
CDNames
}
\description{
This object contains the names of various chemical descriptors.
}
\details{
One descriptor can be associated with multiple features, e.g. the BCUT
descriptor corresponds to the following features: BCUTw.1l, BCUTw.1h,
BCUTc.1l, BCUTc.1h, BCUTp.1l, BCUTp.1h. Some descriptors produce warnings for
certain molecules., e.g. "The AtomType null could not be found" or "Molecule
must have 3D coordinates" and return NA in such cases. Descriptors that
produce only NAs in our test datasets will be excluded. To see which
descriptors produce only NAs, run \code{analyzeCDNames}. The
"LongestAliphaticChain" descriptors sometimes even produces \verb{Error: segfault from C stack overflow} error, e.g. for SMILES \verb{c1ccccc1C(Cl)(Cl)Cl} (==
\code{rcdk::bpdata$SMILES[200]}) when using \verb{OpenJDK Runtime Environment (build 11.0.23+9-post-Ubuntu-1ubuntu122.04.1)}. Therefore, this descriptor is also
excluded.
}
\examples{
str(CDNames)
}
\seealso{
\code{\link[=analyzeCDNames]{analyzeCDNames()}}, \link{CDFeatures}
}
\keyword{internal}
