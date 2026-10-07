% Generated manually for Bioconductor preparation
\name{hotspots}
\alias{hotspots}
\alias{hotspots_snvs}
\alias{hotspots_inframes}
\alias{hotspots_splice}
\alias{hotspots_splice_hg19}
\alias{hotspots_splice_hg38}
\docType{data}
\title{CancerHotspots v2 Annotation Keys}
\format{
Character vectors used as compact lookup keys for mutation hotspot annotation:
\describe{
  \item{hotspots_snvs}{Gene and amino-acid position keys for recurrent SNV hotspots, stored as strings such as \code{"AKT1 17"}.}
  \item{hotspots_inframes}{Gene and amino-acid position keys for recurrent in-frame indel hotspots.}
  \item{hotspots_splice}{Genomic coordinate keys for splice hotspots retained for backwards compatibility.}
  \item{hotspots_splice_hg19}{Genome-build-specific splice hotspot coordinate keys for hg19.}
  \item{hotspots_splice_hg38}{Genome-build-specific splice hotspot coordinate keys for hg38.}
}
}
\source{
Derived from CancerHotspots.org v2 resources. CancerHotspots.org v2 is distributed under the Open Database License (ODbL); downstream use must preserve required attribution and license obligations.
}
\usage{
data(hotspots_snvs)
data(hotspots_inframes)
data(hotspots_splice)
data(hotspots_splice_hg19)
data(hotspots_splice_hg38)
}
\description{
Compact lookup vectors used by Jumble plotting and annotation helpers to flag recurrent cancer mutation hotspots in variant tables. These resources support visualization and annotation overlays and are not a replacement for clinical variant interpretation databases.
}
\details{
The vectors are intentionally small package resources used for offline plotting and annotation during tests, examples, and reports. The splice hotspot vectors use \code{chromosome:position} keys, while SNV and in-frame indel vectors use gene-position keys.
}
\examples{
data(hotspots_snvs)
head(hotspots_snvs)
length(hotspots_snvs)
}
\keyword{datasets}
