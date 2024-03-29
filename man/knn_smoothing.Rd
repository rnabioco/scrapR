% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scraps_priming_region.R
\name{knn_smoothing}
\alias{knn_smoothing}
\title{calculation mean diff and effect size of utr3 priming}
\usage{
knn_smoothing(countlist, knn_matrix, cell_ids)
}
\arguments{
\item{countlist}{countlist}

\item{knn_matrix}{nearest neightbors in matrix format, can be calculated by enabling return.neighbor in Seurat::FindNeighbors}

\item{cell_ids}{cell barcodes to use}
}
\value{
countlist after smoothing
}
\description{
calculation mean diff and effect size of utr3 priming
}
