% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scraps_priming_region.R
\name{fill_mat}
\alias{fill_mat}
\title{Fill matrix with 0 rows or columns as needed}
\usage{
fill_mat(x, vec, row = TRUE)
}
\arguments{
\item{x}{target matrix}

\item{vec}{target row(or column) names to ensure exists and in right order}

\item{row}{if TRUE, target rows, if FALSE, columns}
}
\value{
0 filled and ordered matrix
}
\description{
Fill matrix with 0 rows or columns as needed
}
