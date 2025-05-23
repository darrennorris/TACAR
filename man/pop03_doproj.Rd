% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pop03_doproj.R
\name{pop03_doproj}
\alias{pop03_doproj}
\title{Create data.frame with deterministic projection results.}
\usage{
pop03_doproj(x)
}
\arguments{
\item{x}{data.frame created by "PopPrep.R".}
}
\value{
Creates data.frame with projection results.
}
\description{
Creates data.frame with projection results for use in subsequent
population projection function ("PopProj.R").
}
\details{
This function applies a deterministic population projection moodel.
}
\examples{
\dontrun{
dfproj_vals <- pop03_doproj()
}

}
