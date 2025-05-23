% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/project_threshold.R
\name{project_threshold}
\alias{project_threshold}
\title{Project density dependent population dynamics}
\usage{
project_threshold(
  A,
  vector = NULL,
  time = 100,
  return.vec = TRUE,
  threshold = NULL,
  A_DD = NULL
)
}
\arguments{
\item{A}{A matrix, or list of matrices.}

\item{vector}{A numeric vector describiing the age/stage distribution.}

\item{time}{The number of projection intervals.}

\item{return.vec}{If TRUE, returns the time series of demographic (st)age
vectors.}

\item{threshold}{The population threshold.}

\item{A_DD}{A single high density population matrix.}
}
\value{
A list with projection results.
}
\description{
Project density dependent dynamics of a population matrix
projection model.
}
\details{
A version of popbio::project that will project forward
until reaching a threshold population size, at which point
it will switch to a specified high-density population matrix.
}
\examples{
\dontrun{
# Try out the DD projection function when we use a single matrix:
nohunt_threshold<- project_threshold(pop_mat_100_none, vector=pop_n,
                                    time=100, threshold=1000, A_DD=nohunt_highD)
}

}
\author{
Chrissy Hernandez \href{mailto:chrissy.hernandez@biology.ox.ac.uk}{chrissy.hernandez@biology.ox.ac.uk}
}
