\name{melt.dat}
\alias{melt.dat}
\title{Transforms a matrix of data to a molten data frame for use in rank.test..}
\usage{
  melt.dat(Ymat, indiv.names, var.names,
    time.index = 1:nrow(Ymat), freq = 1)
}
\arguments{
  \item{Ymat}{A matrix of data with variables in columns.
  The variables are assumed grouped by countries.}
}
\value{
  A molten dataframe ready for use in rank test.
}
\description{
  Transforms a matrix of data to a molten data frame for
  use in rank.test..
}
\author{
  Laurent Callot \email{l.callot@vu.nl}
}

