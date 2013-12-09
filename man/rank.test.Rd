\name{rank.test}
\alias{rank.test}
\title{Applies the sequential panel rank test procedure of Callot (2013) to a pcvar model.}
\usage{
  rank.test(Y, W = NULL, lags = 2, ddet = 0, BS = 99,
    Y0 = NULL, bs.method = "resample")
}
\arguments{
  \item{Y}{A molten data set containing the 4 columns.
  indiv: the individual (country, sector,...) name as
  factor levels. Year: the year of observation. Period: the
  period (day, week, month, quarter) index as consecutive
  integers. Variable: the name of the variable. value: the
  value of the variable at the given time for a given
  individual.}

  \item{W}{NULL, A matrix, 3 dimensional array with first
  and second dimension equal to the number of individuals.
  Columns must sum to one, and values on the diagonal must
  be zero. In case the input is an array, its third
  dimension should be equal to the number of observations.
  In case of NULL, equal weights are used.}

  \item{lags}{An integer indicating the number of lags of
  the VAR in level.}

  \item{dett}{An integer indicating the type of
  deterministics to use, following the typology by Johansen
  1988.}

  \item{BS}{The number of bootstrap iterations for the
  test.}

  \item{Y0}{A molten set of global exogenous data to use.
  Its columns must be identical to those of Y, and have the
  same number of observations.}

  \item{bs.method}{'resample' for iid resampling, 'wild'
  for gaussian wild bootstrap.}
}
\value{
  A list.
}
\description{
  Applies the sequential panel rank test procedure of
  Callot (2013) to a pcvar model.
}
\author{
  Laurent Callot \email{l.callot@vu.nl}
}

