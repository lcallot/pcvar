\name{trace.test}
\alias{trace.test}
\title{Computes the trace test for a pcvar model.}
\usage{
  trace.test(Y, W = NULL, lags = 2, ddet = 1, Y0 = NULL)
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

  \item{Y0}{Optional: a molten set of global exogenous data
  to use. Its columns must be identical to those of Y, and
  have the same number of observations.}
}
\value{
  A list with 'est' the individual CVAR and 'trace' the
  trace test statistics for each model and each hypothesis.
}
\description{
  Computes the trace test for a pcvar model.
}
\author{
  Laurent Callot \email{l.callot@vu.nl}
}

