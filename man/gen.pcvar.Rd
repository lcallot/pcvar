\name{gen.pcvar}
\alias{gen.pcvar}
\title{Generates data from a PCVAR model.}
\usage{
  gen.pcvar(obs, N, nvar, W = NULL, Alpha = NULL,
    Beta = NULL, Lambda0, Gammal = NULL, Omega,
    err.dist = "gaussian", burn.smpl = 10, t.df = 1,
    det.type = 1, cdet.load = c(1, 1), res.dep = "iid",
    garch.spec = NULL)
}
\arguments{
  \item{Y}{A molten data set containing the 4 columns.
  indiv: the individual (country, sector,...) name as
  factor levels. Year: the year of observation. Period: the
  period (day, week, month, quarter) index as consecutive
  integers. Variable: the name of the variable. value: the
  value of the variable at the given time for a given
  individual.}

  \item{bs.method}{'resample' for iid resampling, 'wild'
  for gaussian wild bootstrap.}
}
\value{
  A list.
}
\description{
  Generates data from a PCVAR model.
}
\author{
  Laurent Callot \email{l.callot@vu.nl}
}

