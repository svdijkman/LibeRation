\name{nm_advan_template}
\alias{nm_advan_template}
\title{Create an Editable Standard ADVAN Model Template}
\usage{
nm_advan_template(advan = 1L, trans = NULL, n_state = NULL, problem = NULL)
}
\arguments{
\item{advan}{ADVAN number from 1 through 14.}
\item{trans}{Optional TRANS number from 1 through 6. ODE ADVANs use their
template default.}
\item{n_state}{Optional number of states for ADVAN 6, 8, 9, 13, or 14.}
\item{problem}{Optional display name for the model.}
}
\value{
An editable \code{nm_model} object.
}
\description{
Returns the same ADVAN template used by the LibeRation New Model Version
workflow, allowing other LibeR packages to create compatible models without
duplicating template definitions.
}
