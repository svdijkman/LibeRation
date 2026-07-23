\name{liber_support_matrix}
\alias{liber_support_matrix}
\title{Report ecosystem capability and evidence tiers}
\usage{
liber_support_matrix(package = NULL, evidence_tier = NULL)
}
\arguments{
\item{package}{Optional package name or names to retain.}
\item{evidence_tier}{Optional tier or tiers to retain.}
}
\value{
A data frame with capability status, evidence, gate, and recommended use.
}
\description{
Returns the versioned support declaration for all LibeR packages. A
\code{validated} capability has an independent analytic or external-software
comparison, \code{verified} denotes deterministic internal and recovery
testing, and \code{experimental} denotes a research implementation without
sufficient qualification for routine use.
}
