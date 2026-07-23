\name{liber_support_bundle}
\alias{liber_support_bundle}
\title{Create a redacted LibeR diagnostic support bundle}
\usage{
liber_support_bundle(
  path = file.path(
    getwd(),
    paste0("LibeR-support-", format(Sys.time(), "\%Y\%m\%d-\%H\%M\%S"), ".zip")
  ),
  workspace = NULL,
  model = NULL,
  data = NULL,
  fit = NULL,
  logs = NULL,
  include_model_code = FALSE,
  include_logs = FALSE,
  overwrite = FALSE
)
}
\arguments{
\item{path}{Output \code{.zip} path.}
\item{workspace}{Optional \code{nm_workspace} or workspace directory to verify.}
\item{model}{Optional \code{nm_model} to describe structurally.}
\item{data}{Optional data frame to describe without retaining values.}
\item{fit}{Optional \code{nm_fit} to describe without estimates or embedded data.}
\item{logs}{Optional character vector of log file paths or log text.}
\item{include_model_code}{Include model source code.}
\item{include_logs}{Include supplied logs after secret and path redaction.}
\item{overwrite}{Replace an existing archive.}
}
\value{
The normalized archive path, invisibly.
}
\description{
The default archive contains compatibility, compiled-engine, runtime,
structural model/data, and fit metadata. It never contains data values,
parameter estimates, ETAs, environment variables, credentials, workspace
contents, or model code. Model code and redacted logs are opt-in and should
still be inspected before sharing.
}
