\name{nm_model_diagram}
\alias{nm_model_diagram}
\alias{nm_diagram_preview}
\alias{nm_diagram_generate}
\alias{nm_model_diagram_get}
\title{Define and generate a visual compartment model}
\usage{
nm_model_diagram(compartments, flows = NULL, parameters = NULL,
  advan = 6L, residual = c("additive", "proportional", "combined"),
  covariates = NULL, title = "Diagram model")

nm_diagram_preview(diagram)

nm_diagram_generate(diagram,
  input = c("ID", "TIME", "EVID", "AMT", "RATE", "II", "SS", "CMT", "DV", "MDV"),
  output = c("PRED", "IPRED", "CWRES"))

nm_model_diagram_get(model)
}
\arguments{
  \item{compartments}{Data frame or list of compartment rows with identifiers,
  names, type, scale/volume parameters, dose/observation flags, and layout.}
  \item{flows}{Data frame or list of directed flows. Supported types are
  first-order rate, clearance, bidirectional clearance, Michaelis--Menten,
  zero-order, and custom nonlinear expressions. External source/sink is id 0.}
  \item{parameters}{Optional parameter definitions containing initial values,
  bounds, fixed flags, IIV flags, and ETA variances.}
  \item{advan}{General ODE solver: ADVAN6 or ADVAN13.}
  \item{residual}{Residual-error scaffold.}
  \item{covariates}{Input covariate names.}
  \item{title}{Model title.}
  \item{diagram}{An \code{nm_model_diagram}.}
  \item{input}{NONMEM-compatible input columns.}
  \item{output}{Generated output columns retained in model runs.}
  \item{model}{An \code{nm_model} or \code{NMEngine}.}
}
\value{
\code{nm_model_diagram()} returns a serializable diagram;
\code{nm_diagram_preview()} returns generated code and parameter tables;
\code{nm_diagram_generate()} returns an editable \code{nm_model}; and
\code{nm_model_diagram_get()} returns the stored graph or \code{NULL}.
}
\description{
Creates the semantic graph used by the LibeRation visual model builder. Model
code generation is explicit: previewing never changes the current model, and
applying generated code does not prevent subsequent manual editing of
\code{$PK/$PRED} or \code{$DES}. Structural parameters use log-normal ETAs by
default (\code{PARAM = THETA(i) * exp(ETA(i))}).
}
\examples{
diagram <- nm_model_diagram(
  compartments = data.frame(
    id = 1, name = "CENTRAL", kind = "amount",
    volume_parameter = "V", scale_parameter = "V",
    dose = TRUE, observe = TRUE, x = 200, y = 150
  ),
  flows = data.frame(
    id = "elimination", from = 1, to = 0, type = "clearance",
    parameter = "CL", secondary_parameter = "", expression = "", label = ""
  )
)
preview <- nm_diagram_preview(diagram)
model <- nm_diagram_generate(diagram)
}
