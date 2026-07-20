\name{nm_report_design}
\alias{nm_report_block}
\alias{nm_report_design}
\alias{nm_report_design_save}
\alias{nm_report_design_load}
\alias{nm_report_design_render}
\title{Build and render a visual modelling report workflow}
\usage{
nm_report_block(type, title = NULL, source = c("user", "ai"), text = "",
  file = NULL, run_ids = NULL, elements = NULL, options = list(), id = NULL)

nm_report_design(blocks, title = "LibeRation modelling report",
  formats = c("docx", "pdf"), style = list(), id = NULL)

nm_report_design_save(workspace, project, design)

nm_report_design_load(workspace, project, id = NULL)

nm_report_design_render(design, workspace, project, directory = getwd(),
  name = "liberation-report", formats = design$formats, manifest = TRUE)
}
\arguments{
  \item{type}{Block type, including narrative, run, comparison, or page break.}
  \item{title}{Block heading or report title.}
  \item{source}{\code{"user"} or browser-local \code{"ai"} narrative source.}
  \item{text}{Editable narrative text. AI drafting is completed before render.}
  \item{file}{Optional TXT, Markdown, or DOCX narrative source.}
  \item{run_ids}{Saved immutable model-run identifiers.}
  \item{elements}{Evidence to include: summary, parameters, code, GOF,
  VPC, categorical/TTE VPC, NPDE, NPC, covariance, or run information.}
  \item{options}{Named block options such as instructions and template text.}
  \item{id}{Stable block or design identifier.}
  \item{blocks}{Ordered list of report blocks.}
  \item{formats}{One or both of \code{"docx"} and \code{"pdf"}.}
  \item{style}{Named rendering options, including an optional reference DOCX.}
  \item{workspace}{An \code{nm_workspace} or workspace path.}
  \item{project}{Project identifier.}
  \item{design}{An \code{nm_report_design}.}
  \item{directory}{Output directory.}
  \item{name}{Safe output base name.}
  \item{manifest}{Write a machine-readable JSON provenance manifest.}
}
\value{
A report block, report design, saved design metadata, or an
\code{nm_report_bundle} containing generated DOCX/PDF/JSON paths.
}
\description{
Defines the serializable backend for the drag-and-drop report designer. Render
uses saved run evidence only and never contacts an AI provider. DOCX rendering
requires Pandoc; PDF rendering uses the R graphics device.
}
