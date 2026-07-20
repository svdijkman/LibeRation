## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")


## ----workspace----------------------------------------------------------------
library(LibeRation)

workspace <- nm_workspace(tempfile("liber-workspace-"))
project <- nm_project_create(
  workspace,
  name = "Oral PK development",
  description = "Example project created from R"
)
nm_project_list(workspace)


## ----save-version, eval=FALSE-------------------------------------------------
# version_id <- nm_project_save(
#   workspace,
#   project$id,
#   model = model,
#   data = data,
#   label = "Mod001"
# )
# 
# run_id <- nm_project_save_run(
#   workspace,
#   project$id,
#   version = version_id,
#   result = fit,
#   label = "FOCEI base model"
# )
# 
# nm_project_save_diagnostics(
#   workspace,
#   project$id,
#   run_id,
#   list(gof = nm_gof(fit), vpc = nm_vpc(fit, nsim = 200))
# )


## ----gui, eval=FALSE----------------------------------------------------------
# liber_gui()
# 
# # Open an explicit workspace and preload model/data objects.
# liber_gui(
#   workspace = "D:/projects/LibeR/workspace",
#   model = model,
#   data = data
# )


## ----import-control, eval=FALSE-----------------------------------------------
# control <- nm_control_read("run001.ctl", strict = FALSE)
# control$model
# control$compatibility
# 
# fit <- nm_est(control$model, data, method = "FOCEI")


## ----export-control, eval=FALSE-----------------------------------------------
# text <- nm_control_write(
#   model,
#   data = "theo.csv IGNORE=@",
#   estimation = "METHOD=1 INTERACTION MAXEVAL=9999 PRINT=10",
#   covariance = "UNCONDITIONAL"
# )
# cat(text)
# 
# nm_control_write(model, file = "exported.ctl", data = "theo.csv")


## ----report, eval=FALSE-------------------------------------------------------
# report <- nm_report(
#   fit,
#   file = "run001-report.pdf",
#   sections = c("summary", "parameters", "gof", "eta", "vpc"),
#   vpc = vpc,
#   manifest = TRUE
# )


## ----report-design, eval=FALSE------------------------------------------------
# design <- nm_report_design(list(
#   nm_report_block("introduction", text = "Study background"),
#   nm_report_block(
#     "run", run_ids = run_id,
#     elements = c("summary", "parameters", "gof", "vpc", "run_info")
#   ),
#   nm_report_block("discussion", source = "ai", text = "Reviewed local draft"),
#   nm_report_block("conclusion", text = "Conclusions")
# ))
# 
# nm_report_design_save(workspace, project$id, design)
# bundle <- nm_report_design_render(
#   design, workspace, project$id,
#   directory = "reports", formats = c("docx", "pdf")
# )

