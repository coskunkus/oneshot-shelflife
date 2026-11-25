# =============================================================================
# export_shinylive.R
# Turns the Shiny application (shiny/app.R + code/core_oneshot.R) into a
# static website that runs entirely in the browser (webR), written to docs/.
# The docs/ folder can then be served by GitHub Pages: no server, no quota.
# Needs, once: install.packages(c("shinylive", "httpuv"))
# =============================================================================
locate_root <- function(core = "core_oneshot.R") {
  cand <- character(0)
  for (i in rev(seq_len(sys.nframe()))) {
    f <- tryCatch(sys.frame(i)$ofile, error = function(e) NULL)
    if (!is.null(f)) cand <- c(cand, dirname(normalizePath(f)))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      tryCatch(rstudioapi::isAvailable(), error = function(e) FALSE)) {
    p <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
    if (is.character(p) && length(p) == 1 && nzchar(p)) cand <- c(cand, dirname(normalizePath(p)))
  }
  fa <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(fa)) cand <- c(cand, dirname(normalizePath(sub("^--file=", "", fa[1]))))
  cand <- c(cand, getwd())
  for (cd in cand) for (up in c(cd, dirname(cd))) {
    if (file.exists(file.path(up, "code", core))) return(up)
    if (file.exists(file.path(up, core)) && basename(up) == "code") return(dirname(up))
  }
  stop("project root not found")
}
setwd(locate_root())

# Stage the app and the core in an ASCII-only temporary folder (the project
# path contains non-ASCII characters, which some tools do not accept).
app_dir <- file.path(tempdir(), "oneshot_app")
unlink(app_dir, recursive = TRUE); dir.create(app_dir)
file.copy(file.path("shiny", "app.R"), app_dir)
file.copy(file.path("code", "core_oneshot.R"), app_dir)

out_dir <- file.path(tempdir(), "docs")
unlink(out_dir, recursive = TRUE)
shinylive::export(app_dir, out_dir)

# Copy the result into the project folder
unlink("docs", recursive = TRUE)
dir.create("docs")
file.copy(list.files(out_dir, full.names = TRUE, all.files = TRUE, no.. = TRUE), "docs", recursive = TRUE)
cat("Static app written to docs/.\n")
cat("Test it locally with:  httpuv::runStaticServer(\"docs\")\n")
