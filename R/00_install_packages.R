required <- c("survival", "ggplot2", "patchwork", "svglite", "cmprsk")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing)) {
  install.packages(missing, repos = "https://cloud.r-project.org")
}

message("Required R packages are available.")

