suppressPackageStartupMessages({
  library(ggplot2)
  library(svglite)
})

root <- getwd()
if (!dir.exists(file.path(root, "outputs"))) {
  script_arg <- sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))][1])
  script_dir <- dirname(normalizePath(script_arg, winslash = "/", mustWork = FALSE))
  root <- if (basename(script_dir) == "R") dirname(script_dir) else script_dir
}
annual_path <- file.path(root, "outputs", "formal_next_phase", "formal_incidence_annual_fay_feuer.csv")
model_path <- file.path(root, "outputs", "joinpoint", "GCA_D1_D2_D3_Joinpoint.Export.Model.Estimates.txt")
figure_dir <- file.path(root, "outputs", "formal_next_phase", "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

annual <- read.csv(annual_path, check.names = FALSE)
annual$Definition <- c(
  D1_strict_GCA = "D1",
  D2_legacy_GCC_adenocarcinoid = "D2",
  D3_historical_expanded_spectrum = "D3"
)[annual$definition]

models <- read.delim(model_path, check.names = FALSE, na.strings = "NA")
fit_rows <- list()
k <- 1
for (definition in c("D1", "D2", "D3")) {
  m <- models[models$Definition == definition, ]
  years <- 2000:2023
  for (year in years) {
    if (definition == "D1") {
      seg <- if (year <= 2019) 0 else 1
    } else if (definition == "D2") {
      seg <- if (year <= 2018) 0 else 1
    } else {
      seg <- 0
    }
    mr <- m[m$Segment == seg, ][1, ]
    fit_rows[[k]] <- data.frame(
      Definition = definition,
      year = year,
      fitted_rate = exp(mr$`Intercept Estimate` + mr$`Slope Estimate` * year)
    )
    k <- k + 1
  }
}
fits <- do.call(rbind, fit_rows)

definition_labels <- c(
  D1 = "D1: strict 8243/3",
  D2 = "D2: 8243/3 + 8245/3",
  D3 = "D3: 8243/3 + 8244/3 + 8245/3"
)
annual$Definition <- factor(annual$Definition, levels = names(definition_labels), labels = definition_labels)
fits$Definition <- factor(fits$Definition, levels = names(definition_labels), labels = definition_labels)

jp <- data.frame(
  Definition = factor(definition_labels[c("D1", "D2")], levels = definition_labels),
  estimate = c(2019, 2018),
  lower = c(2013, 2016),
  upper = c(2021, 2020)
)
annotations <- data.frame(
  Definition = factor(definition_labels, levels = definition_labels),
  x = 2000.4,
  y = c(0.305, 0.305, 0.305),
  label = c(
    "Joinpoint 2019 (95% CI 2013-2021)\nAPC 6.19% to 16.96%",
    "Joinpoint 2018 (95% CI 2016-2020)\nAPC 2.31% to 15.28%",
    "No joinpoint selected\nAPC 4.56%, 2000-2023"
  )
)

pal <- c(
  "D1: strict 8243/3" = "#24557A",
  "D2: 8243/3 + 8245/3" = "#D8873A",
  "D3: 8243/3 + 8244/3 + 8245/3" = "#8E6C9E"
)

p <- ggplot(annual, aes(year, age_adjusted_rate_per_100k, colour = Definition, fill = Definition)) +
  geom_ribbon(aes(ymin = lcl_95, ymax = ucl_95), alpha = 0.10, colour = NA) +
  geom_rect(
    data = jp,
    aes(xmin = lower, xmax = upper, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE, fill = "grey65", alpha = 0.10
  ) +
  geom_vline(
    data = jp, aes(xintercept = estimate), inherit.aes = FALSE,
    colour = "grey35", linetype = "22", linewidth = 0.35
  ) +
  geom_point(shape = 21, size = 1.25, stroke = 0.35, colour = "white") +
  geom_line(
    data = fits, aes(year, fitted_rate, colour = Definition),
    inherit.aes = FALSE, linewidth = 0.85
  ) +
  geom_text(
    data = annotations, aes(x, y, label = label), inherit.aes = FALSE,
    hjust = 0, vjust = 1, size = 2.15, lineheight = 0.95, colour = "#333333"
  ) +
  facet_wrap(~Definition, ncol = 1) +
  scale_colour_manual(values = pal, guide = "none") +
  scale_fill_manual(values = pal, guide = "none") +
  scale_x_continuous(breaks = seq(2000, 2020, 5), limits = c(1999.5, 2023.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 0.33), breaks = seq(0, 0.3, 0.1), expand = c(0, 0)) +
  labs(x = "Year of diagnosis", y = "Age-adjusted incidence per 100,000") +
  theme_classic(base_size = 7.2, base_family = "Arial") +
  theme(
    axis.line = element_line(linewidth = 0.35),
    axis.ticks = element_line(linewidth = 0.35),
    strip.background = element_rect(fill = "#EEF3F7", colour = NA),
    strip.text = element_text(face = "bold", size = 7.2, hjust = 0),
    panel.spacing = grid::unit(1.8, "mm"),
    plot.margin = margin(3, 4, 3, 3)
  )

stem <- "Figure2_joinpoint_trends_by_case_definition"
w <- 183 / 25.4
h <- 155 / 25.4
tmp_dir <- file.path(tempdir(), "gca_joinpoint_figure")
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

svglite::svglite(file.path(tmp_dir, paste0(stem, ".svg")), width = w, height = h)
print(p)
dev.off()
grDevices::cairo_pdf(file.path(tmp_dir, paste0(stem, ".pdf")), width = w, height = h, family = "Arial")
print(p)
dev.off()
grDevices::png(file.path(tmp_dir, paste0(stem, ".png")), width = 3600, height = 3050, res = 500, type = "cairo")
print(p)
dev.off()

for (ext in c("svg", "pdf", "png")) {
  ok <- file.copy(file.path(tmp_dir, paste0(stem, ".", ext)), file.path(figure_dir, paste0(stem, ".", ext)), overwrite = TRUE)
  if (!ok) stop("Could not copy ", ext, " figure output")
}

write.csv(
  merge(annual, fits, by = c("Definition", "year"), all.x = TRUE),
  file.path(figure_dir, "Figure2_joinpoint_trends_source_data.csv"),
  row.names = FALSE
)

cat("Joinpoint publication figure created in", figure_dir, "\n")
