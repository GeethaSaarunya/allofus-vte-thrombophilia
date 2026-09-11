# ======================================================================
# All of Us VTE Final Paper — Journal of Thrombosis and Haemostasis
# Script 11 (LOCAL R) — Render Figures 1–3 from frozen Script 17 outputs
# Version: 11_R_v2_3_20260911
#
# PURPOSE
# -------
# Render publication-ready Figures 1, 2, and 3 from frozen Script 17
# aggregate outputs, with captions and an audit trail.
#
# IMPORTANT
# ---------
# - Rendering only: no statistical model is fitted in this script.
# - No participant-level data are read or written.
# - Figure 1 is a cohort-flow diagram built from frozen manuscript counts.
# - Figure 2 is rendered from the frozen targeted-variant results table.
# - Figure 3 is rendered from the frozen Figure 3 aggregate table.
# - Frozen values are validated before any figure is produced.
# - v2.1 fixes Figure 1 exclusion-box text wrapping and Figure 2
#   left/right text clipping without changing any scientific values.
# - v2.2 uses the Okabe-Ito color-blind-safe palette for plotted data.
# - v2.3 uses the SAME warm Okabe-Ito two-color combination across Figures 1–3.
# ======================================================================


# ----------------------------------------------------------------------
# 0. Controls
# ----------------------------------------------------------------------

SCRIPT_VERSION <- "11_R_v2_3_20260911"

WORKING_DIRECTORY <- normalizePath(
  Sys.getenv("VTE_REPO_ROOT", unset = getwd()),
  mustWork = FALSE
)

SCRIPT17_ZIP <- Sys.getenv(
  "VTE_RESULTS17_ZIP",
  unset = file.path(
    WORKING_DIRECTORY,
    "data",
    "aggregate",
    "VTE_Results_17_20260906.zip"
  )
)

INPUT_DIR <- file.path(
  WORKING_DIRECTORY,
  "outputs",
  "render_input"
)

OUTPUT_ROOT <- file.path(
  WORKING_DIRECTORY,
  "outputs",
  "figures"
)

FIGURE_DIR <- file.path(OUTPUT_ROOT, "figures", "main")
CAPTION_DIR <- file.path(OUTPUT_ROOT, "captions")
LOG_DIR <- file.path(OUTPUT_ROOT, "logs")

PACKAGE_ZIP <- file.path(
  OUTPUT_ROOT,
  "VTE_JTH_figures.zip"
)

OPEN_OUTPUT_DIRECTORY <- FALSE


# ----------------------------------------------------------------------
# 1. Packages
# ----------------------------------------------------------------------

required_packages <- c(
  "ggplot2",
  "ragg",
  "svglite",
  "zip",
  "grid"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0L) {
  stop(
    "Missing required R package(s): ",
    paste(missing_packages, collapse = ", "),
    "\nInstall once with: install.packages(c(",
    paste(sprintf('"%s"', missing_packages), collapse = ", "),
    "))"
  )
}


# ----------------------------------------------------------------------
# 2. Helpers
# ----------------------------------------------------------------------

find_unique_file <- function(root, filename) {
  x <- list.files(
    root,
    pattern = paste0("^", filename, "$"),
    recursive = TRUE,
    full.names = TRUE
  )
  if (length(x) != 1L) {
    stop("Expected exactly one ", filename, "; found ", length(x), ".")
  }
  x[[1]]
}

check_columns <- function(x, required, label) {
  missing <- setdiff(required, names(x))
  if (length(missing) > 0L) {
    stop(label, " is missing: ", paste(missing, collapse = ", "))
  }
}

check_files <- function(paths) {
  bad <- paths[
    !file.exists(paths) |
      is.na(file.info(paths)$size) |
      file.info(paths)$size <= 0
  ]
  if (length(bad) > 0L) {
    stop("Missing/empty output(s):\n  ", paste(bad, collapse = "\n  "))
  }
}

format_p <- function(x) {
  if (is.na(x)) return("")
  if (x < 0.001) return("P<0.001")
  paste0("P=", sprintf("%.3f", x))
}

save_ggplot_bundle <- function(
    plot,
    stem,
    width,
    height
) {
  paths <- c(
    png = paste0(stem, ".png"),
    pdf = paste0(stem, ".pdf"),
    svg = paste0(stem, ".svg"),
    tiff = paste0(stem, ".tiff")
  )

  ragg::agg_png(
    filename = paths[["png"]],
    width = width,
    height = height,
    units = "in",
    res = 450,
    background = "white"
  )
  print(plot)
  grDevices::dev.off()

  grDevices::pdf(
    file = paths[["pdf"]],
    width = width,
    height = height,
    family = "Helvetica",
    useDingbats = FALSE
  )
  print(plot)
  grDevices::dev.off()

  svglite::svglite(
    file = paths[["svg"]],
    width = width,
    height = height,
    bg = "white"
  )
  print(plot)
  grDevices::dev.off()

  ragg::agg_tiff(
    filename = paths[["tiff"]],
    width = width,
    height = height,
    units = "in",
    res = 450,
    compression = "lzw",
    background = "white"
  )
  print(plot)
  grDevices::dev.off()

  check_files(paths)
  paths
}

save_grid_bundle <- function(
    draw_fun,
    stem,
    width,
    height
) {
  paths <- c(
    png = paste0(stem, ".png"),
    pdf = paste0(stem, ".pdf"),
    svg = paste0(stem, ".svg"),
    tiff = paste0(stem, ".tiff")
  )

  ragg::agg_png(
    filename = paths[["png"]],
    width = width,
    height = height,
    units = "in",
    res = 450,
    background = "white"
  )
  grid::grid.newpage()
  draw_fun()
  grDevices::dev.off()

  grDevices::pdf(
    file = paths[["pdf"]],
    width = width,
    height = height,
    family = "Helvetica",
    useDingbats = FALSE
  )
  grid::grid.newpage()
  draw_fun()
  grDevices::dev.off()

  svglite::svglite(
    file = paths[["svg"]],
    width = width,
    height = height,
    bg = "white"
  )
  grid::grid.newpage()
  draw_fun()
  grDevices::dev.off()

  ragg::agg_tiff(
    filename = paths[["tiff"]],
    width = width,
    height = height,
    units = "in",
    res = 450,
    compression = "lzw",
    background = "white"
  )
  grid::grid.newpage()
  draw_fun()
  grDevices::dev.off()

  check_files(paths)
  paths
}


# ----------------------------------------------------------------------
# 3. Theme and style
# ----------------------------------------------------------------------

DARK <- "#222222"

# Okabe-Ito color-blind-safe palette
OKABE_ITO <- c(
  orange = "#E69F00",
  sky_blue = "#56B4E9",
  bluish_green = "#009E73",
  yellow = "#F0E442",
  blue = "#0072B2",
  vermillion = "#D55E00",
  reddish_purple = "#CC79A7",
  black = "#000000"
)

# ------------------------------------------------------------------
# Warm two-color palette used consistently across ALL THREE figures.
# ------------------------------------------------------------------
PRIMARY_WARM <- OKABE_ITO[["reddish_purple"]]  # #CC79A7
SECONDARY_WARM <- OKABE_ITO[["orange"]]         # #E69F00

# Light fills for Figure 1 boxes; saturated versions are used for
# Figure 2 markers/CIs and Figure 3 bars.
PRIMARY_WARM_LIGHT <- grDevices::adjustcolor(
  PRIMARY_WARM,
  alpha.f = 0.18
)
SECONDARY_WARM_LIGHT <- grDevices::adjustcolor(
  SECONDARY_WARM,
  alpha.f = 0.18
)

# Compatibility aliases used later in the script.
BLUE <- PRIMARY_WARM
CARRIER_COLOR <- PRIMARY_WARM
NONCARRIER_COLOR <- SECONDARY_WARM

journal_theme <- function(base_size = 10) {
  ggplot2::theme_classic(
    base_size = base_size,
    base_family = "sans"
  ) +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(color = DARK),
      axis.text = ggplot2::element_text(color = DARK),
      legend.title = ggplot2::element_blank(),
      legend.position = "bottom",
      plot.margin = ggplot2::margin(12, 18, 12, 18)
    )
}


# ----------------------------------------------------------------------
# 4. Folders and frozen source extraction
# ----------------------------------------------------------------------

dir.create(
  WORKING_DIRECTORY,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(SCRIPT17_ZIP)) {
  stop("Frozen Script 17 ZIP not found:\n  ", SCRIPT17_ZIP)
}

for (d in c(INPUT_DIR, OUTPUT_ROOT)) {
  if (dir.exists(d)) {
    unlink(d, recursive = TRUE, force = TRUE)
  }
}

for (d in c(INPUT_DIR, FIGURE_DIR, CAPTION_DIR, LOG_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

utils::unzip(
  SCRIPT17_ZIP,
  exdir = INPUT_DIR,
  overwrite = TRUE
)

TABLE3_PATH <- find_unique_file(
  INPUT_DIR,
  "17D_Main_Table3_targeted_variants.csv"
)

FIGURE3_DATA_PATH <- find_unique_file(
  INPUT_DIR,
  "17E_Main_Figure3_data.csv"
)

FIGURE3_TEST_PATH <- find_unique_file(
  INPUT_DIR,
  "17E_Main_Figure3_chisquare.csv"
)

table3 <- utils::read.csv(
  TABLE3_PATH,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

figure3_data <- utils::read.csv(
  FIGURE3_DATA_PATH,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

figure3_tests <- utils::read.csv(
  FIGURE3_TEST_PATH,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


# ----------------------------------------------------------------------
# 5. Validation
# ----------------------------------------------------------------------

# Figure 1 frozen counts from the final manuscript.
fig1_counts <- list(
  eligible_ehr = 483707,
  excluded_no_short_read_wgs = 125174,
  analysis_cohort = 358533,
  observed_vte = 7553,
  no_qualifying_vte = 350980
)

if (
  fig1_counts$analysis_cohort +
    fig1_counts$excluded_no_short_read_wgs !=
    fig1_counts$eligible_ehr
) {
  stop("Figure 1 count lock failed: eligible != excluded + analysis.")
}

if (
  fig1_counts$observed_vte +
    fig1_counts$no_qualifying_vte !=
    fig1_counts$analysis_cohort
) {
  stop("Figure 1 count lock failed: analysis != VTE + no VTE.")
}

# Figure 2 source validation.
check_columns(
  table3,
  c(
    "marker_label",
    "odds_ratio",
    "ci_95_low",
    "ci_95_high",
    "adjusted_p",
    "focused_role"
  ),
  "Table 3"
)

table3$odds_ratio <- as.numeric(table3$odds_ratio)
table3$ci_95_low <- as.numeric(table3$ci_95_low)
table3$ci_95_high <- as.numeric(table3$ci_95_high)

if (nrow(table3) != 7L) {
  stop("Figure 2 source must contain exactly 7 targeted variants.")
}

expected_markers <- c(
  "Factor V Leiden (F5 rs6025)",
  "Prothrombin G20210A (F2 rs1799963)",
  "ABO rs8176719 ALT allele",
  "F11 rs2036914 ALT allele",
  "F11 rs2289252 ALT allele",
  "FGG rs2066865 ALT allele",
  "PROCR rs867186 ALT allele"
)

if (!setequal(table3$marker_label, expected_markers)) {
  stop("Figure 2 marker-label validation failed.")
}

# Figure 3 source validation.
check_columns(
  figure3_data,
  c(
    "burden_category",
    "carrier",
    "n",
    "vte_n",
    "vte_percent",
    "wilson_95_low_percent",
    "wilson_95_high_percent"
  ),
  "Figure 3 data"
)

check_columns(
  figure3_tests,
  c(
    "burden_category",
    "pearson_chi2",
    "df",
    "p_value"
  ),
  "Figure 3 chi-square data"
)

figure3_data$burden_category <- as.character(figure3_data$burden_category)
figure3_data$carrier <- as.integer(figure3_data$carrier)
figure3_data$n <- as.numeric(figure3_data$n)
figure3_data$vte_n <- as.numeric(figure3_data$vte_n)
figure3_data$vte_percent <- as.numeric(figure3_data$vte_percent)
figure3_data$wilson_95_low_percent <- as.numeric(figure3_data$wilson_95_low_percent)
figure3_data$wilson_95_high_percent <- as.numeric(figure3_data$wilson_95_high_percent)

figure3_tests$burden_category <- as.character(figure3_tests$burden_category)
figure3_tests$p_value <- as.numeric(figure3_tests$p_value)

if (
  nrow(figure3_data) != 6L ||
  any(!figure3_data$burden_category %in% c("0", "1", "2+")) ||
  any(!figure3_data$carrier %in% c(0L, 1L))
) {
  stop("Figure 3 source failed structural validation.")
}

frozen_cells <- data.frame(
  burden_category = c("0", "0", "1", "1", "2+", "2+"),
  carrier = c(0L, 1L, 0L, 1L, 0L, 1L),
  expected_n = c(208024, 12912, 64496, 4007, 64830, 4264),
  expected_vte_n = c(1475, 215, 1701, 204, 3589, 369),
  expected_percent = c(0.71, 1.67, 2.64, 5.09, 5.54, 8.65)
)

locked <- merge(
  figure3_data,
  frozen_cells,
  by = c("burden_category", "carrier"),
  all.x = TRUE,
  sort = FALSE
)

if (
  any(locked$n != locked$expected_n) ||
  any(locked$vte_n != locked$expected_vte_n) ||
  any(abs(locked$vte_percent - locked$expected_percent) > 0.02)
) {
  stop("Figure 3 frozen cell lock failed.")
}

cat("Frozen Figure 1, Figure 2, and Figure 3 validation: PASS\n")


# ----------------------------------------------------------------------
# 6. Render Figure 1 — cohort flow
# ----------------------------------------------------------------------

draw_figure1 <- function() {
  fmt <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)

  draw_box <- function(
      x, y, w, h, lines,
      line_heights = NULL,
      fill_color = "white",
      border_color = DARK
  ) {
    if (is.null(line_heights)) {
      line_heights <- seq(0.68, 0.18, length.out = length(lines))
    }

    grid::grid.rect(
      x = x, y = y, width = w, height = h,
      gp = grid::gpar(
        col = border_color,
        fill = fill_color,
        lwd = 1.2
      )
    )

    for (i in seq_along(lines)) {
      line <- lines[[i]]
      fontface <- if (!is.null(line$fontface)) line$fontface else "plain"
      fontsize <- if (!is.null(line$fontsize)) line$fontsize else 16
      grid::grid.text(
        label = line$text,
        x = x,
        y = y - h / 2 + h * line_heights[[i]],
        gp = grid::gpar(
          fontsize = fontsize,
          fontface = fontface,
          col = "black"
        )
      )
    }
  }

  grid::pushViewport(
    grid::viewport(x = 0.5, y = 0.5, width = 1, height = 1)
  )

  # Top box
  draw_box(
    x = 0.50, y = 0.86, w = 0.72, h = 0.16,
    lines = list(
      list(text = "Eligible All of Us EHR cohort", fontface = "bold", fontsize = 20),
      list(text = "Adults with pre-enrollment clinical-risk ascertainment", fontsize = 16),
      list(text = paste0("n = ", fmt(fig1_counts$eligible_ehr)), fontface = "bold", fontsize = 18)
    ),
    fill_color = PRIMARY_WARM_LIGHT,
    border_color = PRIMARY_WARM
  )

  # Right exclusion box.
  # The title is intentionally wrapped so every line stays inside the box.
  draw_box(
    x = 0.78, y = 0.58, w = 0.30, h = 0.19,
    lines = list(
      list(text = "Excluded from", fontface = "bold", fontsize = 16),
      list(text = "WGS analysis", fontface = "bold", fontsize = 16),
      list(text = "No qualifying short-read WGS", fontsize = 13),
      list(text = "for targeted genetic analyses", fontsize = 13),
      list(
        text = paste0("n = ", fmt(fig1_counts$excluded_no_short_read_wgs)),
        fontface = "bold",
        fontsize = 17
      )
    ),
    line_heights = c(0.80, 0.66, 0.47, 0.34, 0.14),
    fill_color = SECONDARY_WARM_LIGHT,
    border_color = SECONDARY_WARM
  )

  # Middle box
  draw_box(
    x = 0.50, y = 0.40, w = 0.56, h = 0.16,
    lines = list(
      list(text = "WGS analysis cohort", fontface = "bold", fontsize = 20),
      list(text = "Linked EHR + qualifying short-read WGS", fontsize = 16),
      list(text = paste0("n = ", fmt(fig1_counts$analysis_cohort)), fontface = "bold", fontsize = 18)
    ),
    fill_color = PRIMARY_WARM_LIGHT,
    border_color = PRIMARY_WARM
  )

  # Bottom left box
  draw_box(
    x = 0.25, y = 0.11, w = 0.32, h = 0.16,
    lines = list(
      list(text = "Observed VTE", fontface = "bold", fontsize = 20),
      list(text = paste0("n = ", fmt(fig1_counts$observed_vte)), fontface = "bold", fontsize = 18)
    ),
    line_heights = c(0.62, 0.18),
    fill_color = PRIMARY_WARM_LIGHT,
    border_color = PRIMARY_WARM
  )

  # Bottom right box
  draw_box(
    x = 0.68, y = 0.11, w = 0.32, h = 0.16,
    lines = list(
      list(text = "No qualifying VTE", fontface = "bold", fontsize = 20),
      list(text = "during follow-up", fontsize = 16),
      list(text = paste0("n = ", fmt(fig1_counts$no_qualifying_vte)), fontface = "bold", fontsize = 18)
    ),
    line_heights = c(0.64, 0.42, 0.18),
    fill_color = SECONDARY_WARM_LIGHT,
    border_color = SECONDARY_WARM
  )

  # Arrows
  arrow_gp <- grid::gpar(col = "black", lwd = 1.5)

  grid::grid.lines(
    x = grid::unit(c(0.50, 0.50), "npc"),
    y = grid::unit(c(0.78, 0.48), "npc"),
    gp = arrow_gp,
    arrow = grid::arrow(type = "closed", length = grid::unit(0.12, "inches"))
  )

  # Branch to the exclusion box stops at the box border rather than
  # entering the text area.
  grid::grid.lines(
    x = grid::unit(c(0.50, 0.625), "npc"),
    y = grid::unit(c(0.61, 0.61), "npc"),
    gp = arrow_gp,
    arrow = grid::arrow(type = "closed", length = grid::unit(0.10, "inches"))
  )

  grid::grid.lines(
    x = grid::unit(c(0.50, 0.31), "npc"),
    y = grid::unit(c(0.32, 0.19), "npc"),
    gp = arrow_gp,
    arrow = grid::arrow(type = "closed", length = grid::unit(0.12, "inches"))
  )

  grid::grid.lines(
    x = grid::unit(c(0.50, 0.59), "npc"),
    y = grid::unit(c(0.32, 0.19), "npc"),
    gp = arrow_gp,
    arrow = grid::arrow(type = "closed", length = grid::unit(0.12, "inches"))
  )

  grid::popViewport()
}

figure1_paths <- save_grid_bundle(
  draw_fun = draw_figure1,
  stem = file.path(FIGURE_DIR, "Figure1_CohortFlow_v2_3"),
  width = 8.6,
  height = 10.6
)

FIGURE1_CAPTION <- paste0(
  "Figure 1. Cohort flow diagram. Adults with pre-enrollment clinical-risk ",
  "ascertainment in the eligible All of Us EHR cohort were screened for the ",
  "availability of qualifying short-read whole-genome sequencing (WGS) for ",
  "targeted genetic analyses. The final WGS analysis cohort comprised 358,533 ",
  "participants, including 7,553 participants with observed VTE during follow-up ",
  "and 350,980 participants without a qualifying VTE record during follow-up."
)

writeLines(
  FIGURE1_CAPTION,
  file.path(CAPTION_DIR, "Figure1_caption.txt"),
  useBytes = TRUE
)


# ----------------------------------------------------------------------
# 7. Render Figure 2 — targeted thrombosis-related variants
# ----------------------------------------------------------------------

figure2_data <- table3

label_map <- c(
  "Factor V Leiden (F5 rs6025)" = "F5 rs6025\n(Factor V Leiden)",
  "Prothrombin G20210A (F2 rs1799963)" = "F2 rs1799963\n(prothrombin G20210A)",
  "ABO rs8176719 ALT allele" = "ABO rs8176719",
  "F11 rs2036914 ALT allele" = "F11 rs2036914",
  "F11 rs2289252 ALT allele" = "F11 rs2289252",
  "FGG rs2066865 ALT allele" = "FGG rs2066865",
  "PROCR rs867186 ALT allele" = "PROCR rs867186"
)

desired_order <- c(
  "Factor V Leiden (F5 rs6025)",
  "Prothrombin G20210A (F2 rs1799963)",
  "ABO rs8176719 ALT allele",
  "F11 rs2036914 ALT allele",
  "F11 rs2289252 ALT allele",
  "FGG rs2066865 ALT allele",
  "PROCR rs867186 ALT allele"
)

figure2_data$order <- match(figure2_data$marker_label, desired_order)
figure2_data <- figure2_data[order(figure2_data$order), ]

figure2_data$marker_display <- unname(label_map[figure2_data$marker_label])
figure2_data$right_label <- sprintf(
  "%.2f (%.2f–%.2f)",
  figure2_data$odds_ratio,
  figure2_data$ci_95_low,
  figure2_data$ci_95_high
)
figure2_data$y <- rev(seq_len(nrow(figure2_data)))
figure2_data$is_highlight <- figure2_data$marker_label %in% c(
  "Factor V Leiden (F5 rs6025)",
  "Prothrombin G20210A (F2 rs1799963)"
)

figure2_data$warm_group <- factor(
  ifelse(
    figure2_data$is_highlight,
    "F5/F2 focus",
    "Other targeted variants"
  ),
  levels = c(
    "F5/F2 focus",
    "Other targeted variants"
  )
)

# Text columns intentionally sit outside the forest-plot panel.
# coord_cartesian(clip = "off") and generous plot margins keep them visible.
left_x <- 0.69
right_x <- 3.12

figure2 <- ggplot2::ggplot(
  figure2_data,
  ggplot2::aes(x = odds_ratio, y = y)
) +
  ggplot2::geom_vline(
    xintercept = 1.0,
    linetype = "dashed",
    linewidth = 0.8,
    color = DARK
  ) +
  ggplot2::geom_errorbarh(
    ggplot2::aes(
      xmin = ci_95_low,
      xmax = ci_95_high,
      color = warm_group
    ),
    height = 0,
    linewidth = 1.0
  ) +
  ggplot2::geom_point(
    ggplot2::aes(color = warm_group),
    size = 3.8
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      x = left_x,
      label = marker_display,
      fontface = ifelse(is_highlight, "bold", "plain")
    ),
    hjust = 1,
    size = 4.0,
    lineheight = 0.95,
    color = "black",
    show.legend = FALSE
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      x = right_x,
      label = right_label,
      fontface = ifelse(is_highlight, "bold", "plain")
    ),
    hjust = 0,
    size = 4.0,
    color = "black",
    show.legend = FALSE
  ) +
  ggplot2::annotate(
    "text",
    x = right_x,
    y = max(figure2_data$y) + 1.0,
    label = "Adjusted OR (95% CI)",
    hjust = 0,
    fontface = "bold",
    size = 4.3
  ) +
  ggplot2::scale_color_manual(
    values = c(
      "F5/F2 focus" = PRIMARY_WARM,
      "Other targeted variants" = SECONDARY_WARM
    ),
    guide = "none"
  ) +
  ggplot2::scale_x_continuous(
    breaks = c(0.75, 1.00, 1.50, 2.00, 3.00)
  ) +
  ggplot2::coord_cartesian(
    xlim = c(0.75, 3.00),
    clip = "off"
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0.3, max(figure2_data$y) + 1.2),
    breaks = NULL
  ) +
  ggplot2::labs(
    x = "Adjusted odds ratio",
    y = NULL
  ) +
  ggplot2::theme_classic(
    base_size = 10,
    base_family = "sans"
  ) +
  ggplot2::theme(
    axis.title.x = ggplot2::element_text(size = 10.5, color = DARK),
    axis.text.x = ggplot2::element_text(size = 9.2, color = DARK),
    axis.line.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank(),
    legend.position = "none",
    # Extra room is deliberate: marker names and effect-size text are
    # publication columns outside the plotting panel.
    plot.margin = ggplot2::margin(
      t = 18,
      r = 150,
      b = 16,
      l = 165,
      unit = "pt"
    )
  )

figure2_paths <- save_ggplot_bundle(
  plot = figure2,
  stem = file.path(FIGURE_DIR, "Figure2_TargetedVariants_v2_3"),
  width = 11.0,
  height = 6.8
)

FIGURE2_CAPTION <- paste0(
  "Figure 2. Adjusted associations between targeted thrombosis-related variants ",
  "and observed VTE. Points show adjusted odds ratios and horizontal lines show ",
  "95% confidence intervals from separate adjusted models for each targeted variant. ",
  "Factor V Leiden (F5 rs6025) and prothrombin G20210A (F2 rs1799963) had the two ",
  "largest adjusted associations and were therefore carried forward into the focused ",
  "clinical-risk-stratified analyses."
)

writeLines(
  FIGURE2_CAPTION,
  file.path(CAPTION_DIR, "Figure2_caption.txt"),
  useBytes = TRUE
)


# ----------------------------------------------------------------------
# 8. Render Figure 3 — clinical-risk burden by F5/F2 status
# ----------------------------------------------------------------------

burden_levels <- c("0", "1", "2+")
burden_labels <- c("0 factors", "1 factor", "≥2 factors")
DODGE_OFFSET <- 0.20

figure3_data$burden_index <- match(
  figure3_data$burden_category,
  burden_levels
)

# Carrier first (left), matching manuscript presentation.
figure3_data$x_position <- ifelse(
  figure3_data$carrier == 1L,
  figure3_data$burden_index - DODGE_OFFSET,
  figure3_data$burden_index + DODGE_OFFSET
)

figure3_data$variant_label <- factor(
  figure3_data$carrier,
  levels = c(1L, 0L),
  labels = c(
    "F5 rs6025 or F2 rs1799963 carrier",
    "Neither variant"
  )
)

figure3_data$count_label <- paste0(
  sprintf("%.2f%%", figure3_data$vte_percent),
  "\n",
  format(figure3_data$vte_n, big.mark = ",", scientific = FALSE, trim = TRUE),
  "/",
  format(figure3_data$n, big.mark = ",", scientific = FALSE, trim = TRUE)
)

figure3_data$label_y <- figure3_data$wilson_95_high_percent + 0.24

category_upper <- stats::aggregate(
  wilson_95_high_percent ~ burden_index,
  data = figure3_data,
  FUN = max
)

names(category_upper)[2] <- "max_upper_ci"
category_upper$bracket_y <- category_upper$max_upper_ci + 0.82
category_upper$cap_y <- category_upper$bracket_y - 0.14
category_upper$p_y <- category_upper$bracket_y + 0.25

figure3_tests$burden_index <- match(
  figure3_tests$burden_category,
  burden_levels
)

annotation_data <- merge(
  figure3_tests,
  category_upper,
  by = "burden_index",
  all.x = TRUE,
  sort = FALSE
)

annotation_data$p_label <- vapply(
  annotation_data$p_value,
  format_p,
  FUN.VALUE = character(1)
)

annotation_data$x_left <- annotation_data$burden_index - DODGE_OFFSET
annotation_data$x_right <- annotation_data$burden_index + DODGE_OFFSET

max_y <- max(
  annotation_data$p_y,
  figure3_data$label_y,
  na.rm = TRUE
) + 0.45

figure3 <- ggplot2::ggplot() +
  ggplot2::geom_col(
    data = figure3_data,
    ggplot2::aes(
      x = x_position,
      y = vte_percent,
      fill = variant_label
    ),
    width = 0.37,
    color = DARK,
    linewidth = 0.35
  ) +
  ggplot2::geom_errorbar(
    data = figure3_data,
    ggplot2::aes(
      x = x_position,
      ymin = wilson_95_low_percent,
      ymax = wilson_95_high_percent
    ),
    width = 0.10,
    linewidth = 0.65,
    color = DARK
  ) +
  ggplot2::geom_text(
    data = figure3_data,
    ggplot2::aes(
      x = x_position,
      y = label_y,
      label = count_label
    ),
    vjust = 0,
    lineheight = 0.92,
    size = 3.0,
    color = DARK,
    show.legend = FALSE
  ) +
  ggplot2::geom_segment(
    data = annotation_data,
    ggplot2::aes(
      x = x_left,
      xend = x_right,
      y = bracket_y,
      yend = bracket_y
    ),
    linewidth = 0.48,
    color = DARK,
    inherit.aes = FALSE
  ) +
  ggplot2::geom_segment(
    data = annotation_data,
    ggplot2::aes(
      x = x_left,
      xend = x_left,
      y = cap_y,
      yend = bracket_y
    ),
    linewidth = 0.48,
    color = DARK,
    inherit.aes = FALSE
  ) +
  ggplot2::geom_segment(
    data = annotation_data,
    ggplot2::aes(
      x = x_right,
      xend = x_right,
      y = cap_y,
      yend = bracket_y
    ),
    linewidth = 0.48,
    color = DARK,
    inherit.aes = FALSE
  ) +
  ggplot2::geom_text(
    data = annotation_data,
    ggplot2::aes(
      x = burden_index,
      y = p_y,
      label = p_label
    ),
    fontface = "bold",
    size = 3.15,
    color = DARK,
    inherit.aes = FALSE
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "F5 rs6025 or F2 rs1799963 carrier" = CARRIER_COLOR,
      "Neither variant" = NONCARRIER_COLOR
    ),
    breaks = c(
      "F5 rs6025 or F2 rs1799963 carrier",
      "Neither variant"
    )
  ) +
  ggplot2::scale_x_continuous(
    breaks = 1:3,
    labels = burden_labels,
    limits = c(0.55, 3.45),
    expand = ggplot2::expansion(mult = c(0, 0))
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0, max_y),
    breaks = seq(0, ceiling(max_y), by = 2),
    labels = function(x) paste0(x, "%"),
    expand = ggplot2::expansion(mult = c(0, 0))
  ) +
  ggplot2::labs(
    x = "Five-factor clinical-risk burden",
    y = "Observed VTE proportion"
  ) +
  ggplot2::theme_classic(
    base_size = 10,
    base_family = "sans"
  ) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.title = ggplot2::element_text(color = DARK),
    axis.text = ggplot2::element_text(color = DARK),
    axis.text.x = ggplot2::element_text(size = 9.4),
    legend.title = ggplot2::element_blank(),
    legend.position = "bottom",
    legend.text = ggplot2::element_text(size = 8.8),
    plot.margin = ggplot2::margin(12, 22, 12, 18)
  )

figure3_paths <- save_ggplot_bundle(
  plot = figure3,
  stem = file.path(FIGURE_DIR, "Figure3_ClinicalRisk_F5F2_v2_3"),
  width = 8.5,
  height = 5.9
)

FIGURE3_CAPTION <- paste0(
  "Figure 3. Observed VTE occurrence by five-factor clinical-risk burden and ",
  "F5/F2 inherited thrombophilia. Bars show observed, unadjusted VTE proportions ",
  "among participants carrying either F5 rs6025 (Factor V Leiden) or F2 rs1799963 ",
  "(prothrombin G20210A) and among participants carrying neither variant within ",
  "categories containing 0, 1, or ≥2 of the five captured clinical factors: ",
  "cancer history, major surgery, inpatient hospitalization, infection/sepsis, ",
  "and fracture/major trauma. Error bars show Wilson 95% confidence intervals; ",
  "labels show VTE events/group denominators. Brackets and P values show ",
  "unadjusted Pearson chi-square comparisons within each clinical-risk category. ",
  "These comparisons are not adjusted for covariates or person-time and are not ",
  "formal interaction tests. Formal multiplicative and adjusted probability-",
  "difference interaction analyses are reported in Supplementary Table S6."
)

writeLines(
  FIGURE3_CAPTION,
  file.path(CAPTION_DIR, "Figure3_caption.txt"),
  useBytes = TRUE
)


# ----------------------------------------------------------------------
# 9. Audit, session info, and package
# ----------------------------------------------------------------------

source_audit <- data.frame(
  item = c(
    "Script version",
    "Frozen Script 17 ZIP",
    "Figure 1 source",
    "Figure 2 source",
    "Figure 3 data source",
    "Figure 3 chi-square source",
    "Figure 2 source MD5",
    "Figure 3 data MD5",
    "Figure 3 chi-square MD5",
    "Statistical analysis performed here"
  ),
  value = c(
    SCRIPT_VERSION,
    normalizePath(SCRIPT17_ZIP),
    "Frozen manuscript cohort counts",
    basename(TABLE3_PATH),
    basename(FIGURE3_DATA_PATH),
    basename(FIGURE3_TEST_PATH),
    unname(tools::md5sum(TABLE3_PATH)),
    unname(tools::md5sum(FIGURE3_DATA_PATH)),
    unname(tools::md5sum(FIGURE3_TEST_PATH)),
    "No — rendering only"
  ),
  stringsAsFactors = FALSE
)

SOURCE_AUDIT_PATH <- file.path(
  LOG_DIR,
  "11_v2_3_all_figures_source_audit.csv"
)

SESSION_INFO_PATH <- file.path(
  LOG_DIR,
  "11_v2_3_sessionInfo.txt"
)

utils::write.csv(
  source_audit,
  SOURCE_AUDIT_PATH,
  row.names = FALSE
)

writeLines(
  capture.output(utils::sessionInfo()),
  SESSION_INFO_PATH
)

all_outputs <- c(
  unname(figure1_paths),
  unname(figure2_paths),
  unname(figure3_paths),
  file.path(CAPTION_DIR, "Figure1_caption.txt"),
  file.path(CAPTION_DIR, "Figure2_caption.txt"),
  file.path(CAPTION_DIR, "Figure3_caption.txt"),
  SOURCE_AUDIT_PATH,
  SESSION_INFO_PATH
)

check_files(all_outputs)

if (file.exists(PACKAGE_ZIP)) {
  file.remove(PACKAGE_ZIP)
}

zip::zipr(
  zipfile = PACKAGE_ZIP,
  files = list.files(
    OUTPUT_ROOT,
    recursive = TRUE,
    full.names = FALSE
  ),
  root = OUTPUT_ROOT
)

check_files(PACKAGE_ZIP)


# ----------------------------------------------------------------------
# 10. Completion
# ----------------------------------------------------------------------

cat("\n")
cat(strrep("=", 76), "\n", sep = "")
cat("SCRIPT 11 v2.3 FIGURES 1–3 RENDER COMPLETE\n")
cat(strrep("=", 76), "\n", sep = "")
cat("Figure 1 PNG: ", figure1_paths[["png"]], "\n", sep = "")
cat("Figure 2 PNG: ", figure2_paths[["png"]], "\n", sep = "")
cat("Figure 3 PNG: ", figure3_paths[["png"]], "\n", sep = "")
cat("Package ZIP:  ", PACKAGE_ZIP, "\n", sep = "")
cat("\nNo statistical analysis was performed.\n")
cat("All plotted values came from frozen Script 17 aggregate outputs or\n")
cat("from frozen manuscript counts for Figure 1.\n")

if (
  isTRUE(OPEN_OUTPUT_DIRECTORY) &&
  Sys.info()[["sysname"]] == "Darwin"
) {
  try(
    system2("open", shQuote(OUTPUT_ROOT)),
    silent = TRUE
  )
}
