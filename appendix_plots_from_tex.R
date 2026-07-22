#!/usr/bin/env Rscript
# Generate 10 appendix RMSE plots purely from appendix_tables.tex.
#
# Extracts RMSE for all 6 methods (BRITS, SAITS, CSDI, GAIN, notMIWAE, NIMIWAE)
# from the four per-distribution tables (A1–A4) in appendix_tables.tex and
# plots them using the same ggplot2 style as appendix_tables.R.
#
# Output (10 plots):
#   figs/appendix_rmse_mar.pdf           overview, facet_grid(dist ~ N_label)
#   figs/appendix_rmse_mnar.pdf          overview, facet_grid(dist ~ N_label)
#   figs/appendix_rmse_pp1_mar.pdf       per-dist, facet_wrap(~ N_label)
#   figs/appendix_rmse_pp1_mnar.pdf
#   figs/appendix_rmse_pp2_mar.pdf
#   figs/appendix_rmse_pp2_mnar.pdf
#   figs/appendix_rmse_pp3_mar.pdf
#   figs/appendix_rmse_pp3_mnar.pdf
#   figs/appendix_rmse_pp4_mar.pdf
#   figs/appendix_rmse_pp4_mnar.pdf
#
# Usage:
#   Rscript appendix_plots_from_tex.R

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# ── Configuration ──────────────────────────────────────────────────────────────

tex_path <- "appendix_tables.tex"

pp_labels <- c("1" = "Normal", "2" = "Log-normal", "3" = "t", "4" = "Outlier")

mech_info <- data.frame(
  mech_raw  = c("MAR-5%",  "MAR-15%",  "MAR-30%",
                "MNAR-5%", "MNAR-15%", "MNAR-30%"),
  miss_rate = c(0.05, 0.15, 0.30, 0.05, 0.15, 0.30),
  mech_type = c("MAR", "MAR", "MAR", "MNAR", "MNAR", "MNAR"),
  stringsAsFactors = FALSE
)

N_levels     <- c(100, 200, 500, 1000, 5000, 10000)
method_order <- c("BRITS", "SAITS", "CSDI", "GAIN", "notMIWAE", "NIMIWAE")

# ── Parse appendix_tables.tex ──────────────────────────────────────────────────

parse_tex_tables <- function(path) {
  lines <- readLines(path)

  # Locate each table by \label{tab:rmse_ppN}
  table_starts <- grep("\\\\label\\{tab:rmse_pp([1-4])\\}", lines, perl = TRUE)
  stopifnot(length(table_starts) == 4)

  table_ends <- integer(4)
  for (i in seq_along(table_starts)) {
    end_candidates <- grep("\\\\end\\{table\\}", lines)
    table_ends[i] <- min(end_candidates[end_candidates > table_starts[i]])
  }

  results <- list()

  for (tbl_idx in seq_along(table_starts)) {
    pp_val <- as.integer(sub(".*pp([1-4]).*", "\\1",
                             lines[table_starts[tbl_idx]]))
    tbl_lines <- lines[table_starts[tbl_idx]:table_ends[tbl_idx]]

    current_mech <- NA_character_

    for (line in tbl_lines) {
      if (!grepl("&.*&.*\\\\\\\\", line)) next

      # Extract mechanism from \multirow{6}{*}{MECH-X%} lines
      if (grepl("\\\\multirow", line)) {
        m <- regmatches(line,
                        regexpr("\\{\\*\\}\\{([^}]+)\\}", line, perl = TRUE))
        if (length(m) > 0) {
          current_mech <- gsub("\\{\\*\\}\\{", "", m)
          current_mech <- gsub("\\}$", "", current_mech)
          current_mech <- gsub("\\\\%", "%", current_mech)
        }
      }

      # Skip header line
      if (grepl("BRITS.*SAITS", line)) next

      # Clean LaTeX formatting to extract numeric values
      clean <- line
      clean <- gsub("\\\\textbf\\{([0-9.]+)\\}", "\\1", clean)
      clean <- gsub("\\\\multirow\\{[^}]+\\}\\{\\*\\}\\{[^}]+\\}", "", clean)
      clean <- gsub("\\$", "", clean)
      # Remove thousand-separator commas and braces (but not decimal dots)
      clean <- gsub(",", "", clean)
      clean <- gsub("[{}]", "", clean)
      clean <- trimws(clean)

      # Split by & and parse
      parts <- trimws(strsplit(clean, "&")[[1]])
      parts <- gsub("\\\\\\\\", "", parts)
      parts <- trimws(parts)

      if (length(parts) < 7) next

      # parts[1] is either empty (continuation) or the mechanism text from multirow residue
      # parts[2] is the N value
      n_val <- suppressWarnings(as.numeric(parts[2]))
      if (is.na(n_val)) next

      vals <- suppressWarnings(as.numeric(parts[3:8]))
      if (any(is.na(vals))) next

      for (j in seq_along(method_order)) {
        results[[length(results) + 1]] <- data.frame(
          pp       = pp_val,
          dist     = pp_labels[as.character(pp_val)],
          N        = n_val,
          mech     = current_mech,
          method   = method_order[j],
          rmse     = vals[j],
          stringsAsFactors = FALSE
        )
      }
    }
  }

  do.call(rbind, results)
}

cat("Parsing appendix_tables.tex ...\n")
rmse_df <- parse_tex_tables(tex_path)
expected <- 4 * 6 * 6 * 6
cat(sprintf("Extracted %d rows (expected %d = 4 pp x 6 mech x 6 N x 6 methods)\n",
            nrow(rmse_df), expected))
stopifnot(nrow(rmse_df) == expected)
stopifnot(all(method_order %in% unique(rmse_df$method)))

# ── Prepare plot data ──────────────────────────────────────────────────────────

miss_breaks   <- c(0.05, 0.15, 0.30)
miss_labels   <- c("5%", "15%", "30%")
N_label_levels <- c("N = 100", "N = 200", "N = 500",
                    "N = 1k", "N = 5k", "N = 10k")
dist_levels   <- c("Normal", "Log-normal", "t", "Outlier")

plot_df <- rmse_df %>%
  left_join(mech_info, by = c("mech" = "mech_raw")) %>%
  mutate(
    N_label = factor(case_when(
      N == 100   ~ "N = 100",
      N == 200   ~ "N = 200",
      N == 500   ~ "N = 500",
      N == 1000  ~ "N = 1k",
      N == 5000  ~ "N = 5k",
      N == 10000 ~ "N = 10k"
    ), levels = N_label_levels),
    method = factor(method, levels = method_order),
    dist   = factor(dist,   levels = dist_levels)
  )

# ── Shared aesthetics (matching appendix_tables.R) ─────────────────────────────

theme_mda <- function(legend_pos = "bottom") {
  theme_bw() +
    theme(
      legend.position      = legend_pos,
      legend.text          = element_text(size = 8),
      legend.title         = element_blank(),
      legend.background    = element_blank(),
      legend.key.width     = unit(0.8, "cm"),
      legend.spacing.x     = unit(0.2, "cm"),
      legend.margin        = margin(t = 0, r = 0, b = 0, l = 0),
      strip.text           = element_text(size = 9, face = "bold"),
      axis.title           = element_text(size = 10, face = "bold"),
      axis.text            = element_text(size = 8),
      panel.grid.minor     = element_blank(),
      legend.direction     = "horizontal"
    )
}

method_colors <- c(
  "BRITS"      = "#377EB8",
  "SAITS"      = "#00BFC4",
  "CSDI"       = "#4DAF4A",
  "GAIN"       = "#984EA3",
  "notMIWAE"   = "#000000",
  "NIMIWAE"    = "#F0E442"
)

method_linetypes <- c(
  "BRITS"      = "solid",
  "SAITS"      = "solid",
  "CSDI"       = "solid",
  "GAIN"       = "solid",
  "notMIWAE"   = "dotted",
  "NIMIWAE"    = "solid"
)

scale_method_aes <- list(
  scale_color_manual(values = method_colors, breaks = method_order),
  scale_linetype_manual(values = method_linetypes, breaks = method_order),
  guides(color = guide_legend(nrow = 1), linetype = guide_legend(nrow = 1))
)

dir.create("figs", showWarnings = FALSE)

# ── Overview plots: MAR and MNAR, facet_grid(dist ~ N_label) ───────────────────

cat("\nGenerating overview plots ...\n")
for (mtype in c("MAR", "MNAR")) {
  fig <- plot_df %>% filter(mech_type == mtype)
  stopifnot(nrow(fig) > 0)

  p <- ggplot(fig, aes(x = miss_rate, y = rmse, group = method)) +
    geom_line(aes(linetype = method, color = method), linewidth = 0.6) +
    geom_point(aes(color = method), size = 1.2) +
    scale_x_continuous(breaks = miss_breaks, labels = miss_labels) +
    scale_method_aes +
    facet_grid(dist ~ N_label, scales = "free_y") +
    labs(x = "Missingness Rate", y = "RMSE",
         title = paste0("Imputation RMSE \u2014 ", mtype)) +
    theme_mda()

  fname <- sprintf("figs/appendix_rmse_%s.pdf", tolower(mtype))
  ggsave(fname, p, width = 30, height = 20, units = "cm")
  cat(sprintf("  Saved \u2192 %s\n", fname))
}

# ── Per-distribution plots: 4 dists x 2 mech types = 8 plots ──────────────────

cat("\nGenerating per-distribution plots ...\n")
for (pp_val in 1:4) {
  pp_name <- dist_levels[pp_val]
  for (mtype in c("MAR", "MNAR")) {
    fig <- plot_df %>% filter(pp == pp_val, mech_type == mtype)
    stopifnot(nrow(fig) > 0)

    p <- ggplot(fig, aes(x = miss_rate, y = rmse, group = method)) +
      geom_line(aes(linetype = method, color = method), linewidth = 0.6) +
      geom_point(aes(color = method), size = 1.2) +
      scale_x_continuous(breaks = miss_breaks, labels = miss_labels) +
      scale_method_aes +
      facet_wrap(~ N_label, ncol = 3, scales = "free_y") +
      labs(x = "Missingness Rate", y = "RMSE",
           title = paste0(pp_name, " \u2014 ", mtype)) +
      theme_mda()

    fname <- sprintf("figs/appendix_rmse_pp%d_%s.pdf",
                     pp_val, tolower(mtype))
    ggsave(fname, p, width = 24, height = 14, units = "cm")
    cat(sprintf("  Saved \u2192 %s\n", fname))
  }
}

cat("\nDone. 10 plots saved to figs/\n")
