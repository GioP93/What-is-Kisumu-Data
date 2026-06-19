## =====================================================================
## Cross-insecticide Spearman rank correlation of site LC50s
##
## Run this AFTER all four insecticides have been fitted with
## kisumu_one_insecticide.R. It reads each insecticide's per-site LC50
## posteriors back from the output folders and asks: do sites RANK the
## same way (most to least resistant) across insecticides?
##
## Two outputs:
##   1. Point-estimate heatmap (median LC50 per site), matching the
##      spearman_heatmap.png style.
##   2. Uncertainty-propagated version: Spearman R computed across matched
##      posterior draws, giving a median R and 95% credible interval per
##      pair (so you can see whether a correlation is well-determined).
##
## Only sites NOT excluded by the LC50 rule are used (the excluded curves
## have no reliable LC50). Each pair uses only sites with data for both
## insecticides (pairwise-complete).
## =====================================================================

pacman::p_load(here, tidyverse, reshape2, ggplot2)

# ---------------------------------------------------------------------
# 0/ CONFIG -- relative paths; run AFTER scripts/1_variability_analysis.R
# ---------------------------------------------------------------------
# Reads the per-insecticide LC50 posteriors written into /outputs by the
# variability script, then computes the cross-insecticide Spearman
# correlation of site LC50 rankings.
outputs_root <- here::here("outputs")

# Output subfolder name (created by script 1) -> tidy display label.
insecticide_folders <- c(
  "ddt"                = "DDT",
  "alpha_cypermethrin" = "Alpha-cypermethrin",
  "permethrin"         = "Permethrin",
  "pirimiphos_methyl"  = "Pirimiphos methyl"
)

# Fixed display order
ins_order <- c("DDT", "Alpha-cypermethrin", "Permethrin", "Pirimiphos methyl")

spearman_dir <- file.path(outputs_root, "spearman")
dir.create(spearman_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------
# 1/ LOAD PER-SITE LC50 POSTERIORS FROM EACH INSECTICIDE ----
# ---------------------------------------------------------------------
# For each insecticide folder, read every lc50_<site>.csv (one row per
# posterior draw; columns LC, Strain, Insecticide). Files in the top of
# the folder are the INCLUDED (valid) sites; the excluded/ subfolder is
# ignored here on purpose.
read_one_insecticide <- function(folder, label) {
  dir_path <- file.path(outputs_root, folder)
  if (!dir.exists(dir_path)) {
    warning("Folder not found, skipping: ", dir_path)
    return(NULL)
  }
  files <- list.files(dir_path, pattern = "^lc50_.*\\.csv$", full.names = TRUE)
  if (length(files) == 0) {
    warning("No lc50_*.csv files in ", dir_path)
    return(NULL)
  }
  do.call(rbind, lapply(files, function(f) {
    d <- tryCatch(read.csv(f), error = function(e) NULL)
    if (is.null(d) || !all(c("LC", "Strain") %in% names(d))) return(NULL)
    d$Site  <- sub(" \\| .*$", "", d$Strain)
    d$InsLabel <- label
    # add a within-site draw index so draws can be aligned across insecticides
    d <- d %>% group_by(Site) %>% mutate(draw = row_number()) %>% ungroup()
    d[, c("Site", "InsLabel", "draw", "LC")]
  }))
}

all_post <- do.call(rbind, Map(read_one_insecticide,
                               names(insecticide_folders),
                               insecticide_folders))
if (is.null(all_post) || nrow(all_post) == 0) {
  stop("No LC50 posteriors loaded - check outputs_root and that all four ",
       "insecticides have been fitted.")
}
all_post$InsLabel <- factor(all_post$InsLabel, levels = ins_order)

# Report what actually loaded - so a missing/misnamed folder is obvious.
loaded <- levels(droplevels(all_post$InsLabel))
missing <- setdiff(ins_order, loaded)
message("Insecticides loaded: ", paste(loaded, collapse = ", "))
if (length(missing) > 0) {
  warning("NOT loaded (folder missing or empty): ", paste(missing, collapse = ", "),
          "\n  The correlation will only cover the loaded insecticides.",
          "\n  Check the folder names in `insecticide_folders` against what is",
          "\n  actually in: ", outputs_root)
}

# ---------------------------------------------------------------------
# 2/ POINT-ESTIMATE MATRIX (median LC50 per site x insecticide) ----
# ---------------------------------------------------------------------
med_tab <- all_post %>%
  group_by(Site, InsLabel) %>%
  summarise(LC50 = median(LC), .groups = "drop") %>%
  pivot_wider(names_from = InsLabel, values_from = LC50)

mat <- med_tab %>% select(any_of(ins_order)) %>% as.data.frame()
rownames(mat) <- med_tab$Site

# Spearman on the median LC50s
cor_mat <- cor(mat, method = "spearman", use = "pairwise.complete.obs")
cat("\nSpearman correlation (median LC50):\n"); print(round(cor_mat, 2))

# p-values per pair
pmat <- matrix(NA, ncol(mat), ncol(mat),
               dimnames = list(colnames(mat), colnames(mat)))
for (i in seq_len(ncol(mat))) for (j in seq_len(ncol(mat))) {
  ok <- complete.cases(mat[, i], mat[, j])
  if (sum(ok) >= 3)
    pmat[i, j] <- suppressWarnings(
      cor.test(mat[ok, i], mat[ok, j], method = "spearman")$p.value)
}
cat("\nP-values:\n"); print(round(pmat, 3))

write.csv(round(cor_mat, 3),
          file.path(spearman_dir, "spearman_correlation_matrix.csv"))
write.csv(round(pmat, 3),
          file.path(spearman_dir, "spearman_pvalues.csv"))

# ---------------------------------------------------------------------
# 3/ UNCERTAINTY-PROPAGATED SPEARMAN (across posterior draws) ----
# ---------------------------------------------------------------------
# For each pair of insecticides, and each posterior draw, compute the
# Spearman R across the sites shared by both. The spread of R over draws
# is the uncertainty in the rank correlation. Uses the common number of
# draws available across the loaded curves.
#
# NOTE: each curve was fitted independently, so aligning "draw k" across
# insecticides is just a Monte Carlo sample from the (independent) joint
# posterior - it is NOT paired sampling. This correctly propagates each
# LC50's marginal uncertainty into the rank correlation.
n_draw_common <- all_post %>% count(Site, InsLabel) %>% pull(n) %>% min()
n_use <- min(2000, n_draw_common)   # cap for speed
message("Using ", n_use, " posterior draws per curve for uncertainty propagation.")

draw_spearman <- function(insA, insB) {
  wide <- all_post %>%
    filter(InsLabel %in% c(insA, insB), draw <= n_use) %>%
    select(Site, InsLabel, draw, LC) %>%
    pivot_wider(names_from = InsLabel, values_from = LC)
  
  # If either insecticide produced no rows, its column won't exist -> skip pair
  if (!all(c(insA, insB) %in% names(wide))) return(NULL)
  
  # sites with data for both
  wide <- wide[stats::complete.cases(wide[[insA]], wide[[insB]]), ]
  if (nrow(wide) == 0) return(NULL)
  shared_sites <- unique(wide$Site)
  if (length(shared_sites) < 3) return(NULL)  # need >=3 sites for a rank corr
  
  rs <- vapply(sort(unique(wide$draw)), function(dd) {
    w <- wide[wide$draw == dd, ]
    if (nrow(w) < 3) return(NA_real_)
    xa <- w[[insA]]; xb <- w[[insB]]
    if (length(xa) < 3 || length(xb) < 3) return(NA_real_)
    suppressWarnings(cor(xa, xb, method = "spearman"))
  }, numeric(1))
  rs <- rs[is.finite(rs)]
  if (length(rs) == 0) return(NULL)
  data.frame(
    Insecticide_A = insA, Insecticide_B = insB,
    n_sites = length(shared_sites),
    R_median = median(rs),
    R_lower  = quantile(rs, 0.025, names = FALSE),
    R_upper  = quantile(rs, 0.975, names = FALSE),
    prob_positive = mean(rs > 0)
  )
}

pairs <- combn(ins_order, 2, simplify = FALSE)
spear_unc <- do.call(rbind, lapply(pairs, function(p) draw_spearman(p[1], p[2])))
if (!is.null(spear_unc)) {
  cat("\nUncertainty-propagated Spearman R (median [95% CrI]):\n")
  print(spear_unc, row.names = FALSE)
  write.csv(spear_unc,
            file.path(spearman_dir, "spearman_uncertainty.csv"),
            row.names = FALSE)
}

# ---------------------------------------------------------------------
# 4/ HEATMAP 
# ---------------------------------------------------------------------
cm <- cor_mat
cm[upper.tri(cm)] <- NA                 # keep diagonal + lower triangle
ord <- intersect(ins_order, colnames(cm))

melted <- reshape2::melt(cm, na.rm = TRUE)
melted$Var1 <- factor(melted$Var1, levels = rev(ord))   # rows top->bottom
melted$Var2 <- factor(melted$Var2, levels = ord)        # cols left->right

p <- ggplot(melted, aes(Var2, Var1, fill = value)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", value)), size = 4, color = "black") +
  scale_fill_gradient2(low = "#B2182B", mid = "#FDFBFB", high = "#3B3B9C",
                       midpoint = 0, limits = c(-1, 1),
                       name = "Spearman R",
                       breaks = c(-1, -0.5, 0, 0.5, 1)) +
  coord_fixed() +
  labs(title = "Spearman rank correlation between insecticides",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid  = element_blank(),
    plot.title  = element_text(hjust = 0.5)
  )

ggsave(file.path(spearman_dir, "spearman_heatmap.png"),
       p, width = 8, height = 5.5, dpi = 300)

message("\nSaved to ", spearman_dir, ":")
message("  - spearman_heatmap.png")
message("  - spearman_correlation_matrix.csv")
message("  - spearman_pvalues.csv")
message("  - spearman_uncertainty.csv  (R with 95% CrI per insecticide pair)")