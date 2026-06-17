
## =====================================================================
## Intensity bioassay framework - SINGLE-INSECTICIDE runner
##
## Runs the analysis on ONE insecticide's Excel file (produced by
## split_by_insecticide.R). One dose-response curve per SITE.
## Re-running is cheap: curves whose fit .rds already exists are SKIPPED
## (checkpointing), so you can stop and resume freely.
##
## Input file columns: Site | Insecticide | Dose | Responded | Subjects | Replicate
##   (the file contains one treated insecticide + its matched controls)
## =====================================================================

pacman::p_load(
  rio,
  here,
  tidyverse,
  reshape2,
  rstan,
  loo,
  boot,
  egg,
  gridExtra,
  ggpubr,
  ggsci,
  viridisLite,
  patchwork
)


# ---------------------------------------------------------------------
# 0/ CONFIG -- paths are RELATIVE to the project root ----
# ---------------------------------------------------------------------
# This script expects to be run from the project root (the folder that
# contains /data, /scripts, /stan and /outputs). The `here` package makes
# paths work regardless of the working directory, as long as you open the
# .Rproj file or run from the project root.
#
# It loops over ALL FOUR anonymised insecticide data files in /data and
# writes one output subfolder per insecticide under /outputs.

data_files <- c(
  "permethrin"         = here::here("data", "permethrin.csv"),
  "alpha_cypermethrin" = here::here("data", "alpha_cypermethrin.csv"),
  "ddt"                = here::here("data", "ddt.csv"),
  "pirimiphos_methyl"  = here::here("data", "pirimiphos_methyl.csv")
)

stan_path     <- here::here("stan", "example_model.stan")
outputs_root  <- here::here("outputs")
dir.create(outputs_root, showWarnings = FALSE, recursive = TRUE)

# Toggle: set FALSE to force a refit of every curve (ignore checkpoints)
skip_existing <- TRUE

if (!file.exists(stan_path)) {
  stop("Stan model not found at: ", stan_path,
       "\n  Place the model file 'example_model.stan' in the /stan folder.")
}

# Sites are already anonymised to letter codes (A-P) in the data files;
# the code "K_star" in the data is displayed as "K*".
recode_display <- function(s) ifelse(s == "K_star", "K*", s)

# Fixed colour per site code, consistent across all plots and all runs.
all_codes <- sort(unique(unlist(lapply(data_files, function(f) {
  if (file.exists(f)) unique(read.csv(f)$Site) else character(0)
}))))
site_colours <- setNames(viridisLite::turbo(length(all_codes)),
                         recode_display(all_codes))

# Control label -> treated insecticide is detected case-insensitively:
#   acetone -> PM ; py_control -> permethrin & alpha-cyper ; oc_control -> DDT
control_patterns <- c("acetone", "py_control", "oc_control")

# =====================================================================
# MAIN LOOP -- runs the full analysis for each insecticide file
# =====================================================================
run_one_insecticide <- function(insecticide_file) {

ins_tag    <- tools::file_path_sans_ext(basename(insecticide_file))
output_dir <- file.path(outputs_root, ins_tag)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "summary"), showWarnings = FALSE, recursive = TRUE)
message("\n=========================================================")
message("Running insecticide: ", ins_tag)
message("=========================================================")

# ---------------------------------------------------------------------
# 1/ LOAD & PREP DATA ----
# ---------------------------------------------------------------------
df_raw <- import(insecticide_file)

df_raw <- df_raw %>% mutate(Insecticide = as.character(Insecticide),
                            Site        = as.character(Site))

is_control <- tolower(df_raw$Insecticide) %in% control_patterns

present_controls <- unique(df_raw$Insecticide[is_control])
if (length(present_controls) == 0) {
  warning("No control rows (Acetone / PY_control / OC_Control) found in this file. ",
          "Background mortality will rely entirely on the model's A parameter. ",
          "Insecticide labels present: ", paste(unique(df_raw$Insecticide), collapse = ", "))
} else {
  message("Control(s) present in this file: ", paste(present_controls, collapse = ", "))
}

df_controls <- df_raw %>% filter(is_control)
df_treated  <- df_raw %>% filter(!is_control)

## Observed control mortality, per Site x control type (reported, not fitted)
control_summary <- df_controls %>%
  mutate(Responded = as.integer(Responded), Subjects = as.integer(Subjects)) %>%
  group_by(Site, Insecticide) %>%
  summarise(n_replicates    = n(),
            total_tested    = sum(Subjects),
            total_responded = sum(Responded),
            control_mort_pct = 100 * sum(Responded) / sum(Subjects),
            .groups = "drop")
write.csv(control_summary, file.path(output_dir, "summary/control_mortality.csv"),
          row.names = FALSE)
message("Observed control mortality written to summary/control_mortality.csv")
print(control_summary)

# --- Automated exclusion rule for LC50 estimation -------------------
# A Site x insecticide curve is FLAGGED as excluded from LC50 *presentation*
# if EITHER:
#   (1) it was tested at fewer than `min_doses` distinct concentrations -
#       too few points to identify the 5-parameter logistic (its LC50,
#       slope and asymptotes are not estimable); OR
#   (2) mortality at the LOWEST tested concentration is >= `max_lowdose_mort`%
#       - the assay did not bracket the LC50 (response saturated at/above
#       the lowest dose), so there is no estimable rising portion of the curve.
#
# IMPORTANT: excluded curves are STILL FITTED and all their per-curve outputs
# are retained (written to an `excluded/` subfolder) for audit. They are only
# omitted from the SUMMARY plots and the LC50 summary table. These two criteria
# reproduce exactly the "-" (no LC50) rows in the LC50 tables.
min_doses        <- 4    # require >= 4 distinct doses to attempt a fit
max_lowdose_mort <- 98   # exclude if lowest-dose mortality >= this (%)

# Per-curve diagnostics used by the rule
curve_diag <- df_treated %>%
  mutate(Dose = as.numeric(Dose),
         Responded = as.integer(Responded),
         Subjects  = as.integer(Subjects)) %>%
  group_by(Site, Insecticide, Dose) %>%
  summarise(resp = sum(Responded), subj = sum(Subjects), .groups = "drop_last") %>%
  summarise(
    n_doses        = dplyr::n_distinct(Dose),
    lowest_dose    = min(Dose),
    lowdose_mort   = 100 * resp[which.min(Dose)] / subj[which.min(Dose)],
    .groups = "drop"
  ) %>%
  mutate(
    fail_few_doses = n_doses < min_doses,
    fail_saturated = lowdose_mort >= max_lowdose_mort,
    excluded       = fail_few_doses | fail_saturated,
    reason = dplyr::case_when(
      fail_few_doses & fail_saturated ~ paste0("only ", n_doses,
                                               " doses; lowest-dose mortality ",
                                               round(lowdose_mort), "%"),
      fail_few_doses ~ paste0("only ", n_doses, " distinct doses (<", min_doses, ")"),
      fail_saturated ~ paste0("lowest-dose mortality ", round(lowdose_mort),
                              "% (>=", max_lowdose_mort, "%)"),
      TRUE ~ NA_character_
    )
  )

# Write the exclusion decisions (audit trail for the manuscript)
write.csv(curve_diag, file.path(output_dir, "summary/exclusion_decisions.csv"),
          row.names = FALSE)

excluded_curves <- curve_diag %>% filter(excluded)
if (nrow(excluded_curves) > 0) {
  message(nrow(excluded_curves), " curve(s) flagged excluded from summaries ",
          "(still fitted, outputs kept in excluded/):")
  for (k in seq_len(nrow(excluded_curves))) {
    message("  - ", excluded_curves$Site[k], " [", excluded_curves$Insecticide[k],
            "]: ", excluded_curves$reason[k])
  }
} else {
  message("LC50 exclusion rule: no curves flagged.")
}

# Keys of excluded curves (used to route outputs + filter summaries)
excluded_keys <- curve_diag %>% filter(excluded) %>%
  transmute(key = paste(Site, Insecticide)) %>% pull(key)

# One dose-response curve is defined by Site x Insecticide (treated only).
# We build a combined "Strain" key so the framework code (which loops on
# Strain and uses it for filenames/plots) works unchanged - each curve is
# genuinely one Strain. Site and Insecticide are kept as separate columns
# so the between-lab comparison (block 6) can compare sites WITHIN a compound.
# ALL curves are kept here (excluded ones are fitted too, for audit).
df <- df_treated %>%
  mutate(
    Site        = as.character(Site),
    Insecticide = as.character(Insecticide),
    Dose      = as.numeric(Dose),
    Responded = as.integer(Responded),
    Subjects  = as.integer(Subjects),
    # Curve key: e.g. "O | DDT" (site code | insecticide)
    Strain    = paste(Site, Insecticide, sep = " | "),
    # Mortality_perc is required by the variability / residual sections below
    Mortality_perc = 100 * Responded / Subjects
  )

# Basic sanity checks
stopifnot(all(c("Site","Strain","Insecticide","Dose","Responded","Subjects") %in% names(df)))
if (any(df$Responded > df$Subjects)) {
  bad <- df %>% filter(Responded > Subjects)
  print(bad)
  stop("Some rows still have Responded > Subjects - see printed rows above.")
}

strains_sub <- unique(df$Strain) # one entry per Site x Insecticide = one DR curve
# Filesystem-safe versions (must match safe_name built inside IB_model)
safe_strains <- gsub("[^A-Za-z0-9._-]+", "_", strains_sub)
maxdose     <- max(df$Dose)

# Per-curve excluded flag, aligned to strains_sub order.
# (Strain key is "Site | Insecticide"; excluded_keys use "Site Insecticide".)
strain_key   <- sub(" \\| ", " ", strains_sub)
is_excluded  <- strain_key %in% excluded_keys
names(is_excluded) <- strains_sub

# Excluded curves are still fitted; their outputs go in this subfolder.
excluded_dir <- file.path(output_dir, "excluded")
dir.create(excluded_dir, showWarnings = FALSE, recursive = TRUE)

message("Curves to model (Site | Insecticide; * = excluded from summaries):\n  ",
        paste0(strains_sub, ifelse(is_excluded, " *", "")) %>% paste(collapse = "\n  "))

# ---------------------------------------------------------------------
# 2/ LOAD STAN MODEL + PARALLEL CONFIG ----
# ---------------------------------------------------------------------
# Parallel plan (Windows / PSOCK):
#   - Run several curves at once, each fit using a few chains.
#   - Hard cap of 10 cores total so the PC stays usable.
#   5 curves x 2 chains = 10 cores.
n_workers     <- 5   # curves fitted simultaneously
chains_per    <- 2   # chains per fit  (n_workers * chains_per must be <= 10)
stopifnot(n_workers * chains_per <= 10)

rstan::rstan_options(auto_write = TRUE)

# Compile ONCE in the main session, then save the compiled model so each
# worker can load it instead of recompiling (recompiling per worker is slow).
model <- rstan::stan_model(stan_path)
compiled_model_path <- file.path(output_dir, "compiled_model.rds")
saveRDS(model, compiled_model_path)

# ---------------------------------------------------------------------
# 3/ MODEL FUNCTION (one curve per Site) ----
# ---------------------------------------------------------------------
IB_model <- function(df_single,
                     output_dir, maxdose, skip_existing,
                     compiled_model_path, chains_per) {
  
  # Load the precompiled Stan model once per worker (cached in worker env so
  # repeated curves on the same worker don't reload). Workers are fresh R
  # sessions on Windows, so we cannot rely on a global `model` object.
  if (!exists(".cached_stan_model", envir = .GlobalEnv)) {
    assign(".cached_stan_model", readRDS(compiled_model_path), envir = .GlobalEnv)
  }
  model <- get(".cached_stan_model", envir = .GlobalEnv)
  
  site_name <- unique(df_single$Strain)
  # Filesystem-safe version for filenames (no spaces, pipes, slashes etc.)
  safe_name <- gsub("[^A-Za-z0-9._-]+", "_", site_name)
  
  # --- Checkpointing: skip if this curve's fit already exists ---------
  fit_file <- file.path(output_dir, paste0("fit_", safe_name, ".rds"))
  if (skip_existing && file.exists(fit_file)) {
    message("== SKIP (already fitted): ", site_name, " ==")
    return(invisible(NULL))
  }
  message("\n==== Fitting curve: ", site_name, " ====")
  
  concentration_sim <- seq(0, sqrt(maxdose), length.out = 500)
  N_1 <- length(concentration_sim)
  
  data_stan <- list(
    N             = nrow(df_single),
    mortality     = as.integer(df_single$Responded),
    tested        = as.integer(df_single$Subjects),
    concentration = sqrt(as.vector(df_single$Dose)),
    N_1           = as.integer(N_1),
    concentration_sim = as.vector(concentration_sim)
  )
  
  ## Fit (5000 iters; re-run at 10000 if not converged)
  fit <- sampling(model, data = data_stan, iter = 5000,
                  chains = chains_per, cores = chains_per)
  
  Rhat_df1 <- summary(fit)$summary[, 10]
  modelcheck <- data.frame(
    not_converged        = sum(Rhat_df1 > 1.01, na.rm = TRUE),
    is_NA                = ifelse(length(names(table(is.na(Rhat_df1)))) == 1, "no", "yes"),
    divergent_iterations = rstan::get_num_divergent(fit),
    n_runs               = 5000
  )
  
  if (modelcheck$not_converged > 0) {
    message("Not converged at 5000 iters - re-running at 10000 for ", site_name)
    fit <- sampling(model, data = data_stan, iter = 10000,
                    chains = chains_per, cores = chains_per)
    Rhat_df1 <- summary(fit)$summary[, 10]
    modelcheck <- data.frame(
      not_converged        = sum(Rhat_df1 > 1.01, na.rm = TRUE),
      is_NA                = ifelse(length(names(table(is.na(Rhat_df1)))) == 1, "no", "yes"),
      divergent_iterations = rstan::get_num_divergent(fit),
      n_runs               = 10000
    )
  }
  
  ## LOO
  logLikelihood_lab <- extract_log_lik(fit, "LogLikelihood")
  LOO_lab <- loo(logLikelihood_lab)
  
  ## Simulated mortality curve
  mortality <- rstan::extract(fit, "mean_mortality_sim")[[1]]
  
  mean_mort <- apply(mortality, 2, mean) %>% as.data.frame() %>%
    mutate(concentration = concentration_sim^2) %>% rename(Test_mort_perc = 1) %>% mutate(dat = "sim")
  mean_lower <- apply(mortality, 2, function(x) quantile(x, 0.025)) %>% as.data.frame() %>%
    mutate(concentration = concentration_sim^2) %>% rename(Test_mort_perc = 1) %>% mutate(dat = "sim")
  mean_upper <- apply(mortality, 2, function(x) quantile(x, 0.975)) %>% as.data.frame() %>%
    mutate(concentration = concentration_sim^2) %>% rename(Test_mort_perc = 1) %>% mutate(dat = "sim")
  
  mean_mort <- mean_mort %>% mutate(lower = mean_lower$Test_mort_perc,
                                    upper = mean_upper$Test_mort_perc)
  
  fit_mort <- tibble(
    Dose           = mean_mort$concentration,
    Mortality_perc = mean_mort$Test_mort_perc * 100,
    lower          = mean_mort$lower * 100,
    upper          = mean_mort$upper * 100,
    dat            = mean_mort$dat,
    Strain         = unique(df_single$Strain),
    Insecticide    = unique(df_single$Insecticide)
  )
  
  df1_s <- df_single %>% bind_rows(fit_mort)
  
  plot <- ggplot(df1_s, aes(x = Dose, y = Mortality_perc)) +
    geom_point(data = filter(df1_s, is.na(dat))) +
    geom_ribbon(data = filter(df1_s, !is.na(dat)), aes(ymin = lower, ymax = upper), alpha = 0.2) +
    geom_line(data = filter(df1_s, !is.na(dat)), aes(y = Mortality_perc)) +
    theme_classic() +
    theme(panel.grid.major = element_line(colour = "grey93"),
          panel.grid.minor = element_line(colour = "grey97")) +
    ylab("Mortality (%)") +
    xlab("Insecticide concentration (%, ug/bottle, mg/m2)") +
    ggtitle("Bioassay mortality (5-param logistic model, sqrt scale)",
            paste0(df_single$Insecticide, ", Site: ", site_name)) +
    scale_x_log10()
  
  ## LC values ----
  lcx <- function(y, B, C, E) exp(((log(((-1/(y-1))^(1/E)) - 1)) / B) + C)
  
  Bfit <- rstan::extract(fit)[["B"]]
  Cfit <- rstan::extract(fit)[["C"]]
  Efit <- rstan::extract(fit)[["E"]]
  
  sapply(c(10, 50, 90, 99), function(a) {
    message("Running LC summary ", a, " for ", site_name)
    
    temp <- do.call(rbind, sapply(1:length(Bfit), function(i) {
      data.frame("LC" = lcx(a/100, Bfit[i], Cfit[i], Efit[i]))
    }, simplify = FALSE))
    temp$LC_value <- a
    
    lc <- temp %>% mutate(LC = LC^2,
                          Strain = unique(df_single$Strain),
                          Insecticide = unique(df_single$Insecticide))
    
    lc_summ <- lc %>% group_by(Strain, Insecticide) %>%
      summarise(LC_mean   = mean(LC),
                LC_median = median(LC),
                LC_lower  = quantile(LC, 0.025, na.rm = TRUE),
                LC_upper  = quantile(LC, 0.975, na.rm = TRUE), .groups = "drop")
    
    write.csv(lc,      file = file.path(output_dir, paste0("lc", a, "_", safe_name, ".csv")))
    write.csv(lc_summ, file = file.path(output_dir, paste0("summ_lc", a, "_", safe_name, ".csv")))
    
    if (a == 50) {
      plot_LC <- plot + geom_segment(data = lc_summ,
                                     aes(x = LC_mean, y = -Inf, xend = LC_mean, yend = a),
                                     linetype = "dashed")
      LC50_density_plot <- ggplot(lc, aes(x = LC)) +
        stat_density(position = "identity", aes(alpha = 0.9)) +
        ggtitle(paste0("Distribution of LC50 values: ", site_name)) +
        theme_bw() + guides(alpha = "none")
      
      ggsave(plot_LC, file = file.path(output_dir, paste0("DR-LC50_", safe_name, ".png")),
             width = 4.85, height = 4.39)
      ggsave(LC50_density_plot, file = file.path(output_dir, paste0("LC50-density_", safe_name, ".png")),
             width = 5, height = 4.39)
    }
    NULL
  }, simplify = FALSE)
  
  ## Concentrations at which mosquitoes die ----
  modelfun <- function(x, B, C, E) 1 - (1 / ((1 + exp(B * (log(x) - C)))^E))
  
  lower <- log10(1e-6); upper <- log10(maxdose)
  x <- 10^seq(lower, upper, length.out = 1000); x <- sqrt(x)
  
  getCDF <- do.call(rbind, sapply(1:length(x), function(a) {
    data.frame(sim_y = mean(sapply(1:length(Bfit), function(i) modelfun(x[a], Bfit[i], Cfit[i], Efit[i]))),
               sim_x = x[a])
  }, simplify = FALSE)) %>% mutate(sim_x2 = sim_x^2)
  
  mortx <- function(x, y) {
    f <- approxfun(y, x); u <- runif(20000); data.frame("Mort_C" = f(u))
  }
  
  getPDF <- data.frame(Mortx = mortx(getCDF$sim_x, getCDF$sim_y),
                       Strain = unique(df_single$Strain),
                       Insecticide = unique(df_single$Insecticide))
  
  mortx_plot <- ggplot(getPDF, aes(x = Mort_C)) +
    geom_density(alpha = .5, position = "identity", fill = "black") +
    labs(x = expression("Mortality at" ~ italic(x) ~ "concentration"), y = "Density") +
    theme_classic() +
    ggtitle(paste0("Concentrations at which mosquitoes die (", site_name, ")")) +
    scale_x_log10()
  
  ## Background mortality ----
  BM <- data.frame(mean_BM   = mean(rstan::extract(fit)[["A"]]),
                   median_BM = median(rstan::extract(fit)[["A"]]),
                   lower_BM  = quantile(rstan::extract(fit)[["A"]], 0.025),
                   upper_BM  = quantile(rstan::extract(fit)[["A"]], 0.975),
                   Strain    = unique(df_single$Strain),
                   Insecticide = unique(df_single$Insecticide))
  
  ## Mortality variability (intra-curve) ----
  sim_y <- rstan::extract(fit, "mean_mortality")[[1]] %>% as.data.frame()
  act_x <- (data_stan$concentration)^2
  act_y <- as.integer(df_single$Mortality_perc)
  diff_it <- do.call(rbind, sapply(1:nrow(sim_y), function(i) abs((sim_y[i, ] * 100) - act_y),
                                   simplify = FALSE))
  varmort_it <- apply(diff_it, 1, sum) / length(act_x)
  
  var_mort <- data.frame(Strain = unique(df_single$Strain),
                         Insecticide = unique(df_single$Insecticide),
                         mean_vmort   = mean(varmort_it),
                         median_vmort = median(varmort_it),
                         lower_vmort  = quantile(varmort_it, 0.025),
                         upper_vmort  = quantile(varmort_it, 0.975))
  
  ## Residuals & actual-vs-predicted ----
  mean_mortality_bm <- rstan::extract(fit, "mean_mortality_bm")[[1]]
  
  df1_residuals <- apply(mean_mortality_bm, 2, median) %>%
    melt(value.name = "sim_y") %>%
    bind_cols(Concentration = data_stan$concentration^2,
              act_y = df_single$Mortality_perc / 100) %>%
    mutate(Residuals = act_y - sim_y)
  
  resplot <- ggplot(df1_residuals, aes(x = Concentration, y = Residuals)) +
    geom_point() + geom_hline(aes(yintercept = 0), linetype = "dashed") +
    ggtitle(paste0("Residuals for ", site_name)) + scale_x_log10()
  
  mort_sim <- df_single %>%
    bind_cols(Mortality_sim_perc = apply(mean_mortality_bm, 2, median) * 100)
  
  mod_lm <- mort_sim %>% group_by(Insecticide, Strain) %>%
    do(mod = lm(Mortality_perc ~ Mortality_sim_perc, data = .))
  df_coeff <- mod_lm %>% do(data.frame(
    Strain = .$Strain, Insecticide = .$Insecticide,
    var = names(coef(.$mod)), coef(summary(.$mod)),
    r2 = summary(.$mod)$r.square,
    RMSE = sqrt(mean(.$mod$residuals^2))))
  
  ap_all <- ggplot(mort_sim, aes(x = Mortality_sim_perc, y = Mortality_perc)) +
    geom_point(size = 4) + geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    ylab("Actual mortality (%)") + xlab("Predicted mortality (%)") + theme_bw() +
    ggtitle(paste0("Predictive accuracy (", site_name, ")")) +
    geom_text(data = df_coeff,
              aes(x = 15, y = 75,
                  label = paste0('R2=', round(r2, 2), "\nRMSE=", round(RMSE, 0), "%"),
                  fontface = 3), size = 5) +
    xlim(c(0, 100)) + ylim(c(0, 100))
  
  ## Save objects ----
  saveRDS(fit,     file = fit_file)
  write.csv(modelcheck, file = file.path(output_dir, paste0("modelcheck_", safe_name, ".csv")))
  saveRDS(LOO_lab, file = file.path(output_dir, paste0("LOO_", safe_name, ".rds")))
  write.csv(df1_s, file = file.path(output_dir, paste0("df1_s_", safe_name, ".csv")))
  ggsave(plot, file = file.path(output_dir, paste0("DR_", safe_name, ".png")), width = 4.85, height = 4.39)
  write.csv(getPDF, file = file.path(output_dir, paste0("PDF_", safe_name, ".csv")))
  ggsave(mortx_plot, file = file.path(output_dir, paste0("mortxplot_", safe_name, ".png")), width = 4.85, height = 4.39)
  write.csv(BM, file = file.path(output_dir, paste0("BM_", safe_name, ".csv")))
  write.csv(var_mort, file = file.path(output_dir, paste0("varmort_", safe_name, ".csv")))
  write.csv(df1_residuals, file = file.path(output_dir, paste0("residuals_", safe_name, ".csv")))
  ggsave(resplot, file = file.path(output_dir, paste0("residuals_", safe_name, ".png")), width = 9.5, height = 4.39)
  write.csv(mort_sim, file = file.path(output_dir, paste0("mortsim_", safe_name, ".csv")))
  write.csv(df_coeff, file = file.path(output_dir, paste0("AP-coeff_", safe_name, ".csv")))
  ggsave(ap_all, file = file.path(output_dir, paste0("AP_", safe_name, ".png")), width = 4.85, height = 4.39)
}

# ---------------------------------------------------------------------
# 4/ RUN LOOP OVER SITES ----
# ---------------------------------------------------------------------
# ---------------------------------------------------------------------
# Helper: read per-curve output, skipping (with warning) any curve whose
# file is missing - e.g. not yet fitted in a checkpointed/partial run.
# Reads each curve from its own directory (excluded curves live in excluded/).
#   only_included = TRUE  -> return only curves kept for summaries (default)
#   only_included = FALSE -> return ALL curves (for full-audit read-backs)
# ---------------------------------------------------------------------
read_curve_outputs <- function(prefix, add_strain = FALSE, only_included = TRUE) {
  excl_flag <- unname(is_excluded)            # in strains_sub order
  idx <- seq_along(strains_sub)
  if (only_included) idx <- idx[!excl_flag]
  do.call(rbind, lapply(idx, function(i) {
    dir_i <- if (excl_flag[i]) excluded_dir else output_dir
    f <- file.path(dir_i, paste0(prefix, safe_strains[i], ".csv"))
    if (!file.exists(f)) {
      warning("Missing output for '", strains_sub[i], "' (", prefix,
              ") - curve not yet fitted? Skipping in summary.")
      return(NULL)
    }
    out <- read.csv(f)
    if (add_strain) out$strain <- strains_sub[i]
    out
  }))
}

ldf <- sapply(seq_along(strains_sub),
              function(x) df %>% subset(Strain == strains_sub[x]),
              simplify = FALSE)

# --- Safe wrapper: one failing curve must not kill the whole run -----
safe_IB <- function(df_single, output_dir, maxdose, skip_existing,
                    compiled_model_path, chains_per) {
  tryCatch(
    IB_model(df_single, output_dir, maxdose, skip_existing,
             compiled_model_path, chains_per),
    error = function(e) {
      sn <- tryCatch(unique(df_single$Strain), error = function(...) "<unknown>")
      message("!! ERROR fitting '", sn, "': ", conditionMessage(e),
              " - continuing with remaining curves.")
      data.frame(failed_curve = sn, error = conditionMessage(e))
    }
  )
}

# --- Parallel run over curves (Windows PSOCK cluster) ---------------
# n_workers curves at once, each fit using chains_per cores.
# Total cores = n_workers * chains_per (capped at 10 above).
library(parallel)
cl <- makeCluster(n_workers)
on.exit(try(stopCluster(cl), silent = TRUE), add = TRUE)

# Load packages on every worker
clusterEvalQ(cl, {
  suppressMessages({
    library(tidyverse); library(reshape2); library(rstan)
    library(loo); library(boot); library(ggplot2)
  })
  rstan::rstan_options(auto_write = TRUE)
})

# Per-curve target directory: excluded curves' outputs go to excluded/.
curve_dirs <- ifelse(unname(is_excluded), excluded_dir, output_dir)

# Export the function(s) and constants the workers need
clusterExport(cl, c("IB_model", "safe_IB", "output_dir", "excluded_dir",
                    "maxdose", "skip_existing", "compiled_model_path",
                    "chains_per", "ldf", "curve_dirs"),
              envir = environment())

routput_lab <- parLapply(cl, seq_along(ldf), function(i)
  safe_IB(ldf[[i]], curve_dirs[i], maxdose, skip_existing,
          compiled_model_path, chains_per))

stopCluster(cl)

## Report any curves that errored
fails <- Filter(function(x) is.data.frame(x) && "failed_curve" %in% names(x),
                routput_lab)
if (length(fails) > 0) {
  fail_df <- do.call(rbind, fails)
  write.csv(fail_df, file.path(output_dir, "summary/failed_curves.csv"),
            row.names = FALSE)
  message("\n", nrow(fail_df), " curve(s) failed - see summary/failed_curves.csv")
}

## Convergence check (ALL sites, including excluded - full audit)
conv_ll <- read_curve_outputs("modelcheck_", add_strain = TRUE, only_included = FALSE)
write.csv(conv_ll, file.path(output_dir, "summary/modelcheck_all.csv"))

# ---------------------------------------------------------------------
# 5/ SUMMARY PLOTS ACROSS SITES ----
# ---------------------------------------------------------------------
## DR curves, all sites overlaid
sdf1_s <- read_curve_outputs("df1_s_")
write.csv(sdf1_s, file.path(output_dir, "summary/df_all_sim.csv"))

sdf1_s <- sdf1_s %>% rename(Concentration = Dose)

# Derive Site + letter code from the Strain key for labelling
sdf1_s <- sdf1_s %>%
  mutate(Site = sub(" \\| .*$", "", Strain),
         Code = recode_display(Site))

p1 <- ggplot(sdf1_s, aes(x = Concentration, y = Mortality_perc)) +
  geom_point(data = filter(sdf1_s, is.na(dat)), aes(colour = Code)) +
  geom_line(data = filter(sdf1_s, !is.na(dat)), aes(colour = Code)) +
  geom_ribbon(data = filter(sdf1_s, !is.na(dat)),
              aes(ymin = lower, ymax = upper, fill = Code), alpha = .2) +
  ylab("Mortality (%)") + xlab("Concentration") +
  labs(colour = "Site", fill = "Site") +
  ggtitle("Bioassay mortality", "5PL model (sqrt) - by site") +
  theme_classic() +
  theme(panel.grid.major = element_line(colour = "grey93"),
        panel.grid.minor = element_line(colour = "grey97"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  scale_x_sqrt() +
  scale_fill_manual(values = site_colours) +
  scale_colour_manual(values = site_colours)
ggsave(file.path(output_dir, "summary/DR_all_sites.png"), p1, width = 7, height = 5.52, scale = 0.8)

## Mortality variability across sites
var_mort_all <- read_curve_outputs("varmort_")

# Letter codes on the y-axis, sorted by variability; no colour/legend needed
var_mort_all <- var_mort_all %>%
  mutate(Site = sub(" \\| .*$", "", Strain),
         Code = recode_display(Site)) %>%
  mutate(Code = factor(Code, levels = Code[order(median_vmort)]))

p3 <- ggplot(var_mort_all, aes(y = Code, x = median_vmort, colour = Code)) +
  geom_pointrange(aes(xmin = lower_vmort, xmax = upper_vmort),
                  size = 1, shape = 15) +
  theme_minimal() +
  xlab("Median variability estimate\n(% variability in mortality from best-fit line)") +
  ylab("Site") +
  scale_colour_manual(values = site_colours) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
ggsave(file.path(output_dir, "summary/mortvar_all.png"), p3, width = 5, height = 4)

# ---------------------------------------------------------------------
# 6/ LC50 OVERLAP: between-site vs within-site variability ----
#    This is the part that answers the reviewer question directly.
#    We compare each pair of sites' full posterior LC50 distributions and
#    report (a) the ratio of LC50 medians and (b) the posterior probability
#    that the two LC50s differ, alongside each site's own posterior width
#    (the intra-lab uncertainty already captured by the model).
# ---------------------------------------------------------------------
lc50_post <- read_curve_outputs("lc50_")
lc50_post <- lc50_post %>% dplyr::select(LC, Strain, Insecticide)

## (a) Per-curve posterior summary (intra-lab uncertainty)
##     Recover Site from the combined Strain key ("Site | Insecticide").
lc50_post <- lc50_post %>%
  mutate(Site = sub(" \\| .*$", "", Strain))

lc50_within <- lc50_post %>%
  group_by(Insecticide, Site, Strain) %>%
  summarise(median = median(LC),
            lower  = quantile(LC, 0.025, na.rm = TRUE),
            upper  = quantile(LC, 0.975, na.rm = TRUE),
            CI_width_log10 = log10(upper) - log10(lower),  # width on log scale
            .groups = "drop") %>%
  mutate(Code = recode_display(Site))
write.csv(lc50_within, file.path(output_dir, "summary/lc50_within_site.csv"), row.names = FALSE)

## (b) Pairwise between-site comparison, WITHIN each insecticide.
##     Sites are only compared against other sites tested with the SAME
##     compound - comparing LC50s across different insecticides is meaningless.
n_draw <- 4000  # draws sampled from each posterior for the pairwise comparison

insecticides <- unique(lc50_post$Insecticide)

pairwise <- do.call(rbind, lapply(insecticides, function(ins) {
  sites_ins <- lc50_post %>% filter(Insecticide == ins) %>% pull(Site) %>% unique()
  if (length(sites_ins) < 2) return(NULL)  # need >=2 sites to compare
  
  prs <- combn(sites_ins, 2, simplify = FALSE)
  do.call(rbind, lapply(prs, function(pr) {
    a <- lc50_post %>% filter(Insecticide == ins, Site == pr[1]) %>% pull(LC)
    b <- lc50_post %>% filter(Insecticide == ins, Site == pr[2]) %>% pull(LC)
    a <- sample(a, min(n_draw, length(a)))
    b <- sample(b, min(n_draw, length(b)))
    ratio <- a / b
    data.frame(
      Insecticide     = ins,
      Site_A          = pr[1],
      Site_B          = pr[2],
      ratio_median    = median(ratio),
      ratio_lower     = quantile(ratio, 0.025, na.rm = TRUE),
      ratio_upper     = quantile(ratio, 0.975, na.rm = TRUE),
      # posterior prob that LC50_A > LC50_B (closer to 0.5 = more overlap)
      prob_A_gt_B     = mean(a > b),
      # 2-sided "probability of difference": how far prob is from 0.5
      prob_difference = 2 * abs(mean(a > b) - 0.5)
    )
  }))
}))
write.csv(pairwise, file.path(output_dir, "summary/lc50_between_site_pairwise.csv"), row.names = FALSE)

## Plot: LC50 posterior per site, faceted by insecticide
## (visual between- vs within-site comparison, comparable only within a panel)
## Letter codes on the y-axis, sorted by median LC50 (lowest at bottom);
## sites are on the axis so no colour/legend is needed.
lc50_within_plot <- lc50_within %>%
  arrange(Insecticide, median) %>%
  mutate(Code = factor(Code, levels = unique(Code[order(median)])))

p_lc50 <- ggplot(lc50_within_plot, aes(y = Code, x = median, colour = Code)) +
  geom_pointrange(aes(xmin = lower, xmax = upper), size = 1, shape = 15) +
  facet_wrap(~ Insecticide, scales = "free") +
  scale_x_log10() +
  xlab("LC50 (log scale, 95% CrI = intra-lab uncertainty)") +
  ylab("Site") +
  scale_colour_manual(values = site_colours) +
  theme_minimal() +
  theme(legend.position = "none") +
  ggtitle("LC50 by site, within each insecticide",
          "Overlapping intervals = between-site difference within intra-lab uncertainty")
ggsave(file.path(output_dir, "summary/lc50_between_vs_within.png"), p_lc50,
       width = 8, height = 5)

## Plot: full posterior DISTRIBUTION of LC50 for every valid site,
## one panel per insecticide, overlaid density curves coloured by site.
## Legend ordered by median LC50 (ascending) to match the tables.
dist_plots <- list()
for (ins in insecticides) {
  d <- lc50_post %>% filter(Insecticide == ins)
  if (nrow(d) == 0) next
  
  ord <- d %>% group_by(Site) %>% summarise(m = median(LC), .groups = "drop") %>%
    arrange(m)
  d <- d %>%
    mutate(Code = recode_display(Site),
           Code = factor(Code, levels = recode_display(ord$Site)))
  
  p_dist <- ggplot(d, aes(x = LC, colour = Code, fill = Code)) +
    geom_density(alpha = 0.15, linewidth = 0.6) +
    scale_x_log10() +
    labs(x = "LC50 (log scale)", y = "Posterior density",
         colour = "Site", fill = "Site",
         title = paste0("LC50 posterior distributions - ", ins),
         subtitle = "One density per valid site; ordered by median LC50") +
    theme_classic() +
    theme(panel.grid.major = element_line(colour = "grey93")) +
    scale_colour_manual(values = site_colours) +
    scale_fill_manual(values = site_colours)
  
  dist_plots[[ins]] <- p_dist
  
  fn <- file.path(output_dir,
                  paste0("summary/lc50_distributions_",
                         gsub("[^A-Za-z0-9._-]+", "_", ins), ".png"))
  ggsave(fn, p_dist, width = 7, height = 5)
  message("LC50 distribution plot written: ", basename(fn))
}

# ---------------------------------------------------------------------
# 7/ LC50 SUMMARY TABLE (standalone image only) ----
#    Format matches the manuscript LC50 tables:
#    Lab/Site | LC50 | Lower CI | Upper CI | Resistance/fold change
#    (Significance-letter superscripts and the "<100% mortality with DD"
#     column in the manuscript tables are external/post-hoc and are not
#     reproduced here.)
# ---------------------------------------------------------------------
# Format numbers like the manuscript: small values in E-notation, larger
# values as plain decimals, all to 3 significant figures.
fmt_lc <- function(x) {
  vapply(x, function(v) {
    if (is.na(v)) return(NA_character_)
    if (v != 0 && abs(v) < 0.01) formatC(v, format = "e", digits = 2)
    else formatC(signif(v, 3), format = "g")
  }, character(1))
}

lc50_table_df <- lc50_within %>%
  arrange(Insecticide, median) %>%
  group_by(Insecticide) %>%
  mutate(fold_change = median / min(median)) %>%   # vs most-susceptible site
  ungroup() %>%
  transmute(
    `Lab/Site`                 = Code,
    Insecticide                = Insecticide,
    `LC50`                     = fmt_lc(median),
    `Lower CI`                 = fmt_lc(lower),
    `Upper CI`                 = fmt_lc(upper),
    `Resistance/fold change`   = formatC(fold_change, format = "f", digits = 2)
  )
write.csv(lc50_table_df, file.path(output_dir, "summary/lc50_table.csv"),
          row.names = FALSE)

# --- AUDIT: computed LC50s for the EXCLUDED curves -------------------
# These were still fitted (outputs in excluded/) but are NOT shown in the
# figures or main table. Their LC50 values are retained here for audit only,
# with the exclusion reason attached. Treat these LC50s as UNRELIABLE -
# they come from curves the rule deemed unfit for LC50 estimation.
if (any(is_excluded)) {
  lc50_excl_post <- read_curve_outputs("lc50_", only_included = FALSE) %>%
    dplyr::select(LC, Strain, Insecticide) %>%
    mutate(Site = sub(" \\| .*$", "", Strain)) %>%
    filter(paste(Site, Insecticide) %in% excluded_keys)
  
  if (nrow(lc50_excl_post) > 0) {
    lc50_excl_tbl <- lc50_excl_post %>%
      group_by(Insecticide, Site) %>%
      summarise(LC50_median = median(LC),
                LC50_lower  = quantile(LC, 0.025, na.rm = TRUE),
                LC50_upper  = quantile(LC, 0.975, na.rm = TRUE),
                .groups = "drop") %>%
      mutate(Code = recode_display(Site)) %>%
      left_join(curve_diag %>% select(Site, Insecticide, reason),
                by = c("Site", "Insecticide"))
    write.csv(lc50_excl_tbl,
              file.path(excluded_dir, "lc50_excluded_audit.csv"),
              row.names = FALSE)
    message("Audit: excluded-curve LC50s written to excluded/lc50_excluded_audit.csv")
  }
}

# Render as a graphical table (gridExtra::tableGrob)
tt <- gridExtra::ttheme_minimal(
  base_size = 11,
  core    = list(fg_params = list(hjust = 0, x = 0.05)),
  colhead = list(fg_params = list(fontface = "bold"))
)
lc50_tablegrob <- gridExtra::tableGrob(lc50_table_df, rows = NULL, theme = tt)

# Standalone table image (height scales with number of rows)
n_rows <- nrow(lc50_table_df)
ggsave(file.path(output_dir, "summary/lc50_table.png"),
       lc50_tablegrob, width = 8, height = max(2, 0.32 * n_rows + 1))
message("LC50 summary table written: summary/lc50_table.png (+ .csv)")

# ---------------------------------------------------------------------
# 8/ COMBINED SUMMARY IMAGE (2x2 plots only) ----
#    Panels: (A) DR curves  (B) mortality variability
#            (C) LC50 intervals  (D) LC50 posterior densities
#    The LC50 table is now a SEPARATE figure (summary/lc50_table.png).
#    Individual PNGs above are still written; this is an extra combined file.
# ---------------------------------------------------------------------
# Use the first insecticide's density panel (per-file runs have just one).
p_dist_for_grid <- if (length(dist_plots) >= 1) dist_plots[[1]] else patchwork::plot_spacer()

combined <- (p1 + p3) / (p_lc50 + p_dist_for_grid) +
  patchwork::plot_annotation(
    title = paste0("Summary - ", paste(insecticides, collapse = ", ")),
    tag_levels = "A"
  )
ggsave(file.path(output_dir, "summary/combined_summary.png"),
       combined, width = 16, height = 12, dpi = 200)
message("Combined 2x2 summary image written: summary/combined_summary.png")

message("\nDone. Key outputs for the lab-comparison question:")
message("  - summary/lc50_within_site.csv          (intra-lab LC50 uncertainty per site)")
message("  - summary/lc50_table.csv / .png          (formatted LC50 + CI summary table)")
message("  - summary/lc50_between_site_pairwise.csv (between-site LC50 ratios + prob. of difference)")
message("  - summary/lc50_between_vs_within.png     (visual: do the credible intervals overlap?)")
message("  - summary/lc50_distributions_<ins>.png   (overlaid LC50 posterior densities per insecticide)")
message("  - summary/mortvar_all.csv via mortvar_all.png (intra-curve mortality variability)")
message("  - summary/combined_summary.png           (2x2 summary plots in one image)")
message("  - summary/exclusion_decisions.csv        (audit: why each curve was kept/excluded)")
message("  - excluded/                              (full fitted outputs for excluded curves)")
message("  - excluded/lc50_excluded_audit.csv       (audit: computed LC50s for excluded curves)")

  message("Finished: ", ins_tag)
  invisible(output_dir)
}  # end run_one_insecticide()

# ---------------------------------------------------------------------
# DRIVER: run all four insecticides
# ---------------------------------------------------------------------
for (f in data_files) {
  if (file.exists(f)) {
    run_one_insecticide(f)
  } else {
    warning("Data file not found, skipping: ", f)
  }
}

message("\nAll insecticides processed. Outputs in: ", outputs_root)
message("You can now run scripts/2_spearman_correlation.R")
