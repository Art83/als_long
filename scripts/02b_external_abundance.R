################################################################################
#  02b. External plasma-abundance re-adjustment (Supplementary)
#
#  The adjusted model in 02 uses baseline SomaScan signal (signal_mean) as the
#  abundance covariate, i.e. the same assay as the outcome, so one can
#  argue the adjustment removed a real effect. Here the covariate is an
#  external measure of plasma concentration:
#  the Human Protein Atlas blood concentration (MS based, and
#  immunoassay based where available).
#
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})


cat("================================================================\n")
cat("  02b: external plasma-abundance re-adjustment\n")
cat("================================================================\n\n")

merged <- read.csv("../output/grounding/progression_with_grounding.csv", stringsAsFactors = FALSE)

hpa_path <- "../raw/proteinatlas.tsv"
hpa <- read.delim(hpa_path, check.names = FALSE, quote = "", stringsAsFactors = FALSE)

col_ms <- grep("Conc. blood MS", colnames(hpa), value = TRUE)[1]
col_im <- grep("Conc. blood IM", colnames(hpa), value = TRUE)[1]
if (is.na(col_ms)) stop("No 'Conc. blood MS' column in proteinatlas.tsv. Columns present: ",
                        paste(grep("blood", colnames(hpa), value = TRUE, ignore.case = TRUE), collapse = " | "))
cat("  HPA columns used:", col_ms, if (!is.na(col_im)) paste("|", col_im) else "", "\n")



max_or_na <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)

hpa_blood <- hpa %>%
  transmute(gene_symbol = Gene,
            conc_ms_pgL = suppressWarnings(as.numeric(.data[[col_ms]])),
            conc_im_pgL = if (!is.na(col_im)) suppressWarnings(as.numeric(.data[[col_im]])) else NA_real_) %>%
  filter(!is.na(gene_symbol), gene_symbol != "") %>%
  group_by(gene_symbol) %>%
  summarise(conc_ms_pgL = max_or_na(conc_ms_pgL), conc_im_pgL = max_or_na(conc_im_pgL), .groups = "drop") %>%
  mutate(log10_ms = log10(conc_ms_pgL), log10_im = log10(conc_im_pgL))

dat <- merged %>%
  filter(!is.na(grounding_score)) %>%
  left_join(hpa_blood, by = "gene_symbol") %>%
  mutate(abs_r = abs(cor_r))

cat("  Aptamers in grounding map:", nrow(dat),
    "| with HPA MS blood concentration:", sum(!is.na(dat$log10_ms)),
    sprintf("(%.1f%%)", 100 * mean(!is.na(dat$log10_ms))),
    "| with HPA immunoassay concentration:", sum(!is.na(dat$log10_im)), "\n\n")

# 1. Does SomaScan signal track external concentration? (sanity check on the covariate)
cat("-- 1. SomaScan baseline signal vs external concentration\n")
ok <- !is.na(dat$log10_ms)
ct <- cor.test(dat$signal_mean[ok], dat$log10_ms[ok], method = "spearman", exact = FALSE)
cat(sprintf("  Spearman rho(signal_mean, log10 HPA MS conc) = %.3f (p = %.2g, n = %d)\n", ct$estimate, ct$p.value, sum(ok)))
cat("  Mean signal_mean and HPA MS conc by grounding tier (covered aptamers):\n")
print(dat %>% filter(ok) %>% group_by(grounding_tier) %>%
        summarise(n = n(), mean_signal = round(mean(signal_mean), 2),
                  mean_log10_ms = round(mean(log10_ms), 2), .groups = "drop") %>% as.data.frame(), row.names = FALSE)
cat("\n")

# 2. Models on the covered subset: naive, SomaScan-adjusted, externally adjusted, both
cat("-- 2. |r| ~ grounding score on the HPA-covered subset\n")
sub <- dat %>% filter(ok)
m_naive_sub <- lm(abs_r ~ grounding_score, data = sub)
m_soma_sub <- lm(abs_r ~ grounding_score + signal_mean + signal_cv, data = sub)
m_ext_sub <- lm(abs_r ~ grounding_score + log10_ms + signal_cv, data = sub)
m_both_sub <- lm(abs_r ~ grounding_score + log10_ms + signal_mean + signal_cv, data = sub)

gs_row <- function(m, label, n) {
  s <- summary(m)$coefficients
  data.frame(model = label, n = n,
             beta_grounding = s["grounding_score", "Estimate"],
             p_grounding = s["grounding_score", "Pr(>|t|)"],
             r2 = summary(m)$r.squared)
}
ext_tbl <- bind_rows(
  gs_row(m_naive_sub, "naive (covered subset)", nrow(sub)),
  gs_row(m_soma_sub,  "adjusted: SomaScan signal + CV", nrow(sub)),
  gs_row(m_ext_sub,   "adjusted: HPA MS conc + CV", nrow(sub)),
  gs_row(m_both_sub,  "adjusted: HPA MS conc + SomaScan signal + CV", nrow(sub))
)
print(ext_tbl %>% mutate(across(where(is.numeric), ~ signif(.x, 3))), row.names = FALSE)
cat("  Abundance covariate p values: SomaScan signal in model 2 =",
    signif(summary(m_soma_sub)$coefficients["signal_mean", "Pr(>|t|)"], 3),
    "| HPA MS conc in model 3 =", signif(summary(m_ext_sub)$coefficients["log10_ms", "Pr(>|t|)"], 3), "\n\n")

# 3. Permutation null for the externally adjusted grounding coefficient
cat("-- 3. Permutation test, externally adjusted model\n")
set.seed(42)
B <- 2000
obs_beta <- coef(m_ext_sub)["grounding_score"]
perm <- replicate(B, {
  td <- sub; td$grounding_score <- sample(td$grounding_score)
  coef(lm(abs_r ~ grounding_score + log10_ms + signal_cv, data = td))["grounding_score"]
})
p_perm <- mean(abs(perm) >= abs(obs_beta))
cat(sprintf("  Observed beta = %.6f | permutation p = %.3f | null 95%% range %.6f to %.6f\n\n",
            obs_beta, p_perm, quantile(perm, 0.025), quantile(perm, 0.975)))

# 4. Secondary: confirmed vs transcript-only, externally adjusted
cat("-- 4. Within transcript-positive aptamers, externally adjusted\n")
rna_sub <- sub %>% filter(mn_rna_detected %in% TRUE)
m_rna_ext <- lm(abs_r ~ mn_confirmed + log10_ms + signal_cv, data = rna_sub)
s_rna <- summary(m_rna_ext)$coefficients
cat(sprintf("  n = %d | mn_confirmed beta = %.5f, p = %.3f\n\n", nrow(rna_sub),
            s_rna["mn_confirmedTRUE", "Estimate"], s_rna["mn_confirmedTRUE", "Pr(>|t|)"]))

# 5. Immunoassay-based concentration where available (smaller set, reported for completeness)
if (sum(!is.na(dat$log10_im)) >= 200) {
  cat("-- 5. Same models with HPA immunoassay concentration (smaller covered set)\n")
  sub_im <- dat %>% filter(!is.na(log10_im))
  m_im <- lm(abs_r ~ grounding_score + log10_im + signal_cv, data = sub_im)
  s_im <- summary(m_im)$coefficients
  cat(sprintf("  n = %d | grounding beta = %.6f, p = %.3f | HPA IM conc p = %.3g\n\n", nrow(sub_im),
              s_im["grounding_score", "Estimate"], s_im["grounding_score", "Pr(>|t|)"], s_im["log10_im", "Pr(>|t|)"]))
  ext_tbl <- bind_rows(ext_tbl, gs_row(m_im, "adjusted: HPA IM conc + CV", nrow(sub_im)))
}

ext_tbl$permutation_p_external <- NA_real_
ext_tbl$permutation_p_external[ext_tbl$model == "adjusted: HPA MS conc + CV"] <- p_perm
write.csv(ext_tbl, "../output/summaries/fb_external_abundance.csv", row.names = FALSE)
write.csv(data.frame(n_covered = sum(ok), rho_signal_vs_hpa = unname(ct$estimate), p = ct$p.value),
          "../output/summaries/fb_external_abundance_coverage.csv", row.names = FALSE)

cat("== SCRIPT 02b COMPLETE ==\n")
