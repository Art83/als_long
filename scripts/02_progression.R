################################################################################
#  02. Progression across timepoints
#
#  Main Q: Does motor neuron molecular grounding predict which plasma proteins
#  track ALS disease progression?
#
#  For each SomaScan protein:
#    1. Per-patient BLUP slope (log2 RFU ~ time, random slope per patient)
#    2. Correlate with ALSFRS-R slope (BLUP from lmer)
#    3. Record effect size, p-value, signal quality metrics
#  Then test: do MN-grounded proteins show stronger associations?
#
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lme4)
  library(lmerTest)
  library(ggplot2)
  library(gridExtra)
})

# helper fo r timepoints
source("00_visit_time.R")

cat("================================================================\n")
cat("  02: GROUNDING vs PROGRESSION\n")
cat("================================================================\n\n")

#----- 0. Load data -------------

grounding <- read.csv("../output/grounding/grounding_map.csv", stringsAsFactors = FALSE)
annotation <- read.csv("../raw/somascan/documents/annotation.csv")
dataportal <- read.csv("../raw/als-metadata/answerals/aals_dataportal_datatable.csv", check.names = FALSE)
alsfrs <- read.csv("../raw/als-metadata/answerals/clinical/ALSFRS_R.csv")

df <- read.csv("../raw/somascan/normalized/csv/SS-2342615_v4.1_anmlSMP__AALS-v1.0-participant-samples.csv",
               check.names = FALSE)


misc_columns <- c("Barcode2d", "SubjectUID", "TimePoint")
seq_cols     <- colnames(df)[!colnames(df) %in% misc_columns]


seq2uni <- annotation %>%
  dplyr::select(SeqId, UniProt.ID) %>%
  dplyr::filter(UniProt.ID != "" & !is.na(UniProt.ID)) %>%
  dplyr::rename(uniprot = UniProt.ID)

valid_seq <- intersect(seq_cols, seq2uni$SeqId)
df <- df[, c(misc_columns, valid_seq)]

cat("SomaScan aptamers (mapped):", length(valid_seq), "\n")


# ------- 1. ALSFRS-R trajectory model ----------------------

cat("\n-- 1. ALSFRS-R trajectory model\n\n")

# All ALSFRS-R assessments (numbered visits and follow-up / phone visits) go into the
# trajectory model.
alsfrs_proc <- alsfrs %>%
  dplyr::select(SubjectUID, Visit_Date, alsfrst) %>%
  dplyr::rename(GUID = SubjectUID) %>%
  filter(!is.na(alsfrst), !is.na(Visit_Date)) %>%
  group_by(GUID) %>%
  filter(n() >= 2) %>%
  arrange(Visit_Date, .by_group = TRUE) %>%
  mutate(time_years = Visit_Date / 365.25) %>%
  ungroup()

als_guids <- dataportal %>%
  filter(SUBJECT_GROUP == "ALS") %>%
  pull(GUID) %>%
  unique()

alsfrs_als <- alsfrs_proc %>% filter(GUID %in% als_guids)

m_alsfrs <- lme4::lmer(alsfrst ~ time_years + (time_years | GUID), data = alsfrs_als)

als_slopes <- lme4::ranef(m_alsfrs)$GUID %>%
  tibble::rownames_to_column("GUID") %>%
  dplyr::rename(ALSFRS_int = `(Intercept)`, ALSFRS_slope = time_years)

cat("Patients with ALSFRS slopes:", nrow(als_slopes), "\n\n")


# ------ 2. Longitudinal SomaScan preparation ------------------

cat("-- 2. Longitudinal SomaScan preparation\n\n")

vt <- soma_visit_time("../raw/als-metadata/answerals/clinical/Plasma_Sample.csv",
                      "../raw/als-metadata/answerals/clinical/ALSFRS_R.csv")
write.csv(vt, "../output/summaries/fb_visit_time_table.csv", row.names = FALSE)

soma_long <- df %>%
  filter(SubjectUID %in% als_guids) %>%
  attach_soma_time(vt)
cat("  Time source for ALS draws:"); print(table(soma_long$time_source)); cat("\n")

guids_long <- soma_long %>%
  group_by(GUID) %>%
  summarise(n_visits = n(), .groups = "drop") %>%
  filter(n_visits >= 2) %>%
  pull(GUID)

guids_both <- intersect(guids_long, als_slopes$GUID)
soma_long  <- soma_long %>% filter(GUID %in% guids_both)

cat("ALS patients with >=2 SomaScan visits + ALSFRS:", length(guids_both), "\n")
cat("Total protein-visit observations:", nrow(soma_long), "\n\n")


# ------ 3. Per-protein progression associations ----------------

cat("-- 3. Computing protein-progression associations\n")

slopes_lookup <- setNames(
  als_slopes$ALSFRS_slope[als_slopes$GUID %in% guids_both],
  als_slopes$GUID[als_slopes$GUID %in% guids_both]
)

# stage 1 protein slopes via BLUP (random slope per patient),
# not per-patient OLS. Earlier simulations showed OLS slope reliability ~0.40 vs BLUP ~0.83, so the
# OLS two-stage injects avoidable noise. Singular random-slope fits are kept (BLUPs just
# shrink to the fixed effect, which is the desired behaviour). Full OLS fallback only if
# lmer errors outright.
test_one_protein <- function(seq_id, data, slopes_lookup) {

  vals <- data[, c("GUID", "time_years", seq_id)]
  colnames(vals)[3] <- "value"
  vals$value <- log2(as.numeric(vals$value))
  vals <- vals[is.finite(vals$value), ]

  # Signal quality at baseline (unchanged)
  base_vals <- vals %>%
    group_by(GUID) %>%
    filter(time_years == min(time_years)) %>%
    summarise(base_value = value[1], .groups = "drop")
  signal_mean <- mean(base_vals$base_value, na.rm = TRUE)
  signal_cv <- sd(base_vals$base_value, na.rm = TRUE) / abs(signal_mean)

  # need enough patients with >= 2 observations
  npat <- vals %>% group_by(GUID) %>% filter(n() >= 2) %>%
    summarise(.groups = "drop") %>% nrow()
  if (npat < 20) return(NULL)

  # Stage 1: BLUP per-patient slope
  slope_method <- "blup"
  bl <- tryCatch(
    suppressMessages(suppressWarnings(
      lme4::lmer(value ~ time_years + (time_years | GUID), data = vals,
                 control = lme4::lmerControl(calc.derivs = FALSE)))),
    error = function(e) NULL)

  if (!is.null(bl)) {
    fe <- lme4::fixef(bl)["time_years"]
    re <- lme4::ranef(bl)$GUID
    pslope <- setNames(as.numeric(fe) + re[["time_years"]], rownames(re))
  } else {
    slope_method <- "ols"
    ps <- vals %>% group_by(GUID) %>% filter(n() >= 2) %>%
      summarise(prot_slope = coef(lm(value ~ time_years))[2], .groups = "drop")
    pslope <- setNames(ps$prot_slope, ps$GUID)
  }

  patient_slopes <- data.frame(GUID = names(pslope),
                               prot_slope = as.numeric(pslope),
                               stringsAsFactors = FALSE)
  patient_slopes$alsfrs_slope <- slopes_lookup[patient_slopes$GUID]
  patient_slopes <- patient_slopes %>% filter(!is.na(alsfrs_slope) & is.finite(prot_slope))
  if (nrow(patient_slopes) < 20) return(NULL)

  # Stage 2: association (unchanged)
  fit <- lm(prot_slope ~ alsfrs_slope, data = patient_slopes)
  s <- summary(fit)

  data.frame(
    SeqId = seq_id,
    n_patients = nrow(patient_slopes),
    slope_method = slope_method,
    beta = coef(fit)["alsfrs_slope"],
    se = s$coefficients["alsfrs_slope", "Std. Error"],
    t_stat = s$coefficients["alsfrs_slope", "t value"],
    p_value = s$coefficients["alsfrs_slope", "Pr(>|t|)"],
    r_squared = s$r.squared,
    cor_r = cor(patient_slopes$prot_slope, patient_slopes$alsfrs_slope, use = "complete.obs"),
    signal_mean = signal_mean,
    signal_cv = signal_cv,
    stringsAsFactors = FALSE
  )
}

# parallelise (BLUP per protein is much heavier than OLS).
# mclapply forks, which does not work on Windows: it silently falls back to serial.
# if somebody re-runs this - use Linux. 
# Serial on windows takes ages to converge.
library(parallel)
ncores <- max(1, detectCores() - 1)
cat("  Fitting BLUP protein slopes on", ncores, "cores for", length(valid_seq), "proteins...\n")

if (.Platform$OS.type == "windows") {
  cl <- makeCluster(ncores)
  clusterEvalQ(cl, { suppressPackageStartupMessages({ library(dplyr); library(lme4) }) })
  clusterExport(cl, c("test_one_protein", "soma_long", "slopes_lookup", "valid_seq"),
                envir = environment())
  results_list <- parLapply(cl, seq_along(valid_seq), function(i) {
    tryCatch(test_one_protein(valid_seq[i], soma_long, slopes_lookup), error = function(e) NULL)
  })
  stopCluster(cl)
} else {
  results_list <- mclapply(seq_along(valid_seq), function(i) {
    tryCatch(test_one_protein(valid_seq[i], soma_long, slopes_lookup), error = function(e) NULL)
  }, mc.cores = ncores)
}

progression <- do.call(rbind, results_list[!sapply(results_list, is.null)])
progression$p_fdr <- p.adjust(progression$p_value, method = "fdr")

cat("\n  Stage-1 slope method used:\n"); print(table(progression$slope_method)); cat("\n")

# Map SeqId -> UniProt (clean, no make.names issues)
progression <- progression %>%
  left_join(seq2uni, by = "SeqId")

cat("\n Proteins tested:", nrow(progression), "\n")
cat("Nominal p < 0.05:", sum(progression$p_value < 0.05), "\n")
cat("FDR < 0.05:", sum(progression$p_fdr < 0.05), "\n\n")


# ------- 4. Merge with grounding map ------------------------------------

cat("-- 4. Merging progression with grounding\n\n")

grounding_dedup <- grounding %>%
  group_by(uniprot) %>%
  summarise(
    gene_symbol = gene_symbol[1],
    grounding_score = grounding_score[1],
    grounding_tier = grounding_tier[1],
    mn_confirmed = mn_confirmed[1],
    mn_rna_detected = mn_rna_detected[1],
    ms_measured = if ("ms_measured" %in% names(grounding)) ms_measured[1] else NA,
    mn_rna_rank = mn_rna_rank[1],
    mn_protein_rank = mn_protein_rank[1],
    .groups = "drop"
  )

merged <- progression %>%
  left_join(grounding_dedup, by = "uniprot")

merged$grounding_tier[is.na(merged$grounding_tier)] <- "unmapped"
merged$grounding_score[is.na(merged$grounding_score)] <- NA

cat("Matched to grounding map:", sum(!is.na(merged$grounding_score)), "/", nrow(merged), "\n\n")


# -------- 5. Primary analysis -----------------------------------------------

cat("-- 5. Primary analysis: grounding vs biomarker quality\n\n")

cat("Mean |r| by grounding tier:\n")
tier_summary <- merged %>%
  filter(grounding_tier != "unmapped") %>%
  group_by(grounding_tier) %>%
  summarise(n_proteins  = n(),
            mean_abs_r  = mean(abs(cor_r), na.rm = TRUE),
            median_abs_r = median(abs(cor_r), na.rm = TRUE),
            pct_nom_sig = mean(p_value < 0.05, na.rm = TRUE) * 100,
            pct_fdr_sig = mean(p_fdr < 0.05, na.rm = TRUE) * 100,
            mean_signal = mean(signal_mean, na.rm = TRUE),
            mean_cv = mean(signal_cv, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(desc(mean_abs_r))
print(as.data.frame(tier_summary)); cat("\n")

cat("Mean |r| by grounding score:\n")
score_summary <- merged %>%
  filter(!is.na(grounding_score)) %>%
  group_by(grounding_score) %>%
  summarise(n = n(), mean_abs_r = mean(abs(cor_r), na.rm = TRUE),
            pct_sig = mean(p_value < 0.05, na.rm = TRUE) * 100,
            mean_signal = mean(signal_mean, na.rm = TRUE), .groups = "drop")
print(as.data.frame(score_summary)); cat("\n")

cat("Naive regression: |r| ~ grounding_score\n")
m_naive <- lm(abs(cor_r) ~ grounding_score, data = merged, subset = !is.na(grounding_score))
print(summary(m_naive)); cat("\n")

cat("Controlled: |r| ~ grounding_score + signal_mean + signal_cv\n")
m_controlled <- lm(abs(cor_r) ~ grounding_score + signal_mean + signal_cv,
                   data = merged, subset = !is.na(grounding_score))
print(summary(m_controlled)); cat("\n")


# ---------- 6. Secondary: RNA+Protein vs RNA-only --------------------------------------

cat("-- 6. Secondary: RNA+Protein vs RNA-only (among RNA-detected)\n\n")

rna_detected <- merged %>% filter(mn_rna_detected == TRUE, grounding_tier != "unmapped")

cat("  MN protein confirmed:", sum(rna_detected$mn_confirmed, na.rm = TRUE), "\n")
cat("  RNA-only (MS gap):", sum(!rna_detected$mn_confirmed, na.rm = TRUE), "\n\n")

rna_detected %>%
  group_by(mn_confirmed) %>%
  summarise(n = n(), mean_abs_r = mean(abs(cor_r), na.rm = TRUE),
            pct_sig = mean(p_value < 0.05, na.rm = TRUE) * 100,
            mean_signal = mean(signal_mean, na.rm = TRUE), .groups = "drop") %>%
  as.data.frame() %>% print()

cat("\nWilcoxon:\n")
print(wilcox.test(abs(cor_r) ~ mn_confirmed, data = rna_detected))

cat("\nControlled:\n")
m_rna <- lm(abs(cor_r) ~ mn_confirmed + signal_mean + signal_cv, data = rna_detected)
print(summary(m_rna))


# --------- 7. Permutation test ----------------------------------------------

cat("\n-- 7. Permutation test\n\n")

obs_beta <- coef(m_controlled)["grounding_score"]
set.seed(42)
B <- 2000
test_data <- merged %>% filter(!is.na(grounding_score))

perm_betas <- replicate(B, {
  td <- test_data
  td$grounding_score <- sample(td$grounding_score)
  coef(lm(abs(cor_r) ~ grounding_score + signal_mean + signal_cv, data = td))["grounding_score"]
})

p_perm <- mean(abs(perm_betas) >= abs(obs_beta))
cat("Observed beta:", round(obs_beta, 6), "\n")
cat("Permutation p:", p_perm, "\n")
cat("Null 95% range:", round(quantile(perm_betas, c(0.025, 0.975)), 6), "\n\n")


# ---------- 8. Plots ----------------------------------------------------------

cat("-- 8. Plots\n\n")

tier_order  <- c("no_mn_expression", "rna_only_ms_gap", "ms_only_rna_gap",
                 "confirmed_partial", "confirmed_full")
tier_colors <- c(no_mn_expression = "grey70", rna_only_ms_gap = "grey50",
                 ms_only_rna_gap = "#B2ABD2",
                 confirmed_partial = "#4393C3", confirmed_full = "#2166AC")

tiff("../output/grounding/progression_by_tier.tiff", width = 6, height = 5, units = "in", res = 300)
merged %>%
  filter(grounding_tier != "unmapped") %>%
  mutate(grounding_tier = factor(grounding_tier, levels = tier_order)) %>%
  ggplot(aes(x = grounding_tier, y = abs(cor_r), fill = grounding_tier)) +
  geom_violin(trim = FALSE, alpha = 0.5) + geom_boxplot(width = 0.15, outlier.shape = NA) +
  scale_fill_manual(values = tier_colors) +
  labs(x = "Motor neuron grounding", y = "|Correlation with ALSFRS-R decline|") +
  theme_bw() + theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1))
dev.off()

tiff("../output/grounding/progression_by_score.tiff", width = 5, height = 5, units = "in", res = 300)
merged %>% filter(!is.na(grounding_score)) %>%
  ggplot(aes(x = factor(grounding_score), y = abs(cor_r))) +
  geom_violin(trim = FALSE, alpha = 0.4, fill = "steelblue") +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  labs(x = "Grounding score", y = "|Correlation with ALSFRS-R decline|") + theme_bw()
dev.off()

tiff("../output/grounding/confirmed_vs_rnaonly.tiff", width = 5, height = 5, units = "in", res = 300)
rna_detected %>%
  mutate(group = ifelse(mn_confirmed, "RNA + Protein", "RNA only (MS gap)")) %>%
  ggplot(aes(x = group, y = abs(cor_r), fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.5) + geom_boxplot(width = 0.15, outlier.shape = NA) +
  scale_fill_manual(values = c("RNA + Protein" = "#2166AC", "RNA only (MS gap)" = "grey60")) +
  labs(x = "", y = "|Correlation with ALSFRS-R decline|") +
  theme_bw() + theme(legend.position = "none")
dev.off()

tiff("../output/grounding/signal_quality_check.tiff", width = 8, height = 4, units = "in", res = 300)
p1 <- rna_detected %>% mutate(g = ifelse(mn_confirmed, "Confirmed", "RNA only")) %>%
  ggplot(aes(x = g, y = signal_mean, fill = g)) + geom_boxplot(alpha = 0.5) +
  labs(y = "Mean log2 RFU", x = "") + theme_bw() + theme(legend.position = "none")
p2 <- rna_detected %>% mutate(g = ifelse(mn_confirmed, "Confirmed", "RNA only")) %>%
  ggplot(aes(x = g, y = signal_cv, fill = g)) + geom_boxplot(alpha = 0.5) +
  labs(y = "CV of baseline signal", x = "") + theme_bw() + theme(legend.position = "none")
grid.arrange(p1, p2, ncol = 2)
dev.off()

tiff("../output/grounding/permutation_null.tiff", width = 5, height = 4, units = "in", res = 300)
ggplot(data.frame(beta = perm_betas), aes(x = beta)) +
  geom_histogram(bins = 50, fill = "grey70", color = "grey40") +
  geom_vline(xintercept = obs_beta, color = "red", linewidth = 1) +
  annotate("text", x = obs_beta, y = Inf, label = paste("observed\np =", round(p_perm, 3)),
           vjust = 2, hjust = -0.1, color = "red", size = 4) +
  labs(x = "Beta (grounding_score) under permutation", y = "Count") + theme_bw()
dev.off()


# Numbers for manuscript

primary <- data.frame(
  term = c("naive_grounding_score", "controlled_grounding_score", "permutation_grounding_score"),
  beta = c(coef(m_naive)["grounding_score"], coef(m_controlled)["grounding_score"], obs_beta),
  p_value = c(summary(m_naive)$coefficients["grounding_score", "Pr(>|t|)"],
              summary(m_controlled)$coefficients["grounding_score", "Pr(>|t|)"],
              p_perm)
)
write.csv(primary, "../output/summaries/fb_primary_result.csv", row.names = FALSE)

sec <- rna_detected %>% group_by(mn_confirmed) %>%
  summarise(n = n(), mean_abs_r = mean(abs(cor_r), na.rm = TRUE),
            pct_sig = mean(p_value < 0.05, na.rm = TRUE) * 100, .groups = "drop")
write.csv(sec, "../output/summaries/fb_secondary_rnaonly.csv", row.names = FALSE)

write.csv(data.frame(n_tested = nrow(progression),
                     n_nom = sum(progression$p_value < 0.05),
                     n_fdr = sum(progression$p_fdr < 0.05)),
          "../output/summaries/fb_progression_counts.csv", row.names = FALSE)

merged %>% arrange(p_value) %>% head(15) %>%
  dplyr::select(gene_symbol, uniprot, cor_r, p_value, p_fdr, grounding_tier, slope_method) %>%
  write.csv("../output/summaries/fb_progression_top15.csv", row.names = FALSE)

write.csv(tier_summary, "../output/summaries/fb_tier_summary.csv", row.names = FALSE)


# --- 9. Save ---------------------------------------

write.csv(merged, "../output/grounding/progression_with_grounding.csv", row.names = FALSE)
write.csv(tier_summary, "../output/grounding/tier_summary.csv", row.names = FALSE)

# Also export grounding_dedup for downstream scripts
saveRDS(grounding_dedup, "../output/grounding/grounding_dedup.rds")

cat("\n== SCRIPT 02 COMPLETE ==\n")
