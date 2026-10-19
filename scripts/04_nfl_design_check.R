################################################################################
#  04. NfL DESIGN-SCOPE CHECK + IMMUNOASSAY BENCHMARK
#
#  Purpose: attribute the script-02 NULL for the SomaScan NEFL aptamer correctly.
#
#  Part I (aptamer panel; NfL primary, dynamic top hits as positive controls):
#  Part II (serum NfL by Siemens Atellica immunoassay):
################################################################################


suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lme4)
  library(lmerTest)
  library(ggplot2)
})

source("00_visit_time.R")

cat("================================================================\n")
cat("  04: NfL DESIGN SCOPE CHECK + IMMUNOASSAY BENCHMARK\n")
cat("================================================================\n\n")

# --- 0. Panel (UniProt-keyed) ---------------------------
panel <- c(
  NEFL = "P07196",   # neurofilament light  -- expect flat within-person on the aptamer
  NEFH = "P12036",   # neurofilament heavy
  MMP9 = "P14780",   # significant progression (slope) hit in script 02
  CLEC3B = "P05452",   # top slope-tracking hit  -- positive control, dynamic
  ANTXR2 = "P58335"    # top slope-tracking hit  -- positive control, dynamic
)

# --- 1. Load -------------
annotation <- read.csv("../raw/somascan/documents/annotation.csv")
dataportal <- read.csv("../raw/als-metadata/answerals/aals_dataportal_datatable.csv", check.names = FALSE)
alsfrs <- read.csv("../raw/als-metadata/answerals/clinical/ALSFRS_R.csv")
df <- read.csv("../raw/somascan/normalized/csv/SS-2342615_v4.1_anmlSMP__AALS-v1.0-participant-samples.csv",
                       check.names = FALSE)

seq2uni <- annotation %>%
  dplyr::select(SeqId, UniProt.ID) %>%
  filter(UniProt.ID != "" & !is.na(UniProt.ID)) %>%
  dplyr::rename(uniprot = UniProt.ID)

panel_map <- seq2uni %>%
  filter(uniprot %in% panel) %>%
  mutate(gene = names(panel)[match(uniprot, panel)]) %>%
  filter(SeqId %in% colnames(df))

cat("Panel aptamers found in SomaScan matrix:\n")
print(as.data.frame(panel_map), row.names = FALSE); cat("\n")
if (nrow(panel_map) == 0) stop("No panel aptamers matched. Check UniProt IDs / SeqId column.")

# ---- 2. ALSFRS-R BLUP slopes (identical construction to script 02) ----------
alsfrs_proc <- alsfrs %>%
  dplyr::select(SubjectUID, Visit_Date, alsfrst) %>%
  dplyr::rename(GUID = SubjectUID) %>%
  filter(!is.na(alsfrst), !is.na(Visit_Date)) %>%
  group_by(GUID) %>% filter(n() >= 2) %>%
  arrange(Visit_Date, .by_group = TRUE) %>%
  mutate(time_years = Visit_Date / 365.25) %>%
  ungroup()

als_guids  <- dataportal %>% filter(SUBJECT_GROUP == "ALS") %>% pull(GUID) %>% unique()
alsfrs_als <- alsfrs_proc %>% filter(GUID %in% als_guids)

m_alsfrs   <- lmer(alsfrst ~ time_years + (time_years | GUID), data = alsfrs_als)
als_slopes <- ranef(m_alsfrs)$GUID %>%
  tibble::rownames_to_column("GUID") %>%
  dplyr::rename(ALSFRS_int = `(Intercept)`, ALSFRS_slope = time_years)
# slope on the real (points/year) scale for interpretable plots; adding a constant does
# not change any correlation or p-value below.
als_slopes$ALSFRS_slope <- als_slopes$ALSFRS_slope + fixef(m_alsfrs)["time_years"]

cat("ALS patients with an ALSFRS-R slope:", nrow(als_slopes), "\n\n")

# ---- 3. SomaScan: draw times, longitudinal set (>=2 draws) and baseline draw ---------
vt <- soma_visit_time("../raw/als-metadata/answerals/clinical/Plasma_Sample.csv",
                      "../raw/als-metadata/answerals/clinical/ALSFRS_R.csv" )

soma_long <- df %>% filter(SubjectUID %in% als_guids) %>% attach_soma_time(vt)

guids_long <- soma_long %>% group_by(GUID) %>% summarise(nv = n(), .groups = "drop") %>%
  filter(nv >= 2) %>% pull(GUID)

soma_base <- df %>%
  dplyr::rename(GUID = SubjectUID) %>%
  filter(GUID %in% als_guids) %>%
  group_by(GUID) %>% filter(TimePoint == min(TimePoint)) %>% dplyr::slice(1) %>% ungroup()

cat("Patients in slope design (>=2 visits):", length(guids_long), "\n")
cat("Patients available to baseline design (>=1 visit + ALSFRS slope):",
    length(intersect(soma_base$GUID, als_slopes$GUID)), "\n\n")

# ---- 4. Generic design functions (used for aptamers AND immunoassay) -----------
# long: GUID, time_years, value (already log2). base: GUID, base_value (log2).
run_designs <- function(long, base, label) {

  slope_sd <- resid_sd <- NA_real_; prot_blup <- NULL
  if (nrow(long) > 0) {
    dyn <- tryCatch(suppressMessages(suppressWarnings(
             lmer(value ~ time_years + (time_years | GUID), data = long,
                  control = lmerControl(calc.derivs = FALSE)))), error = function(e) NULL)
    if (!is.null(dyn)) {
      vc <- as.data.frame(lme4::VarCorr(dyn))
      slope_sd <- vc$sdcor[vc$grp == "GUID" & vc$var1 == "time_years" & is.na(vc$var2)]
      resid_sd <- vc$sdcor[vc$grp == "Residual"]
      fe <- lme4::fixef(dyn)["time_years"]; re <- lme4::ranef(dyn)$GUID
      prot_blup <- data.frame(GUID = rownames(re), prot_slope = as.numeric(fe) + re[["time_years"]],
                              stringsAsFactors = FALSE)
    }
  }

  slope_beta <- slope_p <- slope_r <- NA_real_; n_slope <- NA_integer_
  if (!is.null(prot_blup)) {
    ss <- prot_blup %>% inner_join(als_slopes, by = "GUID")
    if (nrow(ss) >= 20) {
      f <- lm(prot_slope ~ ALSFRS_slope, data = ss); s <- summary(f)
      slope_beta <- coef(f)["ALSFRS_slope"]; slope_p <- s$coefficients["ALSFRS_slope", "Pr(>|t|)"]
      slope_r <- cor(ss$prot_slope, ss$ALSFRS_slope); n_slope <- nrow(ss)
    }
  }

  base_beta <- base_p <- base_r <- base_rho <- base_adj_beta <- base_adj_p <- NA_real_; n_base <- NA_integer_
  bb <- base %>% inner_join(als_slopes, by = "GUID")
  if (nrow(bb) >= 20) {
    fb <- lm(ALSFRS_slope ~ base_value, data = bb); sb <- summary(fb)
    base_beta <- coef(fb)["base_value"]; base_p <- sb$coefficients["base_value", "Pr(>|t|)"]
    base_r <- cor(bb$base_value, bb$ALSFRS_slope)
    base_rho <- cor(bb$base_value, bb$ALSFRS_slope, method = "spearman")
    n_base <- nrow(bb)
    fba <- lm(ALSFRS_slope ~ base_value + ALSFRS_int, data = bb); sba <- summary(fba)
    base_adj_beta <- coef(fba)["base_value"]; base_adj_p <- sba$coefficients["base_value", "Pr(>|t|)"]
  }

  data.frame(
    marker = label,
    within_slope_sd = slope_sd, resid_sd = resid_sd, flatness = slope_sd / resid_sd,
    n_slope = n_slope, slope_beta = slope_beta, slope_r = slope_r, slope_p = slope_p,
    n_base = n_base, base_beta = base_beta, base_r = base_r, base_rho = base_rho, base_p = base_p,
    base_adj_beta = base_adj_beta, base_adj_p = base_adj_p,
    stringsAsFactors = FALSE
  )
}

aptamer_long <- function(seq_id) {
  vl <- soma_long %>% filter(GUID %in% guids_long) %>% dplyr::select(GUID, TimePoint, time_years, all_of(seq_id))
  colnames(vl)[4] <- "value"; vl$value <- log2(as.numeric(vl$value)); vl[is.finite(vl$value), ]
}
aptamer_base <- function(seq_id) {
  vb <- soma_base %>% dplyr::select(GUID, all_of(seq_id))
  colnames(vb)[2] <- "base_value"; vb$base_value <- log2(as.numeric(vb$base_value)); vb[is.finite(vb$base_value), ]
}

res <- do.call(rbind, Map(function(s, g) {
  out <- run_designs(aptamer_long(s), aptamer_base(s), paste0(g, "_aptamer"))
  out$gene <- g; out$SeqId <- s; out
}, panel_map$SeqId, panel_map$gene))
rownames(res) <- NULL

cat("-- Part I: design comparison, SomaScan aptamers -------------------\n")
print(res %>%
        dplyr::select(gene, within_slope_sd, flatness, n_slope, slope_r, slope_p,
                      n_base, base_r, base_p, base_adj_p) %>%
        mutate(across(where(is.numeric), ~ signif(.x, 3))), row.names = FALSE)
cat("\n")
write.csv(res, "../output/summaries/fb_nfl_designs.csv", row.names = FALSE)

# ---- 5. NfL aptamer plots ---------------------------------
nfl_seq <- panel_map$SeqId[panel_map$gene == "NEFL"][1]
if (length(nfl_seq) == 1 && !is.na(nfl_seq)) {
  tj <- aptamer_long(nfl_seq)
  tiff("../output/nfl_check/nfl_aptamer_trajectories.tiff", width = 5, height = 4, units = "in", res = 300)
  print(ggplot(tj, aes(time_years, value, group = GUID)) + geom_line(alpha = 0.25) +
          labs(x = "Years from first draw", y = "NfL aptamer (log2 RFU)") + theme_bw())
  dev.off()
  bb <- aptamer_base(nfl_seq) %>% inner_join(als_slopes, by = "GUID")
  tiff("../output/nfl_check/nfl_aptamer_baseline_prognostic.tiff", width = 5, height = 4, units = "in", res = 300)
  print(ggplot(bb, aes(base_value, ALSFRS_slope)) + geom_point(alpha = 0.4) +
          geom_smooth(method = "lm", se = TRUE, color = "red") +
          labs(x = "Baseline NfL aptamer (log2 RFU)", y = "ALSFRS-R slope (points/yr)") + theme_bw())
  dev.off()
}

# ------ 6. Aptamer verdict -------------------------------------
nfl_row <- res %>% filter(gene == "NEFL")
cat("-- NfL aptamer verdict --------------------------------------------\n")
if (nrow(nfl_row)) {
  r <- nfl_row[1, ]
  cat(sprintf("  Within-person slope SD: %.4g   (flatness vs residual: %.3g)\n", r$within_slope_sd, r$flatness))
  cat(sprintf("  Slope design:    n = %s, r = %+.3f, p = %.3g\n", r$n_slope, r$slope_r, r$slope_p))
  cat(sprintf("  Baseline design: n = %s, r = %+.3f, p = %.3g  (adj p = %.3g)\n",
              r$n_base, r$base_r, r$base_p, r$base_adj_p))
}
cat("\n")


# ================================================================
# 7. PART II: SERUM NfL IMMUNOASSAY BENCHMARK
# ================================================================
cat("-- Part II: serum NfL immunoassay benchmark -----------------------\n\n")

nfl_path <- "../raw/NfL/NfL_Seimens_Atellica_Immunoassay_Proteomics_2022_07.csv"
if (!file.exists(nfl_path)) stop("Immunoassay NfL file not found: ", nfl_path)

ia_raw <- read.csv(nfl_path, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "latin1")
ia <- ia_raw %>%
  transmute(GUID = SubjID,
            TimePoint = as.integer(`Sample Timepoint`),
            days_file = as.numeric(`Days from Baseline`),
            days_est = !is.na(`Days from Baseline_Comment`) & `Days from Baseline_Comment` != "",
            raw_value = `Result Value`)

# 'QNS' / 'TNP' are unmeasured -> drop. '>500.0' is above the assay range -> keep at 500
# and flag (a censored high value is still informative for rank-based checks).
ia$censored_high <- grepl("^>", ia$raw_value)
ia$value_pg <- suppressWarnings(as.numeric(sub("^>", "", ia$raw_value)))
n_unmeas <- sum(is.na(ia$value_pg))
cat("  Immunoassay rows:", nrow(ia), "| unmeasured (QNS/TNP) dropped:", n_unmeas,
    "| above-range (>500) kept at 500:", sum(ia$censored_high), "\n")
ia <- ia %>% filter(!is.na(value_pg), !is.na(GUID), !is.na(TimePoint))

# ALS only, then attach visit time: prefer the visit-dated table (same convention as the
# aptamer analyses)
ia <- ia %>% filter(GUID %in% als_guids) %>%
  left_join(vt %>% dplyr::select(GUID, TimePoint, time_years), by = c("GUID", "TimePoint")) %>%
  mutate(time_source = ifelse(!is.na(time_years), "visit_table", ifelse(days_est, "file_estimate", "file_days")),
         time_years = ifelse(is.na(time_years), days_file / 365.25, time_years),
         value = log2(value_pg)) %>%
  filter(is.finite(time_years), is.finite(value))

cat("  ALS immunoassay draws:", nrow(ia), "| subjects:", n_distinct(ia$GUID), "\n")
cat("  Time source:"); print(table(ia$time_source)); cat("\n")
# NOTE: the immunoassay is SERUM; the aptamer is EDTA PLASMA. Same visit, different tube.

ia_guids_long <- ia %>% dplyr::count(GUID) %>% filter(n >= 2) %>% pull(GUID)
ia_long <- ia %>% filter(GUID %in% ia_guids_long) %>% dplyr::select(GUID, TimePoint, time_years, value)
ia_base <- ia %>% group_by(GUID) %>% filter(TimePoint == min(TimePoint)) %>% dplyr::slice(1) %>% ungroup() %>%
  transmute(GUID, base_value = value)

# 7a. Aptamer vs immunoassay concordance at matched draws --------------------
cat("  7a. Aptamer vs immunoassay concordance\n")
conc <- NULL
if (length(nfl_seq) == 1 && !is.na(nfl_seq)) {
  apt_all <- soma_long %>% dplyr::select(GUID, TimePoint, all_of(nfl_seq))
  colnames(apt_all)[3] <- "apt"; apt_all$apt <- log2(as.numeric(apt_all$apt))
  matched <- ia %>% inner_join(apt_all, by = c("GUID", "TimePoint")) %>% filter(is.finite(apt))

  r_p  <- cor.test(matched$apt, matched$value)
  r_s  <- cor.test(matched$apt, matched$value, method = "spearman", exact = FALSE)
  cat(sprintf("    matched draws: %d (%d subjects) | Pearson r = %+.3f (p = %.2g) | Spearman rho = %+.3f (p = %.2g)\n",
              nrow(matched), n_distinct(matched$GUID), r_p$estimate, r_p$p.value, r_s$estimate, r_s$p.value))

  # baseline-only cross-sectional (one draw per person)
  mb <- matched %>% group_by(GUID) %>% filter(TimePoint == min(TimePoint)) %>% dplyr::slice(1) %>% ungroup()
  rb <- cor.test(mb$apt, mb$value)
  cat(sprintf("    baseline draws only: n = %d | Pearson r = %+.3f (p = %.2g)\n", nrow(mb), rb$estimate, rb$p.value))

  # within-person change: last minus first matched draw, per subject with >=2 matched draws
  wp <- matched %>% group_by(GUID) %>% filter(n() >= 2) %>% arrange(TimePoint, .by_group = TRUE) %>%
    summarise(d_apt = apt[n()] - apt[1], d_ia = value[n()] - value[1],
              dt = time_years[n()] - time_years[1], .groups = "drop")
  if (nrow(wp) >= 10) {
    rw <- cor.test(wp$d_apt, wp$d_ia)
    cat(sprintf("    within-person change (last - first): n = %d | r = %+.3f (p = %.2g)\n",
                nrow(wp), rw$estimate, rw$p.value))
  } else rw <- NULL

  cat(sprintf("    spread across subjects at baseline: immunoassay SD(log2) = %.2f | aptamer SD(log2) = %.2f\n",
              sd(mb$value), sd(mb$apt)))

  conc <- data.frame(n_matched = nrow(matched), n_subjects = n_distinct(matched$GUID),
                     pearson_r = r_p$estimate, pearson_p = r_p$p.value,
                     spearman_rho = r_s$estimate, spearman_p = r_s$p.value,
                     baseline_n = nrow(mb), baseline_r = rb$estimate, baseline_p = rb$p.value,
                     within_n = nrow(wp), within_r = if (is.null(rw)) NA else rw$estimate,
                     within_p = if (is.null(rw)) NA else rw$p.value,
                     sd_log2_immunoassay = sd(mb$value), sd_log2_aptamer = sd(mb$apt))
  write.csv(conc, "../output/summaries/fb_nfl_concordance.csv", row.names = FALSE)

  tiff("../output/nfl_check/nfl_aptamer_vs_immunoassay.tiff", width = 5, height = 4.5, units = "in", res = 300)
  print(ggplot(matched, aes(value, apt)) + geom_point(alpha = 0.4) +
          geom_smooth(method = "lm", se = TRUE, color = "red") +
          labs(x = "Serum NfL immunoassay (log2 pg/mL)", y = "Plasma NfL aptamer (log2 RFU)",
               subtitle = sprintf("Pearson r = %.2f, Spearman rho = %.2f, n = %d draws", r_p$estimate, r_s$estimate, nrow(matched))) +
          theme_bw())
  dev.off()
}
cat("\n")

# ---- 7b/7c. The two designs on the immunoassay ------------------------------------
cat("  7b/7c. Baseline-prognostic and slope designs on immunoassay NfL\n")
ia_res_all <- run_designs(ia_long, ia_base, "NEFL_immunoassay_allALS")


apt_base_guids <- intersect(soma_base$GUID, als_slopes$GUID)
ia_res_same <- run_designs(ia_long %>% filter(GUID %in% guids_long),
                           ia_base %>% filter(GUID %in% apt_base_guids),
                           "NEFL_immunoassay_aptamerSubjects")

ia_res <- bind_rows(ia_res_all, ia_res_same)
print(ia_res %>%
        dplyr::select(marker, within_slope_sd, flatness, n_slope, slope_r, slope_p,
                      n_base, base_r, base_rho, base_p, base_adj_p) %>%
        mutate(across(where(is.numeric), ~ signif(.x, 3))), row.names = FALSE)
cat("\n")

# secondary: does baseline immunoassay NfL predict decline non-linearly / in tertiles?
bb_ia <- ia_base %>% inner_join(als_slopes, by = "GUID") %>%
  mutate(tertile = dplyr::ntile(base_value, 3))
if (nrow(bb_ia) >= 30) {
  cat("  ALSFRS-R slope (points/yr) by baseline immunoassay NfL tertile:\n")
  print(bb_ia %>% group_by(tertile) %>%
          summarise(n = n(), median_nfl_pg = round(2^median(base_value), 1),
                    mean_slope = round(mean(ALSFRS_slope), 2), sd_slope = round(sd(ALSFRS_slope), 2),
                    .groups = "drop"), row.names = FALSE)
  kw <- kruskal.test(ALSFRS_slope ~ factor(tertile), data = bb_ia)
  cat("  Kruskal-Wallis across tertiles: p =", format(kw$p.value, digits = 3), "\n\n")
}

write.csv(ia_res, "../output/summaries/fb_nfl_immunoassay.csv", row.names = FALSE)

# plots
tiff("../output/nfl_check/nfl_immunoassay_trajectories.tiff", width = 5, height = 4, units = "in", res = 300)
print(ggplot(ia_long, aes(time_years, value, group = GUID)) + geom_line(alpha = 0.2) +
        labs(x = "Years from baseline", y = "Serum NfL immunoassay (log2 pg/mL)") + theme_bw())
dev.off()
tiff("../output/nfl_check/nfl_immunoassay_baseline_prognostic.tiff", width = 5, height = 4, units = "in", res = 300)
print(ggplot(bb_ia, aes(base_value, ALSFRS_slope)) + geom_point(alpha = 0.4) +
        geom_smooth(method = "lm", se = TRUE, color = "red") +
        labs(x = "Baseline serum NfL immunoassay (log2 pg/mL)", y = "ALSFRS-R slope (points/yr)") + theme_bw())
dev.off()
ss_ia <- NULL
if (nrow(ia_long) > 0) {
  dyn <- suppressMessages(suppressWarnings(lmer(value ~ time_years + (time_years | GUID), data = ia_long,
                                                control = lmerControl(calc.derivs = FALSE))))
  re <- ranef(dyn)$GUID
  ss_ia <- data.frame(GUID = rownames(re), nfl_slope = fixef(dyn)["time_years"] + re[["time_years"]]) %>%
    inner_join(als_slopes, by = "GUID")
  tiff("../output/nfl_check/nfl_immunoassay_slope_on_slope.tiff", width = 5, height = 4, units = "in", res = 300)
  print(ggplot(ss_ia, aes(ALSFRS_slope, nfl_slope)) + geom_point(alpha = 0.4) +
          geom_smooth(method = "lm", se = TRUE, color = "red") +
          labs(x = "ALSFRS-R slope (points/yr)", y = "Serum NfL slope (log2 pg/mL per yr)") + theme_bw())
  dev.off()
}

# ---- 8. Combined verdict ---------------
cat("-- Combined verdict ------------------------------------------------\n")
a <- nfl_row[1, ]; i <- ia_res_all[1, ]
cat(sprintf("  Aptamer:      baseline r = %+.3f (p = %.2g) | slope r = %+.3f (p = %.2g) | within-slope SD = %.3g\n",
            a$base_r, a$base_p, a$slope_r, a$slope_p, a$within_slope_sd))
cat(sprintf("  Immunoassay:  baseline r = %+.3f (p = %.2g) | slope r = %+.3f (p = %.2g) | within-slope SD = %.3g\n",
            i$base_r, i$base_p, i$slope_r, i$slope_p, i$within_slope_sd))
if (!is.null(conc)) cat(sprintf("  Aptamer vs immunoassay: r = %+.3f at %d matched draws\n", conc$pearson_r, conc$n_matched))


cat("\n== SCRIPT 04 COMPLETE ==\n")
