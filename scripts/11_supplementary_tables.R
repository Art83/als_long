################################################################################
#  Supplementary tables
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lme4)
  library(readxl)
})

source("00_visit_time.R")

OUT <- "../output/supplementary_tables/"

cat("================================================================\n")
cat("  Supplementary tables\n")
cat("================================================================\n\n")

num <- function(x) suppressWarnings(as.numeric(x))
txt <- function(x) ifelse(is.na(x), "", trimws(as.character(x)))
med_iqr <- function(x, d = 1) {
  x <- x[is.finite(x)]
  if (!length(x)) return("")
  q <- quantile(x, c(0.25, 0.5, 0.75))
  sprintf(paste0("%.", d, "f [%.", d, "f, %.", d, "f] (n = %d)"), q[2], q[1], q[3], length(x))
}
n_pct <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return("")
  sprintf("%d/%d (%.0f%%)", sum(x), length(x), 100 * mean(x))
}
tier_labels <- c(confirmed_full = "Expressed, promoter-supported",
                 confirmed_partial = "Expressed, promoter-unconfirmed",
                 rna_only_ms_gap = "Transcript only",
                 ms_only_rna_gap = "Protein only",
                 no_mn_expression = "Not expressed")
tier_lab <- function(x) ifelse(x %in% names(tier_labels), unname(tier_labels[x]), "No annotated gene")

tables <- list()
add_table <- function(id, title, parts) {
  tables[[id]] <<- list(title = title, parts = parts)
  for (nm in names(parts)) {
    f <- paste0(OUT, id, if (length(parts) > 1) paste0("_", gsub("[^A-Za-z0-9]+", "_", nm)) else "", ".csv")
    write.csv(parts[[nm]], f, row.names = FALSE, na = "")
  }
  cat(sprintf("  %s: %s (%s)\n", id, title, paste(sapply(parts, nrow), collapse = " + ")))
}


# =============================================================================
# Table S1. Answer ALS participant characteristics
# =============================================================================
portal <- read.csv("../raw/als-metadata/answerals/aals_dataportal_datatable.csv", check.names = FALSE)
soma_file <- "../raw/somascan/normalized/csv/SS-2342615_v4.1_anmlSMP__AALS-v1.0-participant-samples.csv"
soma_ids <- if (requireNamespace("data.table", quietly = TRUE)) {
  as.data.frame(data.table::fread(soma_file, select = c("SubjectUID", "TimePoint")))
} else read.csv(soma_file, check.names = FALSE)[, c("SubjectUID", "TimePoint")]

C9_REPEAT_CUTOFF <- 30
pp <- portal %>% distinct(GUID, .keep_all = TRUE) %>%
  transmute(GUID, SUBJECT_GROUP,
            age = num(Age_at_First_PBMC_Collection),
            male = ifelse(txt(SEX) == "", NA, SEX == "Male"),
            age_onset = num(AGE_AT_SYMPTOM_ONSET),
            duration = num(Age_at_First_PBMC_Collection) - num(AGE_AT_SYMPTOM_ONSET),
            bulbar = ifelse(txt(Site_of_Onset) == "", NA, txt(Site_of_Onset) == "Bulbar"),
            definite = ifelse(txt(REVISED_EL_ESCORIAL_CRITERIA) == "", NA, txt(REVISED_EL_ESCORIAL_CRITERIA) == "Definite"),
            gene_info = txt(HAS_GENOMICS) == "Yes" | txt(ClinReport_Mutations_Details) != "" | txt(HasVariant) != "",
            c9 = (is.finite(num(C9orf72_repeat_length)) & num(C9orf72_repeat_length) >= C9_REPEAT_CUTOFF) |
              grepl("C9", ClinReport_Mutations_Details, ignore.case = TRUE) | grepl("C9", HasVariant, ignore.case = TRUE)) %>%
  mutate(c9 = ifelse(gene_info, c9, NA))

fx <- read.csv("../raw/als-metadata/answerals/clinical/ALSFRS_R.csv") %>%
  filter(!is.na(alsfrst), !is.na(Visit_Date)) %>% arrange(SubjectUID, Visit_Date)
fx_base <- fx %>% group_by(SubjectUID) %>% summarise(alsfrs_base = alsfrst[1], n_alsfrs = n(), .groups = "drop") %>%
  rename(GUID = SubjectUID)

base_guids <- soma_ids %>% filter(SubjectUID %in% pp$GUID[pp$SUBJECT_GROUP %in% c("ALS", "Healthy Control")]) %>%
  distinct(SubjectUID) %>% pull(SubjectUID)
vt <- soma_visit_time("../raw/als-metadata/answerals/clinical/Plasma_Sample.csv",
                      "../raw/als-metadata/answerals/clinical/ALSFRS_R.csv")
draws <- soma_ids %>% attach_soma_time(vt, verbose = FALSE) %>% group_by(GUID) %>%
  summarise(n_draws = n(), span_years = max(time_years) - min(time_years), .groups = "drop")

fs_env <- new.env(); load("../output/validation/fast_slow_tracking.RData", envir = fs_env)
long_guids <- fs_env$pt$GUID

cohort <- pp %>% filter(GUID %in% base_guids) %>%
  left_join(fx_base, by = "GUID") %>% left_join(draws, by = "GUID") %>%
  left_join(fs_env$pt %>% select(GUID, als_slope), by = "GUID")
grp <- list(`ALS, baseline plasma` = cohort %>% filter(SUBJECT_GROUP == "ALS"),
            `Controls, baseline plasma` = cohort %>% filter(SUBJECT_GROUP == "Healthy Control"),
            `ALS, serial plasma (slope design)` = cohort %>% filter(GUID %in% long_guids))
s1 <- data.frame(Characteristic = c("Participants", "Age at first PBMC collection, years", "Male",
                                    "Age at symptom onset, years", "Symptom duration at first PBMC collection, years",
                                    "Bulbar onset", "Revised El Escorial definite", "C9orf72 repeat expansion",
                                    "Baseline ALSFRS-R total", "ALSFRS-R assessments", "ALSFRS-R slope, points per year",
                                    "Dated plasma draws", "Plasma follow-up, years"))
for (g in names(grp)) {
  d <- grp[[g]]; als <- !grepl("Controls", g)
  s1[[g]] <- c(nrow(d), med_iqr(d$age), n_pct(d$male),
               if (als) med_iqr(d$age_onset) else "", if (als) med_iqr(d$duration, 2) else "",
               if (als) n_pct(d$bulbar) else "", if (als) n_pct(d$definite) else "", if (als) n_pct(d$c9) else "",
               if (als) med_iqr(d$alsfrs_base) else "", if (als) med_iqr(d$n_alsfrs, 0) else "",
               if (grepl("serial", g)) med_iqr(d$als_slope, 2) else "",
               if (als) med_iqr(d$n_draws, 0) else "", if (als) med_iqr(d$span_years, 2) else "")
}
add_table("TableS1", "Answer ALS participant characteristics",
          list(Characteristics = s1))


# =============================================================================
# Table S2. Differential plasma proteins, ALS against controls
# =============================================================================
dap <- read.csv("../output/landscape/daps_annotated.csv", stringsAsFactors = FALSE)
s2 <- dap %>% filter(sig %in% TRUE) %>% arrange(adj.P.Val) %>%
  transmute(SeqId, UniProt = uniprot, Gene = gene_symbol, log2FC = logFC, Mean_log2_RFU = AveExpr, t,
            p = P.Value, FDR = adj.P.Val, Direction = ifelse(logFC > 0, "Higher in ALS", "Lower in ALS"),
            iMN_map_class = tier_lab(grounding_tier), Route_to_plasma = route,
            Top_tissue = top_tissue, Normalized_tissue_entropy = entropy_norm)
add_table("TableS2", "Differential plasma proteins between patients with ALS and controls (FDR < 0.05), adjusted for sex and age",
          list(DAPs = s2))


# =============================================================================
# Table S3. Differential plasma proteins against baseline ALSFRS-R and symptom duration
# =============================================================================
st <- read.csv("../output/summaries/fb_dap_stage_association.csv", stringsAsFactors = FALSE)
s3 <- st %>% mutate(covariate = recode(covariate, alsfrs_base = "Baseline ALSFRS-R", duration_years = "Years since symptom onset")) %>%
  arrange(covariate, p) %>%
  transmute(Covariate = covariate, SeqId, UniProt = uniprot, Gene = gene_symbol, r, beta, p, FDR = p_fdr,
            log2FC_ALS_vs_control = dap_logFC)
add_table("TableS3", "Association of the differential plasma proteins with baseline ALSFRS-R and years since symptom onset within ALS",
          list(Association = s3))


# =============================================================================
# Table S4. Progression analysis, all aptamers
# =============================================================================
prog <- read.csv("../output/grounding/progression_with_grounding.csv", stringsAsFactors = FALSE)
s4 <- prog %>% arrange(p_value) %>%
  transmute(SeqId, UniProt = uniprot, Gene = gene_symbol, Patients = n_patients,
            Beta_protein_slope_per_ALSFRS_point = beta, SE = se, r = cor_r, p = p_value, FDR = p_fdr,
            Progression_marker = p_fdr < 0.05, Mean_log2_RFU = signal_mean, CV = signal_cv,
            iMN_map_class = tier_lab(grounding_tier), Grounding_score = grounding_score)
add_table("TableS4", "Association of within-patient change in each aptamer with the rate of functional decline (229 patients)",
          list(Aptamers = s4))


# =============================================================================
# Table S5. Target ALS participant characteristics
# =============================================================================
ta_dir <- "../raw/target_als"
na_codes <- function(x) { x <- suppressWarnings(as.numeric(x)); x[!is.na(x) & x < 0] <- NA; x }
clin <- read.csv(file.path(ta_dir, "nhs_clinical_data.csv"), check.names = FALSE, stringsAsFactors = FALSE) %>%
  transmute(subject_id = `Subject ID`, visit_number = as.integer(`Visit Number`),
            group = ifelse(`Clinical Group` == "ALS", "ALS", "Control"),
            sex = ifelse(`Sex Reported` %in% c("Male", "Female"), `Sex Reported`, NA),
            age = na_codes(`Age In Years`), alsfrs = na_codes(`Alsfrs R Score Total`))
clin$alsfrs[!is.na(clin$alsfrs) & clin$alsfrs > 48] <- NA
subj <- clin %>% arrange(subject_id, visit_number) %>% group_by(subject_id) %>%
  summarise(group = group[1], sex = sex[!is.na(sex)][1], age_base = age[!is.na(age)][1],
            alsfrs_base = alsfrs[!is.na(alsfrs)][1], n_alsfrs = sum(!is.na(alsfrs)), .groups = "drop")

tinfo <- read_excel(file.path(ta_dir, "ms/lc_ms_proteomics_tmt_plasma.xlsx"), sheet = "1-plasma_sample_info") %>% as.data.frame()
names(tinfo) <- make.names(names(tinfo))
tvis <- tinfo %>% filter(technical.replicate != -8888) %>% distinct(subject.id, visit.number) %>%
  group_by(subject_id = subject.id) %>% summarise(n_tmt_visits = n(), .groups = "drop")
olk_ids <- rownames(read.csv(file.path(ta_dir, "olink_matrix.csv"), check.names = FALSE, row.names = 1))

ok_cov <- function(d) d %>% filter(!is.na(group), !is.na(sex), !is.na(age_base))
tmt_s <- subj %>% inner_join(tvis, by = "subject_id") %>% ok_cov()
olk_s <- subj %>% filter(subject_id %in% olk_ids) %>% ok_cov()
ta_env <- new.env(); load("../output/validation/target_als_replication.RData", envir = ta_env)
long_s <- subj %>% filter(subject_id %in% ta_env$comp$subject_id) %>% left_join(tvis, by = "subject_id")

tg <- list(`TMT, ALS` = tmt_s %>% filter(group == "ALS"), `TMT, controls` = tmt_s %>% filter(group == "Control"),
           `Olink, ALS` = olk_s %>% filter(group == "ALS"), `Olink, controls` = olk_s %>% filter(group == "Control"),
           `TMT, ALS with serial samples` = long_s)
s5 <- data.frame(Characteristic = c("Participants", "Age at first visit, years", "Male", "First recorded ALSFRS-R total",
                                    "ALSFRS-R assessments", "TMT plasma visits"))
for (g in names(tg)) {
  d <- tg[[g]]; als <- !grepl("controls", g)
  s5[[g]] <- c(nrow(d), med_iqr(d$age_base), n_pct(d$sex == "Male"),
               if (als) med_iqr(d$alsfrs_base) else "", if (als) med_iqr(d$n_alsfrs, 0) else "",
               if (grepl("TMT", g)) med_iqr(d$n_tmt_visits, 0) else "")
}
add_table("TableS5", "Target ALS participant characteristics by platform",
          list(Characteristics = s5))


# =============================================================================
# Table S6. Replication of the differential plasma proteins in Target ALS
# =============================================================================
rep <- read.csv("../output/summaries/fb_replication_dap_concordance.csv", stringsAsFactors = FALSE)
s6 <- rep %>% arrange(platform, desc(abs(aals_logFC))) %>%
  transmute(Platform = platform, Gene = gene_symbol, Answer_ALS_log2FC = aals_logFC, Target_ALS_log2FC = est,
            Target_ALS_p = P.Value, Target_ALS_FDR = adj.P.Val, Same_direction = sign(aals_logFC) == sign(est))
add_table("TableS6", "Differential plasma proteins from Answer ALS in Target ALS (TMT mass spectrometry and Olink), adjusted for sex and age",
          list(Replication = s6))


# =============================================================================
# Table S7. Replication of the progression markers in Target ALS
# =============================================================================
tmt_prog <- read.csv("../output/replication/target_als_tmt_plasma_progression.csv", stringsAsFactors = FALSE)
aals_prog <- prog %>% filter(p_fdr < 0.05, !is.na(gene_symbol)) %>% group_by(gene_symbol) %>%
  summarise(aals_r = cor_r[which.min(p_value)], .groups = "drop")
s7a <- aals_prog %>% inner_join(tmt_prog, by = "gene_symbol") %>% arrange(desc(abs(aals_r))) %>%
  transmute(Gene = gene_symbol, UniProt = uniprot, Answer_ALS_r = aals_r, Target_ALS_patients = n,
            Target_ALS_r = cor_r, Target_ALS_p = p_value, Same_direction = sign(aals_r) == sign(cor_r))
s7b <- ta_env$per_marker %>% arrange(p_ta) %>%
  transmute(Gene = gene, UniProt = uniprot, SeqId, Direction_from_Answer_ALS = sign_aals,
            Target_ALS_patients = n, Target_ALS_r = r_ta, Target_ALS_p = p_ta, Concordant = concordant)
s7c <- data.frame(Statistic = c("Patients", "Markers", "Composite r with ALSFRS-R slope", "Bootstrap 95% CI lower",
                                "Bootstrap 95% CI upper", "One-sided permutation p", "Empirical p against random non-marker sets",
                                "Leave-one-marker-out r, minimum", "Leave-one-marker-out r, maximum"),
                  Value = c(nrow(ta_env$comp), nrow(ta_env$mk), ta_env$r_obs, quantile(ta_env$r_boot, 0.025, na.rm = TRUE),
                            quantile(ta_env$r_boot, 0.975, na.rm = TRUE), ta_env$p_perm, ta_env$p_null,
                            min(ta_env$loo$r), max(ta_env$loo$r)))
add_table("TableS7", "Progression markers from Answer ALS in Target ALS (TMT mass spectrometry, visit number as time)",
          list(`Gene-level concordance` = s7a, `Composite markers` = s7b, `Composite summary` = s7c))


# =============================================================================
# Table S8. Motor neuron expression map, aptamer to gene
# =============================================================================
gm <- read.csv("../output/grounding/grounding_map.csv", stringsAsFactors = FALSE)
s8 <- gm %>% transmute(SeqId, UniProt = uniprot, Target = target_name, Target_full_name = target_full,
                       Gene = gene_symbol, Ensembl = ensembl_id,
                       Transcript_detected = mn_rna_detected, Mean_VST_expression = mn_rna_mean,
                       Measured_by_MS = ms_measured, Protein_detected = mn_protein_detected,
                       Protein_detected_fraction_of_lines = mn_protein_pct_det,
                       Accessible_promoter = mn_chromatin_open, Grounding_score = grounding_score,
                       iMN_map_class = tier_lab(grounding_tier))
add_table("TableS8", "Motor neuron expression map of the SomaScan panel (7,175 aptamers mapped to an annotated gene)",
          list(Map = s8))


# =============================================================================
# Table S9. Map classes, tracking and threshold sensitivity
# =============================================================================
ts <- read.csv("../output/summaries/fb_tier_summary.csv", stringsAsFactors = FALSE)
s9a <- ts %>% mutate(grounding_tier = factor(grounding_tier, levels = names(tier_labels))) %>% arrange(grounding_tier) %>%
  transmute(iMN_map_class = tier_lab(as.character(grounding_tier)), Aptamers = n_proteins, Mean_abs_r = mean_abs_r,
            Median_abs_r = median_abs_r, Percent_nominal = pct_nom_sig, Percent_FDR = pct_fdr_sig,
            Mean_log2_RFU = mean_signal, Mean_CV = mean_cv)
pr <- read.csv("../output/summaries/fb_primary_result.csv", stringsAsFactors = FALSE)
ext <- read.csv("../output/summaries/fb_external_abundance.csv", stringsAsFactors = FALSE)
s9b <- bind_rows(pr %>% transmute(Model = recode(term, naive_grounding_score = "Unadjusted (7,175 aptamers)",
                                                 controlled_grounding_score = "Adjusted for signal and CV",
                                                 permutation_grounding_score = "Adjusted, permutation (2,000)"),
                                  Beta = beta, p = p_value),
                 ext %>% transmute(Model = paste0("HPA subset: ", model, " (n = ", n, ")"), Beta = beta_grounding, p = p_grounding))
ps <- read.csv("../output/summaries/fb_tier_threshold_sensitivity.csv", check.names = FALSE)
names(ps) <- c("Protein_detection_threshold", tier_lab(names(ps)[-1]))
rs <- read.csv("../output/tables/rna_threshold_sensitivity.csv", check.names = FALSE)
names(rs) <- c("Transcript_rule", "Genes_detected", tier_lab(names(rs)[-(1:2)]))
add_table("TableS9", "Map classes, progression tracking by class, grounding score models and detection threshold sensitivity",
          list(`Tracking by class` = s9a, `Grounding score models` = s9b,
               `Protein threshold sensitivity` = ps, `Transcript rule sensitivity` = rs))


# =============================================================================
# Table S10. Progression markers with map class, route and tissue
# =============================================================================
pa <- read.csv("../output/landscape/progression_annotated.csv", stringsAsFactors = FALSE)
s10 <- pa %>% filter(p_fdr < 0.05) %>% arrange(desc(abs(cor_r))) %>% mutate(Rank = row_number()) %>%
  transmute(Rank, SeqId, UniProt = uniprot, Gene = gene_symbol, r = cor_r, p = p_value, FDR = p_fdr,
            iMN_map_class = tier_lab(grounding_tier), Route_to_plasma = route, Localisation = localisation,
            Specificity = specificity, Top_tissue = top_tissue, Normalized_tissue_entropy = entropy_norm)
add_table("TableS10", "The 72 progression markers ranked by strength of association, with motor neuron map class, route to plasma and tissue distribution",
          list(Markers = s10))


# =============================================================================
# Table S11. Route to plasma
# =============================================================================
ro <- read.csv("../output/summaries/fb_concordr_route_enrichment.csv", stringsAsFactors = FALSE)
tri <- bind_rows(read.csv("../output/landscape/triage_daps.csv", stringsAsFactors = FALSE),
                 read.csv("../output/landscape/triage_progression_fdr.csv", stringsAsFactors = FALSE),
                 read.csv("../output/landscape/triage_background.csv", stringsAsFactors = FALSE))
s11a <- tri %>% mutate(set = recode(set, dap = "Differential plasma proteins", prog_fdr = "Progression markers", background = "Panel")) %>%
  count(Set = set, Route = route) %>% group_by(Set) %>% mutate(Percent = 100 * n / sum(n)) %>% ungroup() %>% rename(Proteins = n)
s11b <- ro %>% transmute(Set = set, Intracellular = n_intracellular, Secreted_or_membrane = n_regulated,
                         OR_intracellular_vs_panel = or, CI_lower = ci_lo, CI_upper = ci_hi, p)
s11c <- tri %>% filter(set %in% c("dap", "prog_fdr")) %>%
  transmute(Set = recode(set, dap = "Differential plasma proteins", prog_fdr = "Progression markers"),
            Gene = gene_symbol, Route = route, Localisation = localisation, Specificity = specificity,
            Primary_tissue = primary_tissue) %>% arrange(Set, Route, Gene)
add_table("TableS11", "Route to plasma of the differential plasma proteins and progression markers against the SomaScan panel",
          list(`Route by set` = s11a, `Intracellular enrichment` = s11b, `Per protein` = s11c))


# =============================================================================
# Table S12. Tissue distribution
# =============================================================================
cp <- read.csv("../output/summaries/fb_concordr_context_presence.csv", stringsAsFactors = FALSE)
s12a <- cp %>% filter(set %in% c("background", "dap", "prog_fdr")) %>%
  mutate(set = recode(set, dap = "Differential plasma proteins", prog_fdr = "Progression markers", background = "Panel")) %>%
  transmute(Set = set, Tissue = context_tissue, Proteins = n, Percent_detected = pct_detected,
            Percent_RNA_protein_mismatch = pct_mismatch, Percent_not_in_atlas = pct_not_in_panel)
et <- read.csv("../output/summaries/fb_tissue_entropy.csv", stringsAsFactors = FALSE)
s12b <- et %>% transmute(Set = recode(set, dap = "Differential plasma proteins", prog_fdr = "Progression markers (FDR)",
                                      prog_nominal = "Progression markers (nominal)", background = "Panel"),
                         Genes = n_genes, Median_entropy = entropy_median, Median_normalized_entropy = entropy_norm_median,
                         Median_top_tissue_probability = p_top_median, Genes_top_probability_over_0.5 = n_p_top_gt_0.5,
                         Median_tissues = n_tissues_median, Wilcoxon_p_vs_panel = wilcox_vs_bg_norm)
gs <- bind_rows(read.csv("../output/landscape/tissue_profile_gene_summary_fdr.csv") %>% mutate(Set = "Progression markers"),
                read.csv("../output/landscape/tissue_profile_gene_summary_dap.csv") %>% mutate(Set = "Differential plasma proteins"))
s12c <- gs %>% transmute(Set, Gene = gene_symbol, Top_tissue = top_tissue, Top_tissue_probability = top_prob,
                         Tissues = n_tissues, Entropy = entropy, Normalized_entropy = entropy_norm) %>% arrange(Set, Gene)
tt <- read.csv("../output/summaries/fb_top_tissue_by_set.csv", stringsAsFactors = FALSE)
s12d <- tt %>% mutate(set = recode(set, dap = "Differential plasma proteins", prog_fdr = "Progression markers", background = "Panel")) %>%
  arrange(set, desc(n)) %>% transmute(Set = set, Top_tissue = top_tissue, Proteins = n, Percent = pct)
add_table("TableS12", "Tissue distribution of the differential plasma proteins and progression markers across 59 human tissues",
          list(`Detection in named tissues` = s12a, `Entropy by set` = s12b, `Per protein` = s12c, `Top tissue by set` = s12d))


# =============================================================================
# Table S13. NEFL aptamer and serum NfL
# =============================================================================
des <- read.csv("../output/summaries/fb_nfl_designs.csv", stringsAsFactors = FALSE)
s13a <- des %>% transmute(Aptamer = gene, SeqId, Slope_design_patients = n_slope, Slope_design_r = slope_r,
                          Slope_design_p = slope_p, Baseline_design_patients = n_base, Baseline_design_r = base_r,
                          Baseline_design_p = base_p, Baseline_design_p_adjusted_for_severity = base_adj_p,
                          Within_patient_slope_SD = within_slope_sd)
conc <- read.csv("../output/summaries/fb_nfl_concordance.csv", stringsAsFactors = FALSE)
s13b <- data.frame(Comparison = c("All matched draws", "Baseline draw only", "Within-patient change"),
                   Draws_or_patients = c(conc$n_matched, conc$baseline_n, conc$within_n),
                   r = c(conc$pearson_r, conc$baseline_r, conc$within_r),
                   p = c(conc$pearson_p, conc$baseline_p, conc$within_p))
ia <- read.csv("../output/summaries/fb_nfl_immunoassay.csv", stringsAsFactors = FALSE)
rel <- read.csv("../output/summaries/fb_nfl_reliability.csv", stringsAsFactors = FALSE)
s13c <- ia %>% transmute(Set = recode(marker, NEFL_immunoassay_allALS = "All ALS with immunoassay",
                                      NEFL_immunoassay_aptamerSubjects = "Patients in the SomaScan designs"),
                         Slope_design_patients = n_slope, Slope_design_r = slope_r, Slope_design_p = slope_p,
                         Baseline_design_patients = n_base, Baseline_design_r = base_r, Baseline_design_p = base_p,
                         Baseline_design_p_adjusted_for_severity = base_adj_p) %>%
  mutate(Reliability_NfL_slope = c(rel$reliability_nfl_slope, NA), Reliability_ALSFRS_slope = c(rel$reliability_alsfrs_slope, NA),
         Reliability_corrected_r = c(rel$disattenuated_r, NA))
nfl_env <- new.env(); load("../output/nfl_check/nfl_figure_data.RData", envir = nfl_env)
bb <- nfl_env$bb_ia
s13d <- bb %>% group_by(Third = c("Low", "Intermediate", "High")[tertile]) %>%
  summarise(Patients = n(), Median_NfL_pg_mL = 2^median(base_value), Mean_ALSFRS_slope = mean(ALSFRS_slope),
            SD_ALSFRS_slope = sd(ALSFRS_slope), .groups = "drop") %>%
  mutate(Third = factor(Third, levels = c("Low", "Intermediate", "High"))) %>% arrange(Third)
tb <- read.csv("../output/summaries/fb_nfl_plateau_by_tertile.csv")
fu <- read.csv("../output/summaries/fb_nfl_plateau_followup_only.csv")
s13e <- bind_rows(tb %>% transmute(Draws = "All draws", Third = c("Low", "Intermediate", "High")[base_tertile], Patients = n,
                                   Median_NfL_pg_mL = median_nfl_pg, r = r_slope_slope, p),
                  fu %>% transmute(Draws = "Follow-up draws only", Third = c("Low", "Intermediate", "High")[base_tertile],
                                   Patients = n, r = r_slope_slope, p))
add_table("TableS13", "NEFL aptamer and serum NfL immunoassay in the slope and baseline designs",
          list(`Aptamer designs` = s13a, `Aptamer against immunoassay` = s13b, `Immunoassay designs` = s13c,
               `ALSFRS-R slope by baseline NfL` = s13d, `NfL slope by baseline NfL` = s13e))


# =============================================================================
# Table S14. Tracking across rates of decline: robustness and outcomes
# =============================================================================
fc_env <- new.env(); load("../output/validation/fast_progressor_clinical.RData", envir = fc_env)
lab_t <- function(x) recode(x, Fast = "Fastest third", Middle = "Middle third", Slow = "Slowest third")
s14a <- bind_rows(fs_env$tert_ins %>% mutate(Analysis = "In-sample composite, all patients"),
                  fs_env$tert_cv %>% mutate(Analysis = "Cross-fitted composite, all patients"),
                  fs_env$tert_1y %>% mutate(Analysis = "In-sample composite, at least 1 year of follow-up"),
                  fs_env$term_tab %>% mutate(Analysis = sprintf("In-sample composite, excluding samples within %d months of death", window_months))) %>%
  transmute(Analysis, Third = lab_t(tertile), Patients = n, Beta = beta, CI_lower = lo, CI_upper = hi)
s14b <- fc_env$cut_tab %>% transmute(Groups = groups, Patients_per_group = n_per_group, Beta_fastest = beta_fastest,
                                     Beta_fastest_CI_lower = beta_fastest_lo, Beta_fastest_CI_upper = beta_fastest_hi,
                                     Beta_slowest = beta_slowest, r_fastest = r_fastest, r_slowest = r_slowest,
                                     p_fastest_vs_slowest = p_fastest_vs_slowest)
s14c <- fs_env$pm %>% left_join(read.csv("../output/grounding/progression_with_grounding.csv", stringsAsFactors = FALSE) %>%
                                  select(SeqId, gene_symbol), by = "SeqId") %>%
  transmute(SeqId, Gene = gene_symbol, Beta_all = b_all, Beta_fastest = b_fast, Beta_slowest = b_slow,
            Ratio_fastest = ratio_fast, Ratio_slowest = ratio_slow)
s14d <- fc_env$track %>% transmute(Outcome = outcome, Third = ifelse(tertile == "All", "All patients", lab_t(tertile)),
                                   Patients = n, r = r, CI_lower = lo, CI_upper = hi)
onset <- read.csv("../output/validation/fast_progressor_onset_site.csv", stringsAsFactors = FALSE)
s14e <- onset %>% transmute(Onset_site_in_fastest_third = onset, Patients = n, r, CI_lower = lo, CI_upper = hi, p)
add_table("TableS14", "Tracking of decline by the progression marker composite across rates of decline",
          list(`By third and sensitivity analysis` = s14a, `Cut-points` = s14b, `Per marker` = s14c,
               `ALSFRS-R subdomains, vital capacity and BMI` = s14d, `Onset site` = s14e))


# =============================================================================
# Table S15. Clinical characteristics by rate of decline
# =============================================================================
fix_names <- function(d) { names(d) <- recode(names(d), Fast = "Fastest third", Middle = "Middle third", Slow = "Slowest third",
                                              p_three_groups = "p, three groups", p_fast_vs_slow = "p, fastest against slowest",
                                              fdr_three_groups = "FDR, three groups"); d }
var_labels <- c(age_onset = "Age at symptom onset, years", age_screening = "Age at screening, years",
                symptom_years = "Symptom duration at screening, years", dfs_pre = "Pre-baseline rate of decline (delta FS), points per month",
                alsfrs_baseline = "Baseline ALSFRS-R total", bulbar_baseline = "Baseline ALSFRS-R bulbar",
                fine_baseline = "Baseline ALSFRS-R fine motor", gross_baseline = "Baseline ALSFRS-R gross motor",
                resp_baseline = "Baseline ALSFRS-R respiratory", vc_baseline = "Baseline vital capacity, % predicted",
                bmi_baseline = "Baseline BMI", cbs_baseline = "Baseline ALS Cognitive Behavioral Screen",
                nfl_baseline = "Baseline serum NfL, pg/mL", male = "Male", bulbar_onset = "Bulbar onset",
                el_escorial_definite = "Revised El Escorial definite", c9 = "C9orf72 repeat expansion",
                any_als_gene = "Any ALS gene variant", riluzole = "Riluzole", edaravone = "Edaravone",
                gastrostomy_at_baseline = "Gastrostomy at baseline", niv_at_baseline = "Non-invasive ventilation at baseline",
                bulbar_slope = "Bulbar subdomain slope, points per year", fine_motor_slope = "Fine motor subdomain slope, points per year",
                gross_motor_slope = "Gross motor subdomain slope, points per year",
                respiratory_slope = "Respiratory subdomain slope, points per year",
                vc_slope = "Vital capacity slope, % predicted per year", bmi_slope = "BMI slope, per year",
                gastrostomy_ever = "Gastrostomy during follow-up", niv_ever = "Non-invasive ventilation during follow-up",
                trach_ever = "Tracheostomy during follow-up", died = "Died during follow-up")
relab <- function(d) d %>% mutate(variable = ifelse(variable %in% names(var_labels), unname(var_labels[variable]), variable)) %>%
  rename(Characteristic = variable) %>% fix_names()
s15a <- relab(fc_env$tab_base)
s15b <- relab(fc_env$tab_course)
s15c <- fc_env$tab_genes %>% rename(Gene = gene, Carriers = total) %>% fix_names()
s15d <- fc_env$strat %>% transmute(Stratifier = recode(stratifier, `Pre-baseline Delta FS` = "Pre-baseline rate of decline (delta FS)"),
                                   Third = recode(stratum, Middle = "Intermediate"), Patients = n,
                                   Median_ALSFRS_slope = median_als_slope, r = r, CI_lower = lo, CI_upper = hi,
                                   Beta = beta, Beta_CI_lower = beta_lo, Beta_CI_upper = beta_hi)
add_table("TableS15", "Clinical characteristics at enrollment and during follow-up by third of decline (229 patients)",
          list(`At enrollment` = s15a, `During follow-up` = s15b, `ALS genes` = s15c,
               `Tracking within enrollment strata` = s15d))


# =============================================================================
# Workbook
# =============================================================================
if (requireNamespace("openxlsx", quietly = TRUE)) {
  wb <- openxlsx::createWorkbook()
  bold <- openxlsx::createStyle(textDecoration = "bold")
  openxlsx::addWorksheet(wb, "Contents")
  openxlsx::writeData(wb, "Contents", data.frame(Table = names(tables), Title = sapply(tables, `[[`, "title")))
  for (id in names(tables)) {
    openxlsx::addWorksheet(wb, id)
    openxlsx::writeData(wb, id, paste0(sub("TableS", "Table S", id), ". ", tables[[id]]$title), startRow = 1)
    openxlsx::addStyle(wb, id, bold, rows = 1, cols = 1)
    row <- 3
    for (nm in names(tables[[id]]$parts)) {
      d <- tables[[id]]$parts[[nm]]
      if (length(tables[[id]]$parts) > 1) {
        openxlsx::writeData(wb, id, nm, startRow = row); openxlsx::addStyle(wb, id, bold, rows = row, cols = 1); row <- row + 1
      }
      openxlsx::writeData(wb, id, d, startRow = row, headerStyle = bold, keepNA = FALSE)
      row <- row + nrow(d) + 3
    }
  }
  openxlsx::saveWorkbook(wb, paste0(OUT, "Supplementary_Tables_S1-S15.xlsx"), overwrite = TRUE)
  cat("\n  Saved", paste0(OUT, "Supplementary_Tables_S1-S15.xlsx"), "\n")
} else {
  cat("\n  openxlsx not installed: CSVs only (install.packages('openxlsx') for the workbook)\n")
}

cat("\n== SCRIPT 11 COMPLETE ==\n")
