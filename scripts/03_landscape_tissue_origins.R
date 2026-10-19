################################################################################
#  03. ALS plasma biomarker set + tissue origin (external tool)
#
#  A) Cross-sectional DAPs (ALS vs Control, limma on baseline SomaScan)
#  B) Full progression marker list
#  C) concordR 0.2.0 gates: route, specificity, concordance, audit call
#  C2) Protein presence vote in named nervous vs muscle contexts
#  C3) Tissue-origin profile (entropy, P(top)) from the 59-tissue presence vote
#  D) Grounding map overlay
#  E) Druggability (DGIdb)
################################################################################


suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(limma)
  library(ggplot2)
  library(concordR)
})


cat("================================================================\n")
cat("  01: ALS biomarjer landscape\n")
cat("================================================================\n\n")

# ── 0. Load data ────────────────────────────────────────────────────────────

df <- read.csv("../raw/somascan/normalized/csv/SS-2342615_v4.1_anmlSMP__AALS-v1.0-participant-samples.csv",
               check.names = FALSE)
dataportal <- read.csv("../raw/als-metadata/answerals/aals_dataportal_datatable.csv", check.names = FALSE)
annotation <- read.csv("../raw/somascan/documents/annotation.csv")
grounding <- read.csv("../output/grounding/grounding_map.csv", stringsAsFactors = FALSE)
progression <- read.csv("../output/grounding/progression_with_grounding.csv", stringsAsFactors = FALSE)

grounding_dedup_path <- "../output/grounding/grounding_dedup.rds"
if (file.exists(grounding_dedup_path)) {
  grounding_dedup <- readRDS(grounding_dedup_path)
} else {
    grounding_dedup <- grounding %>%
    group_by(uniprot) %>%
    summarise(gene_symbol = gene_symbol[1],
              grounding_score = grounding_score[1],
              grounding_tier = grounding_tier[1],
              mn_confirmed = mn_confirmed[1],
              ms_measured = if ("ms_measured" %in% names(grounding)) ms_measured[1] else NA,
              .groups = "drop")
}

misc_columns <- c("Barcode2d", "SubjectUID", "TimePoint")
seq_cols     <- colnames(df)[!colnames(df) %in% misc_columns]
uniprot_map  <- setNames(annotation$UniProt.ID, annotation$SeqId)
colnames(df)[!colnames(df) %in% misc_columns] <- uniprot_map[seq_cols]
df <- df[, !is.na(colnames(df)) & colnames(df) != ""]


# ============================
# A. CROSS-SECTIONAL DAPs
# ============================

cat("-- A. Cross-sectional DAPs (ALS vs Control)\n\n")

age_candidates <- grep("AGE", colnames(dataportal), value = TRUE, ignore.case = TRUE)
cat("  Age-like columns in dataportal:", paste(age_candidates, collapse = ", "), "\n")
dp_cols <- unique(c("GUID", "SUBJECT_GROUP", "SEX", age_candidates))

df_base <- df %>%
  dplyr::rename(GUID = SubjectUID) %>%
  group_by(GUID) %>% filter(TimePoint == min(TimePoint)) %>% dplyr::slice(1) %>% ungroup() %>%
  inner_join(dataportal[, dp_cols], by = "GUID") %>%
  filter(SUBJECT_GROUP %in% c("ALS", "Healthy Control"))

cat("  ALS:", sum(df_base$SUBJECT_GROUP == "ALS"),
    " | Control:", sum(df_base$SUBJECT_GROUP == "Healthy Control"), "\n")


meta_cols <- unique(c("Barcode2d", "GUID", "TimePoint", "SUBJECT_GROUP", "SEX", age_candidates))
prot_cols <- setdiff(colnames(df_base), meta_cols)
cat("  Protein columns:", length(prot_cols),
    "| metadata columns excluded:", length(intersect(meta_cols, colnames(df_base))), "\n")

expr_mat <- df_base[, prot_cols] %>%
  mutate(across(everything(), ~ log2(as.numeric(.)))) %>%
  as.matrix() %>% t()
colnames(expr_mat) <- df_base$GUID

n_nonfinite <- sum(!is.finite(expr_mat))
if (n_nonfinite > 0) { expr_mat[!is.finite(expr_mat)] <- NA; cat("  Non-finite log2 values set to NA:", n_nonfinite, "\n") }

group <- factor(df_base$SUBJECT_GROUP, levels = c("Healthy Control", "ALS"))
sexf <- factor(df_base$SEX)

run_limma_mat <- function(mat, design) {
  design <- as.matrix(design)
  colnames(design)[2] <- "groupALS"          # 2nd column is always the ALS contrast
  fit <- limma::eBayes(limma::lmFit(mat, design))
  tt <- limma::topTable(fit, coef = "groupALS", number = Inf, sort.by = "none")
  tt$uniprot <- rownames(tt)
  tt$sig <- !is.na(tt$adj.P.Val) & tt$adj.P.Val < 0.05
  tt$dir <- ifelse(tt$logFC > 0, "up_in_ALS", "down_in_ALS")
  tt
}

# sex-only model kept for the before/after comparison
daps_sexonly <- run_limma_mat(expr_mat, model.matrix(~ group + sexf))

# pick an age column non-missing for both groups
pick_age <- NA_character_
for (ac in age_candidates) {
  ok <- tapply(!is.na(df_base[[ac]]), df_base$SUBJECT_GROUP, mean)
  if (length(ok) == 2 && all(ok > 0.8)) { pick_age <- ac; break }
}
cat("  Age column used for adjustment:", ifelse(is.na(pick_age), "NONE (age skipped)", pick_age), "\n")

if (!is.na(pick_age)) {
  age_num <- as.numeric(df_base[[pick_age]])
  keep <- !is.na(age_num)
  cat("  Median age ALS vs HC:",
      round(median(age_num[keep & group == "ALS"]), 1), "vs",
      round(median(age_num[keep & group == "Healthy Control"]), 1),
      "| dropped for missing age:", sum(!keep), "\n")
  daps <- run_limma_mat(expr_mat[, keep],
                        model.matrix(~ group[keep] + sexf[keep] + age_num[keep]))
} else {
  daps <- daps_sexonly
}

cat("  Total tested:", nrow(daps), " | DAPs:", sum(daps$sig, na.rm = TRUE), "\n")
cat("  Up:", sum(daps$sig & daps$dir == "up_in_ALS", na.rm = TRUE),
    " | Down:", sum(daps$sig & daps$dir == "down_in_ALS", na.rm = TRUE), "\n\n")

# Disease-stage check on the DAPs. Controls have no ALSFRS-R, so case/control cannot
# be conditioned on stage directly. Instead, within ALS at baseline, we ask if a DAP
# tracks baseline ALSFRS-R total or time since symptom onset.
cat("-- Disease-stage association of DAPs within ALS\n\n")

fx_base <- read.csv("../raw/als-metadata/answerals/clinical/ALSFRS_R.csv") %>%
  filter(!is.na(alsfrst), !is.na(Visit_Date)) %>%
  group_by(SubjectUID) %>% arrange(Visit_Date, .by_group = TRUE) %>%
  summarise(alsfrs_base = alsfrst[1], alsfrs_base_day = Visit_Date[1], .groups = "drop") %>%
  dplyr::rename(GUID = SubjectUID)

stage <- df_base %>%
  filter(SUBJECT_GROUP == "ALS") %>%
  dplyr::select(GUID, all_of(intersect(c("AGE_AT_SYMPTOM_ONSET", "Age_at_First_PBMC_Collection"), colnames(df_base)))) %>%
  left_join(fx_base, by = "GUID") %>%
  mutate(duration_years = if (all(c("AGE_AT_SYMPTOM_ONSET", "Age_at_First_PBMC_Collection") %in% colnames(.)))
           as.numeric(Age_at_First_PBMC_Collection) - as.numeric(AGE_AT_SYMPTOM_ONSET) else NA_real_)

cat("  ALS with baseline ALSFRS-R:", sum(!is.na(stage$alsfrs_base)),
    "| with duration since onset:", sum(!is.na(stage$duration_years)), "\n")

dap_ids <- daps$uniprot[daps$sig %in% TRUE]
stage_mat <- expr_mat[dap_ids, stage$GUID, drop = FALSE]

stage_test <- function(y_name) {
  y <- stage[[y_name]]; ok <- !is.na(y)
  if (sum(ok) < 20) return(NULL)
  res <- t(apply(stage_mat[, ok, drop = FALSE], 1, function(x) {
    f <- summary(lm(x ~ y[ok]))$coefficients
    c(beta = f[2, 1], p = f[2, 4], r = cor(x, y[ok], use = "complete.obs"))
  }))
  out <- data.frame(uniprot = rownames(res), covariate = y_name, res, row.names = NULL)
  out$p_fdr <- p.adjust(out$p, method = "fdr")
  out
}
stage_res <- bind_rows(stage_test("alsfrs_base"), stage_test("duration_years"))

if (nrow(stage_res)) {
  stage_res <- stage_res %>%
    left_join(grounding_dedup %>% dplyr::select(uniprot, gene_symbol), by = "uniprot") %>%
    left_join(daps %>% dplyr::select(uniprot, dap_logFC = logFC), by = "uniprot")
  cat("  DAPs associated with stage within ALS (FDR < 0.05 across the DAP set):\n")
  print(stage_res %>% group_by(covariate) %>%
          summarise(n_dap = n(), n_nominal = sum(p < 0.05), n_fdr = sum(p_fdr < 0.05),
                    median_abs_r = round(median(abs(r)), 3), .groups = "drop"))
  cat("  Strongest stage associations:\n")
  print(stage_res %>% arrange(p) %>% head(10) %>%
          dplyr::select(gene_symbol, covariate, r, p, p_fdr, dap_logFC) %>%
          mutate(across(where(is.numeric), ~ signif(.x, 3))), row.names = FALSE)
  write.csv(stage_res, "../output/summaries/fb_dap_stage_association.csv", row.names = FALSE)
}
cat("\n")


# ==================================
# B. PROGRESSION MARKERS
# ==================================

cat("-- B. Progression markers\n\n")

prog_sig <- progression %>% filter(p_fdr < 0.05) %>% arrange(p_value)
prog_nom <- progression %>% filter(p_value < 0.05) %>% arrange(p_value)

cat("  FDR < 0.05:", nrow(prog_sig), " | Nominal:", nrow(prog_nom), "\n")
prog_sig %>% head(10) %>%
  dplyr::select(uniprot, gene_symbol, cor_r, p_value, p_fdr, grounding_score) %>%
  as.data.frame() %>% print(row.names = FALSE)
cat("\n")


# =================================
# C. concordR 0.2.0 GATES (route, specificity, concordance, audit call)
# =================================


cat("-- C. concordR gates (plasma, systemic)\n\n")
cat("  concordR version:", as.character(packageVersion("concordR")), "\n")
if (packageVersion("concordR") < "0.2.0" || !"details" %in% names(formals(concordR::triage)))
  stop("This script needs concordR >= 0.2.0.")

dap_genes <- daps %>% filter(sig) %>%
  left_join(grounding_dedup, by = "uniprot") %>%
  filter(!is.na(gene_symbol)) %>% pull(gene_symbol) %>% unique()
prog_fdr_genes <- prog_sig %>% filter(!is.na(gene_symbol)) %>% pull(gene_symbol) %>% unique()
prog_nom_genes <- prog_nom %>% filter(!is.na(gene_symbol)) %>% pull(gene_symbol) %>% unique()
all_genes <- grounding_dedup %>% filter(!is.na(gene_symbol)) %>% pull(gene_symbol) %>% unique()

cat("  Gene sets: DAPs", length(dap_genes), "| progression FDR", length(prog_fdr_genes),
    "| progression nominal", length(prog_nom_genes), "| background", length(all_genes), "\n\n")

run_gates <- function(genes, label) {
  tr <- triage(genes, sample_type = "plasma", details = TRUE)
  tr$call <- audit_call(tr)
  tr$set <- label
  tr
}
triage_dap <- run_gates(dap_genes, "dap")
triage_fdr <- run_gates(prog_fdr_genes, "prog_fdr")
triage_nom <- run_gates(prog_nom_genes, "prog_nominal")
triage_bg <- run_gates(all_genes, "background")

regulated <- function(tr) tr$route %in% c("secreted", "membrane")

gate_summary <- function(tr, label) {
  data.frame(
    set = label, n_genes = nrow(tr),
    pct_secreted = round(100 * mean(tr$route == "secreted"), 1),
    pct_membrane = round(100 * mean(tr$route == "membrane"), 1),
    pct_regulated = round(100 * mean(regulated(tr)), 1),
    pct_intracellular = round(100 * mean(tr$route == "intracellular"), 1),
    pct_route_unknown = round(100 * mean(tr$route == "unknown"), 1),
    pct_specific = round(100 * mean(tr$specificity == "specific"), 1),
    pct_broad = round(100 * mean(tr$specificity == "broad"), 1),
    pct_sparse = round(100 * mean(tr$specificity == "sparse"), 1),
    pct_spec_na = round(100 * mean(tr$specificity == "na"), 1),
    pct_concordant = round(100 * mean(tr$concordance == "concordant"), 1),
    pct_not_in_atlas = round(100 * mean(tr$concordance == "not_in_atlas"), 1),
    pct_accept = round(100 * mean(tr$call == "accept"), 1),
    pct_caution = round(100 * mean(tr$call == "caution"), 1),
    pct_decline = round(100 * mean(tr$call == "decline"), 1)
  )
}
gate_tbl <- bind_rows(gate_summary(triage_dap, "dap"), gate_summary(triage_fdr, "prog_fdr"),
                      gate_summary(triage_nom, "prog_nominal"), gate_summary(triage_bg, "background"))
cat("  Gate summary (% of genes in each set):\n"); print(gate_tbl, row.names = FALSE); cat("\n")
write.csv(gate_tbl, "../output/summaries/fb_concordr_gates.csv", row.names = FALSE)

# Enrichment of the intracellular (damage release) route vs background: the replacement
# OR < 1 = markers are LESS damage-release than background.
route_or <- function(tr, bg, label) {
  a <- sum(tr$route == "intracellular"); b <- sum(regulated(tr))
  c <- sum(bg$route == "intracellular"); d <- sum(regulated(bg))
  ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2))
  cat(sprintf("    %s: intracellular vs regulated, OR = %.2f (95%% CI %.2f-%.2f), p = %.2e  [%d intracellular, %d regulated]\n",
              label, ft$estimate, ft$conf.int[1], ft$conf.int[2], ft$p.value, a, b))
  data.frame(set = label, n_intracellular = a, n_regulated = b, or = unname(ft$estimate),
             ci_lo = ft$conf.int[1], ci_hi = ft$conf.int[2], p = ft$p.value)
}
cat("  Damage-release route enrichment vs background (route unknown excluded):\n")
route_tbl <- bind_rows(route_or(triage_dap, triage_bg, "DAPs"),
                       route_or(triage_fdr, triage_bg, "Progression FDR"),
                       route_or(triage_nom, triage_bg, "Progression nominal"))
write.csv(route_tbl, "../output/summaries/fb_concordr_route_enrichment.csv", row.names = FALSE)
cat("\n")

cat("  Intracellular-route progression markers (FDR set), with tissue-specificity:\n")
triage_fdr %>% filter(route == "intracellular") %>%
  dplyr::select(gene_symbol, specificity, primary_tissue, concordance, call, localisation) %>%
  as.data.frame() %>% print(row.names = FALSE)
cat("\n")

cat("  Funnel (first failing gate), progression FDR set:\n")
print(funnel_class(triage_fdr, table = TRUE), row.names = FALSE)
cat("\n")


# =================================================================
# C2. TISSUE PRESENCE IN MN RELEVANT VS PERIPHERAL CONTEXTS
# =================================================================
# The tissue gate in 0.2.0 reads the three-source presence vote (HPA IHC, GTEx MS, PaxDb)
# in a named tissue. Q: is each marker's protein detected in cerebral cortex
# and spinal cord, and in skeletal and heart muscle? A marker absent from the nervous
# tissues but present in muscle is the tissue-level version of the peripheral claim.

cat("-- C2. Tissue presence vote in named contexts\n\n")

hpa_vocab <- unique(load_atlas("hpa")$tissue)
pick_tissue <- function(pattern) { h <- hpa_vocab[grepl(pattern, hpa_vocab, ignore.case = TRUE)]; h[1] }
ctx_wanted <- c(nervous_cortex = pick_tissue("cerebral cortex"),
                nervous_spinal = pick_tissue("spinal"),
                muscle_skel = pick_tissue("skeletal muscle"),
                muscle_heart = pick_tissue("heart"),
                liver = pick_tissue("^liver"),
                adipose = pick_tissue("adipose"))
ctx_wanted <- ctx_wanted[!is.na(ctx_wanted)]
cat("  Context tissues resolved from the hpa_tissue vocabulary:\n"); print(ctx_wanted); cat("\n")

presence_by_ctx <- function(genes, label) {
  tr <- triage(genes, sample_type = "plasma", context_tissue = unname(ctx_wanted))
  tr %>% mutate(set = label,
                detected_any = tissue %in% c("detected", "detected_1platform"),
                absent = tissue == "mismatch")
}
ctx_dap <- presence_by_ctx(dap_genes, "dap")
ctx_fdr <- presence_by_ctx(prog_fdr_genes, "prog_fdr")
ctx_nom <- presence_by_ctx(prog_nom_genes, "prog_nominal")
ctx_bg <- presence_by_ctx(all_genes, "background")

ctx_tbl <- bind_rows(ctx_dap, ctx_fdr, ctx_nom, ctx_bg) %>%
  group_by(set, context_tissue) %>%
  summarise(n = n(),
            pct_detected = round(100 * mean(detected_any), 1),
            pct_mismatch = round(100 * mean(absent), 1),
            pct_not_in_panel = round(100 * mean(tissue == "not_in_panel"), 1),
            .groups = "drop") %>%
  arrange(context_tissue, set)
cat("  Protein presence by context tissue (% of genes):\n"); print(as.data.frame(ctx_tbl), row.names = FALSE); cat("\n")
write.csv(ctx_tbl, "../output/summaries/fb_concordr_context_presence.csv", row.names = FALSE)

# per-gene wide table for the FDR set: detected in cortex? in spinal cord? in muscle?
ctx_fdr_wide <- ctx_fdr %>%
  dplyr::select(gene_symbol, context_tissue, tissue) %>%
  tidyr::pivot_wider(names_from = context_tissue, values_from = tissue)
write.csv(ctx_fdr_wide, "../output/landscape/context_presence_progression_fdr.csv", row.names = FALSE)

# Fisher: detected in skeletal muscle vs cortex, markers vs background
ctx_fisher <- function(ctx_set, ctx_bg, tissue_name, label) {
  a <- sum(ctx_set$detected_any[ctx_set$context_tissue == tissue_name])
  b <- sum(ctx_set$context_tissue == tissue_name) - a
  c <- sum(ctx_bg$detected_any[ctx_bg$context_tissue == tissue_name])
  d <- sum(ctx_bg$context_tissue == tissue_name) - c
  ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2))
  data.frame(set = label, tissue = tissue_name, n_detected = a, n_total = a + b,
             or_vs_background = unname(ft$estimate), p = ft$p.value)
}
ctx_enrich <- bind_rows(lapply(names(ctx_wanted), function(k) {
  bind_rows(ctx_fisher(ctx_dap, ctx_bg, ctx_wanted[[k]], "dap"),
            ctx_fisher(ctx_fdr, ctx_bg, ctx_wanted[[k]], "prog_fdr"))
}))
cat("  Presence enrichment vs background:\n")
print(ctx_enrich %>% mutate(across(where(is.numeric), ~ signif(.x, 3))), row.names = FALSE)
write.csv(ctx_enrich, "../output/summaries/fb_concordr_context_enrichment.csv", row.names = FALSE)
cat("\n")


# =============================================================================
# C3. TISSUE-ORIGIN PROFILE FROM THE 59-TISSUE PRESENCE VOTE (entropy, P(top))
# =============================================================================
# Per gene, per HPA tissue, weight =
# IHC level (0-3)/3 + GTEx detected (0/1) + PaxDb detected (0/1), over tissues with at
# least one voter. Entropy is reported raw and normalised by log2(n tissues covered) so
# genes with different platform coverage are comparable. The claim stays the same:
# progression markers are diffuse, with no tissue, nervous or otherwise, dominating.

cat("-- C3. Tissue-origin profile (59-tissue presence vote)\n\n")

ihc_num <- c("High" = 3, "Medium" = 2, "Low" = 1, "Not detected" = 0)

tissue_profile <- function(genes) {
  qp <- query_presence(genes, missing = "silent")
  if (nrow(qp) == 0) return(NULL)
  qp %>%
    mutate(w_ihc = dplyr::coalesce(unname(ihc_num[as.character(ihc_level)]), 0) / 3,
           w_gtex = ifelse(!is.na(gtex_call)  & gtex_call  == "detected", 1, 0),
           w_paxdb = ifelse(!is.na(paxdb_call) & paxdb_call == "detected", 1, 0),
           covered = !is.na(ihc_level) | !is.na(gtex_call) | !is.na(paxdb_call),
           weight = w_ihc + w_gtex + w_paxdb) %>%
    filter(covered) %>%
    group_by(gene_symbol) %>%
    filter(sum(weight) > 0) %>%
    mutate(prob = weight / sum(weight),
           n_tissues = n(),
           top_tissue = tissue[which.max(prob)],
           top_prob = max(prob),
           entropy = -sum(ifelse(prob > 0, prob * log2(prob), 0)),
           entropy_norm = ifelse(n_tissues > 1, entropy / log2(n_tissues), 0)) %>%
    ungroup()
}
profile_summary <- function(pr) {
  pr %>% group_by(gene_symbol) %>%
    summarise(top_tissue = top_tissue[1], top_prob = top_prob[1], n_tissues = n_tissues[1],
              entropy = entropy[1], entropy_norm = entropy_norm[1], .groups = "drop")
}
agg_profile <- function(pr) {
  pr %>% group_by(tissue) %>%
    summarise(total_weight = sum(weight), n_genes = n_distinct(gene_symbol), .groups = "drop") %>%
    mutate(pct = 100 * total_weight / sum(total_weight)) %>% arrange(desc(pct))
}

prof_fdr <- tissue_profile(prog_fdr_genes)
prof_nom <- tissue_profile(prog_nom_genes)
prof_dap <- tissue_profile(dap_genes)
prof_bg  <- tissue_profile(all_genes)

gs_fdr <- profile_summary(prof_fdr); gs_nom <- profile_summary(prof_nom)
gs_dap <- profile_summary(prof_dap); gs_bg  <- profile_summary(prof_bg)

entropy_row <- function(gs, label) {
  data.frame(set = label, n_genes = nrow(gs),
             entropy_median = round(median(gs$entropy), 2),
             entropy_norm_median = round(median(gs$entropy_norm), 3),
             p_top_median = round(median(gs$top_prob), 3),
             n_p_top_gt_0.5 = sum(gs$top_prob > 0.5),
             n_tissues_median = median(gs$n_tissues),
             wilcox_vs_bg_norm = if (label == "background") NA else
               wilcox.test(gs$entropy_norm, gs_bg$entropy_norm)$p.value)
}
entropy_tbl <- bind_rows(entropy_row(gs_dap, "dap"), entropy_row(gs_fdr, "prog_fdr"),
                         entropy_row(gs_nom, "prog_nominal"), entropy_row(gs_bg, "background"))
cat("  Diffuseness (lead with these, not the ranked tissue list):\n")
print(entropy_tbl, row.names = FALSE); cat("\n")
write.csv(entropy_tbl, "../output/summaries/fb_tissue_entropy.csv", row.names = FALSE)

agg_fdr <- agg_profile(prof_fdr); agg_bg <- agg_profile(prof_bg); agg_dap <- agg_profile(prof_dap)
cat("  Aggregate presence weight, progression FDR set (top 10 of", nrow(agg_fdr), "tissues):\n")
print(head(agg_fdr, 10)); cat("\n")
cat("  Same tissues in background, for reference:\n")
print(agg_bg %>% filter(tissue %in% head(agg_fdr, 10)$tissue)); cat("\n")
cat("  Top tissue counts, progression FDR set:\n"); print(dplyr::count(gs_fdr, top_tissue, sort = TRUE)); cat("\n")
cat("  Top tissue counts, DAPs:\n"); print(dplyr::count(gs_dap, top_tissue, sort = TRUE)); cat("\n")

# muscle share of the DAPs by two definitions, so the "muscle-dominated" claim is stated
muscle_names <- unname(c(ctx_wanted["muscle_skel"], ctx_wanted["muscle_heart"]))
dap_muscle_top <- gs_dap %>% filter(top_tissue %in% muscle_names)
dap_muscle_primary <- triage_dap %>% filter(specificity %in% c("specific", "sparse")) %>%
  filter(vapply(strsplit(as.character(primary_tissue), ", ", fixed = TRUE),
                function(x) any(x %in% muscle_names), logical(1)))
cat("  DAPs with a muscle top tissue in the presence profile:", nrow(dap_muscle_top), "of", nrow(gs_dap), "profiled\n")
cat("  DAPs specific/sparse with muscle in the IHC primary set:", nrow(dap_muscle_primary), "of",
    sum(triage_dap$specificity %in% c("specific", "sparse")), "specific/sparse\n")
cat("    ", paste(sort(unique(c(dap_muscle_top$gene_symbol, dap_muscle_primary$gene_symbol))), collapse = ", "), "\n\n")

cat("  Top tissue counts, background (for the denominator):\n")
bg_top <- dplyr::count(gs_bg, top_tissue, sort = TRUE) %>% mutate(pct = round(100 * n / sum(n), 1))
print(head(bg_top, 10))
cortex_cmp <- bind_rows(
  gs_fdr %>% summarise(set = "prog_fdr", n = n(), cortex_top = sum(top_tissue == "Cerebral cortex")),
  gs_dap %>% summarise(set = "dap",      n = n(), cortex_top = sum(top_tissue == "Cerebral cortex")),
  gs_bg  %>% summarise(set = "background", n = n(), cortex_top = sum(top_tissue == "Cerebral cortex"))
) %>% mutate(pct = round(100 * cortex_top / n, 1))
cat("  Cortex as top presence tissue:\n"); print(cortex_cmp)
ft_cortex <- fisher.test(matrix(c(cortex_cmp$cortex_top[1], cortex_cmp$n[1] - cortex_cmp$cortex_top[1],
                                  cortex_cmp$cortex_top[3], cortex_cmp$n[3] - cortex_cmp$cortex_top[3]), 2))
cat(sprintf("  Cortex-top, progression FDR vs background: OR = %.2f, p = %.3f\n\n", ft_cortex$estimate, ft_cortex$p.value))
write.csv(bind_rows(gs_fdr %>% mutate(set = "prog_fdr"), gs_dap %>% mutate(set = "dap"), gs_bg %>% mutate(set = "background")) %>%
            dplyr::count(set, top_tissue) %>% group_by(set) %>% mutate(pct = round(100 * n / sum(n), 1)) %>% ungroup(),
          "../output/summaries/fb_top_tissue_by_set.csv", row.names = FALSE)


write.csv(prof_fdr, "../output/landscape/tissue_profile_progression_fdr.csv", row.names = FALSE)
write.csv(gs_fdr,   "../output/landscape/tissue_profile_gene_summary_fdr.csv", row.names = FALSE)
write.csv(gs_nom,   "../output/landscape/tissue_profile_gene_summary_nominal.csv", row.names = FALSE)
write.csv(gs_dap,   "../output/landscape/tissue_profile_gene_summary_dap.csv", row.names = FALSE)
write.csv(agg_fdr,  "../output/landscape/tissue_profile_aggregate_fdr.csv", row.names = FALSE)
write.csv(agg_dap,  "../output/landscape/tissue_profile_aggregate_dap.csv", row.names = FALSE)
write.csv(agg_bg,   "../output/landscape/tissue_profile_aggregate_background.csv", row.names = FALSE)


# ============================
# D. GROUNDING MAP OVERLAY
# ============================

cat("-- D. Grounding overlay\n\n")

dap_grounding <- daps %>% filter(sig) %>% left_join(grounding_dedup, by = "uniprot")

cat("DAPs by grounding tier:\n")
dap_grounding %>% dplyr::count(grounding_tier) %>% mutate(pct = round(n / sum(n) * 100, 1)) %>%
  as.data.frame() %>% print(row.names = FALSE)

cat("\nProgression (FDR) by grounding tier:\n")
prog_sig %>% dplyr::count(grounding_tier) %>% mutate(pct = round(n / sum(n) * 100, 1)) %>%
  as.data.frame() %>% print(row.names = FALSE)
cat("\n")

# grounding tier x route: are the MN-grounded markers the intracellular ones?
cat("Progression (FDR) grounding tier x concordR route:\n")
prog_sig %>% left_join(triage_fdr %>% dplyr::select(gene_symbol, route), by = "gene_symbol") %>%
  with(print(table(grounding_tier, route, useNA = "ifany")))
cat("\n")

# Hits by tier against the panel: are progression hits or DAPs enriched in any tier?
tier_levels <- c("confirmed_full", "confirmed_partial", "rna_only_ms_gap",
                 "ms_only_rna_gap", "no_mn_expression")
panel_tiers <- grounding_dedup %>% filter(grounding_tier %in% tier_levels)

tier_enrichment <- function(hit_uniprots, label) {
  hit <- panel_tiers$uniprot %in% hit_uniprots
  tab <- table(tier = factor(panel_tiers$grounding_tier, levels = tier_levels), hit = hit)
  ft  <- fisher.test(tab, simulate.p.value = TRUE, B = 1e5)
  conf <- panel_tiers$grounding_tier %in% c("confirmed_full", "confirmed_partial")
  ft2  <- fisher.test(table(conf, hit))
  cat(sprintf("  %s: %d unique hit proteins mapped to a tier | Fisher across 5 tiers p = %.3f | MN-confirmed (full+partial) vs rest: OR = %.2f (95%% CI %.2f-%.2f), p = %.3f\n",
              label, sum(hit), ft$p.value, ft2$estimate, ft2$conf.int[1], ft2$conf.int[2], ft2$p.value))
  print(round(100 * prop.table(tab, 2), 1))
  data.frame(set = label, n_hit_proteins = sum(hit), fisher_5tier_p = ft$p.value,
             confirmed_or = unname(ft2$estimate), confirmed_ci_lo = ft2$conf.int[1],
             confirmed_ci_hi = ft2$conf.int[2], confirmed_p = ft2$p.value)
}
cat("\nTier enrichment of hits vs panel (unique proteins):\n")
tier_enrich_tbl <- bind_rows(
  tier_enrichment(unique(prog_sig$uniprot), "Progression FDR"),
  tier_enrichment(unique(daps$uniprot[daps$sig %in% TRUE]), "DAPs")
)
write.csv(tier_enrich_tbl, "../output/summaries/fb_tier_enrichment.csv", row.names = FALSE)

cat("\nUnique proteins: progression FDR", n_distinct(prog_sig$uniprot), "of", nrow(prog_sig), "aptamers |",
    "DAPs", n_distinct(daps$uniprot[daps$sig %in% TRUE]), "of", sum(daps$sig, na.rm = TRUE), "aptamers\n\n")

# ===================
# E. DRUGGABILITY
# ===================

cat("-- E. Druggability\n\n")

dgidb_candidates <- "../raw/interactions.tsv"

dgidb_hits <- dgidb_candidates[file.exists(dgidb_candidates)]
dgidb_path <- if (length(dgidb_hits)) dgidb_hits[1] else dgidb_candidates[1]

if (file.exists(dgidb_path)) {
  dgidb <- read.delim(dgidb_path, stringsAsFactors = FALSE)
  drug_targets <- unique(dgidb$gene_name)
  cat("  DGIdb:", nrow(dgidb), "interactions |", length(drug_targets), "targets |", basename(dgidb_path), "\n")

  dap_druggable  <- dap_grounding %>% mutate(druggable = gene_symbol %in% drug_targets)
  prog_druggable <- prog_sig %>% mutate(druggable = gene_symbol %in% drug_targets)
  cat("  DAPs druggable:", sum(dap_druggable$druggable, na.rm = TRUE), "/", nrow(dap_druggable),
      "aptamers (", n_distinct(dap_druggable$uniprot[dap_druggable$druggable]), "unique proteins )\n")
  cat("  Progression druggable:", sum(prog_druggable$druggable, na.rm = TRUE), "/", nrow(prog_druggable),
      "aptamers (", n_distinct(prog_druggable$uniprot[prog_druggable$druggable]), "unique proteins )\n")

  cat("\n  Druggable progression markers with concordR gates:\n")
  drug_prog <- prog_druggable %>%
    filter(druggable) %>%
    left_join(triage_fdr %>% dplyr::select(gene_symbol, route, specificity, primary_tissue, concordance, call),
              by = "gene_symbol") %>%
    left_join(gs_fdr %>% dplyr::select(gene_symbol, top_tissue, top_prob), by = "gene_symbol") %>%
    dplyr::select(gene_symbol, uniprot, cor_r, p_fdr, grounding_tier, route, specificity,
                  primary_tissue, call, top_tissue, top_prob) %>%
    mutate(across(where(is.numeric), ~ signif(.x, 3)))
  print(as.data.frame(drug_prog), row.names = FALSE)
  write.csv(drug_prog, "../output/landscape/druggable_progression_fdr.csv", row.names = FALSE)
  cat("\n  Druggable progression markers by grounding tier:\n"); print(table(drug_prog$grounding_tier))
  cat("  ... by route:\n"); print(table(drug_prog$route))
  cat("  ... by audit call:\n"); print(table(drug_prog$call))
} else {
  cat("  DGIdb not found at:", dgidb_path, "\n")
}
cat("\n")


# ==================
# F. VISUALISATIONS
# ==================

cat("-- F. Plots\n\n")

tier_colors <- c(no_mn_expression = "grey70", rna_only_ms_gap = "grey50",
                 ms_only_rna_gap = "#B2ABD2",
                 confirmed_partial = "#4393C3", confirmed_full = "#2166AC")
set_levels <- c("background", "dap", "prog_fdr")
set_colors <- c(background = "grey60", dap = "#2166AC", prog_fdr = "#B2182B")

# F1. Route composition by set
tiff("../output/landscape/route_by_set.tiff", width = 6, height = 4, units = "in", res = 300)
bind_rows(triage_bg, triage_dap, triage_fdr) %>%
  mutate(set = factor(set, levels = set_levels),
         route = factor(route, levels = c("secreted", "membrane", "intracellular", "unknown"))) %>%
  dplyr::count(set, route) %>% group_by(set) %>% mutate(pct = 100 * n / sum(n)) %>% ungroup() %>%
  ggplot(aes(x = set, y = pct, fill = route)) +
  geom_col(position = "stack", alpha = 0.9) +
  scale_fill_manual(values = c(secreted = "#2166AC", membrane = "#92C5DE", intracellular = "#B2182B", unknown = "grey80")) +
  labs(x = "", y = "% of genes", fill = "Route to plasma") + theme_bw()
dev.off()

# F2. Audit call by set
tiff("../output/landscape/audit_call_by_set.tiff", width = 6, height = 4, units = "in", res = 300)
bind_rows(triage_bg, triage_dap, triage_fdr) %>%
  mutate(set = factor(set, levels = set_levels)) %>%
  dplyr::count(set, call) %>% group_by(set) %>% mutate(pct = 100 * n / sum(n)) %>% ungroup() %>%
  ggplot(aes(x = set, y = pct, fill = call)) +
  geom_col(position = "stack", alpha = 0.9) +
  scale_fill_manual(values = c(accept = "#2166AC", caution = "#FDDBC7", decline = "#B2182B")) +
  labs(x = "", y = "% of genes", fill = "concordR audit call") + theme_bw()
dev.off()

# F3. Presence in named contexts, markers vs background
tiff("../output/landscape/context_presence.tiff", width = 7, height = 4.5, units = "in", res = 300)
ctx_tbl %>% filter(set %in% set_levels) %>%
  mutate(set = factor(set, levels = set_levels)) %>%
  ggplot(aes(x = context_tissue, y = pct_detected, fill = set)) +
  geom_col(position = "dodge", alpha = 0.9) +
  scale_fill_manual(values = set_colors) +
  labs(x = "", y = "% of genes with protein detected", fill = "") +
  theme_bw() + theme(axis.text.x = element_text(angle = 20, hjust = 1), legend.position = "top")
dev.off()

# F4. Volcano with grounding (unchanged)
tiff("../output/landscape/volcano_grounding.tiff", width = 7, height = 6, units = "in", res = 300)
volcano_data <- daps %>%
  left_join(grounding_dedup, by = "uniprot") %>%
  mutate(grounding_tier = replace_na(grounding_tier, "unmapped")) %>%
  filter(is.finite(logFC), is.finite(adj.P.Val))
ggplot(volcano_data, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(data = volcano_data %>% filter(!sig), color = "grey80", alpha = 0.3, size = 0.8) +
  geom_point(data = volcano_data %>% filter(sig & grounding_tier != "confirmed_full"),
             aes(color = grounding_tier), alpha = 0.6, size = 1.5) +
  geom_point(data = volcano_data %>% filter(sig & grounding_tier == "confirmed_full"),
             color = "#2166AC", alpha = 0.8, size = 2.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
  scale_color_manual(values = c(tier_colors, unmapped = "grey50")) +
  labs(x = "log2 FC (ALS vs Control)", y = "-log10 adjusted p-value", color = "MN grounding") +
  theme_bw()
dev.off()

# F5. Grounding DAPs vs all (unchanged)
tiff("../output/landscape/grounding_dap_vs_all.tiff", width = 5, height = 4, units = "in", res = 300)
bind_rows(
  dap_grounding %>% filter(!is.na(grounding_score)) %>% dplyr::count(grounding_score) %>%
    mutate(pct = n / sum(n) * 100, group = "DAPs"),
  grounding_dedup %>% filter(!is.na(grounding_score)) %>% dplyr::count(grounding_score) %>%
    mutate(pct = n / sum(n) * 100, group = "All SomaScan")
) %>%
  ggplot(aes(x = factor(grounding_score), y = pct, fill = group)) +
  geom_col(position = "dodge", alpha = 0.8) +
  scale_fill_manual(values = c("All SomaScan" = "grey60", "DAPs" = "#2166AC")) +
  labs(x = "Grounding score", y = "% of proteins", fill = "") +
  theme_bw() + theme(legend.position = "top")
dev.off()

# F6. Normalised tissue-origin entropy by set
tiff("../output/landscape/entropy_comparison.tiff", width = 5.5, height = 4, units = "in", res = 300)
bind_rows(gs_bg %>% mutate(set = "background"), gs_dap %>% mutate(set = "dap"), gs_fdr %>% mutate(set = "prog_fdr")) %>%
  mutate(set = factor(set, levels = set_levels)) %>%
  ggplot(aes(x = set, y = entropy_norm, fill = set)) +
  geom_violin(trim = FALSE, alpha = 0.5) + geom_boxplot(width = 0.15, outlier.shape = NA) +
  scale_fill_manual(values = set_colors) +
  labs(x = "", y = "Normalised tissue-presence entropy\n(1 = evenly present across covered tissues)") +
  theme_bw() + theme(legend.position = "none")
dev.off()

# F7. Aggregate presence weight, progression FDR set vs background, top 15 tissues
tiff("../output/landscape/tissue_presence_profile_fdr.tiff", width = 7, height = 5.5, units = "in", res = 300)
top15 <- head(agg_fdr, 15)$tissue
bind_rows(agg_fdr %>% mutate(set = "prog_fdr"), agg_bg %>% mutate(set = "background")) %>%
  filter(tissue %in% top15) %>%
  mutate(tissue = factor(tissue, levels = rev(top15)), set = factor(set, levels = c("background", "prog_fdr"))) %>%
  ggplot(aes(x = tissue, y = pct, fill = set)) +
  geom_col(position = "dodge", alpha = 0.9) + coord_flip() +
  scale_fill_manual(values = set_colors[c("background", "prog_fdr")]) +
  labs(x = "", y = "% of total presence weight", fill = "") +
  theme_bw() + theme(legend.position = "top")
dev.off()

# F8. P(top) vs normalised entropy, FDR set
tiff("../output/landscape/entropy_vs_top_prob.tiff", width = 5, height = 5, units = "in", res = 300)
ggplot(gs_fdr, aes(x = top_prob, y = entropy_norm)) +
  geom_point(alpha = 0.5) +
  geom_vline(xintercept = 0.5, linetype = "dashed", alpha = 0.4) +
  labs(x = "P(top tissue)", y = "Normalised entropy") + theme_bw()
dev.off()

cat("  Plots saved to output/landscape/\n\n")


# ================
# G. SAVE
# ================

cat("-- G. Saving\n\n")

known_markers <- c("NEFL","NEFH","NEFM","GFAP","UCHL1","CHIT1","CHI3L1","CHI3L2",
                   "TARDBP","SOD1","FUS","STMN2","PVALB","MMP9","KLK6")
km <- grounding_dedup %>% filter(gene_symbol %in% known_markers) %>%
  dplyr::select(gene_symbol, uniprot, grounding_tier) %>%
  left_join(daps %>% dplyr::select(uniprot, dap_logFC = logFC, dap_p = P.Value,
                                   dap_fdr = adj.P.Val, dap_sig = sig), by = "uniprot") %>%
  left_join(progression %>% dplyr::select(uniprot, prog_cor = cor_r, prog_p = p_value,
                                          prog_fdr = p_fdr), by = "uniprot") %>%
  left_join(triage_bg %>% dplyr::select(gene_symbol, route, specificity, primary_tissue, call), by = "gene_symbol") %>%
  arrange(gene_symbol)
cat("-- Known-marker landscape:\n"); print(as.data.frame(km), row.names = FALSE); cat("\n")
write.csv(km, "../output/summaries/fb_known_marker_landscape.csv", row.names = FALSE)

write.csv(data.frame(
  model = c("sex_only", "sex_plus_age"),
  n_sig = c(sum(daps_sexonly$sig, na.rm = TRUE), sum(daps$sig, na.rm = TRUE)),
  age_column = ifelse(is.na(pick_age), "none", pick_age)
), "../output/summaries/fb_dap_age_effect.csv", row.names = FALSE)

# DAP tissue table: top presence tissue and the IHC primary set, both stated
daps %>% filter(sig) %>%
  left_join(grounding_dedup %>% dplyr::select(uniprot, gene_symbol), by = "uniprot") %>%
  left_join(gs_dap %>% dplyr::select(gene_symbol, top_tissue, top_prob), by = "gene_symbol") %>%
  dplyr::count(top_tissue, sort = TRUE) %>%
  write.csv("../output/summaries/fb_dap_tissue.csv", row.names = FALSE)

write.csv(daps, "../output/landscape/daps_full.csv", row.names = FALSE)

gate_cols <- c("gene_symbol", "concordance", "route", "specificity", "primary_tissue", "call",
               "gene_class", "protein_tau", "tau", "localisation", "mechanism_tier")
dap_out <- daps %>% filter(sig) %>%
  left_join(grounding_dedup, by = "uniprot") %>%
  left_join(triage_dap %>% dplyr::select(any_of(gate_cols)), by = "gene_symbol") %>%
  left_join(gs_dap %>% dplyr::select(gene_symbol, top_tissue, top_prob, entropy_norm), by = "gene_symbol")
write.csv(dap_out, "../output/landscape/daps_annotated.csv", row.names = FALSE)

prog_out <- prog_nom %>%
  left_join(triage_nom %>% dplyr::select(any_of(gate_cols)), by = "gene_symbol") %>%
  left_join(gs_nom %>% dplyr::select(gene_symbol, top_tissue, top_prob, entropy_norm), by = "gene_symbol")
write.csv(prog_out, "../output/landscape/progression_annotated.csv", row.names = FALSE)

write.csv(triage_dap, "../output/landscape/triage_daps.csv", row.names = FALSE)
write.csv(triage_fdr, "../output/landscape/triage_progression_fdr.csv", row.names = FALSE)
write.csv(triage_nom, "../output/landscape/triage_progression_nominal.csv", row.names = FALSE)
write.csv(triage_bg, "../output/landscape/triage_background.csv", row.names = FALSE)

# top 15 FDR progression markers with every layer, for the Results table
prog_sig %>% head(15) %>%
  dplyr::select(gene_symbol, uniprot, cor_r, p_fdr, grounding_tier) %>%
  left_join(triage_fdr %>% dplyr::select(gene_symbol, route, specificity, primary_tissue, call), by = "gene_symbol") %>%
  left_join(gs_fdr %>% dplyr::select(gene_symbol, top_tissue, top_prob, entropy_norm), by = "gene_symbol") %>%
  left_join(ctx_fdr_wide, by = "gene_symbol") %>%
  mutate(across(where(is.numeric), ~ signif(.x, 3))) %>%
  write.csv("../output/summaries/fb_progression_top15_annotated.csv", row.names = FALSE)

cat("  All outputs saved to output/landscape/\n")
cat("\n== SCRIPT 03 COMPLETE ==\n")
