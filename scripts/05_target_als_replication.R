################################################################################
#  05. External replication in Target ALS (Natural History Study)
#
#  Part A  cross-sectional ALS vs control at baseline on each platform, then
#          direction and rank concordance of the Answer ALS DA<- present, and a
#          panel-wide concordance.
#
#  Part B  ALSFRS-R BLUP slopes in Target ALS, per-protein BLUP slopes on each
#          platform in ALS patients with >= 2 visits, slope-on-slope, and sign
#          concordance against the Answer ALS progression markers.
################################################################################

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readxl); library(limma); library(lme4)
})

ta_dir <- "../raw/target_als"

cat("================================================================\n")
cat("  05: Target ALS replication\n")
cat("================================================================\n\n")

na_codes <- function(x) { x <- suppressWarnings(as.numeric(x)); x[!is.na(x) & x < 0] <- NA; x }

# --- 0. Clinical -----------------------
clin <- read.csv(file.path(ta_dir, "nhs_clinical_data.csv"), check.names = FALSE, stringsAsFactors = FALSE)
clin <- clin %>% transmute(
  subject_id = `Subject ID`,
  visit_number = as.integer(`Visit Number`),
  group = ifelse(`Clinical Group` == "ALS", "ALS", "Control"),
  sex = ifelse(`Sex Reported` %in% c("Male", "Female"), `Sex Reported`, NA),
  age = na_codes(`Age In Years`),
  alsfrs = na_codes(`Alsfrs R Score Total`),
  days = if ("Days From Baseline" %in% names(.)) na_codes(`Days From Baseline`) else NA_real_
)
time_axis <- if (any(!is.na(clin$days))) "days" else "visit"
clin <- clin %>% mutate(time_years = if (time_axis == "days") days / 365.25 else as.numeric(visit_number - 1))
subj <- clin %>% arrange(subject_id, visit_number) %>% group_by(subject_id) %>%
  summarise(group = group[1], sex = sex[!is.na(sex)][1], age_base = age[!is.na(age)][1], .groups = "drop")
cat("  Clinical:", nrow(clin), "visits |", nrow(subj), "subjects (", sum(subj$group == "ALS"), "ALS ) | time axis:", time_axis, "\n")
cat("  ALS subjects with >= 2 ALSFRS-R:", sum(table(clin$subject_id[clin$group == "ALS" & !is.na(clin$alsfrs)]) >= 2), "\n\n")

# ---- 0b. Answer ALS results to replicate ---------------------------------
aals_dap  <- read.csv("../output/landscape/daps_annotated.csv", stringsAsFactors = FALSE) %>%
  filter(sig %in% TRUE, !is.na(gene_symbol)) %>%
  group_by(gene_symbol) %>% summarise(aals_logFC = logFC[which.min(adj.P.Val)], .groups = "drop")
aals_prog <- read.csv("../output/grounding/progression_with_grounding.csv", stringsAsFactors = FALSE) %>%
  filter(p_fdr < 0.05, !is.na(gene_symbol)) %>%
  group_by(gene_symbol) %>% summarise(aals_cor_r = cor_r[which.min(p_value)], .groups = "drop")
gdedup <- readRDS("../output/grounding/grounding_dedup.rds") %>% dplyr::select(uniprot, gene_symbol)
aals_all <- read.csv("../output/landscape/daps_full.csv", stringsAsFactors = FALSE) %>%
  filter(!is.na(uniprot)) %>% left_join(gdedup, by = "uniprot") %>% filter(!is.na(gene_symbol)) %>%
  group_by(gene_symbol) %>% summarise(aals_logFC = logFC[which.min(adj.P.Val)], .groups = "drop")
aals_allprog <- read.csv("../output/grounding/progression_with_grounding.csv", stringsAsFactors = FALSE) %>%
  filter(!is.na(gene_symbol)) %>% group_by(gene_symbol) %>% summarise(aals_cor_r = cor_r[which.min(p_value)], .groups = "drop")
cat("  Answer ALS: DAP genes", nrow(aals_dap), "| progression FDR genes", nrow(aals_prog), "\n\n")

# --- 1. TMT plasma: load, log2, average replicates ----------------------
tmt_x <- file.path(ta_dir, "ms/lc_ms_proteomics_tmt_plasma.xlsx")
info <- read_excel(tmt_x, sheet = "1-plasma_sample_info") %>% as.data.frame()
prot <- read_excel(tmt_x, sheet = "2-plasma_protein_info") %>% as.data.frame()
rat <- read_excel(tmt_x, sheet = "4-plasma_protein_TMT_ratios") %>% as.data.frame()
names(info) <- make.names(names(info)); names(prot) <- make.names(names(prot))
names(prot)[1] <- "uniprot"; names(rat)[1] <- "uniprot"
prot$gene_symbol <- sub(";.*$", "", prot$Gene.Symbol)

info <- info %>% filter(technical.replicate != -8888) %>%
  mutate(key = set.label, visit.number = as.integer(visit.number),
         grp = paste(subject.id, visit.number, sep = "|"))
mat <- as.matrix(rat[, -1]); rownames(mat) <- rat$uniprot; mode(mat) <- "numeric"
common <- intersect(colnames(mat), info$key)
mat <- log2(mat[, common]); info <- info[match(common, info$key), ]
cat("  TMT plasma:", nrow(mat), "proteins |", ncol(mat), "channels |", n_distinct(info$subject.id), "subjects\n")

mat_v <- t(rowsum(t(mat), info$grp, na.rm = TRUE) / as.vector(table(info$grp)[sort(unique(info$grp))]))
info_v <- info %>% distinct(grp, .keep_all = TRUE); info_v <- info_v[match(colnames(mat_v), info_v$grp), ]
info_v <- info_v %>% left_join(subj, by = c("subject.id" = "subject_id")) %>%
  left_join(clin %>% dplyr::select(subject_id, visit_number, time_years), by = c("subject.id" = "subject_id", "visit.number" = "visit_number"))
cat("  subject-visits:", ncol(mat_v), "| ALS with >= 2 visits:",
    sum(table(info_v$subject.id[info_v$group == "ALS"]) >= 2), "\n")

# --- 2. Cross-sectional ALS vs control, a function used for both platforms ----------
xs_limma <- function(m, meta, label, set_col = NULL) {
  meta$group <- factor(meta$group, levels = c("Control", "ALS"))
  ok <- !is.na(meta$group) & !is.na(meta$sex) & !is.na(meta$age_base)
  cov_note <- "sex + age"
  if (sum(ok) < 0.8 * nrow(meta)) { ok <- !is.na(meta$group); cov_note <- "none (sex/age incomplete)" }
  meta <- meta[ok, ]; m <- m[, ok, drop = FALSE]
  m <- m[rowSums(!is.na(m)) >= 0.5 * ncol(m), , drop = FALSE]
  f <- if (cov_note == "sex + age") "~ group + sex + age_base" else "~ group"
  if (!is.null(set_col)) { meta$set <- factor(meta[[set_col]]); f <- paste(f, "+ set") }
  design <- model.matrix(as.formula(f), data = meta)
  fit <- eBayes(lmFit(m, design))
  res <- topTable(fit, coef = "groupALS", number = Inf, sort.by = "none") %>% tibble::rownames_to_column("id")
  cat(sprintf("  %s: %d ALS vs %d control | covariates: %s%s | %d proteins | DAPs at FDR 5%%: %d\n",
              label, sum(meta$group == "ALS"), sum(meta$group == "Control"), cov_note,
              if (!is.null(set_col)) " + TMT set" else "", nrow(res), sum(res$adj.P.Val < 0.05, na.rm = TRUE)))
  res
}

cat("-- A. Cross-sectional --\n")
base_v <- info_v %>% group_by(subject.id) %>% slice_min(visit.number, n = 1, with_ties = FALSE) %>% ungroup()
tmt_res <- xs_limma(mat_v[, base_v$grp], base_v, "TMT plasma", set_col = "set") %>%
  dplyr::rename(uniprot = id) %>% left_join(prot %>% dplyr::select(uniprot, gene_symbol), by = "uniprot")

olk <- read.csv(file.path(ta_dir, "olink_matrix.csv"), check.names = FALSE, row.names = 1)
olk_meta <- subj[match(rownames(olk), subj$subject_id), ]
olk_res <- xs_limma(t(as.matrix(olk)), olk_meta, "Olink HT") %>% dplyr::rename(gene_symbol = id)
cat("\n")

# --- 3. Concordance with the Answer ALS DAPs ---------------------
concord <- function(res, ref, ref_col, label, what = "logFC") {
  d <- ref %>% inner_join(res %>% dplyr::select(gene_symbol, est = all_of(what), P.Value, adj.P.Val), by = "gene_symbol")
  if (nrow(d) < 5) { cat("  ", label, ": fewer than 5 present\n"); return(NULL) }
  same <- sign(d$est) == sign(d[[ref_col]])
  bt <- binom.test(sum(same), nrow(d), 0.5)
  sp <- cor.test(d[[ref_col]], d$est, method = "spearman", exact = FALSE)
  cat(sprintf("  %s: %d present | same direction %d (%.0f%%, binomial p = %.2g) | Spearman rho = %.2f (p = %.2g) | nominal p < 0.05 same direction %d | FDR < 0.05 same direction %d\n",
              label, nrow(d), sum(same), 100 * mean(same), bt$p.value, sp$estimate, sp$p.value,
              sum(same & d$P.Value < 0.05), sum(same & d$adj.P.Val < 0.05)))
  d$platform <- label; d
}
panel_wide <- function(res, ref, ref_col, label, what = "logFC") {
  d <- ref %>% inner_join(res %>% dplyr::select(gene_symbol, est = all_of(what)), by = "gene_symbol")
  sp <- cor.test(d[[ref_col]], d$est, method = "spearman", exact = FALSE)
  cat(sprintf("  %s, all shared proteins: n = %d | Spearman rho = %.3f (p = %.2g)\n", label, nrow(d), sp$estimate, sp$p.value))
}

cat("-- A2. Answer ALS DAPs on the replication platforms --\n")
cc <- bind_rows(concord(tmt_res, aals_dap, "aals_logFC", "TMT plasma"),
                concord(olk_res, aals_dap, "aals_logFC", "Olink HT"))
print(cc %>% arrange(platform, desc(abs(aals_logFC))) %>%
        dplyr::select(platform, gene_symbol, aals_logFC, est, P.Value, adj.P.Val) %>%
        mutate(across(where(is.numeric), ~ signif(.x, 3))) %>% as.data.frame(), row.names = FALSE)
cat("\n  Background for comparison:\n")
panel_wide(tmt_res, aals_all, "aals_logFC", "TMT plasma"); panel_wide(olk_res, aals_all, "aals_logFC", "Olink HT")
write.csv(cc, "../output/summaries/fb_replication_dap_concordance.csv", row.names = FALSE)
write.csv(tmt_res, "../output/replication/target_als_tmt_plasma_dap.csv", row.names = FALSE)
write.csv(olk_res, "../output/replication/target_als_olink_dap.csv", row.names = FALSE)
cat("\n")

# --- 4. Longitudinal: ALSFRS-R slopes -------------------
cat("-- B. Longitudinal (time axis:", time_axis, ") --\n")
fx <- clin %>% filter(group == "ALS", !is.na(alsfrs), !is.na(time_years)) %>%
  group_by(subject_id) %>% filter(n() >= 2) %>% ungroup()
m_fx  <- lmer(alsfrs ~ time_years + (time_years | subject_id), data = fx)
fx_sl <- data.frame(subject_id = rownames(ranef(m_fx)$subject_id),
                    ALSFRS_slope = unname(fixef(m_fx)["time_years"]) + ranef(m_fx)$subject_id$time_years)
cat("  ALSFRS-R BLUP slopes:", nrow(fx_sl), "ALS subjects | mean slope",
    round(mean(fx_sl$ALSFRS_slope), 2), "per", ifelse(time_axis == "days", "year", "visit"), "\n")

slope_design <- function(long, ids, label, min_n = 15) {
  # long: subject_id, time_years, value ; ids: vector of protein ids to loop over is implicit in `long` being a list
  one <- function(v) {
    v <- v[is.finite(v$value), ]
    if (n_distinct(v$subject_id) < min_n) return(NULL)
    bl <- tryCatch(suppressMessages(suppressWarnings(
      lmer(value ~ time_years + (time_years | subject_id), data = v, control = lmerControl(calc.derivs = FALSE)))),
      error = function(e) NULL)
    if (is.null(bl)) return(NULL)
    re <- ranef(bl)$subject_id
    ps <- data.frame(subject_id = rownames(re), prot_slope = unname(fixef(bl)["time_years"]) + re$time_years) %>%
      inner_join(fx_sl, by = "subject_id")
    if (nrow(ps) < min_n || sd(ps$prot_slope) == 0) return(NULL)
    ct <- cor.test(ps$prot_slope, ps$ALSFRS_slope)
    data.frame(n = nrow(ps), cor_r = unname(ct$estimate), p_value = ct$p.value)
  }
  out <- bind_rows(lapply(names(long), function(id) { r <- one(long[[id]]); if (!is.null(r)) cbind(id = id, r) }))
  out$p_fdr <- p.adjust(out$p_value, "fdr")
  cat(sprintf("  %s: %d proteins tested in %d ALS | nominal %d | FDR %d\n", label, nrow(out),
              max(out$n), sum(out$p_value < 0.05, na.rm=T), sum(out$p_fdr < 0.05, na.rm=T)))
  out
}

# 4a. TMT plasma
li <- info_v %>% filter(group == "ALS", !is.na(time_years), subject.id %in% fx_sl$subject_id) %>%
  group_by(subject.id) %>% filter(n() >= 2) %>% ungroup()
cat("  TMT plasma: ALS with >= 2 timed visits and an ALSFRS-R slope:", n_distinct(li$subject.id), "\n")
tmt_long <- lapply(setNames(rownames(mat_v), rownames(mat_v)), function(u)
  data.frame(subject_id = li$subject.id, time_years = li$time_years, value = mat_v[u, li$grp]))
tmt_prog <- slope_design(tmt_long, NULL, "TMT plasma") %>% dplyr::rename(uniprot = id) %>%
  left_join(prot %>% dplyr::select(uniprot, gene_symbol), by = "uniprot")
write.csv(tmt_prog, "../output/replication/target_als_tmt_plasma_progression.csv", row.names = FALSE)

# 4b Olink but it seems to be on on baseline only
olk_prog <- NULL
# 4c. Concordance with the Answer ALS progression markers
cat("\n-- B2. Answer ALS progression markers on the replication platforms --\n")
pc <- bind_rows(concord(tmt_prog %>% mutate(adj.P.Val = p_fdr, P.Value = p_value), aals_prog, "aals_cor_r", "TMT plasma", what = "cor_r"),
                if (!is.null(olk_prog)) concord(olk_prog %>% mutate(adj.P.Val = p_fdr, P.Value = p_value), aals_prog, "aals_cor_r", "Olink HT", what = "cor_r"))
print(pc %>% arrange(platform, desc(abs(aals_cor_r))) %>%
        dplyr::select(platform, gene_symbol, aals_cor_r, est, P.Value, adj.P.Val) %>%
        mutate(across(where(is.numeric), ~ signif(.x, 3))) %>% as.data.frame(), row.names = FALSE)
cat("\n  Background for comparison:\n")
panel_wide(tmt_prog, aals_allprog, "aals_cor_r", "TMT plasma", what = "cor_r")

cat("\n== SCRIPT 05 COMPLETE ==\n")
