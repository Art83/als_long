################################################################################
#  01. MAP. Somascan in plasma of patients vs MS, RNA-seq, ATAC in iPSC
#
#  For each SomaScan plasma aptamer, determine whether it has
#  molecular support in iPSC-derived motor neurons across three layers:
#    - RNA-seq (transcription)
#    - Mass spec proteomics (translation)
#    - ATAC-seq (chromatin accessibility at promoter) <- primarily supportive
#
################################################################################


suppressPackageStartupMessages({
  library(dplyr)
  library(GenomicRanges)
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(ggplot2)
})


cat("================================================================\n")
cat("  01: GENE-LEVEL GROUNDING MAP\n")
cat("================================================================\n\n")

#-------- 0. Load data (raw files from Answer ALS)--------------------------

mapper <- read.csv("../raw/transcriptomics/4_matrix/Sample Mapping Information/Sample Mapping File Feb 2024.csv", check.names = FALSE)
meta <- read.csv("../proc/meta.csv")
annotation <- read.csv("../raw/somascan/documents/annotation.csv")

meta <- meta %>% filter(SUBJECT_GROUP %in% c("ALS", "Healthy Control"))


#-------- 1. Gene ID harmonisation ------------------------------------------
ens2sym <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
                                 keys = AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype = "ENSEMBL"),
                                 columns = c("SYMBOL", "UNIPROT"),
                                 keytype = "ENSEMBL")

ens2sym_unique <- ens2sym %>%
  filter(!is.na(SYMBOL)) %>%
  group_by(ENSEMBL) %>%
  summarise(gene_symbol = SYMBOL[1],
            uniprot_from_org = na.omit(UNIPROT)[1],
            .groups = "drop")

soma_anno <- annotation %>%
  dplyr::select(SeqId, UniProt.ID, Target.Name, Target.Full.Name) %>%
  filter(UniProt.ID != "" & !is.na(UniProt.ID)) %>%
  dplyr::rename(uniprot = UniProt.ID, target_name = Target.Name,
                target_full = Target.Full.Name)

uni2sym <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
                                 keys = unique(soma_anno$uniprot),
                                 columns = c("SYMBOL", "ENSEMBL"),
                                 keytype = "UNIPROT") %>%
  filter(!is.na(SYMBOL)) %>%
  group_by(UNIPROT) %>%
  summarise(gene_symbol = SYMBOL[1],
            ensembl_id  = na.omit(ENSEMBL)[1],
            .groups = "drop")

soma_genes <- soma_anno %>%
  left_join(uni2sym, by = c("uniprot" = "UNIPROT")) %>%
  filter(!is.na(gene_symbol))

cat(" SomaScan aptamers with gene mapping:", nrow(soma_genes), "/", nrow(soma_anno), "\n\n")

# -------------- 2. RNA-seq: Motor neuron gene expression ---------------------------

cat("RNA-seq: Motor neuron gene expression\n\n")

rna_raw <- read.csv("../raw/transcriptomics/4_matrix/AnswerALS-290-T-v1-release5_rna-counts-matrix.csv", check.names = FALSE)
rownames(rna_raw) <- rna_raw[, 1] # ENSG, unnamed column
rna_raw <- rna_raw[, -1]

mapper_rna <- mapper %>% filter(Sample_ID %in% colnames(rna_raw))
meta_rna <- meta %>%
  inner_join(mapper_rna[, c("Participant_ID", "Sample_ID", "Batch")], by = "Participant_ID") %>%
  filter(Sample_ID %in% colnames(rna_raw))

rna_raw <- rna_raw[, meta_rna$Sample_ID]
stopifnot(all(colnames(rna_raw) == meta_rna$Sample_ID))
cat("RNA-seq samples:", ncol(rna_raw), "\n")

rna_raw <- rna_raw[rowSums(rna_raw) > 0, ]

dds <- DESeq2::DESeqDataSetFromMatrix(rna_raw, meta_rna, design = ~ 1)
dds <- DESeq2::estimateSizeFactors(dds)
vst_mat <- DESeq2::vst(dds, blind = TRUE)
rna_vst <- SummarizedExperiment::assay(vst_mat)

common_genes <- intersect(rownames(rna_vst), rownames(rna_raw))

rna_summary <- data.frame(
  ensembl_id = common_genes,
  mn_rna_mean = rowMeans(rna_vst[common_genes, ]),
  mn_rna_median = apply(rna_vst[common_genes, ], 1, median),
  mn_rna_pct_det = rowMeans(rna_raw[common_genes, ] > 5),
  stringsAsFactors = FALSE
)

vst_threshold <- quantile(rna_summary$mn_rna_mean, 0.10)
rna_summary$mn_rna_detected <- rna_summary$mn_rna_pct_det >= 0.20 &
  rna_summary$mn_rna_mean > vst_threshold

rna_summary$mn_rna_rank <- NA_real_
rna_summary$mn_rna_rank[rna_summary$mn_rna_detected] <-
  rank(rna_summary$mn_rna_mean[rna_summary$mn_rna_detected]) /
  sum(rna_summary$mn_rna_detected)

cat("Genes in RNA matrix:", nrow(rna_summary), "\n")
cat("Detected in MN:", sum(rna_summary$mn_rna_detected), "\n\n")

# ----------- 3. MS proteomics: Motor neuron protein detection ------------

cat(" MS proteomics: Motor neuron protein detection\n\n")

ms_raw <- read.csv("../raw/proteomics/4_matrix/AnswerALS-436-P_proteomics-protein-matrix_correctedImputed.txt", sep = '\t', check.names = FALSE)

ms_raw <- ms_raw[, !colnames(ms_raw) %in% c("nFragment", "nPeptide")]
ms_raw$uniprot <- sapply(strsplit(ms_raw$Protein, "[|]"), "[", 2)
ms_raw <- ms_raw[!is.na(ms_raw$uniprot), ]
ms_raw <- ms_raw[!duplicated(ms_raw$uniprot), ]
rownames(ms_raw) <- ms_raw$uniprot

sample_cols <- colnames(ms_raw)[!colnames(ms_raw) %in% c("Protein", "uniprot")]
sample_cols_clean <- gsub("_", "-", gsub("(-P).+$", "\\1", sample_cols))
col_map <- data.frame(original = sample_cols, clean = sample_cols_clean, stringsAsFactors = FALSE)

mapper_ms <- mapper %>% filter(Sample_ID %in% col_map$clean)
meta_ms <- meta %>%
  inner_join(mapper_ms[, c("Participant_ID", "Sample_ID")], by = "Participant_ID") %>%
  filter(Sample_ID %in% col_map$clean)

keep_cols <- col_map$original[col_map$clean %in% meta_ms$Sample_ID]
ms_expr <- ms_raw[, keep_cols]
rownames(ms_expr) <- ms_raw$uniprot

cat("MS proteomics samples:", ncol(ms_expr), "\n")
cat("Proteins quantified:", nrow(ms_expr), "\n")


ms_nonimputed_path <- "../raw/proteomics/4_matrix/AnswerALS-447-P_proteomics-protein-matrix_prebatchcorrected.txt"

# "Detected" = quantified (non-zero) in >= this fraction of MN samples.
# In this matrix missing vals are encoded as 0, not NA, and no protein sits below 0.21
# detection (the matrix is already filtered), so 0.50 is a real call, not a floor artefact.
ms_detect_min_frac <- 0.50
ms_detect_sensitivity <- c(0.30, 0.50, 0.70, 0.90)   # reported for a sensitivity check

detection_from_nonimputed <- function(path, keep_clean_ids) {
  if (!file.exists(path)) {
    warning("Non-imputed matrix not found at: ", path)
    return(NULL)
  }
  raw <- read.csv(path, sep = "\t", check.names = FALSE)
  raw <- raw[, !colnames(raw) %in% c("nFragment", "nPeptide")]

  # IDs look like "1/sp|A0AVT1|UBA6_HUMAN"; element 2 after splitting on "|" is the accession.
  # Rows without a pipe are standards/contaminants (e.g. "1/iRT_protein") and are dropped.
  pr    <- as.character(raw$Protein)
  parts <- strsplit(pr, "[|]")
  raw$uniprot <- vapply(parts, function(x) if (length(x) >= 2) x[2] else NA_character_, character(1))
  n_drop <- sum(is.na(raw$uniprot))
  if (n_drop > 0) cat("  Non-imputed: dropped", n_drop, "rows without a UniProt accession\n")
  raw <- raw[!is.na(raw$uniprot) & raw$uniprot != "", ]

  scols <- colnames(raw)[!colnames(raw) %in% c("Protein", "uniprot")]
  sclean <- gsub("_", "-", gsub("(-P).+$", "\\1", scols))
  keep <- scols[sclean %in% keep_clean_ids]
  if (length(keep) == 0) {
    warning("no MN samples matched: check column-name cleaning")
    return(NULL)
  }
  cat("Non-imputed: matched", length(keep), "of", length(scols), "sample columns to MN metadata\n")

  m <- as.matrix(raw[, keep, drop = FALSE])
  storage.mode(m) <- "numeric"
  present <- !is.na(m) & m > 0          # missingness encoded as 0 in this matrix
  pct <- rowMeans(present)

  # collapse any duplicate accessions (max = most permissive; none expected in the protein matrix)
  agg <- tapply(pct, raw$uniprot, max)
  out <- data.frame(uniprot = names(agg), mn_protein_pct_det = as.numeric(agg),
                    stringsAsFactors = FALSE)

  cat("Non-imputed: detection rate quantiles (10/25/50/75/90):",
      paste(sprintf("%.2f", quantile(out$mn_protein_pct_det, c(.1, .25, .5, .75, .9))),
            collapse = " "), "\n")
  for (thr in ms_detect_sensitivity) {
    cat(sprintf(" threshold %.2f -> %d / %d proteins detected (%.1f%%)\n",
                thr, sum(out$mn_protein_pct_det >= thr), nrow(out),
                100 * mean(out$mn_protein_pct_det >= thr)))
  }
  out
}

nonimp <- detection_from_nonimputed(ms_nonimputed_path, meta_ms$Sample_ID)

ms_expr_mat <- as.matrix(ms_expr)
ms_summary <- data.frame(
  uniprot = rownames(ms_expr),
  mn_protein_mean = rowMeans(ms_expr_mat, na.rm = TRUE),
  mn_protein_cv = apply(ms_expr_mat, 1, function(x) sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)),
  stringsAsFactors = FALSE
)


cat("non-imputed:", nrow(nonimp), "| imputed:", nrow(ms_summary),
      "| overlap:", length(intersect(nonimp$uniprot, ms_summary$uniprot)),
      "| non-imputed only:", length(setdiff(nonimp$uniprot, ms_summary$uniprot)),
      "| imputed only:", length(setdiff(ms_summary$uniprot, nonimp$uniprot)), "\n")

ms_summary <- full_join(ms_summary, nonimp, by = "uniprot") %>%
    mutate(
      ms_in_imputed    = !is.na(mn_protein_mean),
      ms_in_nonimputed = !is.na(mn_protein_pct_det),
      mn_protein_detected = ms_in_nonimputed & mn_protein_pct_det >= ms_detect_min_frac
    )

detect_method <- sprintf("non-imputed protein matrix, quantified in >= %.0f%% of MN samples",
                           ms_detect_min_frac * 100)


ms_summary$mn_protein_rank <- rank(ms_summary$mn_protein_mean, na.last = "keep") /
  sum(!is.na(ms_summary$mn_protein_mean))

cat("Detection method:", detect_method, "\n")
n_elig <- sum(ms_summary$ms_in_nonimputed)
n_det  <- sum(ms_summary$mn_protein_detected)
cat("Proteins detected in MN:", n_det, "/", n_elig, "with detection evidence",
    sprintf("(%.1f%%)", 100 * n_det / n_elig), "\n")
cat("ms_summary rows (union of both matrices):", nrow(ms_summary),
    "| no detection evidence (imputed matrix only):", nrow(ms_summary) - n_elig, "\n\n")

# ---- 4. ATAC-seq: Promoter accessibility -------------------------

cat(" ATAC-seq: Promoter accessibility\n\n")

atac_raw <- read.csv("../raw/epigenomics/epigenomics/4_matrix/AnswerALS-620-E-v1-release6_DiffBind-raw-counts-minOverlap-0.1.csv", check.names = FALSE)
rownames(atac_raw) <- atac_raw[, 1]
atac_raw <- atac_raw[, -1]

mapper_atac <- mapper %>% filter(Sample_ID %in% colnames(atac_raw))
meta_atac <- meta %>%
  inner_join(mapper_atac[, c("Participant_ID", "Sample_ID", "Batch")], by = "Participant_ID") %>%
  filter(Sample_ID %in% colnames(atac_raw))

atac_raw <- atac_raw[, meta_atac$Sample_ID]
cat("ATAC-seq samples:", ncol(atac_raw), "\n")
cat("Total peaks:", nrow(atac_raw), "\n")

peak_ids <- rownames(atac_raw)
parts <- do.call(rbind, strsplit(peak_ids, "-"))
gr <- GRanges(seqnames = parts[, 1],
              ranges = IRanges(as.integer(parts[, 2]), as.integer(parts[, 3])))

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
anno <- annotatePeak(gr, TxDb = txdb, tssRegion = c(-2000, 2000), annoDb = "org.Hs.eg.db", verbose = FALSE)
anno_df <- as.data.frame(anno)
anno_df$peak_id <- peak_ids

promoter_peaks <- anno_df %>%
  filter(grepl("Promoter", annotation), !is.na(SYMBOL))

cat("Promoter peaks:", nrow(promoter_peaks), "\n")
cat("Unique genes with promoter peaks:", n_distinct(promoter_peaks$SYMBOL), "\n")

atac_counts <- atac_raw[promoter_peaks$peak_id, ]
promoter_peaks$detected <- rowMeans(atac_counts > 5) >= 0.20
promoter_peaks$mean_count <- rowMeans(atac_counts)

atac_summary <- promoter_peaks %>%
  group_by(SYMBOL) %>%
  summarise(mn_chromatin_open = any(detected),
            mn_atac_max_score = max(mean_count),
            mn_atac_n_peaks   = n(),
            .groups = "drop") %>%
  dplyr::rename(gene_symbol = SYMBOL)

cat("Genes with open chromatin:", sum(atac_summary$mn_chromatin_open), "\n\n")

#------- 5. Assemble grounding map -----------------------------------

cat("-- 5. Assembling grounding map\n\n")

grounding <- soma_genes %>%
  left_join(rna_summary %>% dplyr::select(ensembl_id, mn_rna_detected, mn_rna_mean, mn_rna_rank), by = "ensembl_id") %>%
  left_join(ms_summary %>% dplyr::select(uniprot, mn_protein_detected, mn_protein_mean,
                                         mn_protein_rank, mn_protein_pct_det, ms_in_nonimputed),
            by = "uniprot") %>%
  left_join(atac_summary %>% dplyr::select(gene_symbol, mn_chromatin_open, mn_atac_max_score), by = "gene_symbol") %>%
  mutate(
    mn_rna_detected = tidyr::replace_na(mn_rna_detected, FALSE),
    mn_protein_detected = tidyr::replace_na(mn_protein_detected, FALSE),
    mn_chromatin_open = tidyr::replace_na(mn_chromatin_open, FALSE),
    ms_measured = tidyr::replace_na(ms_in_nonimputed, FALSE),
    grounding_score = as.integer(mn_rna_detected) + as.integer(mn_protein_detected) + as.integer(mn_chromatin_open),
    mn_confirmed = mn_rna_detected & mn_protein_detected,
    grounding_tier = case_when(
      mn_confirmed & mn_chromatin_open       ~ "confirmed_full",
      mn_confirmed & !mn_chromatin_open      ~ "confirmed_partial",
      mn_rna_detected & !mn_protein_detected ~ "rna_only_ms_gap",
      !mn_rna_detected & mn_protein_detected ~ "ms_only_rna_gap",
      TRUE ~ "no_mn_expression"
    )
  )

# --- 6. Summary -------------------------------------------

cat("Summary\n\n")
cat("Grounding score distribution:\n"); print(table(grounding$grounding_score)); cat("\n")
cat("Grounding tier distribution:\n");  print(table(grounding$grounding_tier));  cat("\n")
cat("RNA x Protein in MN:\n"); print(table(RNA = grounding$mn_rna_detected, Protein = grounding$mn_protein_detected)); cat("\n")

ms_uniprots <- ms_summary$uniprot
cat("SomaScan UniProts:", n_distinct(grounding$uniprot), " - MS proteins:", length(ms_uniprots),
    " - Overlap:", length(intersect(unique(grounding$uniprot), ms_uniprots)), "\n\n")

# --- MS coverage of the plasma panel (the limit on the protein layer) ---
cat("MS coverage of SomaScan panel:\n")
cat("  aptamers measured by MS:", sum(grounding$ms_measured), "/", nrow(grounding),
    sprintf("(%.1f%%)", 100 * mean(grounding$ms_measured)), "\n")
cat("  unique proteins measured by MS:",
    n_distinct(grounding$uniprot[grounding$ms_measured]), "/", n_distinct(grounding$uniprot),
    sprintf("(%.1f%%)", 100 * n_distinct(grounding$uniprot[grounding$ms_measured]) /
              n_distinct(grounding$uniprot)), "\n")
cat("  rna_only_ms_gap split (TRUE = measured but not detected):\n")
print(table(ms_measured = grounding$ms_measured[grounding$grounding_tier == "rna_only_ms_gap"]))
cat("\n")

cat("  genes with promoter peak but below threshold:",
    sum(!atac_summary$mn_chromatin_open), "\n")
cat("  aptamer-level: no peak annotated:", sum(is.na(grounding$mn_atac_max_score)),
    "| peak present but closed:",
    sum(!is.na(grounding$mn_atac_max_score) & !grounding$mn_chromatin_open), "\n")

cat("  rna_only_ms_gap breakdown:\n")
rog <- grounding$grounding_tier == "rna_only_ms_gap"
cat("    quantified, below threshold :", sum(rog & grounding$ms_measured), "\n")
cat("    imputed matrix only         :",
    sum(rog & !grounding$ms_measured & !is.na(grounding$mn_protein_mean)), "\n")
cat("    never quantified            :",
    sum(rog & !grounding$ms_measured & is.na(grounding$mn_protein_mean)), "\n")
# --- Tier sensitivity to the detection threshold ---
if (!is.null(nonimp)) {
  cat("Tier counts across detection thresholds:\n")
  sens <- lapply(ms_detect_sensitivity, function(thr) {
    det <- ms_summary$uniprot[ms_summary$ms_in_nonimputed &
                                ms_summary$mn_protein_pct_det >= thr]
    pd  <- grounding$uniprot %in% det
    tier <- dplyr::case_when(
      grounding$mn_rna_detected & pd & grounding$mn_chromatin_open  ~ "confirmed_full",
      grounding$mn_rna_detected & pd & !grounding$mn_chromatin_open ~ "confirmed_partial",
      grounding$mn_rna_detected & !pd                               ~ "rna_only_ms_gap",
      !grounding$mn_rna_detected & pd                               ~ "ms_only_rna_gap",
      TRUE                                                          ~ "no_mn_expression"
    )
    data.frame(threshold = thr, as.list(table(tier)), check.names = FALSE)
  })
  sens_tbl <- dplyr::bind_rows(sens)
  print(sens_tbl, row.names = FALSE)
  write.csv(sens_tbl, "../output/summaries/fb_tier_threshold_sensitivity.csv",
            row.names = FALSE)
  cat("\n")
}


#  Tier sensitivity to the RNA detection threshold
rna_sens_thresholds <- c(0.10, 0.20, 0.30, 0.50)

recount_tiers <- function(det_ens, label) {
  rd <- grounding$ensembl_id %in% det_ens
  tier <- dplyr::case_when(
    rd & grounding$mn_protein_detected & grounding$mn_chromatin_open ~ "confirmed_full",
    rd & grounding$mn_protein_detected & !grounding$mn_chromatin_open ~ "confirmed_partial",
    rd & !grounding$mn_protein_detected                               ~ "rna_only_ms_gap",
    !rd & grounding$mn_protein_detected                               ~ "ms_only_rna_gap",
    TRUE                                                              ~ "no_mn_expression"
  )
  data.frame(rule = label, n_genes_detected = length(det_ens),
             as.list(table(tier)), check.names = FALSE)
}

cat("Tier counts across RNA detection thresholds:\n")
rna_sens <- lapply(rna_sens_thresholds, function(thr) {
  det <- rna_summary$ensembl_id[rna_summary$mn_rna_pct_det >= thr &
                                  rna_summary$mn_rna_mean > vst_threshold]
  recount_tiers(det, sprintf("pct>=%.2f + vst floor", thr))
})
# same thresholds without the 10th-percentile vst floor
rna_sens <- c(rna_sens, lapply(rna_sens_thresholds, function(thr) {
  det <- rna_summary$ensembl_id[rna_summary$mn_rna_pct_det >= thr]
  recount_tiers(det, sprintf("pct>=%.2f only", thr))
}))

# named method comparison: edgeR::filterByExpr

keep_fbe <- edgeR::filterByExpr(as.matrix(rna_raw))
fbe_genes <- rownames(rna_raw)[keep_fbe]
rna_sens <- c(rna_sens, list(recount_tiers(fbe_genes, "edgeR::filterByExpr")))
  
ours <- rna_summary$ensembl_id[rna_summary$mn_rna_detected]
cat("  filterByExpr keeps", length(fbe_genes), "genes; current rule keeps", length(ours), "\n")
cat("  gene-level Jaccard:",
      round(100 * length(intersect(fbe_genes, ours)) / length(union(fbe_genes, ours)), 1), "%\n")
cat("  aptamer-level concordance:",
      round(100 * mean((grounding$ensembl_id %in% fbe_genes) == grounding$mn_rna_detected), 1), "%\n")


rna_sens_tbl <- dplyr::bind_rows(rna_sens)
print(rna_sens_tbl, row.names = FALSE)
write.csv(rna_sens_tbl, "../output/tables/rna_threshold_sensitivity.csv",
          row.names = FALSE)
cat("\n")

# ---- 7. Plots -----------------------------------------------------

cat("-- 7. Plots\n\n")

tier_colors <- c(no_mn_expression = "grey70", rna_only_ms_gap = "grey50",
                 ms_only_rna_gap = "#B2ABD2",
                 confirmed_partial = "#4393C3", confirmed_full = "#2166AC")

tiff("../output/grounding/grounding_score_distribution.tiff", width = 5, height = 4, units = "in", res = 300)
grounding %>% dplyr::count(grounding_score) %>% mutate(pct = n / sum(n) * 100) %>%
  ggplot(aes(x = factor(grounding_score), y = n)) +
  geom_col(fill = c("grey80", "grey60", "grey40", "grey20")) +
  geom_text(aes(label = sprintf("%d\n(%.0f%%)", n, pct)), vjust = -0.3, size = 3) +
  labs(x = "Grounding score", y = "SomaScan aptamers") + theme_bw()
dev.off()

tiff("../output/grounding/grounding_tier_distribution.tiff", width = 6, height = 4, units = "in", res = 300)
grounding %>% dplyr::count(grounding_tier) %>%
  mutate(grounding_tier = factor(grounding_tier, levels = names(tier_colors))) %>%
  ggplot(aes(x = grounding_tier, y = n, fill = grounding_tier)) +
  geom_col() + scale_fill_manual(values = tier_colors) +
  labs(x = "MN grounding tier", y = "SomaScan targets") +
  theme_bw() + theme(legend.position = "none", axis.text.x = element_text(angle = 15, hjust = 1))
dev.off()

tiff("../output/grounding/detection_pattern.tiff", width = 6, height = 5, units = "in", res = 300)
grounding %>%
  group_by(mn_rna_detected, mn_protein_detected, mn_chromatin_open) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(label = paste0(ifelse(mn_rna_detected, "RNA+", "RNA-"), " / ",
                        ifelse(mn_protein_detected, "Prot+", "Prot-"), " / ",
                        ifelse(mn_chromatin_open, "ATAC+", "ATAC-"))) %>%
  ggplot(aes(x = reorder(label, n), y = n)) +
  geom_col(fill = "steelblue") + coord_flip() +
  labs(x = "Detection pattern", y = "SomaScan targets") + theme_bw()
dev.off()

tiff("../output/grounding/rna_vs_protein_rank.tiff", width = 5, height = 5, units = "in", res = 300)
bd <- grounding %>% filter(mn_rna_detected & mn_protein_detected)
r_val <- cor(bd$mn_rna_rank, bd$mn_protein_rank, use = "complete.obs")
ggplot(bd, aes(x = mn_rna_rank, y = mn_protein_rank)) +
  geom_point(alpha = 0.4, size = 1) + geom_smooth(method = "lm", se = TRUE, color = "red") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = 0.5) +
  labs(x = "RNA rank in MN", y = "Protein rank in MN") +
  annotate("text", x = 0.2, y = 0.9, label = paste0("r = ", round(r_val, 2)), size = 4) +
  theme_bw()
dev.off()


known_markers <- c("NEFL","NEFH","NEFM","GFAP","UCHL1","CHIT1","CHI3L1","CHI3L2",
                   "TARDBP","SOD1","FUS","STMN2","PVALB","MMP9","KLK6")

marker_grounding <- grounding %>%
  filter(gene_symbol %in% known_markers) %>%
  dplyr::select(gene_symbol, uniprot, mn_rna_detected, ms_measured, mn_protein_pct_det,
                mn_protein_detected, mn_chromatin_open, grounding_score, grounding_tier) %>%
  arrange(gene_symbol)
cat("-- Known-marker grounding:\n"); print(as.data.frame(marker_grounding), row.names = FALSE); cat("\n")
write.csv(marker_grounding, "../output/summaries/fb_known_marker_grounding.csv", row.names = FALSE)

tier_dist <- grounding %>% dplyr::count(grounding_tier) %>% mutate(pct = round(n / sum(n) * 100, 1))
write.csv(tier_dist, "../output/summaries/fb_grounding_tiers.csv", row.names = FALSE)


write.csv(data.frame(method = detect_method,
                     n_detected = n_det,
                     n_eligible = n_elig,
                     n_union = nrow(ms_summary),
                     pct = round(100 * n_det / n_elig, 1)),
         "../output/summaries/fb_ms_detection.csv", row.names = FALSE)


# --- 8. Save --------------------------------------------

write.csv(grounding, "../output/grounding/grounding_map.csv", row.names = FALSE)
saveRDS(rna_summary, "../output/grounding/rna_summary.rds")
saveRDS(ms_summary, "../output/grounding/ms_summary.rds")
saveRDS(atac_summary, "../output/grounding/atac_summary.rds")

cat("\n== SCRIPT 01 COMPLETE ==\n")



