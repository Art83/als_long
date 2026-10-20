# Plasma proteomic signal in ALS does not reflect what iPSC-derived motor neurons express

Analysis code for the Answer ALS plasma proteomics manuscript. 

**The question** is whether the plasma proteins that track clinical decline in ALS are produced in iPSC-dervied motor neurons. **Specifically, whether a motor neuron model can reproduce what is measured in the plasma of patients**

The code: 
1. Builds a gene-level map of motor neuron expression for every SomaScan 7K plasma target from iPSC-derived motor neuron RNA-seq, mass spectrometry and ATAC-seq of the same resource.
2. Tests whether that map predicts which plasma proteins track ALSFRS-R decline.
3. Characterises the route to plasma and tissue distribution of the proteins that do.
4. Benchmarks the SomaScan NEFL aptamer against a serum NfL immunoassay.
5. Validate the findings in Target ALS datasets

**The answer** is that iPSC-derived motor neurons, on their own, do not model ALS progression as it is measurable in patients, and that any attempt to connect a motor neuron phenotype to the clinical trajectory will have to go through the periphery. If the clinically meaningful signal (that tracks decline in patients) is made in periphery and reaches blood by secretion and shedding, then a motor neuron culture is disconnected from that signal by construction.

Specifically:
1. Motor neuron expression explains 0.16% of the variance in how well a plasma protein tracks decline.
2. The proteins that cross the significance threshold are depleted of motor neuron proteins. 
3. The progression markers are secreted and membrane proteins with no explicit tissue signature. 
4. The cross-sectional signature is muscle-weighted and reflects severity rather than disease-specific signature. 
5. NfL is prognostic at baseline but its within-person change does not track the rate of decline.

## Layout for reproducibility

```
als/
  scripts/        even scripts, run from this folder
  raw/            Answer ALS release files (not in the repository due to user agreement, see Data)
  proc/           meta.csv, participant metadata restricted to ALS and controls
  output/         everything the scripts write (created on first run)
    grounding/    grounding map, per-aptamer progression results, rds summaries
    landscape/    annotated hit tables, triage tables, tissue profiles, plots
    nfl_check/    NfL plots
    summaries/    small CSVs with the numbers quoted in the manuscript
    tables/       threshold sensitivity tables
    replication/  result of external validation in Target ALS cohort
    [X]_console.txt  console log of each script to facilitate the assessment without data files
```

All paths in the scripts are relative to `scripts/`.

## Data

Raw files come from the Answer ALS data portal (https://dataportal.answerals.org) or Target ALS portal (https://targetals.org). To facilitate the exact files tracking, the original names were kept. 
The scripts expect the following under `raw/`:

| Path                                                 | File                                                                                                                                                                            | Used by    |
|------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|------------|
| `somascan/normalized/csv/`                           | `SS-2342615_v4.1_anmlSMP__AALS-v1.0-participant-samples.csv`                                                                                                                    | 02, 03, 04 |
| `somascan/documents/`                                | `annotation.csv` (SomaScan 7K v4.1 annotation, SeqId to UniProt)                                                                                                                | 01 to 04   |
| `als-metadata/answerals/`                            | `aals_dataportal_datatable.csv`                                                                                                                                                 | 02, 03, 04 |
| `als-metadata/answerals/clinical/`                   | `ALSFRS_R.csv`, `Plasma_Sample.csv`                                                                                                                                             | 02, 03, 04 |
| `transcriptomics/4_matrix/`                          | `AnswerALS-290-T-v1-release5_rna-counts-matrix.csv`                                                                                                                             | 01         |
| `transcriptomics/4_matrix/Sample Mapping Information/` | `Sample Mapping File Feb 2024.csv`                                                                                                                                              | 01         |
| `proteomics/4_matrix/`                               | `AnswerALS-447-P_proteomics-protein-matrix_prebatchcorrected.txt` (detection calls) and `AnswerALS-436-P_proteomics-protein-matrix_correctedImputed.txt` (abundance ranks only) | 01         |
| `epigenomics/epigenomics/4_matrix/`                  | `AnswerALS-620-E-v1-release6_DiffBind-raw-counts-minOverlap-0.1.csv`                                                                                                            | 01         |
| `NfL/`                                               | `NfL_Seimens_Atellica_Immunoassay_Proteomics_2022_07.csv` (serum NfL, pg/mL)                                                                                                    | 04         |
| `raw/`                                               | `interactions.tsv` (DGIdb interactions export, https://dgidb.org/downloads)                                                                                                     | 03         |
| `raw/`                                               | `proteinatlas.tsv` (HPA export, https://www.proteinatlas.org/download)                                                                                                                    | 02b        |
| `target_als/ms`                                      | `lc_ms_proteomics_tmt_plasma.csv`                                                                                                                                               | 05         | 
| `target_als/`                                        | `nhs_clinical_data.csv`                                                                                                                                                         | 05         |
| `target_als/`                                        | `olink_matrix.csv`                                                                                                                                                              | 05         |

`proc/meta.csv` is the data portal participant table restricted to `SUBJECT_GROUP` in ALS and Healthy Control, with `Participant_ID` matching the sample mapping file.

## Software

R 4.3 or later with dplyr, tidyr, ggplot2, gridExtra, lme4, lmerTest, limma, DESeq2, edgeR, GenomicRanges, ChIPseeker, TxDb.Hsapiens.UCSC.hg38.knownGene, org.Hs.eg.db, AnnotationDbi and parallel.

Script 03 needs concordR 0.2.0 or later (gate columns from `triage()`, `audit_call()`, `funnel_class()`, `query_presence()`, `load_atlas()`), which is a separate package.

## Running

From `scripts/`, in order, either in powershell or bash terminal. In Win some scripts run slower. Each script reads the outputs of the ones before it.

```
Rscript 01_grounding_map.R            > ../output/01_console.txt 2>&1
Rscript 02_progression.R              > ../output/02_console.txt 2>&1
Rscript 02b_external_abundance.R      > ../output/02b_console.txt 2>&1
Rscript 03_landscape_tissue_origins.R > ../output/03_console.txt 2>&1
Rscript 04_nfl_design_check.R         > ../output/04_console.txt 2>&1
Rscript 05_target_als_replication.R   > ../output/05_console.txt 2>&1
```

`00_visit_time.R` is sourced by 02 and 04 and is not run on its own. Script 02 fits a random slope model per aptamer and runs on all available cores minus one (socket cluster on Windows, fork elsewhere), and takes the longest. The permutation test in 02 is seeded (`set.seed(42)`).

## Scripts

**00_visit_time.R.** Draw timing. The SomaScan `TimePoint` column is the study visit number, not the rank of a participant's ALSFRS-R assessments. `soma_visit_time()` builds a GUID by visit table from the plasma collection date in `Plasma_Sample.csv`, with the same-visit ALSFRS-R date as fallback, and `attach_soma_time()` joins it to the SomaScan rows.

**01_grounding_map.R.** The grounding map. Every SomaScan aptamer with a UniProt accession is mapped to a gene symbol and Ensembl ID through org.Hs.eg.db. Transcript detection comes from DESeq2 variance-stabilised RNA-seq of the motor neuron lines (five or more counts in at least 20% of lines and mean VST above the 10th percentile). Protein detection comes from the pre-batch correction MS matrix. The imputed matrix supplies abundance ranks only. Promoter accessibility comes from ATAC peaks within 2 kb of a TSS (ChIPseeker, hg38), detected at more than five counts in at least 20% of lines. Threshold sensitivity for both the protein and transcript rules, including edgeR filterByExpr, is written to `output/summaries/` and `output/tables/`.

**02_progression.R.** Does grounding predict progression tracking? ALSFRS-R trajectories are modelled with a random intercept and slope per patient and each patient's decline rate is the BLUP slope. For each aptamer, per-patient slopes of log2 RFU on time are estimated by the same model in patients with two or more dated draws, and the protein slope is regressed on the ALSFRS-R slope. The absolute correlation is then regressed on grounding score, naive and adjusted for baseline signal and its CV, with a permutation null, and the transcript-positive aptamers are split by protein confirmation.

**02b_external_abundance.R.** Plasma abundance re-adjustment using external dataset. Supp script for script 02. 

**03_landscape_tissue_origins.R.** What the hits are? Cross-sectional differential abundance at baseline (limma, sex and age adjusted) and within-ALS association of each DAP with baseline ALSFRS-R and disease duration. concordR gates on the progression and differential sets against the whole panel as background, route to plasma, tissue specificity, concordance and audit call, presence votes in named tissues, and a 59-tissue presence profile summarised as normalised entropy and the probability of the top tissue. Grounding tier overlay with a Fisher test of hits by tier against the panel, counted as unique proteins.

**04_nfl_design_check.R.** NfL. A pre-specified panel is run through the slope design and the conventional baseline-prognostic design. The serum NfL immunoassay is then benchmarked against the aptamer at matched draws, and run through both designs itself, with tertiles of baseline NfL against ALSFRS-R slope.

**05_target_als_replication.** Does the findings replicate in external cohort? Cross-sectional ALS vs control at baseline on each platform and ALSFRS-R BLUP slopes, per-protein BLUP slopes on each platform.

## Note on annotations

The SomaScan matrix is keyed by SeqId in every script. UniProt is carried as a column and is never used as a column or row name, because 1,125 of the 7,550 aptamers share an accession with another aptamer and R's name repair silently mangles duplicates (Q13508 becomes Q13508.1, which then fails every join). Counts in the manuscript are stated as aptamers or as unique proteins and the scripts print both.

## Tier labels in the code map to the manuscript as follows.

| Code | Manuscript | Definition |
|---|---|---|
| `confirmed_full` | MN-confirmed | transcript, protein and accessible promoter |
| `confirmed_partial` | MN-confirmed-partial | transcript and protein, no annotated accessible promoter |
| `rna_only_ms_gap` | MN-transcript | transcript, protein not detected or not measured |
| `ms_only_rna_gap` | MN-protein | protein without transcript (10 aptamers) |
| `no_mn_expression` | MN-absent | neither |

The MN-transcript tier records the limit of MS coverage rather than biology. Only 30% of panel proteins were quantified in the motor neuron proteome, and the scripts report the split between quantified below threshold and never quantified. Transcript absence is the robust negative call.

The slope design uses the 229 patients with two or more dated SomaScan blood draws and an ALSFRS-R slope. The baseline design in 04 uses every patient with a baseline draw and a slope (242 for the aptamers, 608 for the immunoassay). Spinal cord is not in the HPA tissue vocabulary and cannot be tested as a context tissue. The concordR concordance gate returns `caution` or `not_in_atlas` for almost every secreted protein and is reported but not interpreted.

Files in `output/summaries/` prefixed `fb_` hold the numbers quoted in the Results and are small enough to read without the raw data.
