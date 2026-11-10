# Runs the whole analysis in order, each step in a fresh R session.
# From the repository root:   Rscript run_all.R        (or Rscript run_all.R 6 to start at step 6)
# Needs the input data in resources/ (see resources/README.md), the chip_lod code in external/chip_lod (see README.md)
# and pandoc for rmarkdown (in RStudio's installation, or set RSTUDIO_PANDOC).

steps <- c(
  "scripts/merge_mutect2_output.R",                                       # resources/chip_calls.csv
  "analysis/1_variant_calls/01_per_sample_LoD_and_QC.Rmd",                # sample QC, bam-readcount regions (Figure S2 needs later steps: skipped on a first run)
  "analysis/1_variant_calls/02_SHCS_CHIP_trajectories.Rmd",               # SHCS trajectories with the LoD, Table S1 data
  "analysis/1_variant_calls/03_nonCHIP_variants_QC.Rmd",                  # non-CHIP variants as controls
  "analysis/1_variant_calls/04_SardiNIA_CHIP_trajectories.Rmd",           # SardiNIA trajectories with the same LoD
  "analysis/2_trajectory_models/01_spline_fits.Rmd",                      # REML splines of both cohorts
  "analysis/2_trajectory_models/02_ages_at_LoD2_and_2pct.Rmd",            # extrapolated curves, ages at LoD/2 and 2% VAF
  "analysis/2_trajectory_models/03_consensus_clusters.Rmd",               # consensus DTW clusters, Table S6
  "analysis/3_clinical_data/01_SHCS_metadata.Rmd",                        # demographics, HIV and ART history
  "analysis/3_clinical_data/02_SHCS_blood_counts.Rmd",                    # interpolated blood counts
  "analysis/4_sensitivity_analyses/01_spline_sensitivity_and_simulation.Rmd",  # Figures S3-S5, Tables S2-S4
  "analysis/5_cohort_comparison/01_SHCS_SardiNIA_matching.Rmd",           # matched and unmatched sets
  "analysis/5_cohort_comparison/02_cohort_characteristics_sex_smoking.Rmd",
  "analysis/6_figures_and_tables/Figure1_FigureS1.Rmd",
  "analysis/6_figures_and_tables/Figure2_FiguresS7_S8.Rmd",
  "analysis/6_figures_and_tables/Figure3_FiguresS9_S10.Rmd",
  "analysis/6_figures_and_tables/Figure4_FiguresS11_S18.Rmd",
  "analysis/6_figures_and_tables/FigureS6_cooccurring.Rmd",
  "analysis/1_variant_calls/01_per_sample_LoD_and_QC.Rmd",                # again, for Figure S2 (section 8)
  "analysis/6_figures_and_tables/Tables.Rmd"                              # Table 1, Tables S1-S12
)

stopifnot(file.exists("CHIPinHIV_longitudinal.Rproj"), file.exists("external/chip_lod/r/R/lod.R"))
for (d in c("results", "figures", "output")) dir.create(d, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
first <- if (length(args)) as.integer(args[1]) else 1L

for (i in seq(first, length(steps))) {
  f <- steps[i]
  message(sprintf("[%d/%d] %s", i, length(steps), f))
  cmd <- if (grepl("\\.R$", f)) shQuote(f) else c("-e", shQuote(sprintf("rmarkdown::render('%s', quiet = TRUE)", f)))
  status <- system2(file.path(R.home("bin"), "Rscript"), cmd)
  if (status != 0) stop("step ", i, " (", f, ") failed")
}
