# CHIP trajectories in people with HIV

Analysis code for:

> Timonina V, van der Graaf A, Calmy A, Cavassini M, Capoferri G, Günthard HF, Walti LN, Schmid P, Tarr PE, Bick AG, Fellay J and the
> Swiss HIV Cohort Study. **Longitudinal modelling of clonal hematopoiesis reveals altered early clonal dynamics in
> people with HIV.** *medRxiv* 2026.04.08.26350407 (2026). <https://doi.org/10.64898/2026.04.08.26350407>

If you use this code, please cite the paper ([medRxiv](https://www.medrxiv.org/content/10.64898/2026.04.08.26350407v2)):

```bibtex
@article{Timonina2026CHIPtrajectories,
  title   = {Longitudinal modelling of clonal hematopoiesis reveals altered early clonal dynamics in people with HIV},
  author  = {Timonina, Valeriia and van der Graaf, Adriaan and Calmy, Alexandra and Cavassini, Matthias and Capoferri, Gioele and
             G{\"u}nthard, Huldrych F. and Walti, Laura N. and Schmid, Patrick and Tarr, Philip E. and
             Bick, Alexander G. and Fellay, Jacques and {Swiss HIV Cohort Study}},
  journal = {medRxiv},
  year    = {2026},
  doi     = {https://doi.org/10.64898/2026.04.08.26350407}
}
```

The R notebooks reconstruct the variant allele frequency (VAF) trajectories of CHIP clones from serial targeted
sequencing of people with HIV in the Swiss HIV Cohort Study (SHCS), with a per-sample limit of detection (LoD). They
fit a penalised regression spline (REML) to each trajectory and derive the VAF speed and acceleration. The trajectories
are clustered with consensus dynamic time warping, related to blood parameters, and compared with the SardiNIA general
population cohort (Fabre et al., *Nature* 2022).

## Structure

```
R/                          shared functions: spline fits and trajectory statistics, figure helpers, gene colours
scripts/                    merge of the per-sample Mutect2 VCFs
analysis/
  1_variant_calls/          limit of detection, sample QC, SHCS and SardiNIA CHIP trajectories
  2_trajectory_models/      REML splines of both cohorts, ages at LoD/2 and 2% VAF, consensus clusters
  3_clinical_data/          SHCS participant metadata and blood counts
  4_sensitivity_analyses/   spline settings, read-count bootstrap, simulation study
  5_cohort_comparison/      matching of SHCS and SardiNIA trajectories, cohort characteristics, sex and smoking
  6_figures_and_tables/     main and supplementary figures, Table 1 and Tables S1-S12
assets/                     the two drawn schemes (Figure 1a, Figure 2g)
resources/                  input data (not included, see resources/README.md) and intermediate tables
run_all.R                   runs all steps in order
```

The notebooks write intermediate tables to `resources/` and `results/`, working figures to `figures/`, and the
figures and tables of the paper to `output/`:

```
output/Figure1 ... Figure4/           Figure1.png|svg ... Figure4.png|svg, and their panels in panels/
output/Supplementary_Figures/         FigS1_... to FigS18_... (png and svg), panels in panels/
output/Tables/                        Table1.txt, Supplementary_Tables.xlsx (Tables S1-S12, one sheet per table)
figures/paper/additional/             figures of additional analyses that are not display items of the paper
results/paper/                        the numbers behind the figures, tables and text
```

| Paper item | Notebook (in `analysis/`) |
|---|---|
| Figure 1, Figure S1 | `6_figures_and_tables/Figure1_FigureS1.Rmd` |
| Figure 2, Figures S7-S8 | `6_figures_and_tables/Figure2_FiguresS7_S8.Rmd` |
| Figure 3, Figures S9-S10 | `6_figures_and_tables/Figure3_FiguresS9_S10.Rmd` |
| Figure 4, Figures S11-S18 | `6_figures_and_tables/Figure4_FiguresS11_S18.Rmd` |
| Figure S2 | `1_variant_calls/01_per_sample_LoD_and_QC.Rmd` (section 8; panels a-d and the assembled figure) |
| Figures S3-S5, Tables S2-S4 | `4_sensitivity_analyses/01_spline_sensitivity_and_simulation.Rmd` |
| Figure S6, Table S5 | `6_figures_and_tables/FigureS6_cooccurring.Rmd` |
| Table S6 | `2_trajectory_models/03_consensus_clusters.Rmd` |
| Tables S7, S8 | `6_figures_and_tables/Tables.Rmd` (computed there) |
| Tables S9-S12 | the Figure 3 and Figure 4 notebooks and `5_cohort_comparison/02_cohort_characteristics_sex_smoking.Rmd` |
| Table 1, Tables S1-S12 (workbook) | `6_figures_and_tables/Tables.Rmd` |

## Data

The SHCS data are available to researchers upon project submission to and approval by the SHCS Scientific Board. The
SardiNIA data are published with Fabre et al. (2022). `resources/README.md` lists the input files the notebooks expect.
Participants are named by the ids used in the paper (`SHCS_1`, `SHCS_2`, ...); the link to SHCS ids is part of the
data and is not distributed.

The intermediate tables in `resources/` and `results/` and the rendered notebooks (`*.html`) hold participant-level
data with the SHCS ids. They are excluded by `.gitignore` and must not be committed or shared; the figures and tables in
`output/` use the paper ids only.

## Running

R 4.3 or later with the packages tidyverse, mgcv, gratia, lme4, lmerTest, broom, broom.mixed, dtwclust, proxy, clue,
cluster, patchwork, gridExtra, ggpubr, ggrepel, ggExtra, scales, Hmisc, vcfR, paletteer, png, svglite, writexl, knitr,
rmarkdown and here, and pandoc for rmarkdown (part of RStudio; outside RStudio set `RSTUDIO_PANDOC`):

```r
install.packages(c("tidyverse", "mgcv", "gratia", "lme4", "lmerTest", "broom", "broom.mixed", "dtwclust", "proxy", "clue",
                   "cluster", "patchwork", "gridExtra", "ggpubr", "ggrepel", "ggExtra", "scales", "Hmisc", "vcfR",
                   "paletteer", "png", "svglite", "writexl", "knitr", "rmarkdown", "here"))
```

The code was run with R 4.5.0, tidyverse 2.0.0 (ggplot2 4.0.3, dplyr 1.2.1), mgcv 1.9.4, lme4 2.0.6, lmerTest 3.2.1,
dtwclust 6.0.0 and patchwork 1.3.2; every notebook ends with its `sessionInfo()`. The limit of detection uses
[chip_lod](https://github.com/adriaan-vd-graaf/chip_lod) at commit `91ebf7a`:

```bash
git clone https://github.com/adriaan-vd-graaf/chip_lod external/chip_lod && git -C external/chip_lod checkout 91ebf7a
```

Then, from the repository root:

```bash
Rscript run_all.R
```

The whole run takes about half an hour on a laptop. `run_all.R` knits `1_variant_calls/01_per_sample_LoD_and_QC.Rmd`
twice: its last section (Figure S2) needs the trajectory tables of both cohorts, which later steps write. To knit a
single notebook, open `CHIPinHIV_longitudinal.Rproj` in RStudio: every notebook uses the repository root as its working
directory.

Alignment, Mutect2 variant calling, annotation and curation of the CHIP variants (with the filters given in the
Methods of the paper) and the bam-readcount pile-ups at the regions of all Mutect2 calls were run before these
notebooks; their outputs are inputs here (`resources/README.md`). A few sections compare with tables of earlier
versions of the analysis and run only when those files are present.
