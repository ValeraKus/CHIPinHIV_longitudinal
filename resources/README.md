# Input data

The data are not part of this repository. The SHCS data (sequencing results and clinical database) are available to
researchers upon project submission to and approval by the Scientific Board of the Swiss HIV Cohort Study; the
SardiNIA read counts and clone fitness estimates are published with Fabre et al., *Nature* 606, 335–342 (2022).

The notebooks expect these files under `resources/`:

| Path | Content |
|---|---|
| `mutect2_Calls/mutect2_output/*.vcf.gz` | Mutect2 VCF of every sample of the longitudinal targeted sequencing |
| `CHIP_calls_longitudinal/CHIP_calls_all_long.csv` | Annotated CHIP calls of the longitudinal sequencing |
| `CHIP_screening/CHIP_candidate_variants.csv` | CHIP variant definition (curation criteria) |
| `CHIP_screening/SHCS_Chip_all_vars.txt`, `CHIP_screening/SHCS_Chip_per_person.txt` | CHIP calls of the earlier cross-sectional screening |
| `chip_vars_for_bam_readcount.txt`, `nonchip_vars_for_bam_readcount.txt` | Region lists (chromosome, start, end) at which bam-readcount was run: all CHIP and all non-CHIP variants called by Mutect2 |
| `All.chip.bam_readcount.txt`, `All.nonchip.bam_readcount.txt` | bam-readcount pile-up of every sample at these regions |
| `samples_with_dates.csv` | Sample manifest: sample, participant and collection date |
| `SHCS_change_ids.txt` | SHCS participant id and the participant id used in the paper (`SHCS_1`, ...) |
| `excluded_individuals.txt` | SHCS ids of the two sequenced individuals without a CHIP call, one per line |
| `pileup_sample_names.txt` | The sample that is twice in the pile-up: tab-separated `pileup_name`, `action` (`use` / `drop`), `sample` |
| `SHCS_database/2024-12/pat.csv` | SHCS database export (December 2024) |
| `SHCS_database/2025-09/{pat,lab,med_drug,med_drug_code}.csv` | SHCS database export (September 2025) |
| `SHCS_database/2025-12/admin.csv` | SHCS database export (December 2025) |
| `spVL/spVL_calc.txt`, `spVL/spVL_all_patients.txt` | Set-point viral loads (participant metadata only) |
| `SardiNIA/ALLvariants_exclSynonymous_Xadj.txt` | SardiNIA read counts per variant and phase (Fabre et al. 2022) |
| `SardiNIA/smoke.txt` | Smoking status of the SardiNIA participants (Fabre et al. 2022) |
| `SardiNIA/Sardinia_growth_coefficients_and_age_at_onset_estimates_for_all_clones.csv` | Per-clone fitness estimates of Fabre et al. 2022 |

The notebooks also write their intermediate tables here (for example `chip_calls_long_filtered.txt`,
`SHCS_Metadata_univ.txt`, `SHCS_CHIP_blood_trajectories_merged.txt` and the curves in `LoD/`).
