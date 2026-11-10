# Merges the per-sample Mutect2 VCFs of the longitudinal sequencing (resources/mutect2_Calls/mutect2_output/*.vcf.gz)
# with the annotated CHIP variants of the longitudinal assay (matched on position +/- 1 bp) and writes
# resources/chip_calls.csv, the input of analysis/1_variant_calls/02_SHCS_CHIP_trajectories.Rmd, which gives each
# record the annotation of its own allele.
# Run from the repository root: Rscript scripts/merge_mutect2_output.R
setwd("resources/mutect2_Calls/")
library(tidyr)
library(tidyverse)

files=list.files("./mutect2_output/", pattern = "vcf.gz$")

chip_calls_long <- read.csv("../CHIP_calls_longitudinal/CHIP_calls_all_long.csv", fileEncoding = "UTF-8-BOM") %>% 
  select(Gene.refGene, Chr, Start, End, Ref, Alt, Func.refGene,GeneDetail.refGene,
         ExonicFunc.refGene, AAChange.refGene,   cosmic70 , Accession,transcriptOI,NonsynOI) %>% distinct()

pot_vars <- rbind(chip_calls_long, chip_calls_long) %>% distinct()

result_vcf <- c()

for (i in files){

file=paste0("./mutect2_output/", i)
sample=gsub(".vcf.gz|./mutect2_output/", "", file)
vcf <- read.table(gzfile(file), header = F, sep = "\t")
vcf$pos2 <- vcf$V2 + 1
vcf$pos3 <- vcf$V2 - 1


vcf_merged <- merge(vcf, pot_vars, by.x = c("V1", "V2"), by.y = c("Chr", "Start"))
vcf_merged2 <- merge(vcf, pot_vars, by.x = c("V1", "pos2"), by.y = c("Chr", "Start"))
vcf_merged3 <- merge(vcf, pot_vars, by.x = c("V1", "pos3"), by.y = c("Chr", "Start"))

vcf_merged <- rbind(vcf_merged, vcf_merged2, vcf_merged3)
if (nrow(vcf_merged) > 0){
  

vcf_merged <- vcf_merged[,c("Gene.refGene","V1", "V2", "End", "V4", "V5", "Func.refGene","GeneDetail.refGene",
                            "ExonicFunc.refGene", "AAChange.refGene",   "cosmic70" , "V7", "V8", "V9", "V10",
                            "Accession","transcriptOI","NonsynOI")]

vcf_merged <- vcf_merged[!duplicated(vcf_merged),]
vcf_merged$Sample=sample

vcf_merged <- vcf_merged %>%
  separate(col = V10, str_split_fixed("GT:AD:AF:DP:F1R2:F2R1:FAD:SB", ":", 8), sep = ":", remove = T, fill = "right")

result_vcf <- rbind(result_vcf, vcf_merged)
}
}

write.csv(result_vcf, "../chip_calls.csv")
