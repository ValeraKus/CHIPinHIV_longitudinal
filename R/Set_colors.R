library(paletteer)

cols <- paletteer_d("colorBlindness::paletteMartin")
# Gene colors

gene_colors=c("All genes"=cols[12], "Combined"=cols[12], "DNMT3A"=cols[3],
              "TET2"=cols[4], "ASXL1"=cols[6], "PPM1D"=cols[7], "TP53"=cols[5], "JAK2"=cols[8], "GNAS"=cols[11],
              "SF3B1"=cols[9], "SRSF2"=cols[10], "KRAS"=cols[1], "Other"=cols[13]) #fixed

hiv_colors=c("HIV"="red3","SHCS"="red3", "non_HIV"="grey44", "Fabre"="grey44", "Sard"="grey44", "Sardinia"="grey44")

cluster_colors=c("1"="#65CCFF", "cluster_1"="#65CCFF", "cl1"="#65CCFF", "2"="#FF7F00", "cluster_2"="#FF7F00", "cl2"="#FF7F00")
