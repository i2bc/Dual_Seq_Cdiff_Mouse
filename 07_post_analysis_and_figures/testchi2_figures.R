library(ggplot2)
library(reshape2)
library(dplyr)
library(ComplexHeatmap)
library(cowplot)
library(eulerr)
library(tidyverse)

contingencie_table <- function(updo, t){
# NA suppression for padj and l2FC:
if (updo == "up") {
   t_ud1 <- filter(t, padj1 < 0.05 & l2FC1 > 0.5)
   t_ud2 <- filter(t, padj2 < 0.05 & l2FC2 > 0.5)
   t_ud12 <- filter(t, padj1 < 0.05 & l2FC1 > 0.5 & padj2 < 0.05 & l2FC2 > 0.5)
} else { # down case
   t_ud1 <- filter(t, padj1 < 0.05 & l2FC1 < -0.5)
   t_ud2 <- filter(t, padj2 < 0.05 & l2FC2 < -0.5)
   t_ud12 <- filter(t, padj1 < 0.05 & l2FC1 < -0.5 & padj2 < 0.05 & l2FC2 < -0.5)
}
nb_ud12 <- dim(t_ud12)[1]
nb_ud1 <- dim(t_ud1)[1]-dim(t_ud12)[1]
nb_ud2 <- dim(t_ud2)[1]-dim(t_ud12)[1]
nb_noUD <- dim(t)[1]-nb_ud12-nb_ud1-nb_ud2
contingence <- as.table( rbind( c( nb_ud12, nb_ud1 ), c( nb_ud2, nb_noUD )))    
return(contingence)
}

euler_plot <- function(ec, updo, feature) {
   fit <- euler(c(A=ec[5,2],B=ec[3,2],"A&B"=ec[7,2],C=ec[2,2],"A&C"=ec[6,2],"B&C"=ec[4,2],"A&B&C"=ec[8,2]))
   e_plot <- plot(fit, 
       fills = c("lightcoral", "mediumslateblue", "palegreen3"), 
       quantities = TRUE, 
       labels = c("kreis","fletcher","pruss"), 
       edges = FALSE, 
       main=paste(updo, feature))
   ggsave(paste(feature, updo, "euler.png", sep="_"), e_plot)
}

run_chi2test <- function(kp4ct, feat, m1, m2) {
   colnames(kp4ct) <- c("Id1","mA1","mB1","l2FC1","padj1","Id2","mA2","mB2","l2FC2","padj2")
   for (ud in c("up","do")){
      ct <- contingencie_table(ud,kp4ct)
      dimnames(ct) <- list(c(paste(ud,m1,sep="_"),paste("no",ud,m1,sep="_")), c(paste(ud,m2,sep="_"),paste("no",ud,m2,sep="_")))
      test<-chisq.test(ct)
      print(paste(ud,feat,m1,"vs",m2,sep=" "))
      print(ct)
      print(paste("chisq.test p.value:",test$p.value))
   }
}

filenames <- "kreis_fletcher_pruss_mAmBl2FCpadj.txt" 
kfp <- read.csv(filenames, sep="\t", header=TRUE)
colnames(kfp) <- c("Idk","mAk","mBk","l2FCk","padjk","Idf","mAf","mBf","l2FCf","padjf","Idp","mAp","mBp","l2FCp","padjp")
# NAfisrt: NA suppression of Log2FC & padj :
kfp_noNA <- filter(kfp, !((is.na(padjk)|is.na(padjf)|is.na(padjp)|is.na(l2FCk)|is.na(l2FCf)|is.na(l2FCp))))

# comparison kreis - fletcher:
# NAfisrt: kf <- select(kfp_noNA, "Idk","mAk","mBk","l2FCk","padjk","Idf","mAf","mBf","l2FCf","padjf")
kf <- kfp %>% select("Idk","mAk","mBk","l2FCk","padjk","Idf","mAf","mBf","l2FCf","padjf") %>% filter(!((is.na(padjk)|is.na(padjf)|is.na(l2FCk)|is.na(l2FCf))))
kf_sRNA=kf[!grepl ("CD630_0|CD630_1|CD630_2|CD630_3",kf$Idk),]
run_chi2test(kf_sRNA, "sRNAs", "kreis", "fletcher")
kf_CDS=kf[grepl ("CD630_0|CD630_1|CD630_2|CD630_3",kf$Idk),]
run_chi2test(kf_CDS, "CDS", "kreis", "fletcher")
# run_chi2test(kf_noNA, "genes", "kreis", "fletcher")

# comparison kreis - pruss:
# NAfisrt: kp <- select(kfp_noNA, "Idk","mAk","mBk","l2FCk","padjk","Idp","mAp","mBp","l2FCp","padjp")
kp <- kfp %>% select("Idk","mAk","mBk","l2FCk","padjk","Idp","mAp","mBp","l2FCp","padjp") %>% filter(!((is.na(padjk)|is.na(padjp)|is.na(l2FCk)|is.na(l2FCp))))
kp_sRNA=kp[!grepl ("CD630_0|CD630_1|CD630_2|CD630_3",kp$Idk),]
run_chi2test(kp_sRNA, "sRNAs", "kreis", "pruss")
kp_CDS=kp[grepl ("CD630_0|CD630_1|CD630_2|CD630_3",kp$Idk),]
run_chi2test(kp_CDS, "CDS", "kreis", "pruss")
# run_chi2test(kp_noNA, "genes", "kreis", "pruss")

# comparison fletcher - pruss:
# NAfisrt: fp <- select(kfp_noNA, "Idf","mAf","mBf","l2FCf","padjf","Idp","mAp","mBp","l2FCp","padjp")
fp <- kfp %>% select("Idf","mAf","mBf","l2FCf","padjf","Idp","mAp","mBp","l2FCp","padjp") %>% filter(!((is.na(padjf)|is.na(padjp)|is.na(l2FCf)|is.na(l2FCp))))
fp_sRNA=fp[!grepl ("CD630_0|CD630_1|CD630_2|CD630_3",fp$Idf),]
run_chi2test(fp_sRNA, "sRNAs", "fletcher", "pruss")
fp_CDS=fp[grepl ("CD630_0|CD630_1|CD630_2|CD630_3",fp$Idf),]
run_chi2test(fp_CDS, "CDS", "fletcher", "pruss")
# run_chi2test(fp_noNA, "genes", "fletcher", "pruss")

# Euler (due to small quantity for sRNA, only for all genes) up & down:
ec_up <- kfp_noNA %>% 
   select("Idk", "l2FCk","padjk","l2FCf","padjf","l2FCp","padjp") %>% 
   mutate(k = ifelse(padjk < 0.05 & l2FCk > 0.5, 1, 0), f = ifelse(padjf < 0.05 & l2FCf > 0.5, 1, 0), p = ifelse(padjp < 0.05 & l2FCp > 0.5, 1, 0)) %>% 
   unite(eulercounts,c(k,f,p)) %>% 
   count(eulercounts)
euler_plot(ec_up, "up", "genes")
ec_do <- kfp_noNA %>% select("Idk", "l2FCk","padjk","l2FCf","padjf","l2FCp","padjp") %>% 
   mutate(k = ifelse(padjk < 0.05 & l2FCk < -0.5, 1, 0), f = ifelse(padjf < 0.05 & l2FCf < -0.5, 1, 0), p = ifelse(padjp < 0.05 & l2FCp < -0.5, 1, 0)) %>%   
   unite(eulercounts,c(k,f,p)) %>% 
   count(eulercounts)
euler_plot(ec_do, "down", "genes")

# heatmap
# select DE (Differentially Expressed) sRNA:
kfp_noNA_sRNA <- kfp_noNA[!grepl ("CD630_0|CD630_1|CD630_2|CD630_3",kfp_noNA$Idk),]
kfp_sRNA_DE <- filter(kfp_noNA_sRNA, (((padjk<0.05)&((l2FCk>0.5)|(l2FCk<-0.5)))|((padjf<0.05)&((l2FCf>0.5)|(l2FCf<-0.5)))|((padjp<0.05)&((l2FCp>0.5)|(l2FCp<-05)))))
rownames(kfp_sRNA_DE) <- kfp_sRNA_DE[,c("Idk")]

# prepare for heatmap graph: 
#    i) set log2FC to 0 when non DE 
#   ii) suppress sRNA when non DE in less than 2 experiments
kfp_4hm <- as.matrix(kfp_sRNA_DE %>% 
  mutate(l2FCkNA = ifelse(padjk > 0.05, 0.0, l2FCk)) %>%  mutate(l2FCfNA = ifelse(padjf > 0.05, 0.0, l2FCf)) %>%
  mutate(l2FCpNA = ifelse(padjp > 0.05, 0.0, l2FCp)) %>%  select(l2FCkNA, l2FCfNA, l2FCpNA) %>% 
  filter(!((l2FCkNA==0&l2FCfNA==0)|(l2FCkNA==0&l2FCpNA==0)|(l2FCfNA==0&l2FCpNA==0))))
# heatmap graph:
colnames(kfp_4hm) <- c("Kreis","Fletcher","Pruss")
png("kfp_sRNA_heatmap.png", width=300, height = 1000)
Heatmap(kfp_4hm, 
        name = "Log2FC", 
        column_names_rot = 45, 
        cluster_columns = FALSE, 
        show_row_dend = FALSE,
        row_names_gp = gpar(fontsize = 7))
dev.off()

# select 25 best DE (Differentially Expressed) genes:
kfp_noNA_CDS <- kfp_noNA[grepl ("CD630_0|CD630_1|CD630_2|CD630_3",kfp_noNA$Idk),]
# prepare for heatmap graph: 
#    1) select all DE from each experiment
#    2) set log2FC to 0 when non DE 
#    3) suppress CDS when non DE in less than 2 experiments
kfp_CDS_DE <- kfp_noNA_CDS %>% 
  filter(( (padjk<0.05)&((l2FCk>0.5)|(l2FCk<-0.5)) )|( (padjf<0.05)&((l2FCf>0.5)|(l2FCf<-0.5)) )|( (padjp<0.05)&((l2FCp>0.5)|(l2FCp<-05)) )) %>%
  mutate(meanFCkfp = (l2FCk+l2FCf+l2FCp)/3) %>% 
  mutate(l2FCkNA = ifelse(padjk > 0.05, 0.0, l2FCk)) %>% 
  mutate(l2FCfNA = ifelse(padjf > 0.05, 0.0, l2FCf)) %>%
  mutate(l2FCpNA = ifelse(padjp > 0.05, 0.0, l2FCp)) %>%
  filter(!((l2FCkNA==0&l2FCfNA==0)|(l2FCkNA==0&l2FCpNA==0)|(l2FCfNA==0&l2FCpNA==0)))
rownames(kfp_CDS_DE) <- kfp_CDS_DE[,c("Idk")]
#    4) sort by the log2FC mean of the 3 experiments and keep only n genes up and down:
# kfp_CDS_DE_4hm <- as.matrix(
#  union(kfp_CDS_DE %>% slice_min(meanFCkfp, n = 50), kfp_CDS_DE %>% slice_max(meanFCkfp, n = 50)) %>% 
#  select(l2FCkNA, l2FCfNA, l2FCpNA))
#    4) keep all DE CDS:
kfp_CDS_DE_4hm <- as.matrix(kfp_CDS_DE %>% select(l2FCkNA, l2FCfNA, l2FCpNA))
# heatmap graph:
colnames(kfp_CDS_DE_4hm) <- c("Kreis","Fletcher","Pruss")
png("kfp_CDS_heatmap.png", width=300, height = 1000)
Heatmap(kfp_CDS_DE_4hm, 
        name = "Log2FC", 
        column_names_rot = 45, 
        cluster_columns = FALSE, 
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        row_names_gp = gpar(fontsize = 7))
dev.off()
