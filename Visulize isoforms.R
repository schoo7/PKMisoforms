library(ggplot2)
library(reshape2)
library(preprocessCore)
setwd("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/PKM_bioinformatics")
# study the PKM isoforms
# start from human PCa cell lines
cell=read.csv("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/Siyuan/Raw PCa cell lines study/siyuan analysis/Combined_transcript_tpm.csv")
tx=read.table("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/Siyuan/Raw PCa cell lines study/siyuan analysis/salmon_tx2gene.tsv",
              sep = '\t', header = F)
tx=as.character(subset(tx,V3=="PKM")$V1)
rownames(cell)=cell$X
cell=cell[,-1]
cell=log2(cell+1)
cell=normalize.quantiles(as.matrix(cell),keep.names = T)
cell=cell[tx,]
cell=as.data.frame(t(cell))
rownames(cell)=gsub("X22RV1","22RV1",rownames(cell))
cell$group=gsub("_.*","",rownames(cell))
table(cell$group)
cell$PKM1=cell$NM_182470.2+cell$NM_182471.2
cell$PKM2=cell$NM_002654.4
cell$PKM1_like=cell$NM_182470.2+cell$NM_182471.2+cell$NM_001206796.1+cell$NM_001206799.1
cell$PKM2_like=cell$NM_002654.4+cell$NM_001206797.1+cell$NM_001206798.1+cell$XM_005254443.1+cell$XM_005254445.2+cell$XM_006720570.1
cell=melt(cell,id="group")
colnames(cell)=c("Cell_line","Transcript","log2TPM")
cell$Cell_line=factor(cell$Cell_line,levels = c("VCaP","LNCaP","C42","22RV1","PC3","DU145","H660"))
cell$Transcript=factor(cell$Transcript,levels = c("PKM1","PKM2","PKM1_like","PKM2_like","NM_182470.2","NM_182471.2","NM_001206796.1","NM_001206799.1","NM_002654.4",
                                                  "NM_001206797.1","NM_001206798.1","XM_005254443.1","XM_005254445.2","XM_006720570.1"))
ggplot(cell[grep("PKM",cell$Transcript),],aes(x=Cell_line,y=log2TPM,group=Cell_line,color=Cell_line))+geom_boxplot(outlier.colour = NA)+geom_jitter(alpha=0.3)+
  facet_grid(cell[grep("PKM",cell$Transcript),]$Transcript)+
  theme_classic()+theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black",face = "bold", size = 12),
                        axis.text = element_text(color = "black",face = "bold", size = 12),
                        text = element_text(color = "black",face = "bold", size = 12))
ggsave("Cell_PKM1_2.tiff",width = 10,height = 14,units = "cm",dpi=600,compression = "lzw")

ggplot(cell[-grep("PKM",cell$Transcript),],aes(x=Cell_line,y=log2TPM,group=Cell_line,color=Cell_line))+geom_boxplot(outlier.colour = NA)+geom_jitter(alpha=0.3)+
  facet_grid(cell[-grep("PKM",cell$Transcript),]$Transcript)+
  theme_classic()+theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black",face = "bold", size = 12),
                        axis.text = element_text(color = "black",face = "bold", size = 12),
                        text = element_text(color = "black",face = "bold", size = 12))
ggsave("Cell_PKM isoforms.tiff",width = 10,height = 35,units = "cm",dpi=600,compression = "lzw")      


# now study PKM isoforms in human samples
patient=read.csv("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/PCa_patient_atlas/PCa_patients_atlas_transcript_tpm.csv")
rownames(patient)=patient$X
patient=patient[,-c(1:2)]
patient=log2(patient+1)
patient=normalize.quantiles(as.matrix(patient),keep.names = T)
patient=patient[tx,]
patient=as.data.frame(t(patient))
patient_group=read.csv("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/HuPSA/HuPSA_integration/Assigned_sample_groups.csv")
patient_group=patient_group[,c(1,22)]
colnames(patient_group)=c("ID","class")
patient$ID=rownames(patient)
patient=merge(patient,patient_group,by="ID")
patient=subset(patient,patient$class!="contaminate")
patient=subset(patient,patient$class!="KRT7_AdPCa_mix")
patient=subset(patient,patient$class!="NEPCa_AdPCa_mix")
patient=subset(patient,patient$class!="Progenitor_like_AdPCa_mix")
patient=subset(patient,patient$class!="unclassified")
patient=subset(patient,patient$class!="Progenitor_like")
patient=subset(patient,patient$class!="KRT7")
patient$class=gsub("_.*","",patient$class)
rownames(patient)=patient$ID
patient=patient[,-1]


patient$PKM1=patient$NM_182470.2+patient$NM_182471.2
patient$PKM2=patient$NM_002654.4
patient$PKM1_like=patient$NM_182470.2+patient$NM_182471.2+patient$NM_001206796.1+patient$NM_001206799.1
patient$PKM2_like=patient$NM_002654.4+patient$NM_001206797.1+patient$NM_001206798.1+patient$XM_005254443.1+patient$XM_005254445.2+patient$XM_006720570.1
patient=melt(patient,id="class")
colnames(patient)=c("group","Transcript","log2TPM")
patient$group=factor(patient$group,levels = c("Normal","AdPCa","NEPCa"))
patient$Transcript=factor(patient$Transcript,levels = c("PKM1","PKM2","PKM1_like","PKM2_like","NM_182470.2","NM_182471.2","NM_001206796.1","NM_001206799.1","NM_002654.4",
                                                  "NM_001206797.1","NM_001206798.1","XM_005254443.1","XM_005254445.2","XM_006720570.1"))
ggplot(patient[grep("PKM",patient$Transcript),],aes(x=group,y=log2TPM,group=group,color=group))+geom_boxplot(outlier.colour = NA)+geom_jitter(alpha=0.2)+
  facet_grid(patient[grep("PKM",patient$Transcript),]$Transcript)+
  theme_classic()+theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black",face = "bold", size = 12),
                        axis.text = element_text(color = "black",face = "bold", size = 12),
                        text = element_text(color = "black",face = "bold", size = 12))
ggsave("Patient_PKM1_2.tiff",width = 8,height = 14,units = "cm",dpi=600,compression = "lzw")

ggplot(patient[-grep("PKM",patient$Transcript),],aes(x=group,y=log2TPM,group=group,color=group))+geom_boxplot(outlier.colour = NA)+geom_jitter(alpha=0.2)+
  facet_grid(patient[-grep("PKM",patient$Transcript),]$Transcript)+
  theme_classic()+theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black",face = "bold", size = 12),
                        axis.text = element_text(color = "black",face = "bold", size = 12),
                        text = element_text(color = "black",face = "bold", size = 12))
ggsave("Patient_PKM isoforms.tiff",width = 8,height = 35,units = "cm",dpi=600,compression = "lzw")   
# calculate P-values
# using T-test
pkm1_value=subset(patient,Transcript=="PKM1")
t.test(pkm1_value[which(pkm1_value$group=="NEPCa"),"log2TPM"],pkm1_value[which(pkm1_value$group=="AdPCa"),"log2TPM"])
pkm1_like_value=subset(patient,Transcript=="PKM1_like")
t.test(pkm1_like_value[which(pkm1_like_value$group=="NEPCa"),"log2TPM"],pkm1_like_value[which(pkm1_like_value$group=="AdPCa"),"log2TPM"])
