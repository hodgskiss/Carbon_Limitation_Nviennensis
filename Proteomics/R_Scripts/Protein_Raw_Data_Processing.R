library(ggplot2)
library(dplyr)
library(factoextra)
#Loads the necessary packages. 


protein_groups <- read.table("proteinGroups.txt", header=TRUE,  dec=".",sep="\t")
#Loads proteomic data.

protein_groups_edit <- protein_groups[,c(2,240:283)]
#Selects only the proteins names and LFQ intensites.

col_names <- colnames(protein_groups_edit)
new_col_names <- gsub("LFQ.intensity.","",col_names)
colnames(protein_groups_edit) <- new_col_names
#Delete the "LFQ.intensity." portion of each name.

prot <- protein_groups_edit[-c(grep("REV_",protein_groups_edit$Majority.protein.IDs),grep("CON_",protein_groups_edit$Majority.protein.IDs)),]
#Remove reverse and contaminant proteins.

row_names <- prot[,1]
row.names(prot) <- row_names
prot[,1]<-NULL
#Replace row names.

LFQ_new <- prot
#Save formatted LFQ data.

write.table(LFQ_new, "LFQ_new.txt", col.names=NA, sep = "\t")
#Write LFQ data to its own table. 
