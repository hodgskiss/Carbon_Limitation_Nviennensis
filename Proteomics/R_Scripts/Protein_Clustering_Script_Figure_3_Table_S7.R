#This script clusters proteins as seen in Figure 3 of the manuscript.
#It also constructs the bulk of information found in Table_S7.

library(ggplot2)
library(gplots)
library(tidyverse)
library(expss)
library(viridis)
library(RColorBrewer)
library(wesanderson)
#Load necessary libraries.

avg_list <- read.table("Protein_Averages.txt", header=TRUE, row.names = 1, dec=".",sep="\t")
#Loads the resulting averages (after removing outliers) of proteins in each condition.

avg_scale <- scale(avg_list)
#Scales the data in each column (i.e. for each protein across conditions).

avg_scale <- t(avg_scale)
#Transpose the scaled data.

all_p_vals <- read.table("all_p_vals", header=TRUE, row.names = 1, dec=".",sep=",")
#Load p-value information for all proteins.

labels <- read.csv(file="Gene_Codes_List.csv", header=TRUE, sep=",")
#Load gene annotation data for all proteins.

carbon <- c("A","D","G","C","F","B","E")
#Can be used to arrange columns in a specific order when creating a heatmap later. 
#This is in order of most to least carbon concentration. 

###CREATE THE PROTEIN CLUSTERING HEATMAP IN FIGURE 3###

x11()
hm <- heatmap.2(avg_scale[,carbon], Colv = FALSE,density.info = "none",labRow = FALSE, dendrogram = "row",srtCol = 0  ,trace=c("none"))
#Create a basic heatmap using scaled protein data.
#Proteins are hierarchically clustered.
#Conditions are ordered based on carbon concentration.

hc <- as.hclust( hm$rowDendrogram )
#Passes the row dendrogram into the hclust fuction.

groups <- cutree( hc, k=7 )
#Separates the dendrogram into a number of clusters.
#The order of proteins when written out is the same as "avg_list" and "avg_scale". 

clusterCols <- brewer.pal(n=length(unique(groups)),"YlGnBu")
#Specfies a clustering color palette for the proteins.

myClusterSideBar <- clusterCols[groups]
#Identifies each cluster sith colored side bars.

x11()
hm <- heatmap.2(avg_scale[,carbon], col = rev(brewer.pal(9, "RdYlBu")), Colv = FALSE,density.info = "none",labRow = FALSE, dendrogram = "none",srtCol = 0  ,trace=c("none"),RowSideColors=myClusterSideBar,key=TRUE)

###EXTRACT INFORMTAION FOR PROTEINS IN EACH CLUSTER###

groups_ordered <- groups[hc$order]
#Organized clusters by appearance on heatmap from bottom to top.

groups_ordered_df <- as.data.frame(groups_ordered)
#Convert to a dataframe.

groups_ordered_df <- mutate(groups_ordered_df, Protein.Code=row.names(groups_ordered_df))
#Makes a data frame and adds the Protein.Code column.
#For unknown reasons, the addition of the Protein.Code column will only work if there is a different name for the data frame.
#Not sure why this should matter. 
#It might be a conflict of the name of the data frame and the name of a column within the data frame. 


###CREATE THE BULK OF INFORMATION FOR TABLE_S7###

values <- hm$carpet
#Extracts the scaled values created for each protein to make the heatmap. 

values <- t(values)
#Transpose data.

values <- as.data.frame(values)
#Converts to a data frame.

values <- mutate(values, Protein.Code=row.names(values))
#Adds the Protein.Code as a column.

values <- mutate(values, Cluster= vlookup(values$Protein.Code,groups_ordered_df,lookup_column = "Protein.Code",result_column = "groups_ordered"))
#Adds the cluster number.

values <- mutate(values, p_adj=vlookup(values$Protein.Code, all_p_vals,lookup_column="Protein",result_column = "p_adj"))
values <- mutate(values, Method=vlookup(values$Protein.Code,all_p_vals,lookup_column = "Protein",result_column = "Method"))
#Adds the adjusted p-value and the method used to obtain the p-value (based on normality and variance tests).

values <- separate(values,Protein.Code, into=c("Protein.Code",NA),sep = "\\.",extra="drop", fill="right")
#Remove multiple protein codes.

values <- mutate(values, Gene_Name= vlookup(values$Protein.Code,labels,lookup_column = "Protein.Code",result_column = "Gene.Name"))
values <- mutate(values, Locus_Tag=vlookup(values$Protein.Code,labels,lookup_column = "Protein.Code",result_column = "Locus.Tag"))
values <- mutate(values, Protein_Product= vlookup(values$Protein.Code,labels,lookup_column = "Protein.Code",result_column = "Protein.Product"))
#Adds other helpful information from the labels data frame. 

values <- arrange(values, Cluster)
#Organizes the values data frame by cluster number. 

values <- values[,c(8,9,13,12,10,11,14,1,6,4,2,7,5,3)] 
#Reorders columns in values data frame.
#IMPORTANT: This organization puts the samples in alphabetical order at the end of the data frame.
#So, NOT by carbon concentration, i.e. NOT the same order as the heatmap. This is easier for data manipulation downstream.
#Be cautious of organizing like this; mistakes are very easy to make and will alter later interpretation. 

#The following loop will make a data frame for each cluster from the values data frame. 
for(i in 1:7) {
  name <- paste("Cluster", i, sep = "_")
  assign(name, filter(values,Cluster==i))
}

means <- aggregate(values[, 8:14], list(values$Cluster), mean)
means <- means[,c(1,2,5,4,8,7,3,6)]
#Calculates the mean scaled heatmap value for each condition in each cluster. 

heatmap_cluster_order <- unique(groups_ordered_df$groups_ordered)
#Lists the clusters from bottom of heatmap to top.

#Pulling from values table

values <- values[,c(1:7,8,11,14,10,13,9,12)] 
values <- mutate(values, sig_check=if_else(p_adj<=0.05,"yes","no"))
#Identify if proteins show a significant difference between conditions.

cond_sub_max <- values[,8:14]
cond_sub_max <- cond_sub_max %>% mutate(Max = names(.)[max.col(.)])
cond_sub_max <- cond_sub_max %>% mutate(Prot=row.names(.))
#Identify the condition with the maximum value for each protein.

cond_sub_min <- values[,8:14]
cond_sub_min <- cond_sub_min %>% mutate(Min = names(.)[max.col(-.)])
cond_sub_min <- cond_sub_min %>% mutate(Prot=row.names(.))
#Identify the condition with the minimum value for each protein.

values <- mutate(values, Max= vlookup(values$Protein.Code,cond_sub_max,lookup_column = "Prot",result_column = "Max"))
values <- mutate(values, Min= vlookup(values$Protein.Code,cond_sub_min,lookup_column = "Prot",result_column = "Min"))
#Add max and min condition to "values" data frame.

values <- values[,c(1,8:14,5:6,3,15,2,16,17,4,7)]
#Organize columns in data frame.

write.table(values,"Table_S7.txt", col.names=NA, sep = "\t")
#Save extracted information to a text file.

