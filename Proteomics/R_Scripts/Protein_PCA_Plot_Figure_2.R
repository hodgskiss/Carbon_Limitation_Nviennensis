#The following script processes proteomic data and creates a protein PCA plot.
#This is the script correpsonding to Figure 2 of the manuscript.

library(ggplot2)
library(ggrepel)
library(dplyr)
library(reshape2)
library(factoextra)
library(naniar)
#Load necessary work packages. 

prot <- read.table("LFQ_new.txt", header=TRUE, row.names = 1, dec=".",sep="\t")
#Loads proteomic data.

prot <- as.data.frame(t(prot))
#Transposes the data so that the protein IDs are column headers. 

Condition <- as.factor(c(rep("A", times=9),"Bstat","B",rep("Bstat",times=4),rep("B",times=4),rep("C",times=5),rep("D",times=5),rep("E",times=5),rep("F",times=5),rep("G",times=5)))
#Creates a column indicating the different conditions for each row. 

prot$Condition <- Condition
#Attaches descriptive column to the data frame.  This helps for subsetting the data later. 

prot <- prot[c(prot$Condition!="Bstat"),]
#Removes conditions that were in stationary phase.

prot$Condition <- NULL
#Removes condition column.

Condition <- as.factor(c(rep("A", times=9),rep("B",times=5),rep("C",times=5),rep("D",times=5),rep("E",times=5),rep("F",times=5),rep("G",times=5)))
#Creates a new condition vector without stationary phase cultures.

prot$Condition <-Condition
#Adds new condition vector to data frame.

protnames <- colnames(prot)
#Creates a vector with all of the column names, i.e. protein IDs.

m <-length(protnames)-1
#Counts the number of protein IDs. -1 is to not count the "Condition" column.

Condition_Names <- c("A","B","C","D","E","F","G")
#Makes a vector with each condition label.
#This has to be done manually to get the order right.  the command "levels" arranges it alphabetically. 

n <- length(Condition_Names)
#Counts the number of conditions. 

k <- melt(prot, id.vars="Condition")
#Melts the data frame into 3 columns. The first is the condition, the second is the protein, and the third is the value.
#This will make it easier to subset. 

e_list = data.frame(row.names = row.names(prot))
#Creates an empty data frame with the sample names as row.names. 

###PROTEIN CLEANING for PRINCIPAL COMPONENT ANALYSIS (PCA)###
#The following will identify outliers from the protein dataset.
#Outliers will be replaced by the average of the non-outliers.
#This will prevent outliers from distorting the PCA analysis.

for(j in c(1:m)){
#Starts a for loop that is to be done for the number of proteins "m". 

f_list = list()
#Creates and empty vector to store data later. 

for(i in c(1:n)) {
  #Starts a for loop for the number of conditions "n".
  #This initiates a process that will be done on each protein within each condition. 
  
  f <- subset(k, variable==protnames[j]&Condition==Condition_Names[i])
  #Creates a subset of the data frame based on the protein and the condition.
  #"j" is the number in the protnames vector signifying the protein.
  #"i" is the number in the condition vector signifying the condtion. 
  q <- boxplot.stats(f$value)$out
  #Checks and stores outlier information based on a boxplot. 
  p <- length(q)
  #Determines the number of outliers. 
  
  if(length(q) > 0){
    #Starts an if loop that is conditional on the presence of outliers, i.e. length of q > 0.
    
    for(l in c(1:p)) local({
    #Starts a loop for the number of outliers, p. 
     
      f <<- replace_with_na(f, replace=list(value=q))
      #For each outlier, this replaces the value with "NA" and re-saves the new vector inro f. 
    
    })
    
    avg <- mean(f$value, na.rm = TRUE)
    #Calculates the average value after the outliers are removed. 
    f[is.na(f)] <- avg
    #Replaces the NA values in vector f with the average of the non-outliers. 
    
  }
  
   f_list <- rbind(f_list,f)
   #Rebuilds the first column from the original data frame with the adjusted values. 
   #This is dependent on the condition labels in "Condition_Names" being in the same order as the original data.frame. 
  
 
}

e_list <- cbind(e_list, f_list[,3])
#Takes the newly constructed values and attaches it to the empty data frame.
#This rebuilds the whole data.frame with the corrected values.  


}


e_list$Condition <- Condition
#Adds the condition column to the reconstructed vector.

names(e_list) <- protnames
#Replaces the column headers with the original Protein IDs.
#IMPORTANT: The protnames vector needs to be the same as the original data.frame headers.  Otherwise the proteins will be out of order. 

prot_edit <- e_list
#Copies the new data frame to a new data object. 

prot_edit$Condition <- NULL
#Removes the condition column.

prot_edit <- prot_edit[,c(colSums(prot_edit) !=0)]
#Removes proteins that have only 0 values.

prot_edit$Condition <- Condition
#Adds the condition vector to the data frame. 

prot_edit_pca <- prot_edit
#Creates a data.frame to be used for pca plots. 

prot_edit_pca$Condition <- NULL
#Removes the Condition column to make the data.frame compatible for PCA analysis.

pca <- prcomp(prot_edit_pca, scale. = TRUE)
#Performs the PCA calculations.

pca_df <- as.data.frame(pca$x)
#Saves the data needed for plotting to a data frame.

pca_df$Condition <- Condition
#Adds the condition vector to the data frame.

sum <- summary(pca)
#Summarizes the PCA analysis.

var <- round(sum$importance[2,]*100,2)
#Identifies the percentage of variance for each principal component.
#Rounds the variance to two decimal places.

rotation <- pca$rotation
#Selects information for the differenct loading vectors of the PCA analysis.

write.table(rotation,"Protein_PCA_Rotations.txt", col.names=NA, sep = "\t")
#Saves the rotation data to a text file.

percentage_label <- paste( colnames(pca_df), "(", paste(as.character(var), "%", ")", sep="") )
#Creates levels for each principal component for plotting.


###Plot PCA using ggplot.###
plot <- ggplot(pca_df, aes(x=PC1,y=PC2,color=Condition))
plot <- plot + geom_point()
plot <- plot + ggtitle("Proteins LFQ Normalized-No Outliers, No Bstat, New") + xlab(percentage_label[1]) + ylab(percentage_label[2])
plot <- plot + theme(plot.title=element_text(hjust=0.5))
plot <- plot + geom_label_repel(aes(label=rownames(pca_df)),point.padding = 0.7)
plot <- plot + stat_ellipse()
plot

###Plot PCA using fviz.###
fplot <- fviz_pca_ind(pca, habillage = pca_df$Condition, addEllipses=TRUE, ellipse.level=0.95, palette=c("#157E41","#F47842","#EFE809","#71AD46","#d70e27","#E2B700","#C2E543"))
fplot <- fplot +ggtitle("Proteins LFQ Normalized-No Outliers, No Bstat, New") + xlab(percentage_label[1]) + ylab(percentage_label[2]) + theme(plot.title=element_text(hjust=0.5))
fplot <- fplot+guides(fill = guide_legend(title = "Condition", override.aes = aes(label = "")))

x11()
fplot

