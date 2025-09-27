#This script processes the proteomic data, removes outliers, and calculates the 
#standard deviations for each protein in each condition.

library(ggplot2)
library(dplyr)
library(reshape2)
library(factoextra)

prot <- read.table("LFQ_new.txt", header=TRUE, row.names = 1, dec=".",sep="\t")
#Loads proteomic data.

prot <- as.data.frame(t(prot))
#Transposes the data so that the protein IDs are column headers.

Condition <- as.factor(c(rep("A", times=9),"Bstat","B",rep("Bstat",times=4),rep("B",times=4),rep("C",times=5),rep("D",times=5),rep("E",times=5),rep("F",times=5),rep("G",times=5)))
#Creates a column indicating the different conditions for each row. 

prot$Condition <- Condition
#Attaches descriptive column to the data frame.  This helps for subsetting the data later. 

prot <- prot[c(prot$Condition!="Bstat"),]
#Remove data for stationary cultures.

prot$Condtion <- NULL
#Remove condition column.

Condition <- as.factor(c(rep("A", times=9),rep("B",times=5),rep("C",times=5),rep("D",times=5),rep("E",times=5),rep("F",times=5),rep("G",times=5)))
#Create a new condition vector without stationary cultures.

prot$Condition <-Condition
#Attach new vector to data frame.

protnames <- colnames(prot)
#Creates a vector with all of the column names, i.e. Protein IDs.

m <-length(protnames)-1
#Counts the number of Protein IDs. -1 is to not count the "Condition" column.

Condition_Names <- as.factor(c("A","B","C","D","E","F","G"))
#Makes a vector with each condition label.
#This has to be done manually to get the order right.  the command "levels" arranges it alphabetically. 

n <- length(Condition_Names)
#Counts the number of conditions. 

k <- melt(prot, id.vars="Condition")
#Melts the data frame into 3 columns. The first is the condition, the second is the metabolite, and the third is the value.
#This will make it easier to subset. 

e_list = data.frame(row.names = row.names(prot))
#Creates an empty data frame with the sample names as row.names. 

h_list=data.frame(row.names = Condition_Names)
#Creates an empty data frame for storing data. 

###CALCULATING PROTEIN STANDARD DEVIATIONS###
#This is the same process as for the PCA plot except that outliers
#are NOT replaced with averages.  They are simply removed and 
#marked as NA.

for(j in c(1:m)){
#Starts a for loop that is to be done for the number of proteins "m". 

  f_list = vector()
  #Creates and empty vector to store data later.
  #This NEEDS to be recreated here each time for the loop to work.
  
  g_list = vector()
  #Creates and empty vector to store data later.
  #This NEEDS to be recreated here each time for the loop to work.

for(i in c(1:n)) {
  #Starts a for loop for the number of conditions "n".
  #This initiates a process that will be done on each protein within each condition. 
  
  f <- subset(k, variable==protnames[j]&Condition==Condition_Names[i])
  #Creates a subset of the data frame based on the protein and the condition.
  #"j" is the number in the protnames vector signifying the protein.
  #"i" is the number in the condition vector signifying the condtion. 
  
  q <- boxplot.stats(f$value)$out
  #Checks and stores outlier information based on a boxplot. 

  f$value[f$value %in% q] <- NA      
  #For each outlier, this replaces the value with "NA" and re-saves the new vector inro f. 
  
  sd <- sd(f$value, na.rm = TRUE)
  #Calculates the standard deviation value after the outliers are removed.
  
  g_list <- rbind(g_list,sd)
  #Builds a vector containing the standard deviations, after outlier removal, for each protein.
  
  f_list <- rbind(f_list,f)
  #Rebuilds the first column from the original data frame with the adjusted values. 
  #This is dependent on the condition labels in "Condition_Names" being in the same order as the original data.frame. 
  
}

e_list <- cbind(e_list, f_list[3])
#Takes the newly constructed values and attaches it to the empty data frame.
#This rebuilds the whole data.frame with the corrected values.  

h_list <- cbind(h_list, g_list)
#Adds the condition column to the reconstructed data frame.

}

e_list$Condition <- Condition
#Attaches condition vector to data frame.

h_list <- h_list
#Adds the condition column to the reconstructed vector. 

names(e_list) <- protnames
names(h_list) <- protnames[1:m]
#Replaces the column headers with the original protein names.
#IMPORTANT: The protnames vector needs to be the same as the original data.frame headers.  Otherwise the proteins will be out of order. 
#"protnames" contains the name "Condition" for its last entry. 
#This is compatible with the e_list as the condition vector was added.
#"protnames" needs to be only length "m" (i.e., 1490) to work for h_list.

prot_edit <- e_list
#Copies the new data frame to a new data object. 

prot_edit$Condition <- NULL
#Removes condition column.

prot_edit <- prot_edit[,c(colSums(prot_edit, na.rm = TRUE) !=0)]
#Removes proteins that only have 0 or NA.
#This creates the same data frame when calculating the averages.

sd_list <- h_list
#Saves protein standard deviations to a new data frame.

sd_list <- sd_list[,c(colSums(sd_list) !=0)]
#Removes proteins with only 0.

write.table(sd_list, 'Protein_Standard_Deviations.txt', col.names=NA, sep = "\t")




