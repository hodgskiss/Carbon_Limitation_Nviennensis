#This script uses either one-way analysis of variance (ANOVA) or 
#Kruskal-Wallis tests to determine proteins that change across
#the different conditions.

library(tidyverse)
library(rstatix)
library(car)
library(FSA)
#Load necessary packages.

labels <- read.csv(file="Gene_Codes_List.csv", header=TRUE, sep=",")
#Load gene annotation information. 

dat <- read.table("Proteins_LFQ_Normalized_No_Outliers_no_Bstat.txt", header=TRUE, row.names = 1, dec=".",sep="\t",na.strings = "NA")
#Load data that has been outlier corrected and removed the stationary B condition. 

condition_info <- read.table("condition_info.txt", header=TRUE, row.names = 1, dec=".",sep="\t")
#Load information about each condition. 

Condition <- condition_info$Condition
#The label for each biological condition of each sample.  Each condition has at least 5 replicates.

ROS <- condition_info$ROS
#Indicates whether pyruvate or catalase was used to counter ROS.

Carbon <- condition_info$Carbon
#Indicates the carbon concentration in each sample. 

Stress <- condition_info$Stress
#Indicates how stressed the cultures were based on growth curves. 

data <- dat%>%mutate(ROS)%>%mutate(Carbon)%>%mutate(Condition)
#Adds conidition information to data frame.

avg_list <- read.table("Protein_Averages.txt", header=TRUE, row.names = 1, dec=".",sep="\t")
#Loads the calculated averages for each protein in each condition. 

avg_list_t <- as.data.frame(t(avg_list))
#Transposes the average list.

###PROTEINS UNIQUE TO ONE CONDITION###
#Before proceeding with ANOVA or Kruskal-Wallis, proteins found
#in only one condition are identified and removed from the dataset.

only_A <- avg_list_t%>%filter(A!=0&B==0&C==0&D==0&E==0&F==0&G==0)
only_B <- avg_list_t%>%filter(A==0&B!=0&C==0&D==0&E==0&F==0&G==0)
only_C <- avg_list_t%>%filter(A==0&B==0&C!=0&D==0&E==0&F==0&G==0)
only_D <- avg_list_t%>%filter(A==0&B==0&C==0&D!=0&E==0&F==0&G==0)
only_E <- avg_list_t%>%filter(A==0&B==0&C==0&D==0&E!=0&F==0&G==0)
only_F <- avg_list_t%>%filter(A==0&B==0&C==0&D==0&E==0&F!=0&G==0)
only_G <- avg_list_t%>%filter(A==0&B==0&C==0&D==0&E==0&F==0&G!=0)
#Identifies and saves proteins found in only one condition. 
#For these proteins, calculating ANOVA or Kruskal-Wallis would not
#ve helpful.


#The following subsets proteins found in only one condition. 
A_set <- as.data.frame(row.names(only_A))
colnames(A_set) <- c("Protein")
A_set$p_value <- NA
A_set$p_adj <- NA
A_set$Method <- c("only in A")

B_set <- as.data.frame(row.names(only_B))
colnames(B_set) <- c("Protein")
B_set$p_value <- NA
B_set$p_adj <- NA
B_set$Method <- c("only in B")

C_set <- as.data.frame(row.names(only_C))
colnames(C_set) <- c("Protein")
C_set$p_value <- NA
C_set$p_adj <- NA
C_set$Method <- c("only in C")

D_set <- as.data.frame(row.names(only_D))
colnames(D_set) <- c("Protein")
D_set$p_value <- NA
D_set$p_adj <- NA
D_set$Method <- c("only in D")

E_set <- as.data.frame(row.names(only_E))
colnames(E_set) <- c("Protein")
E_set$p_value <- NA
E_set$p_adj <- NA
E_set$Method <- c("only in E")

F_set <- as.data.frame(row.names(only_F))
colnames(F_set) <- c("Protein")
F_set$p_value <- NA
F_set$p_adj <- NA
F_set$Method <- c("only in F")

G_set <- as.data.frame(row.names(only_G))
colnames(G_set) <- c("Protein")
G_set$p_value <- NA
G_set$p_adj <- NA
G_set$Method <- c("only in G")

only_sets <- rbind(A_set,B_set,C_set,D_set,E_set,F_set,G_set)
#Combines the information for proteins found in one condition.

only_one_condition <- rbind(only_A,only_B,only_C,only_D,only_E,only_F,only_G)
#Extracts the abundance information for proteins only found in one condition.

only_one_condition_proteins <- noquote(row.names(only_one_condition))
#Creates a vector of proteins only found in one condition.

data_trimmed <- data[, !names(data) %in% only_one_condition_proteins]
#Removes data only found in one condition. 

###TESTING FOR ANOVA ASSUMPTIONS###
#Each protein is tested for normality (Shapiro Test) and homogeneity of 
#variance (Levene Test).  Proteins that meet both assumptions will be tested
#with ANOVA.  Those that fail will be tested with Kruskal-Wallis.

#Shapiro Test Loop
#Identify proteins that are normally distributed. 

names <- colnames(data_trimmed)
names <- as.vector(names)
#Collect identifiers of proteins to be tested.

names <- names[-c(1199:1201)]
#Remove condition information.

Conditions <- c("A","B","C","D","E","F","G")
#Specify the name of each condition.

m <- length(names)
#Determine number of proteins to be tested.

n <- length(Conditions)
#Determine number of conditions to be accounted for. 

SW_p_val <- vector()
#Create an empty vector to store data.

for(i in 1:m){
  
  j <- names[i]
  #Extract the protein to be tested.
  
  filler <- vector()
  #Create an empty vector.
  
  for(f in 1:n){
    
    k <- Conditions[f]
    #Extract the condition to be tested.
    
    summation <- data_trimmed%>%filter(Condition==k)%>%select(all_of(j))%>%data.matrix()%>%as.numeric%>%sum(.,na.rm=TRUE)
    #Extract the protein information for the seleccted protein and condition.
    if(summation == 0){
      
      row = c(j,k,"na","na")
      #Exclude proteins that have zero or NA values in a specified condition.
      
      
    } else{
      
      t <- data_trimmed%>%filter(Condition==k)%>%select(all_of(j))%>%data.matrix()%>%as.numeric%>%shapiro_test()
      #Perform the shapiro test and tore results.
      
      row = c(j,k,as.numeric(t[1,2]),as.numeric(t[1,3]))
      #Extract results.
    }
    
    
    filler <- rbind(filler, row)
    #Save results.
    
  }
  
  SW_p_val <- rbind(SW_p_val,filler)
  #Add results to previously calculated results.
  
}

SW_p_val <- as.data.frame(SW_p_val)
#Compile Shapiro test results for all proteins in a dataframe.

colnames(SW_p_val) <- c("Protein","Condition","statistic","p_value")
#Rename columnbs of Shapiro Test results.

row.names(SW_p_val) <- NULL
#Remove row names.

SW_p_val$p_value <- as.numeric(as.character(SW_p_val$p_value))
#Make sure p-values are numeric.

SW_p_val <- SW_p_val%>%mutate(p_adj=p.adjust(p_value,"BH"))
#Calculate adjusted p-values using the Benjamini-Hochberg method.

write.csv(SW_p_val, "./SW_p_val.csv")
#Write results to csv file.

normality_check <- SW_p_val%>%filter(p_adj<0.05)
#Identify proteins that did not pass the normality check.
#This was done using a cut-off of 0.05. Proteins below this are not
#considered to be distributed normally.
#This is checked for every condition with each protein.

non_normal_proteins <- normality_check$Protein
#Stores information of proteins that did not pass the normality check. 

write.csv(normality_check, "./protein_normality_check.csv")
#Saves non-normal proteins to a csv file.

#Levene Loop Test 
#Testing proteins for homogeneity of variance.

names <- colnames(data_trimmed)
names <- as.vector(names)
#Collect identifiers of proteins to be tested.

names <- names[-c(1199:1201)]
#Remove condition information.

m <- length(names)
#Identify number of proteins

Lev_p_val <- vector()
#Create an empty vector to store data.

for(i in 1:m){
  
  j <- names[i]
  #Identify the protein.
  
  p <- leveneTest( data_trimmed[,i] ~ Condition, data = data_trimmed)
  #Perform the Levene test across all conditions.
  
  row2 = c(j,as.numeric(p[1,2]),as.numeric(p[1,3]))
  #Extract relevant information.
  
  Lev_p_val <- rbind(Lev_p_val,row2)
  #Store with previously calculated proteins.
  
  
}


Lev_p_val <- as.data.frame(Lev_p_val)
#Convert results to a data frame.

colnames(Lev_p_val) <- c("Protein","F_value","p_value")
#Rename columns.

row.names(Lev_p_val) <- NULL
#Remove row names.

Lev_p_val$p_value <- as.numeric(as.character(Lev_p_val$p_value))
#Make sure p-values are numeric.

Lev_p_val <- Lev_p_val%>%mutate(p_adj=p.adjust(p_value,"BH"))
#Adjust p-values using the Benjamini-Hochberg method.

write.csv(Lev_p_val, "./Lev_p_val.csv")
#Save Levene test results to a csv file.

homogeneity_check <- Lev_p_val%>%filter(p_adj<0.05)
#Identify proteins that did not pass the variance check.
#This was done using a p-value cut-off of 0.05.

non_homogenous_proteins <- homogeneity_check$Protein
#Stores proteins that did not pass the variance check.

look_var <- data[,non_homogenous_proteins]
#Subsets proteins that did not pass the variance check.

look_var$Condition <- Condition
#Adds condition information.

look_norm <- data[,non_normal_proteins]
#Subsets proteins that did not pass the normality check.

fail_proteins <- unique(c(as.character(non_homogenous_proteins),as.character(non_normal_proteins)))
#Combines the proteins that did not pass the normality and/or homogenetiy of variance checks. 

look_fail_proteins <- data[,c(fail_proteins,"Condition")]
#Extracts protein data and combines the proteins that did not pass the
#normality and/or homogenetiy of variance checks. 

non_ANOVA <- data_trimmed[,c(noquote(fail_proteins),"Condition")]
#Separates proteins that cannot be analyzed via ANOVA.

for_ANOVA <- data_trimmed[, !names(data_trimmed) %in% fail_proteins]
#Separates proteins that will be analyzed by ANOVA.

write.csv(homogeneity_check, "./protein_homogeneity_check.csv")
#Saves proteins that failed the homogeneity check to a csv file. 

###ANOVA Analysis###

n <- length(colnames(for_ANOVA))-3
#Specifies the number of times to run the analysis based on the number of proteins.

m <- colnames(for_ANOVA)
#Specifies the number of times to run the analysis based on the number of proteins.

p_vals <-vector()
#Creates an empty vector to store p-value information.


for(i in 1:n){
  
  analysis <- aov( for_ANOVA[,i] ~ Condition, data = for_ANOVA)
  #Runs a one-way ANOVA anlysis based on the Condition column.
  
  summary <- summary(analysis)
  #Get a summary of the results.
  
  p <- summary[[1]][1,5]
  #Extracts the p-value.
  
  name <- colnames(for_ANOVA)[i]
  #Extracts the name of the protein.
  
  row <- c(name,p)
  
  p_vals <- rbind(p_vals, row)
  #Builds a new table.
  
  
}

p_vals <- as.data.frame(p_vals)
#Makes the table a data frame.

colnames(p_vals) <- c("Protein","p_value")
#Rename the columns.

row.names(p_vals) <- NULL
#Remove the row names.

p_vals$p_value <- as.numeric(as.character(p_vals$p_value))
#Changes the p-values to a numeric value.  As a data frame these are read as factors.
#This is important for using the numbers as values.
#It changes the numbers from factors to numbers. 

p_vals <- p_vals%>%mutate(p_adj=p.adjust(p_value,"BH"))
#Calculate adjusted p-values based on the Benjamini-Hochberg method.

p_vals <- p_vals%>%arrange(p_adj)
#Order the proteins from lowest to highest based on adjusted p-value.

p_vals$Method <- c("ANOVA")
#Add information about the method used.

###Kruskal-Wallis ANALYSIS###

n2 <- length(colnames(look_fail_proteins))-1
#Specifies the number of times to run the analysis based on protiein number.

m2 <- colnames(look_fail_proteins)
#Specifies the number of times to run the analysis based on protein number.

p_vals2 <-vector()
#Creates an empty vector to store p-value information.


for(i in 1:n2){
  
  analysis_KW <- kruskal.test( look_fail_proteins[,i] ~ Condition, data = look_fail_proteins)
  
  #Runs a one-way Kruskal-Wallis anlysis based on the Condition column.

  p <- analysis_KW[[3]]
  #Extracts the p-value.
  
  name <- colnames(look_fail_proteins)[i]
  #Extracts the name of the protein.
  
  row <- c(name,p)
  #Extracts the necessary information.
  
  p_vals2 <- rbind(p_vals2, row)
  #Builds a new table.
  
  
}

p_vals2 <- as.data.frame(p_vals2)
#Makes the table a data frame.

colnames(p_vals2) <- c("Protein","p_value")
#Renames the columns.

row.names(p_vals2) <- NULL
#Removes the row names.

p_vals2$p_value <- as.numeric(as.character(p_vals2$p_value))
#Changes the p-values to a numeric value.  As a data frame these are read as factors.
#This is important for using the numbers as values.
#It changes the numbers from factors to numbers.

p_vals2 <- p_vals2%>%mutate(p_adj=p.adjust(p_value,"BH"))
#Calculate adjusted p-values based on the Benjamini-Hochberg method.

p_vals2 <- p_vals2%>%arrange(p_adj)
#Order the proteins from lowest to highest based on adjusted p-value.

p_vals2$Method <- c("Kruskal-Wallis")
#Add information about the method used.

###COMBINE RESULTS OF ANOVA, KRUSKAL-WALLIS, AND ONLY ONE CONDITION TESTS###

all_p_vals <- rbind(p_vals,p_vals2,only_sets)
#Combine the two datasets along with proteins only found in one condition.

all_p_vals <- all_p_vals%>%arrange(p_adj)
#Order the proteins from lowest to highest based on adjusted p-value.

write.csv(all_p_vals,"all_p_vals")
#Save the final results to a csv file.
