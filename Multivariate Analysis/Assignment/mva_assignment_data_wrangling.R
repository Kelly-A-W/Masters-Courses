# Initialise

# Clear environment
rm(list=ls())
# Clear all plots
if(!is.null(dev.list())) dev.off()

# Install packages
library(readxl)

sheet1 <- read_excel("dBCG_Final data file_.xlsx", sheet = 1)
sheet1 <- as.data.frame(sheet1)
head(sheet1)

sheet2 <- read_excel("dBCG_Final data file_.xlsx", sheet = 2)
sheet2 <- as.data.frame(sheet2)
head(sheet2)

sheet3 <- read_excel("dBCG_Final data file_.xlsx", sheet = 3)
sheet3 <- as.data.frame(sheet3)
head(sheet3)

sheet4 <- read_excel("dBCG_Final data file_.xlsx", sheet = 4)
sheet4 <- as.data.frame(sheet4)
# columns 36, 69 and 102 are NAs, so remove
sheet4 <- sheet4[,-c(36, 69, 102)]
head(sheet4)

sheet5 <- read_excel("dBCG_Final data file_.xlsx", sheet = 5)
sheet5 <- as.data.frame(sheet5)
# columns 9 and 15 are NAs, so remove
sheet5 <- sheet5[,-c(9, 15)]
head(sheet5)

sheet6 <- read_excel("dBCG_Final data file_.xlsx", sheet = 6)
sheet6 <- as.data.frame(sheet6)
head(sheet6)

sheet7 <- read_excel("dBCG_Final data file_.xlsx", sheet = 7)
sheet7 <- as.data.frame(sheet7)
head(sheet7)


#----------------
# SHEET 1
#----------------
# sheet1 is demographic data
names(sheet1) <- c("PID", "Group", "Status", "Race", "Sex", "Age")
head(sheet1)
sheet1$Group <- as.factor(sheet1$Group)
sheet1$Status <- as.factor(sheet1$Status)
sheet1$Race <- as.factor(sheet1$Race)
sheet1$Sex <- as.factor(sheet1$Sex)
# Group refers to age group in weeks
# Age refers to the child's exact age measured in weeks
str(sheet1)


#----------------
# SHEET 2
#----------------
# columns 4-16 are BCG-UNS frequencies
# columns 17-21 are BCG frequencies
names(sheet2) <- c("PID", "Status", "Age", "BCGUNS_CD33+Myeloid_IL6+", "BCGUNS_CD33+Myeloid_IL1B+", 
                   "BCGUNS_CD33+Myeloid_TNF+",  "BCGUNS_CD4+T_Gr+", "BCGUNS_CD4+T_G+", 
                   "BCGUNS_CD4+T_IL1B+", "BCGUNS_CD4+T_IL6+", "BCGUNS_CD4+T_TNF+",
                   "BCGUNS_CD56+CD16- NK_G+", "BCGUNS_CD56++ NK_G+", "BCGUNS_CD56+CD16+ NK_G+",
                   "BCGUNS_CD56-CD16+ NK_G+", "BCGUNS_Total NK_G+", "BCG_CD56+CD16- NK_Perf+", 
                   "BCG_CD56++ NK_Perf+", "BCG_CD56+CD16+ NK_Perf+",
                   "BCG_CD56-CD16+ NK _Perf+", "BCG_Total NK_Perf+")
# we don't need "Status" and "Age" since they are already included in sheet1
sheet2 <- sheet2[,-c(2,3)]
str(sheet2)


#----------------
# SHEET 3
#----------------
# sheet3 contains the perforin MFI values with BCG stim
names(sheet3) <- c("PID", "Status", "Age", "BCG_CD56+CD16-NK_Perf+ MFI",
                   "BCG_CD56++ NK_Perf+ MFI", "BCG_CD56+CD16+NK_Perf+ MFI", 
                   "BCG_CD56-CD16+NK_Perf+ MFI", "BCG_Rainbow NK_Perf+ MFI")
# we don't need "Status" and "Age" since they are already included in sheet1
sheet3 <- sheet3[,-c(2,3)]
str(sheet3)


#------------------
# SHEET 4 & SHEET 5
#------------------
# sheet4 only has measurements for babies in group  5W
# sheet5 only has measurements for babies in group 9W
# since these variables are measured on only one group and not the other, and there is reason to believe that
# age group has a significant influence on these results, we will ignore these two sheets since we want to 
# look at a dataset containing both groups togther
# will maybe return to these sheets if looking only at 5W old babies or 9W old babies


#------------------
# SHEET 6
#------------------

head(sheet6)
# each patient has three rows, corresponding to the BCG, BCG-UNS and UNS measurements respectively
# we only need the BCG-UNS rows
sheet6 <- sheet6[sheet6$Stim=="BCG-UNS",]
# we can remove "Age", "Status" and "Stim" columns
sheet6 <- sheet6[,-c(2,3,4)]
names(sheet6) <- c("PID", "BCGUNS_IFNG", "BCGUNS_IL1B", "BCGUNS_IL10", "BCGUNS_IL13", "BCGUNS_IL17A", "BCGUNS_IL18", 
                   "BCGUNS_IL2", "BCGUNS_IL21", "BCGUNS_IL23", "BCGUNS_IL4", "BCGUNS_IL6" , "BCGUNS_IL9", "BCGUNS_TNFA",
                   "BCGUNS_GMCSF", "BCGUNS_IL12p70", "BCGUNS_IL5", "BCGUNS_IL22", "BCGUNS_IL27")
str(sheet6)



#------------------
# SHEET 7
#------------------
head(sheet7)
str(sheet7)
# IL5 is a character but needs to be a numeric
sheet7$IL5 <- as.numeric(sheet7$IL5)
# NA should be 0.02
# so set it as such
sheet7$IL5[68] <- 0.02
# while some of the variables may look the same as in sheet6, they are NOT the same
# these where measured in terms of PBMC, so will use that to differeniate the variables
# for each PID we have 5 rows: UNS, M.tb, LPS, HKSA, HKCA
# can remove HKSA because not of interest
sheet7 <- sheet7[sheet7$Stim !="HKSA",]
# now we need to subtract UNS rows from the other rows, setting values that end up being zero to 0.19

dummy <- sheet7[sheet7$Stim == "UNS",]
# remove UNS data
sheet7 <- sheet7[sheet7$Stim !="UNS",] 

# match each dummy pid with sheet7 pid
dummy2 <- matrix(NA, nrow = nrow(sheet7), ncol = ncol(sheet7))
dummy2 <- as.data.frame(dummy2)
for(i in 1:nrow(sheet7)){
  n = which(sheet7$PID[i]==dummy$PID)
  dummy2[i,] <- dummy[n,]
}

new <- sheet7[,5:52]-dummy2[,5:52]
# convention is to set all values <= 0 to 0.19
for(i in 1:nrow(new)){
  for(j in 1:ncol(new)) {
    if(new[i,j]<=0){
      new[i,j] <- 0.19
    }
    
  }
}
head(new)

old_sheet7 <- sheet7
sheet7 <- cbind(old_sheet7$PID, old_sheet7$Stim, new)

head(sheet7)


# now we need decompose into smaller datasets, one for each PID
sheet7_a <- sheet7[sheet7$`old_sheet7$Stim` == "M.tb",]
sheet7_b <- sheet7[sheet7$`old_sheet7$Stim` == "LPS",]
sheet7_c <- sheet7[sheet7$`old_sheet7$Stim` == "HKCA",]

nrow(sheet7) == nrow(sheet7_a) + nrow(sheet7_b) + nrow(sheet7_c) 

# remove Stim column and rename variables
sheet7_a <- sheet7_a[,-2]
names(sheet7_a) <- c("PID", "M.tbUNS_PBMC_SCD40L", "M.tbUNS_PBMC_EGF", "M.tbUNS_PBMC_CCL11", "M.tbUNS_PBMC_FGF2",
                     "M.tbUNS_PBMC_FLT3L", "M.tbUNS_PBMC_CX3CL1", "M.tbUNS_PBMC_GCSF", "M.tbUNS_PBMC_GMCSF", "M.tbUNS_PBMC_GROA",
                     "M.tbUNS_PBMC_IFNA2", "M.tbUNS_PBMC_IFNG", "M.tbUNS_PBMC_IL1A", "M.tbUNS_PBMC_IL1B", "M.tbUNS_PBMC_IL1RA",
                     "M.tbUNS_PBMC_IL2", "M.tbUNS_PBMC_IL3", "M.tbUNS_PBMC_IL4", "M.tbUNS_PBMC_IL5", "M.tbUNS_PBMC_IL6", "M.tbUNS_PBMC_IL7", 
                     "M.tbUNS_PBMC_IL8", "M.tbUNS_PBMC_IL9", "M.tbUNS_PBMC_IL10", "M.tbUNS_PBMC_IL12P40", "M.tbUNS_PBMC_IL12P70",
                     "M.tbUNS_PBMC_IL13", "M.tbUNS_PBMC_IL15", "M.tbUNS_PBMC_IL17A", "M.tbUNS_PBMC_IL17E/IL-25", "M.tbUNS_PBMC_IL17F",
                     "M.tbUNS_PBMC_IL18", "M.tbUNS_PBMC_IL22", "M.tbUNS_PBMC_IL27", "M.tbUNS_PBMC_IP10", "M.tbUNS_PBMC_CCL2",            
                     "M.tbUNS_PBMC_CCL7", "M.tbUNS_PBMC_MCSF", "M.tbUNS_PBMC_CCL22", "M.tbUNS_PBMC_CXCL9", "M.tbUNS_PBMC_CCL3",
                     "M.tbUNS_PBMC_CCL4", "M.tbUNS_PBMC_PDGFAA", "M.tbUNS_PBMC_PDGFAB", "M.tbUNS_PBMC_CCL5", "M.tbUNS_PBMC_TGFA",
                     "M.tbUNS_PBMC_TNFA", "M.tbUNS_PBMC_TNFB", "M.tbUNS_PBMC_VEGF" )

sheet7_b <- sheet7_b[,-2]
names(sheet7_b) <- c("PID", "LPSUNS_PBMC_SCD40L", "LPSUNS_PBMC_EGF", "LPSUNS_PBMC_CCL11", "LPSUNS_PBMC_FGF2",
                     "LPSUNS_PBMC_FLT3L", "LPSUNS_PBMC_CX3CL1", "LPSUNS_PBMC_GCSF", "LPSUNS_PBMC_GMCSF", "LPSUNS_PBMC_GROA",
                     "LPSUNS_PBMC_IFNA2", "LPSUNS_PBMC_IFNG", "LPSUNS_PBMC_IL1A", "LPSUNS_PBMC_IL1B", "LPSUNS_PBMC_IL1RA",
                     "LPSUNS_PBMC_IL2", "LPSUNS_PBMC_IL3", "LPSUNS_PBMC_IL4", "LPSUNS_PBMC_IL5", "LPSUNS_PBMC_IL6", "LPSUNS_PBMC_IL7", 
                     "LPSUNS_PBMC_IL8", "LPSUNS_PBMC_IL9", "LPSUNS_PBMC_IL10", "LPSUNS_PBMC_IL12P40", "LPSUNS_PBMC_IL12P70",
                     "LPSUNS_PBMC_IL13", "LPSUNS_PBMC_IL15", "LPSUNS_PBMC_IL17A", "LPSUNS_PBMC_IL17E/IL-25", "LPSUNS_PBMC_IL17F",
                     "LPSUNS_PBMC_IL18", "LPSUNS_PBMC_IL22", "LPSUNS_PBMC_IL27", "LPSUNS_PBMC_IP10", "LPSUNS_PBMC_CCL2",            
                     "LPSUNS_PBMC_CCL7", "LPSUNS_PBMC_MCSF", "LPSUNS_PBMC_CCL22", "LPSUNS_PBMC_CXCL9", "LPSUNS_PBMC_CCL3",
                     "LPSUNS_PBMC_CCL4", "LPSUNS_PBMC_PDGFAA", "LPSUNS_PBMC_PDGFAB", "LPSUNS_PBMC_CCL5", "LPSUNS_PBMC_TGFA",
                     "LPSUNS_PBMC_TNFA", "LPSUNS_PBMC_TNFB", "LPSUNS_PBMC_VEGF" )

sheet7_c <- sheet7_c[,-2]
names(sheet7_c) <- c("PID", "HKCAUNS_PBMC_SCD40L", "HKCAUNS_PBMC_EGF", "HKCAUNS_PBMC_CCL11", "HKCAUNS_PBMC_FGF2",
                     "MHKCAUNS_PBMC_FLT3L", "HKCAUNS_PBMC_CX3CL1", "HKCAUNS_PBMC_GCSF", "HKCAUNS_PBMC_GMCSF", "HKCAUNS_PBMC_GROA",
                     "HKCAUNS_PBMC_IFNA2", "HKCAUNS_PBMC_IFNG", "HKCAUNS_PBMC_IL1A", "HKCAUNS_PBMC_IL1B", "HKCAUNS_PBMC_IL1RA",
                     "HKCAUNS_PBMC_IL2", "HKCAUNS_PBMC_IL3", "HKCAUNS_PBMC_IL4", "HKCAUNS_PBMC_IL5", "HKCAUNS_PBMC_IL6", "HKCAUNS_PBMC_IL7", 
                     "HKCAUNS_PBMC_IL8", "HKCAUNS_PBMC_IL9", "HKCAUNS_PBMC_IL10", "HKCAUNS_PBMC_IL12P40", "HKCAUNS_PBMC_IL12P70",
                     "HKCAUNS_PBMC_IL13", "HKCAUNS_PBMC_IL15", "HKCAUNS_PBMC_IL17A", "HKCAUNS_PBMC_IL17E/IL-25", "HKCAUNS_PBMC_IL17F",
                     "HKCAUNS_PBMC_IL18", "HKCAUNS_PBMC_IL22", "HKCAUNS_PBMC_IL27", "HKCAUNS_PBMC_IP10", "HKCAUNS_PBMC_CCL2",            
                     "HKCAUNS_PBMC_CCL7", "HKCAUNS_PBMC_MCSF", "HKCAUNS_PBMC_CCL22", "HKCAUNS_PBMC_CXCL9", "HKCAUNS_PBMC_CCL3",
                     "HKCAUNS_PBMC_CCL4", "HKCAUNS_PBMC_PDGFAA", "HKCAUNS_PBMC_PDGFAB", "HKCAUNS_PBMC_CCL5", "HKCAUNS_PBMC_TGFA",
                     "HKCAUNS_PBMC_TNFA", "HKCAUNS_PBMC_TNFB", "HKCAUNS_PBMC_VEGF" )


dim(sheet7_a) # 27 obs
dim(sheet7_b) # 21 obs
dim(sheet7_c) # 12 obs

final_sheet7 <- merge(sheet7_a, sheet7_b, by="PID", all.x = T, all.y = T)
final_sheet7 <- merge(final_sheet7, sheet7_c, by="PID", all.x = T, all.y = T)

dim(final_sheet7) # only 27 observations
dim(na.omit(final_sheet7)) # only 12 observation without NAs!!!


#-------------------
# Merge all
#-------------------

final <- merge(sheet1, sheet2, by="PID", all.x = T, all.y = T)
final <- merge(final, sheet3, by="PID", all.x = T, all.y = T)
final <- merge(final, sheet6, by="PID", all.x = T, all.y = T)
head(final)  # without sheet 7 data

# there are a couple of rows which only have demographic data and NA for everything else
# Need to remove these
final[rowSums(is.na(final)) > 0,]
# can see that all rows containing NAs have NAs for all observations except demographic obs
# thus we can just use na.omit to get rid of these "empty" observations
final <- na.omit(final)

dim(final) # 64 observations 
# 47 variables 



# including sheet 7
full <- merge(final, final_sheet7, by="PID", all.x = T, all.y = T)
dim(full)
dim(na.omit(full)) # only 12 full observations, out of 65 !!!!! YIKES
# also 191 variables !!! Too many !

# Maybe if we just included M.tb stim data
full <- merge(final, sheet7_a, by="PID", all.x = T, all.y = T)
dim(full) # 95 variables now
dim(na.omit(full)) # only 26 observations not containing NAs! 
full[rowSums(is.na(full)) > 0,-47:-3] # removing all the columns from final sheet
# can see that all the rows that have NA have NAs for all the sheet7_a variables


#-------------------
# Filter Data
#-------------------
# remove unnecessary demographic information
final <- final[,-c(4,5,6)]

# remove all variables that have constant values for >= 2/3 of observations for both
# vaccinated and unvaccinated groups

unvac <- final[final$Status=="Unvaccinated",]
vac <- final[final$Status=="Vaccinated",]

nrow(unvac)  #22
nrow(vac)    #42

lim1 <- nrow(unvac)/2
lim2 <- nrow(vac)/2

k <- c()

for(i in 4:ncol(final)){
  n <- length(unique(unvac[,i]))
  m <- length(unique(vac[,i]))
  if(n<lim1 & m<lim2){
    k <- append(k, i)
  }
}

final <- final[,-k]

# remove five variables recommended because all have values very close to 100
final <- final[,-c(15,16,17,18,19)]

#-------------------
# Save data
#-------------------


write.csv(final,"mva_assignment_final_data.csv", row.names = FALSE)

write.csv(full, "mva_assignment_FULL_data.csv", row.names = FALSE)



