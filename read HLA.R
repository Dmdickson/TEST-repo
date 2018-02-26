
#  Purpose: to read into R the data files for the IDEA correlates project: HLA data

###############################################
# Original Author :Dorothy Dickson            #
#        Edited by                            #
# Last updated: 19 January 2018               #  
#                                             #
###############################################

# HERE IS A TEST CHANGE to update the data file name

setwd("C:/Users/dmdickso/IDEA")
path="C:/Users/dmdickso/IDEA/"
###### change the above path once if needed to fit the current location of datafiles   ########


hlafile="UVM_GATES_CHALLENGE_HLA.csv"
pathhla<-paste0(path,hlafile, collapse='')

library(tidyr)
library(knitr)
library(plyr)
library(dplyr)
library(reshape2)

hla<-read.csv(pathhla, header=T, stringsAsFactors = T,
                    col.names=c("subject",
                                "altid",
                                "liai",
                                "hla",
                                "a1",
                                "a2",
                                "b1",
                                "b2",
                                "c1",
                                "c2",
                                "dpa1",
                                "dpa2",
                                "dpb11",
                                "dpb12",
                                "dqa11",
                                "dqa12",
                                "dqb11",
                                "dqb12",
                                "drb11",
                                "drb12",
                                "drb31",
                                "drb32",
                                "drb41",
                                "drb42",
                                "drb51",
                                "drb52"))

### tabulate allele frequencies
### default options ignore NA missing
#HLAtable<-table(hla$a1)
############################################################################# 
#flip the dataset HLA to make occurring alleles the variables and 0/1 response per subject depending on whether the allele was present for HLA_A
#HLA_a1<-dcast(hla, subject~a1,length)
#HLA_a2<-dcast(hla, subject~a2,length)

### alleles in common within A_allele1 and A_allele2
#common<-intersect (names(HLA_a1), names(HLA_a2))
## drop subject and Var.2 from the list of common variables names (interested in allele names only)
# common<-common [-c(1:2)]
# HLA_A<-merge(HLA_a1, HLA_a2, by=c("subject"))
# for (variable in common){
#   # Create a summed variable
#   HLA_A[[variable]] <- HLA_A %>% select(starts_with(paste0(variable,"."))) %>% rowSums()
#   # Delete columns with .x and .y suffixes
#   HLA_A <- HLA_A %>% select(-one_of(c(paste0(variable,".x"), paste0(variable,".y"))))
# }
# ### drop var.2.x and var.2.y
# drops <- c("Var.2.x","Var.2.y")
# HLA_A<-HLA_A[ , !(names(HLA_A) %in% drops)]

#####################################################################################
####### two functions = one to accommmodate blank alleles for HLA_A, _B, _C
#######                 one for alleles that have no subjects with blanks for allele1 and allele2
#### function for no blanks
flipandsum<-function (dframe,var1,var2){
  df1<-dcast(dframe, subject~var1,length)
  df2<-dcast(dframe, subject~var2,length)
  
  common<-intersect (names(df1), names(df2))
  common<-common [-c(1)]
  df<-merge(df1, df2, by=c("subject"))
  for (variable in common){
    df[[variable]] <- df %>% select(starts_with(paste0(variable,"."))) %>% rowSums()
    df <- df %>% select(-one_of(c(paste0(variable,".x"), paste0(variable,".y"))))
  }
  
  return (df)
  }

#### function to accommodate blanks
flipandsum2<-function (dframe,var1,var2){
  df1<-dcast(dframe, subject~var1,length)
  df2<-dcast(dframe, subject~var2,length)
  
  common<-intersect (names(df1), names(df2))
  ##### this line different in the two functions, need to get rid of var.2 if there are blanks
  common<-common [-c(1:2)]
  df<-merge(df1, df2, by=c("subject"))
  for (variable in common){
    df[[variable]] <- df %>% select(starts_with(paste0(variable,"."))) %>% rowSums()
    df <- df %>% select(-one_of(c(paste0(variable,".x"), paste0(variable,".y"))))
  }
  drops <- c("Var.2.x","Var.2.y")
  df<-df[ , !(names(df) %in% drops)]
  return (df)
}
#################################################################################

HLA_A<-flipandsum2(hla, hla$a1, hla$a2)
HLA_B<-flipandsum2(hla, hla$b1, hla$b2)
HLA_C<-flipandsum(hla, hla$c1, hla$c2)
HLA_DPB1<-flipandsum(hla, hla$dpb11, hla$dpb12)
HLA_DRB1<-flipandsum(hla, hla$drb11, hla$drb12)
HLA_DQA1<-flipandsum2(hla, hla$dqa11, hla$dqa12)
HLA_DQB1<-flipandsum2(hla, hla$dqb11, hla$dqb12)
HLA_DRB3<-flipandsum2(hla, hla$drb31, hla$drb32)
HLA_DRB4<-flipandsum2(hla, hla$drb41, hla$drb42)
HLA_DRB5<-flipandsum2(hla, hla$drb51, hla$drb52)

# merge all separate HLA-type dataframes then remove individuals
HLAcomplete<-Reduce(function(x,y) merge(x,y, by="subject"), list(HLA_A,HLA_B, HLA_C, HLA_DPB1, HLA_DQA1, HLA_DQB1, HLA_DRB1, HLA_DRB3, HLA_DRB4, HLA_DRB5))
rm(HLA_A, HLA_B, HLA_C,HLA_DPB1, HLA_DRB1, HLA_DQA1, HLA_DQB1, HLA_DRB3,HLA_DRB4, HLA_DRB5)





