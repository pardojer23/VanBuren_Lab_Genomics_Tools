# Title     : DESeq2
# Objective : conduct differential expressiom analysis
# Created by: Jeremy Pardo
# Created on: 2019-04-02
args = commandArgs(TRUE)
outputDir = args[1]
conditionString = args[2]
# set working directory
setwd(outputDir)
# load required packages
# load required packages
library(readr)
library(tximport)
library(DESeq2)
library(dplyr)
#read in sample table
sampleTable <-read_delim("SampleTable.txt",delim="\t")
# load tximport data
load("txi.RData")
#
DEseq_obj <- DESeqDataSetFromTximport(txi=txi,colData = sampleTable,design = ~ Condition)
dds_obj <- DESeq(DEseq_obj)
##
##
get_sig_df = function(df,cont){
  df = mutate(df,GeneID = rownames(df),Contrast = cont)
  df = df[which(df$padj<0.05),]
  return(df)
}
print(conditionString)
compare_list = unlist(strsplit(conditionString,","))
df_list = list()
print(compare_list)
for (i in 1:length(compare_list)){
  print(compare_list[i])
  condition_list = unlist(strsplit(compare_list[i],"-"))
  contrast = paste0(condition_list[1],"v",condition_list[2])
  print(contrast)
  df_list[[i]] = assign(paste0(contrast,"_sig"),get_sig_df(assign(contrast,as.data.frame(results(dds_obj,contrast=c("Condition",condition_list[2],condition_list[1]),alpha=0.05,pAdjustMethod = "fdr"))),contrast))
}
DE_Gene_df = bind_rows(df_list)
write_delim(DE_Gene_df,"DE_Gene_df.txt",delim="\t")
##
save.image("DESeq2.RData")

