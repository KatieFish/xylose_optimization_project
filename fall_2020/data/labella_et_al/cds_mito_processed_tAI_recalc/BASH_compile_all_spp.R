#!/usr/bin/env Rscript

#user input 1 is the raw data table and user input 2 is the name of the chemical

user_input=commandArgs(trailingOnly = TRUE)


df<- read.delim(user_input[1], header=FALSE, stringsAsFactors=FALSE)
sp<- user_input[2]

df$taxa<-sp

ndf<-read.delim("running_table.txt")

ndf1<-rbind(ndf, df)

write.table(ndf1, "running_table.txt", sep="\t", row.names=FALSE)