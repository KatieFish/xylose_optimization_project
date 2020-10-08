
#10-5-20

#read in ortholog cds matrices
xyl1<- read.delim("~/xylose_optimization_project/data/orthologs/cds/xyl1_cds_MSA.tsv",
                  stringsAsFactors=FALSE)
xyl2<- read.delim("~/xylose_optimization_project/data/orthologs/cds/xyl2_cds_MSA.tsv", 
                  stringsAsFactors=FALSE)
xyl3<- read.delim("~/xylose_optimization_project/data/orthologs/cds/xyl3_cds_MSA.tsv", 
                  stringsAsFactors=FALSE)
tal1<- read.delim("~/xylose_optimization_project/data/orthologs/cds/tal1_cds_MSA.tsv", 
                  stringsAsFactors=FALSE)
tkl1<- read.delim("~/xylose_optimization_project/data/orthologs/cds/tkl1_cds_MSA.tsv", 
                  stringsAsFactors=FALSE)

###NOTE - as of 10/5/20 - the ortholog search for TKL genes looks like it was overly conservative, 
#or perhaps not sensitive enough. I am switching to using the orthoMCL cluster containing Scer tkl1
#from Shen et al. 2018 Cell. For now we can use the old file, and I'll push the new file
#ASAP. 


#remove identical sequences using the unique() command
#then, note which taxa remain duplicated with which() command
xyl1<-unique(xyl1)
xyl1_dups<-xyl1[which(duplicated(xyl1$V1)), 1]
xyl2<-unique(xyl2)
xyl2_dups<-xyl2[which(duplicated(xyl2$V1)), 1]
xyl3<-unique(xyl3)
xyl3_dups<-xyl3[which(duplicated(xyl3$V1)), 1]
tal1<-unique(tal1)
tal1_dups<-tal1[which(duplicated(tal1$V1)), 1]
tkl1<-unique(tkl1)
tkl1_dups<-tkl1[which(duplicated(tkl1$V1)), 1]


#make a master list of all taxa for all genes
all_taxa<-append(xyl1$V1, xyl2$V1)
all_taxa<-append(all_taxa, xyl3$V1)
all_taxa<-append(all_taxa, tal1$V1)
all_taxa<-append(all_taxa, tkl1$V1)
all_taxa<-unique(all_taxa)

#make a matrix of taxa x gene ID
##ex. 
######## xyl1 ##### xyl2 ##### xyl3 #####
#taxa1    1           0         1
######
#taxa2    0           1         1

xylose_utilization_gene_presence<- data.frame(all_taxa)
xylose_utilization_gene_presence$xyl1<-0
xylose_utilization_gene_presence$xyl2<-0
xylose_utilization_gene_presence$xyl3<-0
xylose_utilization_gene_presence$tkl1<-0
xylose_utilization_gene_presence$tal1<-0

############ mess around space below

for (i in 1:10){
  print(i+2)
}

if (1==1){
  print("yay")
}
############
# fill in xyl1 genes in the xylose_utilization_gene_presence table

xylose_utilization_gene_presence[,1]<- as.character(xylose_utilization_gene_presence[,1])

for (i in 1:nrow(xylose_utilization_gene_presence)){
  if(xylose_utilization_gene_presence[i,1] %in% xyl1[,1]){
    xylose_utilization_gene_presence[i, 2] <- 1
  }
}












# ignore for now
wi_vals<-read.delim("~/xylose_optimization_project/data/labella_et_al/wi_values.txt",
                    stringsAsFactors=FALSE)
