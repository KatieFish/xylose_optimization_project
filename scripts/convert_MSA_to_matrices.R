
#Code adjusted by KJF 8-24-21

#No need to fix names - Rishitha has already done that in the input file.

#No need to replace gaps with NAs - no longer worried about codon occupancy.


require("seqinr")

alnmnt<-read.alignment("~/Xyl_project_Fall_2021/xyl3/xyl3_cds.fasta", format="fasta")
#also reading in fasta
fs<- read.fasta("~/Xyl_project_Fall_2021/xyl3/xyl3_cds.fasta")

bp<-max(getLength(fs))

##conversion of alnmnt into matrix

alnmnt_matrix<-matrix(nrow=length(alnmnt[[2]]), ncol=bp+3)
alnmnt_matrix[,1]<-unlist(alnmnt[[2]])
for (i in 1:length(alnmnt[[2]])){
  j<-(getLength(fs)[i])
  alnmnt_matrix[i,(2:(j+1))]<-unlist(strsplit(alnmnt$seq[[i]], split=""))
}

########
###
#Fix names - one column for genus species and one column for copy #
##


names<-alnmnt_matrix[,1]
taxon<-names
copy_no.<-names
for(i in 1:length(names)){
  x<-strsplit(names[i], split = "_")
  taxon[i]<-paste(x[[1]][1:length(x[[1]])-1], collapse = " ")
  copy_no.[i]<-as.integer(x[[1]][length(x[[1]])]) 
}
  
alnmnt_matrix[,bp+2]<-taxon
alnmnt_matrix[,bp+3]<-copy_no.

alnmnt_df<-data.frame(alnmnt_matrix)
alnmnt_df<-alnmnt_df[c((bp+2), (bp+3), 2:(bp+1))]

colnames(alnmnt_df)[1:2]<-c("taxon", "copy.no")

write.table(alnmnt_df, "~/Xyl_project_Fall_2021/xyl3/xyl3_cds_matrix_forR.txt",
            sep="\t", quote=FALSE, row.names=FALSE)


