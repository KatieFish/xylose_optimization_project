
require("seqinr")

alnmnt<-read.alignment("~/xylose_optimization_project/spring_2021/subgroup_consensus/xyl3/xyl3_subgroup.fasta", format="fasta")
#also reading in fasta
fs<- read.fasta("~/xylose_optimization_project/spring_2021/subgroup_consensus/xyl3/xyl3_subgroup.fasta")

bp<-max(getLength(fs))

##conversion of alnmnt into matrix

alnmnt_matrix<-matrix(nrow=length(alnmnt[[2]]), ncol=bp+1)
alnmnt_matrix[,1]<-unlist(alnmnt[[2]])
for (i in 1:length(alnmnt[[2]])){
  j<-(getLength(fs)[i])
  alnmnt_matrix[i,(2:(j+1))]<-unlist(strsplit(alnmnt$seq[[i]], split=""))
}

#######fix taxa names
x<-alnmnt_matrix[,1]
taxa_IDs<-x
for (i in 1:length(x)){
  if (grepl("yHMP", x[i]) | grepl("yHAB", x[i])){
   name<-paste(strsplit(x[i], "_")[[1]][2], strsplit(x[i], "_")[[1]][3], sep=" ")
  }
  else{
    name<-paste(strsplit(x[i], "_")[[1]][1], strsplit(x[i], "_")[[1]][2], sep=" ")
  }
taxa_IDs[i]<-name
}
####
alnmnt_matrix[,1]<-taxa_IDs

write.table(alnmnt_matrix, "~/xylose_optimization_project/spring_2021/subgroup_consensus/xyl3/xyl3.subgroup.matrix.txt",
            sep="\t", quote=FALSE, row.names=FALSE)


