
#Code adjusted by KJF 4-21-21 

#names of sequences from ORFfinder different from orginal maker annotations, so I needed
#to adjust the way I cleaned up the names. I will need to re-assure that these correspond 
#to the names used in the Labella et al data we use. 

#Since the new approach is based on a codon-aware alignment, I replaced all
#gaps with NA values. 

#I made all taxon names lower case to match Opulente data and Labella data. 



require("seqinr")

alnmnt<-read.alignment("~/xylose_optimization_project/spring_2021/subgroup_consensus/xyl3/subgroup/xyl3_subgroup-codon_aligned_cds-cleaned.fasta", format="fasta")
#also reading in fasta
fs<- read.fasta("~/xylose_optimization_project/spring_2021/subgroup_consensus/xyl3/subgroup/xyl3_subgroup-codon_aligned_cds-cleaned.fasta")

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
  if(grepl("\\d",strsplit(x[i], "_")[[1]][3])){
    name<-paste(strsplit(x[i], "_")[[1]][1], strsplit(x[i], "_")[[1]][2], sep=" ")
  } else if (grepl("_var_", x[i])){
   name<-paste(strsplit(x[i], "_")[[1]][1], strsplit(x[i], "_")[[1]][2],
               strsplit(x[i], "_")[[1]][3], strsplit(x[i], "_")[[1]][4], sep=" ")
  } else {
    xchar<-which(strsplit(x[i], "_")[[1]]=="x")
    name<-paste(strsplit(x[i], "_")[[1]][1:(xchar-1)], collapse=" ")
  }   
taxa_IDs[i]<-name
}
####
alnmnt_matrix[,1]<-taxa_IDs

########
alnmnt_matrix[ alnmnt_matrix == "-" ] <- NA

alnmnt_matrix[,1]<-tolower(alnmnt_matrix[,1])

write.table(alnmnt_matrix, "~/xylose_optimization_project/spring_2021/subgroup_consensus/xyl3/subgroup/xyl3_subgroup-codon_aligned_cds-cleaned.txt",
            sep="\t", quote=FALSE, row.names=FALSE)


