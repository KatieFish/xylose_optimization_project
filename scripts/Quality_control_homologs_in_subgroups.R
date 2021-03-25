
sequence_mat<-read.delim("~/xylose_optimization_project/spring_2021/subgroup_consensus/xyl1/xyl1.subgroup.matrix.txt",
                         header=FALSE, stringsAsFactors=FALSE)

#how many sequences do not begin with methionine? 
length(which(!sequence_mat[,2] == "m"))

#let's get rid of these for the time being
sequence_mat<-sequence_mat[which(sequence_mat[,2]=="m"),]

#how long are the sequences?
seqlengths<-data.frame(sequence_mat[,1])
seqlengths$aa_length<-NA
for(i in 1:nrow(sequence_mat)){
  ndf<-as.list(sequence_mat[i, 2:ncol(sequence_mat)])
  if("*" %in% ndf){
  length<-as.numeric(which(ndf=="*")-1)
  }
  else{
    length<-(length(ndf)-length(which(is.na(ndf))))
  }
  seqlengths$aa_length[i]<-length  
}

hist(seqlengths$aa_length, breaks=75, ylab="protein length (aa)", col="lightblue")

#based on this histogram, should we eliminate any sequences for being too small or too large?



