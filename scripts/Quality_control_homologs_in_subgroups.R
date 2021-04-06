
sequence_mat<-read.delim("~/xylose_optimization_project/spring_2021/subgroup_consensus/xyl1/xyl1.subgroup.matrix.txt",
                         header=TRUE, stringsAsFactors=FALSE)

#fix names in sequence mat
taxon_name_switches <- read.delim("~/xylose_optimization_project/spring_2021/data/taxon_name_switches.txt", stringsAsFactors=FALSE)





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
#for xyl1, I'm setting cuttoff at 275 - 400 aa

hist(seqlengths$aa_length, breaks=75, ylab="protein length (aa)", col="lightblue")
abline(v=275, lty=2)
abline(v=400, lty=2)

#retaining sequence between thresholds
sequence_mat<-sequence_mat[which(seqlengths$aa_length>275 & seqlengths$aa_length<400),]

#importing growth data
growth_data<-read.delim("~/xylose_optimization_project/spring_2021/data/growth_rate_master_df.txt", stringsAsFactors=FALSE, header=TRUE)

#changing sequence matrix taxa names to lower case
sequence_mat$V1<-tolower(sequence_mat$V1)
#pulling out species with nonzero growth
growers<-growth_data$all_taxa[which(growth_data$Growth.Rate>0)]
#determining which species with nonzero growth are not represented in the xyl1 dataset
missing_taxa<-growers[which(!growers %in% sequence_mat$V1)]
#are they ones that we got rid of during quality control? 
which(missing_taxa %in% seqlengths$sequence_mat...1.)
#none were in the sequences we got rid of for size reasons. 
#some of these might be mispellings

for(i in 1:length(missing_taxa)){
  genus<-strsplit(missing_taxa[i], split = " ")[[1]][1]
  species<-strsplit(missing_taxa[i], split = " ")[[1]][2]
  if(length(which(grepl(genus, sequence_mat$V1)))>0){
    print(missing_taxa[i])
    print(sequence_mat$V1[which(grepl(genus, sequence_mat$V1))])
  }
  if(length(which(grepl(species, sequence_mat$V1)))>0){
    print(missing_taxa[i])
    print(sequence_mat$V1[which(grepl(species, sequence_mat$V1))])
  }
}
