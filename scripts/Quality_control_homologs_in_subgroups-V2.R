
xyl1<-read.delim("~/xylose_optimization_project/spring_2021/subgroup_consensus/xyl1/xyl1.subgroup.matrix.txt",
                         header=TRUE, stringsAsFactors=FALSE)

xyl2<-read.delim("~/xylose_optimization_project/spring_2021/subgroup_consensus/xyl2/xyl2.subgroup.matrix.txt",
                 header=TRUE, stringsAsFactors=FALSE)

xyl3<-read.delim("~/xylose_optimization_project/spring_2021/subgroup_consensus/xyl3/xyl3.subgroup.matrix.txt",
                 header=TRUE, stringsAsFactors=FALSE)



#fix names in sequence mat
taxon_name_switches <- read.delim("~/xylose_optimization_project/spring_2021/data/taxon_name_switches.txt", stringsAsFactors=FALSE)

#xyl1
xyl1->sequence_mat
sequence_mat$V1<-tolower(sequence_mat$V1)

for(i in 1:nrow(sequence_mat)){
  if(sequence_mat$V1[i] %in% taxon_name_switches$sequences_name){
    sequence_mat$V1[i]<-taxon_name_switches$growth_data_name[
     which(taxon_name_switches$sequences_name==sequence_mat$V1[i])]
  }
}
sequence_mat->xyl1

#xyl2
xyl2->sequence_mat
sequence_mat$V1<-tolower(sequence_mat$V1)

for(i in 1:nrow(sequence_mat)){
  if(sequence_mat$V1[i] %in% taxon_name_switches$sequences_name){
    sequence_mat$V1[i]<-taxon_name_switches$growth_data_name[
      which(taxon_name_switches$sequences_name==sequence_mat$V1[i])]
  }
}
sequence_mat->xyl2

#xyl3
xyl3->sequence_mat
sequence_mat$V1<-tolower(sequence_mat$V1)

for(i in 1:nrow(sequence_mat)){
  if(sequence_mat$V1[i] %in% taxon_name_switches$sequences_name){
    sequence_mat$V1[i]<-taxon_name_switches$growth_data_name[
      which(taxon_name_switches$sequences_name==sequence_mat$V1[i])]
  }
}
sequence_mat->xyl3



#importing growth data
growth_data<-read.delim("~/xylose_optimization_project/spring_2021/data/growth_rate_master_df.txt", stringsAsFactors=FALSE, header=TRUE)

#changing sequence matrix taxa names to lower case
xyl1$V1<-tolower(xyl1$V1)
xyl2$V1<-tolower(xyl2$V1)
xyl3$V1<-tolower(xyl3$V1)
#pulling out species with nonzero growth
growers<-growth_data$all_taxa[which(growth_data$Growth.Rate>0)]
#determining which species with nonzero growth are not represented in the xyl1 dataset
xyl1_missing_taxa<-growers[which(!growers %in% xyl1$V1)]
xyl2_missing_taxa<-growers[which(!growers %in% xyl2$V1)]
xyl3_missing_taxa<-growers[which(!growers %in% xyl3$V1)]

#only missing is candida sojae from xyl1 list. I am blasting the sojae genome with 
#the tropicalis homolog. 

#There is not really a great match here, probably because cadida sojae's genome
#is incomplete. #Best hit - identity at 118/316 (37.34% identity at a.a. level)

#####################################################



