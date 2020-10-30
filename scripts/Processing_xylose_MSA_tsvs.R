
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


# fill in xyl1 genes in the xylose_utilization_gene_presence table

xylose_utilization_gene_presence[,1]<- as.character(xylose_utilization_gene_presence[,1])

for (i in 1:nrow(xylose_utilization_gene_presence)){
  if(xylose_utilization_gene_presence[i,1] %in% xyl1[,1]){
    xylose_utilization_gene_presence[i, 2] <- 1
  }
}
for (i in 1:nrow(xylose_utilization_gene_presence)){
  if(xylose_utilization_gene_presence[i,1] %in% xyl2[,1]){
    xylose_utilization_gene_presence[i, 3] <- 1
  }
}
for (i in 1:nrow(xylose_utilization_gene_presence)){
  if(xylose_utilization_gene_presence[i,1] %in% xyl3[,1]){
    xylose_utilization_gene_presence[i, 4] <- 1
  }
}
for (i in 1:nrow(xylose_utilization_gene_presence)){
  if(xylose_utilization_gene_presence[i,1] %in% tkl1[,1]){
    xylose_utilization_gene_presence[i, 5] <- 1
  }
}
for (i in 1:nrow(xylose_utilization_gene_presence)){
  if(xylose_utilization_gene_presence[i,1] %in% tal1[,1]){
    xylose_utilization_gene_presence[i, 6] <- 1
  }
}






# ignore for now
wi_vals<-read.delim("~/xylose_optimization_project/data/labella_et_al/wi_values.txt",
                    stringsAsFactors=FALSE)


# are there any spp. in our lists that are not in abbey's wi values?
which(!xylose_utilization_gene_presence[,1] %in% wi_vals[,1])
xylose_utilization_gene_presence[which(!xylose_utilization_gene_presence[,1] %in% wi_vals[,1]), 1]


# are there any spp. in abbey's list that are not in ours? 
which(!wi_vals[,1] %in% xylose_utilization_gene_presence[,1])->x
wi_vals[x,1]

### start a new section of code! ;) In this section we will
#write code to calculate the mean wi value for each gene (estAI)
require(EnvStats)

stAI_dataframe<-xylose_utilization_gene_presence
stAI_dataframe$xyl1<-NA
stAI_dataframe$xyl2<-NA
stAI_dataframe$xyl3<-NA
stAI_dataframe$tal1<-NA
stAI_dataframe$tkl1<-NA
#make new columns for 2nd copy of genes
stAI_dataframe$xyl1.2<-NA
stAI_dataframe$xyl2.2<-NA
stAI_dataframe$xyl3.2<-NA
stAI_dataframe$tal1.2<-NA
stAI_dataframe$tkl1.2<-NA


#iterate by row of the xylose_ut_df and check for values of 1
for(i in 1:nrow(xylose_utilization_gene_presence)){ 
  #check if xyl1 (column 2)==1
  if(xylose_utilization_gene_presence[i,2]==1){
    # create x variable that refers to row number in xyl1 of species
    x<- which(xyl1[,1] == xylose_utilization_gene_presence[i,1])
    # set up if statement to check if x has 2 numbers
    y<- which(wi_vals[,1] == xylose_utilization_gene_presence[i,1])
    #identify the gene sequence for the respective gene:
    gene_x<-as.character(xyl1[x[1],(5:ncol(xyl1))])
    #below if statement checks to see if NAs are at end of gene and then removes them
    if(is.na(xyl1[x[1], ncol(xyl1)])){
    z<-which(is.na(gene_x))
    gene_x<-gene_x[1:(z[1]-1)]
    }
    #create a vector to store the wi values
    gene_x_wi_vals<-numeric()
    for(j in seq(1,length(gene_x), 3)){
      codon<-paste(gene_x[j], gene_x[j+1], gene_x[j+2], sep="")
      codon_column<-which(grepl(codon, colnames(wi_vals), ignore.case=TRUE))
      codon_wi_val<-as.numeric(wi_vals[y,codon_column])
      gene_x_wi_vals<-append(gene_x_wi_vals, codon_wi_val)
    }
    stAI<-geoMean(gene_x_wi_vals)
    stAI_dataframe[i,2]<-stAI
    # add an if statement for if there's a second copy
    if(length(x)>1){
      gene_x<-as.character(xyl1[x[2],(5:ncol(xyl1))])
      if(is.na(xyl1[x[2], ncol(xyl1)])){
        z<-which(is.na(gene_x))
        gene_x<-gene_x[1:(z[1]-1)]
      }
      #create a vector to store the wi values
      gene_x_wi_vals<-numeric()
      for(j in seq(1,length(gene_x), 3)){
        codon<-paste(gene_x[j], gene_x[j+1], gene_x[j+2], sep="")
        codon_column<-which(grepl(codon, colnames(wi_vals), ignore.case=TRUE))
        codon_wi_val<-as.numeric(wi_vals[y,codon_column])
        gene_x_wi_vals<-append(gene_x_wi_vals, codon_wi_val)
      }
      stAI<-geoMean(gene_x_wi_vals)
      stAI_dataframe[i,7]<-stAI
    }   
  }
}



### copy 2 of loop
for(i in 1:nrow(xylose_utilization_gene_presence)){ 
  #check if xyl1 (column 2)==1
  if(xylose_utilization_gene_presence[i,2]==1){
    # create x variable that refers to row number in xyl1 of species
    x<- which(xyl1[,1] == xylose_utilization_gene_presence[i,1])
    # set up if statement to check if x has 2 numbers
    y<- which(wi_vals[,1] == xylose_utilization_gene_presence[i,1])
    #identify the gene sequence for the respective gene:
    gene_x<-as.character(xyl1[x[1],(5:ncol(xyl1))])
    #below if statement checks to see if NAs are at end of gene and then removes them
    if(is.na(xyl1[x[1], ncol(xyl1)])){
      z<-which(is.na(gene_x))
      gene_x<-gene_x[1:(z[1]-1)]
    }
    #create a vector to store the wi values
    gene_x_wi_vals<-numeric()
    for(j in seq(1,length(gene_x), 3)){
      codon<-paste(gene_x[j], gene_x[j+1], gene_x[j+2], sep="")
      codon_column<-which(grepl(codon, colnames(wi_vals), ignore.case=TRUE))
      codon_wi_val<-as.numeric(wi_vals[y,codon_column])
      gene_x_wi_vals<-append(gene_x_wi_vals, codon_wi_val)
    }
    stAI<-geoMean(gene_x_wi_vals)
    stAI_dataframe[i,2]<-stAI
    # add an if statement for if there's a second copy
    if(length(x)>1){
      gene_x<-as.character(xyl1[x[2],(5:ncol(xyl1))])
      if(is.na(xyl1[x[2], ncol(xyl1)])){
        z<-which(is.na(gene_x))
        gene_x<-gene_x[1:(z[1]-1)]
      }
      #create a vector to store the wi values
      gene_x_wi_vals<-numeric()
      for(j in seq(1,length(gene_x), 3)){
        codon<-paste(gene_x[j], gene_x[j+1], gene_x[j+2], sep="")
        codon_column<-which(grepl(codon, colnames(wi_vals), ignore.case=TRUE))
        codon_wi_val<-as.numeric(wi_vals[y,codon_column])
        gene_x_wi_vals<-append(gene_x_wi_vals, codon_wi_val)
      }
      stAI<-geoMean(gene_x_wi_vals)
      stAI_dataframe[i,7]<-stAI
    }   
  }
}


######## Taylor & Rishitha finished all orthos for stAI#####

#import genome-wide stAI values

try<-read.delim("/Users/katiefisher/xylose_optimization_project/data/labella_et_al/cds_mito_processed_tAI_recalc/running_table.txt")
try$taxa<-as.char

##have to convert names to names in use in our dfs
x<-try[,3]
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




