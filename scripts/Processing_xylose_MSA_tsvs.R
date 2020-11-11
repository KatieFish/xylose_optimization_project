
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
#xyl1_dups<-xyl1[which(duplicated(xyl1$V1)), 1]
xyl2<-unique(xyl2)
#xyl2_dups<-xyl2[which(duplicated(xyl2$V1)), 1]
xyl3<-unique(xyl3)
#xyl3_dups<-xyl3[which(duplicated(xyl3$V1)), 1]
tal1<-unique(tal1)
#tal1_dups<-tal1[which(duplicated(tal1$V1)), 1]
tkl1<-unique(tkl1)
#tkl1_dups<-tkl1[which(duplicated(tkl1$V1)), 1]


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

write.table(xylose_utilization_gene_presence, 
            "~/xylose_optimization_project/data/XYLpthwy_gene_presence_absence_matrix.txt", sep="\t", quote=FALSE, row.names=FALSE)





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
#start 
for(i in 1:nrow(xylose_utilization_gene_presence)){
  #check if gene is present for species
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
  if(xylose_utilization_gene_presence[i,3]==1){
    # create x variable that refers to row number in xyl1 of species
    x<- which(xyl2[,1] == xylose_utilization_gene_presence[i,1])
    # set up if statement to check if x has 2 numbers
    y<- which(wi_vals[,1] == xylose_utilization_gene_presence[i,1])
    #identify the gene sequence for the respective gene:
    gene_x<-as.character(xyl2[x[1],(5:ncol(xyl2))])
    #below if statement checks to see if NAs are at end of gene and then removes them
    if(is.na(xyl2[x[1], ncol(xyl2)])){
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
    stAI_dataframe[i,3]<-stAI
    # add an if statement for if there's a second copy
    if(length(x)>1){
      gene_x<-as.character(xyl2[x[2],(5:ncol(xyl2))])
      if(is.na(xyl2[x[2], ncol(xyl2)])){
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
      stAI_dataframe[i,8]<-stAI
    }   
  }
}
#copy 3 of loop - xyl3
for(i in 1:nrow(xylose_utilization_gene_presence)){ 
  #check if xyl1 (column 2)==1
  if(xylose_utilization_gene_presence[i,4]==1){
    # create x variable that refers to row number in xyl1 of species
    x<- which(xyl3[,1] == xylose_utilization_gene_presence[i,1])
    # set up if statement to check if x has 2 numbers
    y<- which(wi_vals[,1] == xylose_utilization_gene_presence[i,1])
    #identify the gene sequence for the respective gene:
    gene_x<-as.character(xyl3[x[1],(5:ncol(xyl3))])
    #below if statement checks to see if NAs are at end of gene and then removes them
    if(is.na(xyl3[x[1], ncol(xyl3)])){
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
    stAI_dataframe[i,4]<-stAI
    # add an if statement for if there's a second copy
    if(length(x)>1){
      gene_x<-as.character(xyl3[x[2],(5:ncol(xyl3))])
      if(is.na(xyl3[x[2], ncol(xyl3)])){
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
      stAI_dataframe[i,9]<-stAI
    }   
  }
}

#copy 4 of loop - tkl1
for(i in 1:nrow(xylose_utilization_gene_presence)){ 
  #check if xyl1 (column 2)==1
  if(xylose_utilization_gene_presence[i,5]==1){
    # create x variable that refers to row number in xyl1 of species
    x<- which(tkl1[,1] == xylose_utilization_gene_presence[i,1])
    # set up if statement to check if x has 2 numbers
    y<- which(wi_vals[,1] == xylose_utilization_gene_presence[i,1])
    #identify the gene sequence for the respective gene:
    gene_x<-as.character(tkl1[x[1],(5:ncol(tkl1))])
    #below if statement checks to see if NAs are at end of gene and then removes them
    if(is.na(tkl1[x[1], ncol(tkl1)])){
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
    stAI_dataframe[i,5]<-stAI
    # add an if statement for if there's a second copy
    if(length(x)>1){
      gene_x<-as.character(tkl1[x[2],(5:ncol(tkl1))])
      if(is.na(tkl1[x[2], ncol(tkl1)])){
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
      stAI_dataframe[i,10]<-stAI
    }   
  }
}

#copy 5 of loop - tal1
for(i in 1:nrow(xylose_utilization_gene_presence)){ 
  #check if xyl1 (column 2)==1
  if(xylose_utilization_gene_presence[i,6]==1){
    # create x variable that refers to row number in xyl1 of species
    x<- which(tal1[,1] == xylose_utilization_gene_presence[i,1])
    # set up if statement to check if x has 2 numbers
    y<- which(wi_vals[,1] == xylose_utilization_gene_presence[i,1])
    #identify the gene sequence for the respective gene:
    gene_x<-as.character(tal1[x[1],(5:ncol(tal1))])
    #below if statement checks to see if NAs are at end of gene and then removes them
    if(is.na(tal1[x[1], ncol(tal1)])){
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
    stAI_dataframe[i,6]<-stAI
    # add an if statement for if there's a second copy
    if(length(x)>1){
      gene_x<-as.character(tal1[x[2],(5:ncol(tal1))])
      if(is.na(tal1[x[2], ncol(tal1)])){
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
      stAI_dataframe[i,11]<-stAI
    }   
  }
}
rm(all_taxa, codon, codon_column, codon_wi_val, gene_x,
   gene_x_wi_vals,i, j, stAI, x, y, z)

#import genome-wide stAI values

genome_wide_tAI<-read.delim("~/xylose_optimization_project/data/labella_et_al/cds_mito_processed_tAI_recalc/genome_wide_tAI_all_spp.txt")
#after unzipping before import had to correct taxa names in file on command line using sed to replace 
#underscores with spaces. 
#have to fix taxa names
x<-as.character(genome_wide_tAI$taxa)
taxa_IDs<-x
for (i in 1:length(x)){
  if (grepl("yHMP", x[i]) | grepl("yHAB", x[i])){
    name<-paste(strsplit(x[i], " ")[[1]][2], strsplit(x[i], " ")[[1]][3], sep=" ")
  }
  else{
    name<-paste(strsplit(x[i], " ")[[1]][1], strsplit(x[i], " ")[[1]][2], sep=" ")
  }
  taxa_IDs[i]<-name
}
####
genome_wide_tAI$taxa<-taxa_IDs
###NAMES FIXED ABOVE
#Hash out and import in the future. 
write.table(genome_wide_tAI, "~/xylose_optimization_project/data/labella_et_al/cds_mito_processed_tAI_recalc/genome_wide_tAI_all_spp.txt", sep="\t", quote=FALSE, row.names=FALSE)

#Create empty datafrme to populate with estAI values
estAI_df<-stAI_dataframe
for(i in 2:ncol(estAI_df)){
  for(j in 1:nrow(estAI_df)){
    estAI_df[j,i]<-NA
  }
}

genome_wide_tAI$taxa<-as.character(genome_wide_tAI$taxa)
stAI_dataframe$all_taxa<-as.character(stAI_dataframe$all_taxa)

# iterate through stAI df and populate estAI df 
for(i in 2:ncol(stAI_dataframe)){
  for(j in 1:nrow(stAI_dataframe)){
  # grab the species tAI vals for all genes 
  if(!is.na(stAI_dataframe[j,i])){
  species_tAI_vals<-genome_wide_tAI[which(genome_wide_tAI$taxa==stAI_dataframe$all_taxa[j]), 2]
  ecdf_func<-ecdf(species_tAI_vals)
  estAI_df[j,i]<-ecdf_func(stAI_dataframe[j,i])
  }
  }
}

rm(ecdf_func, i, j, name, species_tAI_vals, taxa_IDs, x)

####writing tAI, estAI results as tables to import
write.table(stAI_dataframe, "~/xylose_optimization_project/data/spp_by_gene_tAI_vals.txt", sep="\t", quote=FALSE, row.name=FALSE)
write.table(estAI_df, "~/xylose_optimization_project/data/spp_by_gene_estAI_vals.txt", sep="\t", quote=FALSE, row.name=FALSE)



###Abbe's paper uses the empirical distribution function to find
#the percent of genes with lower tAI vals. Does this differ from 
#our division function? 

# create ecdf distribution function
# ecdf_func<-ecdf(species_tAI_vals)
# ecdf_func(stAI_dataframe[j,i])
# I did a couple by hand and the numbers are the exact same. 
# I guess there's no harm in doing EXACTLY what Abbe did, 
# I'll adjust the code


#Compare histograms of estAI vals between genes
quartz()
par(mar=c(2,2,2,2))
par(mfrow=c(2,3))
hist(estAI_df$xyl1)
hist(estAI_df$xyl2)
hist(estAI_df$xyl3)
hist(estAI_df$tkl1)
hist(estAI_df$tal1)

#strategy for finding good xylose species - find species for which
#all genes are in top 10% of respective distributions
#we can't directly compare tAI values (raw data)
#we have to compare estAI values (normalized data)

#first just making a df of just the highest estAI 
#for those spp with more than one. 
max_estAI_df<-estAI_df
for (i in 2:6){
  for(j in 1:nrow(estAI_df)){
    if(!is.na(estAI_df[j,i]) & !is.na(estAI_df[(j+5), i])){
      max_estAI_df[j,i]<-max(estAI_df[j,i],estAI_df[(j+5),i]) 
    }
  }
}
max_estAI_df<-max_estAI_df[1:6]
#can import from now on
max_estAI_df$all_taxa<-tolower(max_estAI_df$all_taxa)
write.table(max_estAI_df, "~/xylose_optimization_project/data/spp_by_gene_maximum_paralog_estAI_vals.txt", sep="\t", quote=FALSE, row.names=FALSE)


#now finding spp in top ten percent for all
toptens<- as.character(unique(max_estAI_df$all_taxa))
for (i in 2:ncol(max_estAI_df)){
  vec<-sort(max_estAI_df[,i], decreasing = TRUE)
  spp<-max_estAI_df[which(max_estAI_df[,i]>=
                            quantile(vec, 0.90)), 1]
  keep<-which(toptens %in% spp)
  toptens<-toptens[keep]
}
##3 spp in 90th percentile for xyl1, xyl2, xyl3: 
#spathaspora gorwiae
#kluyveromyces aestuarii 
#sugiyamaella lignohabitans
##1 sp in 90th percentile for xyl1, xyl2, xyl3, tkl1:
#spathaspora gorwiae
#0 spp in top 10% for all 5 genes

#wider net of top 25% (75th percentile)
top25s<- as.character(unique(max_estAI_df$all_taxa))
for (i in 2:ncol(max_estAI_df)){
  vec<-sort(max_estAI_df[,i], decreasing = TRUE)
  spp<-max_estAI_df[which(max_estAI_df[,i]>=
                            quantile(vec, 0.75)), 1]
  keep<-which(top25s %in% spp)
  top25s<-top25s[keep]
}

#in top 25 for all 5 genes: 
# spathaspora gorwiae 
# spathaspora hagerdaliae
# spathaspora girioi
# kodamaea ohmeri




###EXAMINGING CORRELLATIONS BTW XYLOSE ESTAI AND XYLOSE GROWTH DATA
#data from Dana Opulente 

#have to merge spp. to keys. 
key<-read.delim("~/xylose_optimization_project/data/Spp_indices.txt", strip.white = TRUE)
key$Species<-tolower(key$Species)
growth_data<-read.delim("~/xylose_optimization_project/data/Xylose_growth_data_DO.txt", strip.white = TRUE)
x<-merge(growth_data, key[c(1,3)], by="PU.")

which(!x$all_taxa %in% max_estAI_df$all_taxa)->not_in_data

#Ok - the spp. we have that are NOT in the tree mostly seem to be due to renaming. 
#looking up possible new spp. names using second name of sp. 
for (i in 1:length(not_in_data)){
  x[(not_in_data[i]), 3]->old_sp
  spName<-strsplit(old_sp, " ")[[1]][2]
  possiblechange<-max_estAI_df$all_taxa[which(grepl(spName, max_estAI_df$all_taxa))]
  if (length(possiblechange)==1){
    x[(not_in_data[i]), 3]<-possiblechange
  }
}

which(!x$all_taxa %in% max_estAI_df$all_taxa)->not_in_data
x[(not_in_data[1]), 3]<-"hanseniaspora vinae"
x[(not_in_data[2]), 3]<-"spencermartinsiella europaea"
x[(not_in_data[3]), 3]<-"martiniozyma abiesophila"
x[(not_in_data[4]), 3]<-"nakaseomyces castellii"
#x[(not_in_data[5]), 3]<-albicans - which we don't have
x[(not_in_data[6]), 3]<-"suhomyces pyralidae"
x[(not_in_data[7]), 3]<-"metschnikowia lockheadii"
x[(not_in_data[8]), 3]<-"metschnikowia dekortum"
#x[(not_in_data[9]), 3]<-"metschnikowia gruessii"-> not sure who this is. 
colnames(x)[3]<-"all_taxa"

growth_rates_df<-merge(max_estAI_df, x, by="all_taxa")
write.table(growth_rates_df, "xylose_optimization_project/data/growth_rate_master_df.txt", sep="\t", quote=FALSE, row.names=FALSE)


#######work from above table from now on
install.packages("ade4")
require(ade4)
library(ade4)
growth_rates<-read.delim("~/xylose_optimization_project/data/growth_rate_master_df.txt", stringsAsFactors = FALSE)

plot(x=growth_rates$xyl1, y=growth_rates$Growth.Rate)
filtered_growth_rates<-growth_rates[which(growth_rates$Growth.Rate>0), ]
plot(x=filtered_growth_rates$xyl1, y=filtered_growth_rates$Growth.Rate)
cor.test(x=filtered_growth_rates$xyl1, y=filtered_growth_rates$Growth.Rate, method = "pearson")
cor.test(x=filtered_growth_rates$xyl2, y=filtered_growth_rates$Growth.Rate, method = "pearson")
cor.test(x=filtered_growth_rates$xyl3, y=filtered_growth_rates$Growth.Rate, method = "pearson")
cor.test(x=filtered_growth_rates$tkl1, y=filtered_growth_rates$Growth.Rate, method = "pearson")
cor.test(x=filtered_growth_rates$tal1, y=filtered_growth_rates$Growth.Rate, method = "pearson")

########
library(ape)
require(stringr)
tree<-read.tree("~/xylose_optimization_project/data/iTol_files/332_Newick_tree.txt")
key<-read.delim("~/xylose_optimization_project/data/Spp_indices.txt", strip.white = TRUE, stringsAsFactors = FALSE)
tree$tip.label
tree$tip.label<-str_replace(string = tree$tip.label, pattern = "_", replacement = " ")
tree$tip.label<-tolower(tree$tip.label)
tree$tip.label
growth_data<-read.delim("~/xylose_optimization_project/data/growth_rate_master_df.txt", stringsAsFactors = FALSE)
x<-growth_data

x$all_taxa[which(!x$all_taxa %in% tree$tip.label)]

#Ok - the spp. we have that are NOT in the tree mostly seem to be due to renaming. 
#looking up possible new spp. names using second name of sp. 
which(!x$all_taxa %in% tree$tip.label)->not_in_data
x[(not_in_data[1]), 1]<-"blastobotrys raffinosifermentans"
#x[(not_in_data[2]), 1]<-"spencermartinsiella europaea"
x[(not_in_data[3]), 1]<-"hanseniaspora vineae"
x[(not_in_data[4]), 1]<-"lachancea fantastica_nom_nud"
x[(not_in_data[5]), 1]<-"magnusiomyces tetraspermus"
x[(not_in_data[6]), 1]<-"metschnikowia dekortorum"
x[(not_in_data[7]), 1]<-"metschnikowia lochheadii"
x[(not_in_data[8]), 1]<-"metschnikowia matae_var._matae"
x[(not_in_data[9]), 1]<-"metschnikowia matae_var._maris"
x[(not_in_data[10]), 1]<-"candida castellii"
x[(not_in_data[11]), 1]<-"ogataea philodendri"
x[(not_in_data[12]), 1]<-"ogataea populialbae"
#pulling candida azyma out b/c it's not in the tree
x<-x[-24, ]
#now need NA rows for those in the tree but not in the data
tree$tip.label[which(!tree$tip.label %in% x$all_taxa)]->all_taxa
tobind<-data.frame(all_taxa)
tobind$"xyl1"<-NA      
tobind$"xyl2"<-NA
tobind$"xyl3"<-NA   
tobind$"tkl1"<-NA
tobind$"tal1"<-NA
tobind$"PU."<-NA
tobind$"Growth.Rate"<-NA
x<-rbind(x, tobind)
labs<-data.frame(tree$tip.label, c(1:length(tree$tip.label)))
colnames(labs)<-c("all_taxa", "phylo_order")
x<-merge(x, labs, by="all_taxa")
x<-x[order(x$phylo_order), ]

####################################################################
##IGNORE ABOVE -> fixed names in data to match 332 tree
###############
##Import the below written dataframe for PIC analysis
###################################################################
write.table(x, "~/xylose_optimization_project/data/growth_rate_master_df_treematchednames.txt", sep="\t", quote=FALSE, row.names = FALSE)
#growth data and estAI vals
growth_data<-read.delim("~/xylose_optimization_project/data/growth_rate_master_df_treematchednames.txt")
#332 tree
library(ape)
require(stringr)
require(ggpubr)
tree<-read.tree("~/xylose_optimization_project/data/iTol_files/332_Newick_tree.txt")
tree$tip.label
tree$tip.label<-str_replace(string = tree$tip.label, pattern = "_", replacement = " ")
tree$tip.label<-tolower(tree$tip.label)
#
require(ade4)
require(adephylo)
require(ape)
#for each orthogram I need to drop the tree taxa that are missing data
#they're in the same order - so should be easy
removes<-which(is.na(growth_data$xyl1))
xyl1_tree<-drop.tip(tree, tree$tip.label[removes])
removes<-which(is.na(growth_data$xyl2))
xyl2_tree<-drop.tip(tree, tree$tip.label[removes])
removes<-which(is.na(growth_data$xyl3))
xyl3_tree<-drop.tip(tree, tree$tip.label[removes])
removes<-which(is.na(growth_data$tkl1))
tkl1_tree<-drop.tip(tree, tree$tip.label[removes])
removes<-which(is.na(growth_data$tal1))
tal1_tree<-drop.tip(tree, tree$tip.label[removes])
removes<-which(is.na(growth_data$Growth.Rate))
growth_tree<-drop.tip(tree, tree$tip.label[removes])


#phylo <- ape::read.tree(text = tree)
xyl1<-as.numeric(na.omit(growth_data$xyl1))
orthogram(xyl1, tre = xyl1_tree)
#xyl1=diffuse dependence
xyl2<-as.numeric(na.omit(growth_data$xyl2))
orthogram(xyl2, tre=xyl2_tree)
#xyl2= no phylogenetic dependence
xyl3<-as.numeric(na.omit(growth_data$xyl3))
orthogram(xyl3, tre=xyl3_tree)
#xyl3=diffuse dependence and specific node importance
tal1<-as.numeric(na.omit(growth_data$tal1))
orthogram(tal1, tre=tal1_tree)
#tal1 = no phylo dependence
tkl1<-as.numeric(na.omit(growth_data$tkl1))
orthogram(tkl1, tre=tkl1_tree)
#tkl1 = no phylo dependence
growth<-as.numeric(na.omit(growth_data$Growth.Rate))
orthogram(growth, growth_tree)
#growth = diffuse phylo dependence

##PIC below: 
PIC_df<-growth_data
#PIC growth data x XYL1
removes<-which(is.na(growth_data$xyl1) | is.na(growth_data$Growth.Rate))
xyl1PICtree<-drop.tip(tree, tree$tip.label[removes])
df<-growth_data[-removes, ]
PIC.xyl1<-pic(df$xyl1, xyl1PICtree)
PIC.growth<-pic(df$Growth.Rate, xyl1PICtree)
cor.test(PIC.xyl1, PIC.growth)


