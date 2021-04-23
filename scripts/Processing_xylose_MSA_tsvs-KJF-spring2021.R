
#04-20-2021

require(stringr)
#read in ortholog cds matrices
xyl1<- read.delim("~/xylose_optimization_project/spring_2021/subgroup_consensus/xyl1/subgroup/xyl1_subgroup-codon_aligned_cds-cleaned.txt",
                  stringsAsFactors=FALSE)
xyl2<- read.delim("~/xylose_optimization_project/spring_2021/subgroup_consensus/xyl2/subgroup/xyl2_subgroup-codon_aligned_cds-cleaned.txt", 
                  stringsAsFactors=FALSE)
xyl3<- read.delim("~/xylose_optimization_project/spring_2021/subgroup_consensus/xyl3/subgroup/xyl3_subgroup-codon_aligned_cds-cleaned.txt", 
                  stringsAsFactors=FALSE)


#remove identical sequences using the unique() command
#updated 4-21-21 KJF
#I need want to find not just identical sequences, but sequences that
#are identical fragments of other sequences belonging to the SAME TAXON

xyl1_taxa<-unique(xyl1[,1])
remove<-vector()
for(i in 1:length(unique(xyl1_taxa))){
  ndf<-xyl1[which(xyl1[1]==xyl1_taxa[i]),]
  if(nrow(ndf)>1){
      for(j in 1:nrow(ndf)){
        seq<-as.character(ndf[j, 2:ncol(ndf)])
        start<-min(which(!is.na(seq)))
        end<-max(which(!is.na(seq)))
        for(k in 1:nrow(ndf)){
           if(k!=j){
            if(str_detect(paste(ndf[k, 2:ncol(ndf)], collapse=""),
              paste(seq[start:end], collapse=""))){
              seq_length<-length(seq)
              match_length<-max(which(!is.na(ndf[k, 2:ncol(ndf)])))-
                min(which(!is.na(ndf[k, 2:ncol(ndf)])))
              if(match_length>seq_length){
                remove<-append(remove, row.names(ndf)[j])
              } else {
                remove<-append(remove, row.names(ndf)[k])
              }
          }
        }
      }
    }
  }
  }
#loop generates a vector called remove which stores the rows of the dataframe 
#that contain sequences that are identical fragments or complete matches to another
#sequence in the dataframe of that same taxon.

remove<-unique(as.numeric(remove))
xyl1<-xyl1[-remove, ]

############Repeating remove loop for xyl2 and xyl3

xyl2_taxa<-unique(xyl2[,1])
remove<-vector()
for(i in 1:length(unique(xyl2_taxa))){
  ndf<-xyl2[which(xyl2[1]==xyl2_taxa[i]),]
  if(nrow(ndf)>1){
    for(j in 1:nrow(ndf)){
      seq<-as.character(ndf[j, 2:ncol(ndf)])
      start<-min(which(!is.na(seq)))
      end<-max(which(!is.na(seq)))
      for(k in 1:nrow(ndf)){
        if(k!=j){
          if(str_detect(paste(ndf[k, 2:ncol(ndf)], collapse=""),
                        paste(seq[start:end], collapse=""))){
            seq_length<-length(seq)
            match_length<-max(which(!is.na(ndf[k, 2:ncol(ndf)])))-
              min(which(!is.na(ndf[k, 2:ncol(ndf)])))
            if(match_length>seq_length){
              remove<-append(remove, row.names(ndf)[j])
            } else {
              remove<-append(remove, row.names(ndf)[k])
            }
          }
        }
      }
    }
  }
}
#loop generates a vector called remove which stores the rows of the dataframe 
#that contain sequences that are identical fragments or complete matches to another
#sequence in the dataframe of that same taxon.

remove<-unique(as.numeric(remove))
xyl2<-xyl2[-remove, ]

#####
xyl3_taxa<-unique(xyl3[,1])
remove<-vector()
for(i in 1:length(unique(xyl3_taxa))){
  ndf<-xyl3[which(xyl3[1]==xyl3_taxa[i]),]
  if(nrow(ndf)>1){
    for(j in 1:nrow(ndf)){
      seq<-as.character(ndf[j, 2:ncol(ndf)])
      start<-min(which(!is.na(seq)))
      end<-max(which(!is.na(seq)))
      for(k in 1:nrow(ndf)){
        if(k!=j){
          if(str_detect(paste(ndf[k, 2:ncol(ndf)], collapse=""),
                        paste(seq[start:end], collapse=""))){
            seq_length<-length(seq)
            match_length<-max(which(!is.na(ndf[k, 2:ncol(ndf)])))-
              min(which(!is.na(ndf[k, 2:ncol(ndf)])))
            if(match_length>seq_length){
              remove<-append(remove, row.names(ndf)[j])
            } else {
              remove<-append(remove, row.names(ndf)[k])
            }
          }
        }
      }
    }
  }
}
#loop generates a vector called remove which stores the rows of the dataframe 
#that contain sequences that are identical fragments or complete matches to another
#sequence in the dataframe of that same taxon.

remove<-unique(as.numeric(remove))
xyl3<-xyl3[-remove, ]

rm(start, end, i, j, k, seq, seq_length, 
   match_length,remove, xyl1_taxa, xyl2_taxa, xyl3_taxa, ndf)

##################################################################

#make a master list of all taxa for all genes
all_taxa<-append(xyl1$V1, xyl2$V1)
all_taxa<-append(all_taxa, xyl3$V1)
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

write.table(xylose_utilization_gene_presence, 
            "~/xylose_optimization_project/spring_2021/data/XYLpthwy_gene_presence_absence_matrix.txt", sep="\t", quote=FALSE, row.names=FALSE)

wi_vals<-read.delim("~/xylose_optimization_project/fall_2020/data/labella_et_al/wi_values.txt",
                    stringsAsFactors=FALSE)


# are there any spp. in our lists that are not in abbe's wi values?
which(!xylose_utilization_gene_presence[,1] %in% wi_vals[,1])
taxa_to_remove<-xylose_utilization_gene_presence[which(!xylose_utilization_gene_presence[,1] %in% wi_vals[,1]), 1]


##04-21-21
#for now we are going to drop the 4 taxa that are not in labella et al wi val dataset.
#we can go back later and see if we can rescue any of these with name changes. 

xylose_utilization_gene_presence<-xylose_utilization_gene_presence[
  -which(xylose_utilization_gene_presence$all_taxa %in% taxa_to_remove),]
xyl1<-xyl1[-which(xyl1$V1 %in% taxa_to_remove),]
xyl2<-xyl2[-which(xyl2$V1 %in% taxa_to_remove),]
xyl3<-xyl3[-which(xyl3$V1 %in% taxa_to_remove),]




### start a new section of code! ;) In this section we will
#write code to calculate the mean wi value for each gene (estAI)
require(EnvStats)

##########looking for max # of paralogs
max((max(data.frame(table(xyl1$V1))[2])),
    (max(data.frame(table(xyl2$V1))[2])),
    (max(data.frame(table(xyl3$V1))[2])))
#max paralogs is 6 - belonging to alloascoidea hylecoeti - probably 
#not real but we don't have any real filter to get rid of these yet. 

stAI_dataframe<-xylose_utilization_gene_presence
#make new columns for paralogs
stAI_dataframe$xyl1.1<-NA
stAI_dataframe$xyl1.2<-NA
stAI_dataframe$xyl1.3<-NA
stAI_dataframe$xyl1.4<-NA
stAI_dataframe$xyl1.5<-NA
stAI_dataframe$xyl1.6<-NA
stAI_dataframe$xyl2.1<-NA
stAI_dataframe$xyl2.2<-NA
stAI_dataframe$xyl2.3<-NA
stAI_dataframe$xyl2.4<-NA
stAI_dataframe$xyl2.5<-NA
stAI_dataframe$xyl2.6<-NA
stAI_dataframe$xyl3.1<-NA
stAI_dataframe$xyl3.2<-NA
stAI_dataframe$xyl3.3<-NA
stAI_dataframe$xyl3.4<-NA
stAI_dataframe$xyl3.5<-NA
stAI_dataframe$xyl3.6<-NA

##Loop to caluculate stAI values for each paralog of XYL1
for(i in 1:nrow(xylose_utilization_gene_presence)){
  #check if gene is present for species
  if(xylose_utilization_gene_presence[i,2]==1){
    # create x variable that refers to row number in xyl1 of species
    x<- which(xyl1[,1] == xylose_utilization_gene_presence[i,1])
    y<- which(wi_vals[,1] == xylose_utilization_gene_presence[i,1])
    for(l in 1:length(x)){
    gene_x<-as.character(xyl1[x[l],((min(which(!is.na(as.character(xyl1[x[l], 2:ncol(xyl1)]))))+1):
                                    (max(which(!is.na(as.character(xyl1[x[l], 2:ncol(xyl1)]))))+1))])
    gene_x<-gene_x[which(!is.na(gene_x))]
    #create a vector to store the wi values
    gene_x_wi_vals<-numeric()
    for(j in seq(1,length(gene_x), 3)){
      codon<-paste(gene_x[j], gene_x[j+1], gene_x[j+2], sep="")
      codon_column<-which(grepl(codon, colnames(wi_vals), ignore.case=TRUE))
      codon_wi_val<-as.numeric(wi_vals[y,codon_column])
      gene_x_wi_vals<-append(gene_x_wi_vals, codon_wi_val)
    }
    stAI<-geoMean(gene_x_wi_vals)
    stAI_dataframe[i,(l+4)]<-stAI 
  }
}
}

##Loop to caluculate stAI values for each paralog of XYL2
for(i in 1:nrow(xylose_utilization_gene_presence)){
  #check if gene is present for species
  if(xylose_utilization_gene_presence[i,3]==1){
    # create x variable that refers to row number in xyl2 of species
    x<- which(xyl2[,1] == xylose_utilization_gene_presence[i,1])
    y<- which(wi_vals[,1] == xylose_utilization_gene_presence[i,1])
    for(l in 1:length(x)){
      gene_x<-as.character(xyl2[x[l],((min(which(!is.na(as.character(xyl2[x[l], 2:ncol(xyl2)]))))+1):
                                        (max(which(!is.na(as.character(xyl2[x[l], 2:ncol(xyl2)]))))+1))])
      gene_x<-gene_x[which(!is.na(gene_x))]
      #create a vector to store the wi values
      gene_x_wi_vals<-numeric()
      for(j in seq(1,length(gene_x), 3)){
        codon<-paste(gene_x[j], gene_x[j+1], gene_x[j+2], sep="")
        codon_column<-which(grepl(codon, colnames(wi_vals), ignore.case=TRUE))
        codon_wi_val<-as.numeric(wi_vals[y,codon_column])
        gene_x_wi_vals<-append(gene_x_wi_vals, codon_wi_val)
      }
      stAI<-geoMean(gene_x_wi_vals)
      stAI_dataframe[i,(l+10)]<-stAI 
    }
  }
}

##Loop to caluculate stAI values for each paralog of XYL3
for(i in 1:nrow(xylose_utilization_gene_presence)){
  #check if gene is present for species
  if(xylose_utilization_gene_presence[i,4]==1){
    # create x variable that refers to row number in xyl3 of species
    x<- which(xyl3[,1] == xylose_utilization_gene_presence[i,1])
    y<- which(wi_vals[,1] == xylose_utilization_gene_presence[i,1])
    for(l in 1:length(x)){
      gene_x<-as.character(xyl3[x[l],((min(which(!is.na(as.character(xyl3[x[l], 2:ncol(xyl3)]))))+1):
                                        (max(which(!is.na(as.character(xyl3[x[l], 2:ncol(xyl3)]))))+1))])
      gene_x<-gene_x[which(!is.na(gene_x))]
      #create a vector to store the wi values
      gene_x_wi_vals<-numeric()
      for(j in seq(1,length(gene_x), 3)){
        codon<-paste(gene_x[j], gene_x[j+1], gene_x[j+2], sep="")
        codon_column<-which(grepl(codon, colnames(wi_vals), ignore.case=TRUE))
        codon_wi_val<-as.numeric(wi_vals[y,codon_column])
        gene_x_wi_vals<-append(gene_x_wi_vals, codon_wi_val)
      }
      stAI<-geoMean(gene_x_wi_vals)
      stAI_dataframe[i,(l+16)]<-stAI 
    }
  }
}

rm(all_taxa, codon, codon_column, codon_wi_val, gene_x, gene_x_wi_vals, i, j, l, stAI, x)

stAI_dataframe<-stAI_dataframe[c(1, 5:ncol(stAI_dataframe))]

#Create max stAI_dataframe

max_stAI_dataframe<-stAI_dataframe[1]
max_stAI_dataframe$max_xyl1_stAI<-NA
max_stAI_dataframe$max_xyl2_stAI<-NA
max_stAI_dataframe$max_xyl3_stAI<-NA
#populate
for(i in 1:nrow(stAI_dataframe)){
  if(length(which(!is.na(as.numeric(stAI_dataframe[i, 2:7]))))>0){
    max_stAI_dataframe$max_xyl1_stAI[i]<-max(na.omit(as.numeric(stAI_dataframe[i,2:7])))
  }
  if(length(which(!is.na(as.numeric(stAI_dataframe[i, 8:13]))))>0){
    max_stAI_dataframe$max_xyl2_stAI[i]<-max(na.omit(as.numeric(stAI_dataframe[i,8:13])))
  } 
  if(length(which(!is.na(as.numeric(stAI_dataframe[i, 14:19]))))>0){
    max_stAI_dataframe$max_xyl3_stAI[i]<-max(na.omit(as.numeric(stAI_dataframe[i,14:19])))
  }  
}


#table contains stAI value for each gene for each species
write.table(stAI_dataframe, "~/xylose_optimization_project/spring_2021/data/spp_by_gene_stAI_vals.txt", 
            row.names=FALSE, quote=FALSE, sep="\t")

#table contains the max stAI value for each species
write.table(max_stAI_dataframe, "~/xylose_optimization_project/spring_2021/data/spp_by_gene_maximum_paralog_stAI_vals.txt",
            row.names=FALSE, quote=FALSE, sep="\t")


#import genome-wide stAI values

genome_wide_tAI<-read.delim("~/xylose_optimization_project/fall_2020/data/labella_et_al/cds_mito_processed_tAI_recalc/genome_wide_tAI_all_spp.txt")
#after unzipping before import had to correct taxa names in file on command line using sed to replace 
#underscores with spaces. 
#have to fix taxa names
#x<-as.character(genome_wide_tAI$taxa)
#taxa_IDs<-x
#for (i in 1:length(x)){
#  if (grepl("yHMP", x[i]) | grepl("yHAB", x[i])){
#    name<-paste(strsplit(x[i], " ")[[1]][2], strsplit(x[i], " ")[[1]][3], sep=" ")
 # } else if (grepl("maris", x[i])){
 #   name<-paste(strsplit(x[i], " ")[[1]][1], strsplit(x[i], " ")[[1]][2], strsplit(x[i], " ")[[1]][3], sep=" ")
#  }
#  else{
#    name<-paste(strsplit(x[i], " ")[[1]][1], strsplit(x[i], " ")[[1]][2], sep=" ")
#  }
#  taxa_IDs[i]<-name
#}
####
#genome_wide_tAI$taxa<-taxa_IDs
###NAMES FIXED ABOVE
#Hash out and import in the future. 
#write.table(genome_wide_tAI, "~/xylose_optimization_project/fall_2020/data/labella_et_al/cds_mito_processed_tAI_recalc/genome_wide_tAI_all_spp.txt", sep="\t", quote=FALSE, row.names=FALSE)

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

####writing estAI results as tables to import
write.table(estAI_df, "~/xylose_optimization_project/spring_2021/data/spp_by_gene_estAI_vals.txt", sep="\t", quote=FALSE, row.name=FALSE)

#create a max estAI_dataframe
max_estAI_df<-estAI_df[1]
max_estAI_df$max_xyl1_estAI<-NA
max_estAI_df$max_xyl2_estAI<-NA
max_estAI_df$max_xyl3_estAI<-NA
#populate
for(i in 1:nrow(estAI_df)){
  if(length(which(!is.na(as.numeric(estAI_df[i, 2:7]))))>0){
  max_estAI_df$max_xyl1_estAI[i]<-max(na.omit(as.numeric(estAI_df[i,2:7])))
  }
  if(length(which(!is.na(as.numeric(estAI_df[i, 8:13]))))>0){
  max_estAI_df$max_xyl2_estAI[i]<-max(na.omit(as.numeric(estAI_df[i,8:13])))
  } 
  if(length(which(!is.na(as.numeric(estAI_df[i, 14:19]))))>0){
  max_estAI_df$max_xyl3_estAI[i]<-max(na.omit(as.numeric(estAI_df[i,14:19])))
  }  
}

#first just making a df of just the highest estAI 
write.table(max_estAI_df, "~/xylose_optimization_project/spring_2021/data/spp_by_gene_maximum_paralog_estAI_vals.txt", sep="\t", quote=FALSE, row.names=FALSE)


#Compare histograms of estAI vals between genes
quartz()
par(mar=c(2,2,2,2))
par(mfrow=c(2,3))
hist(max_estAI_df$max_xyl1_estAI)
hist(max_estAI_df$max_xyl2_estAI)
hist(max_estAI_df$max_xyl3_estAI)


#strategy for finding good xylose species - find species for which
#all genes are in top 10% of respective distributions
#we can't directly compare tAI values (raw data)
#we have to compare estAI values (normalized data)



#now finding spp in top ten percent for all
toptens<- as.character(unique(max_estAI_df$all_taxa))
for (i in 2:ncol(max_estAI_df)){
  vec<-sort(max_estAI_df[,i], decreasing = TRUE)
  spp<-max_estAI_df[which(max_estAI_df[,i]>=
                            quantile(vec, 0.90)), 1]
  keep<-which(toptens %in% spp)
  toptens<-toptens[keep]
}
## 9 spp in 90th percentile for xyl1, xyl2, xyl3: 
"spathaspora passalidarum" "spathaspora gorwiae"      "scheffersomyces stipitis"
"lipomyces doorenjongii"   "spathaspora arborariae"   "spathaspora girioi"      
"spathaspora hagerdaliae"  "ogataea nitratoaversa"    "scheffersomyces lignosus"



###EXAMINGING CORRELLATIONS BTW XYLOSE ESTAI AND XYLOSE GROWTH DATA
#data from Dana Opulente 
growth_data<-read.delim("~/xylose_optimization_project/spring_2021/data/growth_rate_master_df.txt", strip.white = TRUE, stringsAsFactors = FALSE)
#growth_data->x
#taxon_name_switches <- read.delim("~/xylose_optimization_project/spring_2021/data/growth_data_taxon_name_switches.txt", stringsAsFactors=FALSE)
#x$all_taxa<-as.character(x$all_taxa)
#taxon_name_switches$growth_data_name<-as.character(taxon_name_switches$growth_data_name)

#for(i in 1:nrow(x)){
#if(x$all_taxa[i] %in% taxon_name_switches$growth_data_name){
#  x$all_taxa[i]<-taxon_name_switches$sequences_name[which(
#    taxon_name_switches$growth_data_name==x$all_taxa[i])]
#  }
#}

#x<-merge(growth_data, max_estAI_df, by="all_taxa", all=TRUE)
#write.table(x, "xylose_optimization_project/spring_2021/data/growth_rate_master_df.txt", sep="\t", quote=FALSE, row.names=FALSE)

###NAMES FIXED 04-21-21 - JUST IMPORT IN FUTURE. 

#######work from above table from now on
install.packages("ade4")
require(ade4)
growth_rates<-read.delim("~/xylose_optimization_project/spring_2021/data/growth_rate_master_df.txt", stringsAsFactors = FALSE)

########
require(ape)
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


####################################################################
##IGNORE ABOVE -> fixed names in data to match 332 tree
###############
##Import the below written dataframe for PIC analysis
###################################################################
#write.table(x, "~/xylose_optimization_project/data/growth_rate_master_df_treematchednames.txt", sep="\t", quote=FALSE, row.names = FALSE)
#growth data and estAI vals
growth_data<-read.delim("~/xylose_optimization_project/data/growth_rate_master_df_treematchednames.txt")
#S values supplied by abbe's supp data
s_values<-read.delim("~/xylose_optimization_project/data/labella_et_al/s_values.txt", stringsAsFactors = FALSE)
#only taking spp with s values greater than .5
s_values<-s_values[which(s_values$species.name %in% growth_data$all_taxa), ]
colnames(s_values)<-c("all_taxa", "s_value")
growth_data<-merge(growth_data, s_values, by="all_taxa")
growth_data<-growth_data[which(growth_data$s_value > .5),]

#332 tree
library(ape)
require(stringr)
require(ggpubr)
tree<-read.tree("~/xylose_optimization_project/data/iTol_files/332_Newick_tree.txt")
tree$tip.label<-str_replace(string = tree$tip.label, pattern = "_", replacement = " ")
tree$tip.label<-tolower(tree$tip.label)
#
require(ade4)
require(adephylo)
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
#PIC growth data x XYL1
removes<-which(is.na(growth_data$xyl1) | is.na(growth_data$Growth.Rate))
df<-growth_data[-removes, ]
removes<-which(!tree$tip.label %in% df$all_taxa)
xyl1PICtree<-drop.tip(tree, tree$tip.label[removes])
#calculate and compare PIC values
PIC.xyl1<-pic(df$xyl1, xyl1PICtree)
PIC.growth<-pic(df$Growth.Rate, xyl1PICtree)
cor.test(PIC.xyl1, PIC.growth)
cordf<-data.frame(PIC.xyl1, PIC.growth)
ggscatter(data=cordf, x="PIC.xyl1", y="PIC.growth",  
          col="red", size=2, add="reg.line",add.params = list(color = "blue", fill = "gray"),
          cor.coeff.args = list(method = "pearson", label.sep = "\n"), conf.int=TRUE, cor.coef=TRUE, cor.method="pearson",
          xlab="PIC XYL1", ylab="PIC xylose growth rate")


#PIC growth data x XYL2
removes<-which(is.na(growth_data$xyl2) | is.na(growth_data$Growth.Rate))
df<-growth_data[-removes, ]
#remove optical outliers
removes<-which(!tree$tip.label %in% df$all_taxa)
xyl2PICtree<-drop.tip(tree, tree$tip.label[removes])
#calculate and compare PIC values
PIC.xyl2<-pic(df$xyl2, xyl2PICtree)
PIC.growth<-pic(df$Growth.Rate, xyl2PICtree)
cor.test(PIC.xyl2, PIC.growth)
cordf<-data.frame(PIC.xyl2, PIC.growth)
ggscatter(data=cordf, x="PIC.xyl2", y="PIC.growth",  
          col="red", size=2, add="reg.line",add.params = list(color = "blue", fill = "gray"),
          cor.coeff.args = list(method = "pearson", label.sep = "\n"), conf.int=TRUE, cor.coef=TRUE, cor.method="pearson",
          xlab="PIC xyl2", ylab="PIC xylose growth rate")


#PIC growth data x XYL3
removes<-which(is.na(growth_data$xyl3) | is.na(growth_data$Growth.Rate))
df<-growth_data[-removes, ]
#remove optical outliers
removes<-which(!tree$tip.label %in% df$all_taxa)
xyl3PICtree<-drop.tip(tree, tree$tip.label[removes])
#calculate and compare PIC values
PIC.xyl3<-pic(df$xyl3, xyl3PICtree)
PIC.growth<-pic(df$Growth.Rate, xyl3PICtree)
cor.test(PIC.xyl3, PIC.growth)
cordf<-data.frame(PIC.xyl3, PIC.growth)
ggscatter(data=cordf, x="PIC.xyl3", y="PIC.growth",  
          col="red", size=2, add="reg.line",add.params = list(color = "blue", fill = "gray"),
          cor.coeff.args = list(method = "pearson", label.sep = "\n"), conf.int=TRUE, cor.coef=TRUE, cor.method="pearson",
          xlab="PIC xyl3", ylab="PIC xylose growth rate")



###### Comparison of growers and non-growers
growers<-growth_data[which(growth_data$Growth.Rate>0), ]
nongrowers<-growth_data[which(growth_data$Growth.Rate <= 0), ]

growth_data<-growth_data[which(!is.na(growth_data$Growth.Rate)), ]
growth_data$growth.binary<-NA
for(i in 1:nrow(growth_data)){
  if(growth_data$Growth.Rate[i]>0){
    growth_data$growth.binary[i]<-1
  }
  if(growth_data$Growth.Rate[i]==0){
    growth_data$growth.binary[i]<-0
  }
}

boxplot(growth_data$xyl1 ~ as.factor(growth_data$growth.binary),
        col="lightblue", xlab="codon optimization")




getwd()

length(which(!s_values$all_taxa %in% tree$tip.label))
