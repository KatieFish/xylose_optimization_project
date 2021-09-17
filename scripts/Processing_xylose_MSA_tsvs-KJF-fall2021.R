
#09-16-2021

require(stringr)
#read in ortholog cds matrices
xyl1<- read.delim("~/Xyl_project_Fall_2021/xyl1/xyl1_cds_matrix_forR.txt",
                  stringsAsFactors=FALSE)
xyl2<- read.delim("~/Xyl_project_Fall_2021/xyl2/xyl2_cds_matrix_forR.txt", 
                  stringsAsFactors=FALSE)
xyl3<- read.delim("~/Xyl_project_Fall_2021/xyl3/xyl3_cds_matrix_forR.txt", 
                  stringsAsFactors=FALSE)

xyl1[,1]<-tolower(xyl1[,1])
xyl2[,1]<-tolower(xyl2[,1])
xyl3[,1]<-tolower(xyl3[,1])
##################################################################
xyl1$gene<-"xyl1"
xyl2$gene<-"xyl2"
xyl3$gene<-"xyl3"

###################################################################
#make a master list of all taxa present for all 3 genes
all_taxa<-append(unique(xyl1$taxon), unique(xyl2$taxon))
all_taxa<-unique(append(all_taxa, unique(xyl3$taxon)))

wi_vals<-read.delim("~/xylose_optimization_project/fall_2020/data/labella_et_al/wi_values.txt",
                    stringsAsFactors=FALSE)

# are there any spp. in our lists that are not in abbe's wi values?

#Note - from abbe's paper - 
#Nadsonia fulvescens var. fulvescens does not have wi vals, so 
#the nadsonia fulvescens in the wi file must be nadsonia fulvescens var. elongata.
#middelhovenomyces tepae is missing with no alternative name
#kazachstania humatica was not part of the 332 paper and maybe JW added the genome?
#spencermartinsiella europaea also not in abbe's data

#getting rid of middelhovenomyces tepae, kazachstania humatica, 
#spencermartinsiella europaea, and nadsonia fulvesens from our data.

#removing the NA rows from abbe's wi dataset (botryozyma nematodophila, martiniozyma abiesophila)

taxa_to_remove<-all_taxa[which(!all_taxa %in% wi_vals[,1])]

#write code to calculate the mean wi value for each gene (estAI)
require(EnvStats)
######seperately for each gene#########

#####start with xyl1##############
xyl1<-xyl1[which(!xyl1$taxon %in% taxa_to_remove),]

xyl1.stAI_dataframe<-xyl1[c(1,2,ncol(xyl1))]
xyl1.stAI_dataframe$stAI<-NA

##Loop to caluculate stAI values for each entry in the fasta file
for(i in 1:nrow(xyl1.stAI_dataframe)){  
  gene_x<-as.character(xyl1[i,3:(ncol(xyl1)-1)])
    gene_x<-gene_x[which(!is.na(gene_x))]
    y<- which(wi_vals[,1] == xyl1.stAI_dataframe$taxon[i])
    #create a vector to store the wi values
    gene_x_wi_vals<-numeric()
    #iterate through gene, starting at 4th nt (start ignored)
    for(j in seq(4,length(gene_x), 3)){
      codon<-paste(gene_x[j], gene_x[j+1], gene_x[j+2], sep="")
      codon_column<-which(grepl(codon, colnames(wi_vals), ignore.case=TRUE))
      codon_wi_val<-as.numeric(wi_vals[y,codon_column])
      gene_x_wi_vals<-append(gene_x_wi_vals, codon_wi_val)
    }
    stAI<-geoMean(gene_x_wi_vals)
    xyl1.stAI_dataframe$stAI[i]<-stAI 
  }

#####start with xyl2##############
xyl2<-xyl2[which(!xyl2$taxon %in% taxa_to_remove),]

xyl2.stAI_dataframe<-xyl2[c(1,2,ncol(xyl2))]
xyl2.stAI_dataframe$stAI<-NA

##Loop to caluculate stAI values for each entry in the fasta file
for(i in 1:nrow(xyl2.stAI_dataframe)){
  gene_x<-as.character(xyl2[i,3:(ncol(xyl2)-1)])
  gene_x<-gene_x[which(!is.na(gene_x))]
  y<- which(wi_vals[,1] == xyl2.stAI_dataframe$taxon[i])
  #create a vector to store the wi values
  gene_x_wi_vals<-numeric()
  for(j in seq(4,length(gene_x), 3)){
    codon<-paste(gene_x[j], gene_x[j+1], gene_x[j+2], sep="")
    codon_column<-which(grepl(codon, colnames(wi_vals), ignore.case=TRUE))
    codon_wi_val<-as.numeric(wi_vals[y,codon_column])
    gene_x_wi_vals<-append(gene_x_wi_vals, codon_wi_val)
  }
  stAI<-geoMean(gene_x_wi_vals)
  xyl2.stAI_dataframe$stAI[i]<-stAI 
}

#####start with xyl3##############
xyl3<-xyl3[which(!xyl3$taxon %in% taxa_to_remove),]

xyl3.stAI_dataframe<-xyl3[c(1,2,ncol(xyl3))]
xyl3.stAI_dataframe$stAI<-NA

##Loop to caluculate stAI values for each entry in the fasta file
for(i in 1:nrow(xyl3.stAI_dataframe)){
  gene_x<-as.character(xyl3[i,3:(ncol(xyl3)-1)])
  gene_x<-gene_x[which(!is.na(gene_x))]
  y<- which(wi_vals[,1] == xyl3.stAI_dataframe$taxon[i])
  #create a vector to store the wi values
  gene_x_wi_vals<-numeric()
  for(j in seq(4,length(gene_x), 3)){
    codon<-paste(gene_x[j], gene_x[j+1], gene_x[j+2], sep="")
    codon_column<-which(grepl(codon, colnames(wi_vals), ignore.case=TRUE))
    codon_wi_val<-as.numeric(wi_vals[y,codon_column])
    gene_x_wi_vals<-append(gene_x_wi_vals, codon_wi_val)
  }
  stAI<-geoMean(gene_x_wi_vals)
  xyl3.stAI_dataframe$stAI[i]<-stAI 
}

rm(all_taxa,taxa_to_remove, codon, codon_column, codon_wi_val, gene_x, gene_x_wi_vals, i, j, stAI, y)


##########
#Master dataset 
xyl_master_df<-rbind(xyl1.stAI_dataframe, xyl2.stAI_dataframe, xyl3.stAI_dataframe)
write.table(xyl_master_df, "~/Xyl_project_Fall_2021/stAI_dataframe-09-16-2021.txt", 
            row.names=FALSE, quote=FALSE, sep="\t")
###############



#Now we need to find estAI values
####################

#import genome-wide stAI values
genome_wide_tAI<-read.delim("~/xylose_optimization_project/fall_2020/data/labella_et_al/cds_mito_processed_tAI_recalc/genome_wide_tAI_all_spp.txt",
                            stringsAsFactors=FALSE)
genome_wide_tAI$taxa<-as.character(genome_wide_tAI$taxa)

xyl_master_df$estAI<-NA

# iterate through stAI df and populate estAI df 
  for(j in 1:nrow(xyl_master_df)){
  # grab the species tAI vals for all genes 
  species_tAI_vals<-genome_wide_tAI[which(genome_wide_tAI$taxa==xyl_master_df$taxon[j]), 2]
  ecdf_func<-ecdf(species_tAI_vals)
  xyl_master_df$estAI[j]<-ecdf_func(xyl_master_df$stAI[j])
  }

rm(ecdf_func, j, species_tAI_vals)

####writing estAI results as tables to import
write.table(xyl_master_df, "Xyl_project_Fall_2021/estAI_dataframe-09-16-2021.txt", sep="\t", quote=FALSE, row.name=FALSE)
########################


#create a max estAI_dataframe

max_estAI_df<-data.frame(aggregate(estAI ~ gene + taxon, data=xyl_master_df, FUN=max))

write.table(max_estAI_df, "~/Xyl_project_Fall_2021/maximum_estAI_vals.txt", sep="\t", quote=FALSE, row.names=FALSE)


#Compare histograms of estAI vals between genes
#RISHITHA



#now finding spp in top ten percent for all
toptens<- as.character(unique(max_estAI_df$taxon))
genes<-c("xyl1", "xyl2", "xyl3")

for (i in 1:length(genes)){
  rows<-subset(max_estAI_df, subset=max_estAI_df$gene==genes[i])
  vec<-sort(rows$estAI, decreasing = TRUE)
  ##top ten percent by stting quantile to .9##
  spp<-rows$taxon[which(rows$estAI>=
                            quantile(vec, 0.90))]
  keep<-which(toptens %in% spp)
  toptens<-toptens[keep]
}

rm(spp, rows, i, keep, spp, vec, genes)

write.table(toptens, "~/Xyl_project_Fall_2021/spp_in_90th_percentile_codon_opt-all_genes.txt",
            sep="\t", row.names=FALSE, quote=FALSE)
## 8 spp in 90th percentile for xyl1, xyl2, xyl3: 

#ogataea nitratoaversa,scheffersomyces lignosus,scheffersomyces stipitis,
#spathaspora arborariae,spathaspora girioi,spathaspora gorwiae,
#spathaspora hagerdaliae,spathaspora passalidarum

#all but one spp belongs to the genus scheff OR spath**
# WHAT IS THE PROBABILITY OF ANY GIVEN SPP BEING IN THE TOP TEN % IN ALL 3 GENES?
#  -> 0.1^3 = 0.001 
# WHAT IS THE PROBABILIITY OF THESE BEING ALL OF THE SAME 2 GENERA?
# 


