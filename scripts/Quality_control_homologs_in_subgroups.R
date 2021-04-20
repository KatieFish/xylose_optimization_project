
xyl1<-read.delim("~/xylose_optimization_project/spring_2021/subgroup_consensus/xyl1/xyl1.subgroup.matrix.txt",
                         header=TRUE, stringsAsFactors=FALSE)

#fix names in sequence mat
taxon_name_switches <- read.delim("~/xylose_optimization_project/spring_2021/data/taxon_name_switches.txt", stringsAsFactors=FALSE)

sequence_mat$V1<-tolower(sequence_mat$V1)

for(i in 1:nrow(sequence_mat)){
  if(sequence_mat$V1[i] %in% taxon_name_switches$sequences_name){
    sequence_mat$V1[i]<-taxon_name_switches$growth_data_name[
     which(taxon_name_switches$sequences_name==sequence_mat$V1[i])]
  }
}


#let's get rid of these for the time being
#sequence_mat<-sequence_mat[which(sequence_mat[,2]=="m"),]

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


#without any filter - only 1 spp missing. 
    # "candida sojae" 
#with only size filter - 
    #  "candida sojae","lipomyces arxii","lipomyces starkeyi"
    # if I extend the filter to 230+ I would recover seuqence(s) for both lipo spp. 
#getting rid of alternative start codons - 
  # 30 missing taxa: 
# "babjeviella inositovora"        "barnettozyma hawaiiensis"      
# "blastobotrys raffinofermentans" "candida arabinofermentans"     
# "candida ascalaphidarum"         "candida azyma"                 
# "candida freyschussii"           "candida gorgasii"              
# "candida incommunis"             "candida mycetangii"            
# "candida orthopsilosis"          "candida sojae"                 
# "clavispora lusitaniae"          "cyberlindnera jadinii"         
# "cyberlindnera suaveolens"       "deakozyma indianensis"         
#"dipodascus albidus"             "kuraishia molischiana"         
# "lipomyces japonicus"            "metschnikowia borealis"        
# "metschnikowia cerradonensis"    "metschnikowia continentalis"   
# "middelhovenomyces tepae"        "priceomyces haplophilus"       
# "scheffersomyces stipitis"       "suhomyces canberraensis"       
# "tortispora caseinolytica"       "wickerhamomyces bovis"         
# "wickerhamomyces canadensis"     "zygoascus ofunaensis" 


# I think we'll apply the size filter and then try to manually flesh out
#the alt start codon seqs. 

length(which(!sequence_mat[,2] == "m"))
#that would leave 127 sequences. I would need a pipeline here.

#1. Blast each sequence
#2. If top hit >=99% identical, keep full length record and nt fasta. 
