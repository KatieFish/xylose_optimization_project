###Code written to merge D.O. quantitative growth data 
#with qualitative data from Opulente et al BMC Bio 2019

#growth data df put together by RLN from original data sent by Dana (also contains estAI vals)
growth_data<-read.delim("~/Xyl_project_Fall_2021/data tables/growth_data_df-10-19-21.txt",
                        header=TRUE, stringsAsFactors=FALSE, strip.white=TRUE)
growth_data<-growth_data[1:2]

#presence_absence dataframe from processing xyl homologs
presence_absence<-read.delim("~/Xyl_project_Fall_2021/data tables/pathway_presence_absence.txt",
                             stringsAsFactors =FALSE, strip.white=TRUE)

#7 taxa in growth data not in presence_absence - taxa that do not have any pathway genes
growth_data$all_taxa[which(!growth_data$all_taxa %in% presence_absence$taxon)]

#qualitative data curated by Dana from The Yeasts in Opulente et al. BMC Bio 2018
qual_data<-read.delim("~/Xyl_project_Fall_2021/data tables/DO_BMC_growth_qualitative_data.txt",
                      header=TRUE, stringsAsFactors=FALSE, strip.white=TRUE)

qual_data$taxa<-tolower(qual_data$taxa)

#taxa not in qualitative data
not_in_both<-growth_data$all_taxa[which(!growth_data$all_taxa %in% qual_data$taxa)]
as.data.frame(not_in_both)->x
sink("~/Xyl_project_Fall_2021/data tables/taxa_missing_in_qual.txt")
print(x, row.names=FALSE)
sink()

#8 different spellings made consistent 
#s. bayanus var bayanus changed to eubayanus
#s. bayanus var uvarum changed to uvarum
#s. arboricolus changed to arboricola
#reduced to 76 spp. of the 333 in our growth data matrix

#merging quantitative unpublished data and qualitative published data
colnames(qual_data)[1]<-"all_taxa"
qual_data<-qual_data[which(qual_data$all_taxa %in% growth_data$all_taxa),]
quant_qual_growth<-merge(growth_data, qual_data, by="all_taxa", all=TRUE)

#do any disagree?
quant_qual_growth$all_taxa[which(quant_qual_growth$Growth.Rate>0 & quant_qual_growth$opulente.2018.growth.binary<1)]
#no taxa that grow in Dana's assay but not in the qualitative data
quant_qual_growth[which(quant_qual_growth$Growth.Rate==0 & quant_qual_growth$opulente.2018.growth.binary==1),]
#24 spp. Dana has identfied with growth rates of 0 that are predicted to grow in qualitative data
#where they disagree I'm going to go with Dana's data
quant_qual_growth$growth.binary<-NA
quant_qual_growth$growth.binary[which(quant_qual_growth$Growth.Rate>0)]<-1
quant_qual_growth$growth.binary[which(quant_qual_growth$Growth.Rate==0)]<-0
for(i in 1:nrow(quant_qual_growth)){
  if(is.na(quant_qual_growth$Growth.Rate[i])){
    quant_qual_growth$growth.binary[i]<-quant_qual_growth$opulente.2018.growth.binary[i]
  }
}
###
length(which(is.na(quant_qual_growth$Growth.Rate) & !is.na(quant_qual_growth$growth.binary)))
#of 81 missing taxa from Dana's data:
#was able to fill in data for 60 spp. 
length(which(is.na(quant_qual_growth$Growth.Rate) & is.na(quant_qual_growth$growth.binary)))
#still missing 21 spp. from both datasets. 

##compiling all data into tab file to import into analysis
quant_qual_growth<-quant_qual_growth[c(1, 2, 5)]
#write.table(quant_qual_growth, "~/Xyl_project_Fall_2021/Pathway_growthrate_growthbinary.txt", sep="\t", quote=FALSE, row.names=FALSE)

#fixing file to use as iTol file.

#name swap issues 
name_swaps<-read.delim("~/xylose_optimization_project/spring_2021/data/tree_to_seq_data_nameswaps.txt")
for(i in 1:nrow(quant_qual_growth)){
  if(quant_qual_growth$all_taxa[i] %in% name_swaps$data_name){
    quant_qual_growth$all_taxa[i]<- name_swaps$tree_name[which(name_swaps$data_name==quant_qual_growth$all_taxa[i])]
   }
 }
write.table(quant_qual_growth, "~/Xyl_project_Fall_2021/iTOLKey.txt", sep="\t", quote=FALSE, row.names=FALSE)


