#This matches up tip labels from the 332 dataset to our 
#binary xylose gene presence matrix. 

#files written to be used in iTol. 

#332 newick tree from Shen et al. 
#gene_pres_abs_matrix
gene_mat<-read.delim("~/xylose_optimization_project/data/XYLpthwy_gene_presence_absence_matrix.txt")
gene_mat$all_taxa<-tolower(gene_mat$all_taxa)

library(ape)
require(stringr)
tree<-read.tree("~/xylose_optimization_project/data/332_Newick_tree.txt")
tree$tip.label<-tolower(tree$tip.label)
tree$tip.label<-str_replace(string = tree$tip.label, pattern = "_", replacement = " ")
tree$tip.label[234]<-"lachancea fantastica"
tree_taxa<-tree$tip.label
which(!gene_mat$all_taxa %in% tree_taxa)-> not_in_data
gene_mat[not_in_data, 1]
#Ok - the spp. we have that are NOT in the tree mostly seem to be due to renaming. 
#looking up possible new spp. names using second name of sp. 
for (i in 1:length(not_in_data)){
  gene_mat[(not_in_data[i]), 1]->old_sp
  spName<-strsplit(old_sp, " ")[[1]][2]
  possiblechange<-tree$tip.label[which(grepl(spName, tree$tip.label))]
  if (length(possiblechange)==1){
    gene_mat[(not_in_data[i]), 1]<-possiblechange
  }
}
which(!gene_mat$all_taxa %in% tree_taxa)-> not_in_data
gene_mat[not_in_data, 1]
#ok - now there are less. Dealing with these one by one. 
gene_mat[(not_in_data[1]), 1]<-"metschnikowia dekortorum"
gene_mat[(not_in_data[2]), 1]<-"ogataea populialbae"
gene_mat[(not_in_data[3]), 1]<-"ogataea philodendri"
# cannot find
#gene_mat[(not_in_data[4]), 1]<-
gene_mat[(not_in_data[5]), 1]<-"blastobotrys raffinosifermentans"
gene_mat[(not_in_data[6]), 1]<-"candida castellii"
gene_mat[(not_in_data[7]), 1]<-"hanseniaspora vineae"
gene_mat[(not_in_data[8]), 1]<-"nadsonia fulvescens_var._fulvescens"
gene_mat[(not_in_data[9]), 1]<-"magnusiomyces tetraspermus"
gene_mat[(not_in_data[11]), 1]<-"metschnikowia lochheadii"
gene_mat[(not_in_data[12]), 1]<-"metschnikowia matae_var._matae"
gene_mat<-gene_mat[-c(not_in_data[4], not_in_data[10], not_in_data[13]),]
which(!tree$tip.label %in% gene_mat$all_taxa)-> not_in_data
tree$tip.label[not_in_data]
gene_mat[(nrow(gene_mat)+1), ]<-c("candida albicans", rep(NA, 5))               
gene_mat[(nrow(gene_mat)+1), ]<-c("metschnikowia matae_var._maris", rep(NA, 5))               
gene_mat[(nrow(gene_mat)+1), ]<-c("wickerhamomyces sp._yb_2243", rep(NA, 5))               
gene_mat[(nrow(gene_mat)+1), ]<-c("nadsonia fulvescens_var._elongata", rep(NA, 5))               

write.table(gene_mat, "~/xylose_optimization_project/data/XYLpthwy_gene_presence_absence_matrix.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.tree(tree, "~/xylose_optimization_project/data/332_Newick_tree_lowcaselabels.nwk")

### same name matching for growth data: 

growth_data<-read.delim("~/xylose_optimization_project/data/growth_rate_master_df_treematchednames.txt")
library(ape)
require(stringr)
tree<-read.tree("~/xylose_optimization_project/data/iTol_files/332_Newick_tree.txt")
tree$tip.label<-tolower(tree$tip.label)
tree$tip.label<-str_replace(string = tree$tip.label, pattern = "_", replacement = " ")
tree_taxa<-tree$tip.label
which(!growth_data$all_taxa %in% tree_taxa)-> not_in_data
growth_data[not_in_data, 1]
