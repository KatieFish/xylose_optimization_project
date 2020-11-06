##Testing for phylogenetic signal in xylose gene codon optimization
#of xyl genes. 

#332 newick tree from Shen et al. 
library(ape)
tree<-read.tree("~/xylose_optimization_project/data/332_Newick_tree.txt")
tree$tip.label<-tolower(tree$tip.label)
tree$tip.label<-str_replace(string = tree$tip.label, pattern = "_", replacement = " ")
tree$tip.label[234]<-"lachancea fantastica"
tree_taxa<-tree$tip.label
which(!max_estAI_df$all_taxa %in% tree_taxa)-> not_in_data
max_estAI_df[not_in_data, 1]
#Ok - the spp. we have that are NOT in the tree mostly seem to be due to renaming. 
#looking up possible new spp. names using second name of sp. 
for (i in 1:length(not_in_data)){
  max_estAI_df[(not_in_data[i]), 1]->old_sp
  spName<-strsplit(old_sp, " ")[[1]][2]
  possiblechange<-tree$tip.label[which(grepl(spName, tree$tip.label))]
  if (length(possiblechange)==1){
    max_estAI_df[(not_in_data[i]), 1]<-possiblechange
  }
}
#ok - now there are less (11). Dealing with these one by one. 
max_estAI_df[(not_in_data[1]), 1]<-"metschnikowia dekortorum"
max_estAI_df[(not_in_data[2]), 1]<-"ogataea populialbae"
max_estAI_df[(not_in_data[3]), 1]<-"ogataea philodendri"
max_estAI_df[(not_in_data[5]), 1]<-"blastobotrys raffinosifermentans"
max_estAI_df[(not_in_data[6]), 1]<-"candida castellii"
#stopped here
max_estAI_df[(not_in_data[7]), 1]<-"ogataea populialbae"
max_estAI_df[(not_in_data[8]), 1]<-"ogataea philodendri"
max_estAI_df[(not_in_data[9]), 1]<-"blastobotrys raffinosifermentans"


