##Changes needed - change all xyl1 to xyl2/3




library(ape)
library(stringr)
xyl1.fastTree<-read.tree("~/xylose_optimization_project/spring_2021/hmmer_searches/all_hmmer_result_fastTrees/orf_finder_xyl1.pep_fastTree.tree")
allclades<-read.delim("~/xylose_optimization_project/spring_2021/subgroup_consensus/xyl1/xyl1_ko.txt", stringsAsFactors = FALSE, header = FALSE, na.strings=c("","NA"))

data.frame(xyl1.fastTree$tip.label)->tree_tips
tree_tips$xyl1.fastTree.tip.label<-gsub("\\|", "_", tree_tips$xyl1.fastTree.tip.label)
tree_tips$xyl1.fastTree.tip.label<-gsub("\\+", "_", tree_tips$xyl1.fastTree.tip.label)
tree_tips$xyl1.fastTree.tip.label<-gsub("-", "_", tree_tips$xyl1.fastTree.tip.label)
tree_tips$xyl1.fastTree.tip.label<-gsub("__", "_", tree_tips$xyl1.fastTree.tip.label)
allclades$V1<-gsub("\\|", "_", allclades$V1)
allclades$V1<-gsub("\\+", "_", allclades$V1)
allclades$V1<-gsub(":", "_", allclades$V1)
allclades$V1<-gsub("-", "_", allclades$V1)
allclades$V1<-gsub("__", "_", allclades$V1)


#blastKOALA annoyingly truncates long annotations. I need to run a grep loop to
#reassign full annotations. 

for(i in 1:nrow(allclades)){
  treetipsmatch<-which(grepl(allclades$V1[i], tree_tips[,1]))
  if(length(treetipsmatch)==1){
    allclades$V1[i]<-tree_tips[treetipsmatch,1]
  }
  else(print(i))
}

#rename tree tips
xyl1.fastTree$tip.label<-tree_tips$xyl1.fastTree.tip.label
#name coloumns of all clades
colnames(allclades)<-c("annotation", "KO")


allclades$TYPE<-"branch"
allclades$WHAT<-NA
allclades$COLOR<-"#ff0000"
allclades$WIDTH_OR_SIZE_FACTOR<-1
allclades$STYLE<-"normal"

KO_number<-length(unique(allclades$KO))
#18 colors - add more hex codes if you have more than 18 KOs! 
colors<-c("#FF0000","#800000","#FFFF00","#808000","#00FF00","#008000","#00FFFF",
          "#008080","#0000FF","#000080","#FF00FF","#800080", "#C0C0C0", "#808080", 
          "#6495ed","#eee8aa", "#696969", "#663399")

NA.col<-"#000000"          
          
for (i in 1:KO_number){
  allclades$COLOR[which(allclades$KO==(unique(allclades$KO))[i])]<-colors[i]
}
allclades$COLOR[which(is.na(allclades$KO))]<-NA.col

write.table(allclades[c(1,2)], "~/xylose_optimization_project/spring_2021/subgroup_consensus/xyl1/xyl1_KO_annotations.txt", sep="\t", quote=FALSE, row.names=FALSE)

write.table(allclades[c(1,3:7)], "~/xylose_optimization_project/spring_2021/subgroup_consensus/xyl1/xyl1_KO_annotations_iTOLfile.txt", sep="\t", quote=FALSE, row.names=FALSE)

key<-unique(allclades[c(2,5)]) 

write.table(key, "~/xylose_optimization_project/spring_2021/subgroup_consensus/xyl1/xyl1_KO_annotations_iTOLkey.txt", sep="\t", quote=FALSE, row.names=FALSE)

write.tree(xyl1.fastTree, "~/xylose_optimization_project/spring_2021/subgroup_consensus/xyl1/xyl1.hits.TreeFileForItol.tree")

#########
#read in clades
clade1<-read.delim("~/xylose_optimization_project/spring_2021/subgroup_consensus/xyl1/clades/clade1.txt", header=FALSE)
clade1$clade<-1
clade2<-read.delim("~/xylose_optimization_project/spring_2021/subgroup_consensus/xyl1/clades/clade2.txt", header=FALSE)
clade2$clade<-2
clade3<-read.delim("~/xylose_optimization_project/spring_2021/subgroup_consensus/xyl1/clades/clade3.txt", header=FALSE)
clade3$clade<-3
clade4<-read.delim("~/xylose_optimization_project/spring_2021/subgroup_consensus/xyl1/clades/clade4.txt", header=FALSE)
clade4$clade<-4
clade5<-read.delim("~/xylose_optimization_project/spring_2021/subgroup_consensus/xyl1/clades/clade5.txt", header=FALSE)
clade5$clade<-5
clade6<-read.delim("~/xylose_optimization_project/spring_2021/subgroup_consensus/xyl1/clades/clade6.txt", header=FALSE)
clade6$clade<-6
clade7<-read.delim("~/xylose_optimization_project/spring_2021/subgroup_consensus/xyl1/clades/clade7.txt", header=FALSE)
clade7$clade<-7
clade8<-read.delim("~/xylose_optimization_project/spring_2021/subgroup_consensus/xyl1/clades/clade8.txt", header=FALSE)
clade8$clade<-8
clade9<-read.delim("~/xylose_optimization_project/spring_2021/subgroup_consensus/xyl1/clades/clade9.txt", header=FALSE)
clade9$clade<-9
clade10<-read.delim("~/xylose_optimization_project/spring_2021/subgroup_consensus/xyl1/clades/clade10.txt", header=FALSE)
clade10$clade<-10
clade11<-read.delim("~/xylose_optimization_project/spring_2021/subgroup_consensus/xyl1/clades/clade11.txt", header=FALSE)
clade11$clade<-11
clade12<-read.delim("~/xylose_optimization_project/spring_2021/subgroup_consensus/xyl1/clades/clade12.txt", header=FALSE)
clade12$clade<-12
clade13<-read.delim("~/xylose_optimization_project/spring_2021/subgroup_consensus/xyl1/clades/clade13.txt", header=FALSE)
clade13$clade<-13
clade14<-read.delim("~/xylose_optimization_project/spring_2021/subgroup_consensus/xyl1/clades/clade14.txt", header=FALSE)
clade14$clade<-14
clade15<-read.delim("~/xylose_optimization_project/spring_2021/subgroup_consensus/xyl1/clades/clade15.txt", header=FALSE)
clade15$clade<-15
clade16<-read.delim("~/xylose_optimization_project/spring_2021/subgroup_consensus/xyl1/clades/clade16.txt", header=FALSE)
clade16$clade<-16

catclades<-rbind(clade1, clade2, clade3, clade4, clade5, clade6,
                 clade7, clade8, clade9, clade10, clade11, clade12, 
                 clade13, clade14, clade15, clade16)


