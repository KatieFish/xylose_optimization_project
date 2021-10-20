data<-subset(estAI_vals, subset = estAI_vals$gene=="xyl1")
data<-unique(data[c(1,6)])
data<-merge(data, growth_data[c(1,2)], by="all_taxa")
#removing outlier growth rates 
data<-data[which(data$Growth.Rate<2),]
#removing 0
data<-data[which(data$Growth.Rate>0),]


tips<-data.frame(tree$tip.label)
colnames(tips)<-"all_taxa"
tips$phyloorder<-c(1:nrow(tips))
x<-merge(tips, data[c(1,2)], by="all_taxa", all=TRUE)
x<-x[1:332,]
x<-x[order(x$phyloorder),]

x$multiple_copies->copies
names(copies)<-x$all_taxa

pdf(height=60, width=20, file="look.pdf")
plot(tree)
tiplabels(copies)
dev.off()


pairs<-read.delim("~/Xyl_project_Fall_2021/spp_pairs_singl_multi_copy.txt", stringsAsFactors = FALSE)
pairs$pair<-c(1:nrow(pairs))
pairs<-pairs[c(3,1,2)]
x<-melt(pairs, id.vars = "pair")
x$gr<-NA
for(i in 1:nrow(x)){
  x$gr[i]<-data$Growth.Rate[which(data$all_taxa==x$value[i])]  
}

ggplot(x, aes(x=variable, y=gr))+
  geom_point(aes(col=factor(pair)), size=3)+
  geom_line(aes(group = pair))+
  geom_label(aes(label=value),position = position_nudge(y = -0.01) )+
  theme(legend.position = "none")+
  ylab("Xylose growth rate")+
  xlab("XYL1 copy number")+
  theme_bw()

quartz.save("copyno_gr_sisterspp.pdf", type="pdf")




