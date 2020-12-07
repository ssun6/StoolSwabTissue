#Author: Shan Sun

#BioLockJ configuration: Ali Sorgen
#Date: 12-02-20
#Description: Perform permanova test and generate PCoA plots.

#Libraries
library(stringr)
library(vegan)

rm(list = ls())

pipeRoot = dirname(dirname(getwd()))
moduleDir <- dirname(getwd())
input = file.path(pipeRoot,"input/")
output = file.path(moduleDir,"output/")
pathToMeta<-paste0(pipeRoot,"/",str_subset(dir(pipeRoot), "Metadata"),"/output/")

sample_type=c("stool","swab","tissue")

#taxa
bug_tab=read.table(file=paste0(input, "/counts/van_re_metaphlan2.txt"),header=T,row.names=1,sep="\t")

map7=read.table(file=paste0(pathToMeta, "metadata_all.csv"),header=T,sep=",")

map7=map7[match(colnames(bug_tab),map7$sample_id),]

col11=c("red","blue","green","orange","purple","pink","black","lightblue","lightgreen","hotpink","cyan","orchid","tan","grey","gold")
tax_level=sapply(rownames(bug_tab),function(i){length(strsplit(i,"\\|")[[1]])})
tax_levels=c("kingdom","phylum","class","order","family","genus","species","strain")
permanova_mat=matrix(nrow=8,ncol=4)

for (n in 6:6){
  bug_tab2=bug_tab[tax_level==n,]
  bug_tab3=t(t(bug_tab2)/colSums(bug_tab2))*mean(colSums(bug_tab2))
  
  print(dim(bug_tab3))
  fname=paste(output,tax_levels[n],".csv",sep="")
  write.csv(bug_tab3,fname)
  
  tab_s=bug_tab3
  tab_s1=log10(tab_s+1) #the average sequencing depth based on pathway data
  
  shannon=vegan::diversity(tab_s,index = "shannon", MARGIN = 2, base = exp(1))
  wilcox.test(shannon[map7$sample_type!="tissue"]~map7$sample_type[map7$sample_type!="tissue"])
  wilcox.test(shannon[map7$sample_type!="stool"]~map7$sample_type[map7$sample_type!="stool"])
  wilcox.test(shannon[map7$sample_type!="swab"]~map7$sample_type[map7$sample_type!="swab"])
  
  ## Figure 1a
  pdf(paste0(output, "tax_alpha_box_l",n,".pdf"),width=5,height=10)
  par(mar=c(5,5,5,5),mfrow=c(1,1),bty="n")
  type_new=factor(map7$sample_type,levels=c("stool","swab","tissue"))
  boxplot(shannon~type_new, lwd = 2, xlab = "", ylab = 'Shannon Index',main="",cex.lab=2,cex.main=2,cex.axis=2,border=col11,las=2)
  stripchart(shannon~type_new, vertical = TRUE,  method = "jitter", add = TRUE, pch = 16, col = col11)
  dev.off()
  
  #all three samples
  gen_pcoa=capscale(t(tab_s1)~1,distance="bray")
  mds=list()
  mds[[1]]=round((gen_pcoa$CA$eig/sum(gen_pcoa$CA$eig))[1]*100,2)
  mds[[2]]=round((gen_pcoa$CA$eig/sum(gen_pcoa$CA$eig))[2]*100,2)
  mds[[3]]=round((gen_pcoa$CA$eig/sum(gen_pcoa$CA$eig))[3]*100,2)
  mds[[4]]=round((gen_pcoa$CA$eig/sum(gen_pcoa$CA$eig))[4]*100,2)
  
  ado1=adonis(t(tab_s1)~factor(map7$sample_type),permutations=999)
  permanova_mat[n,1]=ado1$aov.tab[1,5]
  permanova_mat[n,2]=ado1$aov.tab[1,6]
  
  permanova_meta=matrix(nrow=5,ncol=6)
  x1=1
  
  for (y2 in c(8,9,11,13,23)){
    for(y1 in 1:3){
      tab_s1_t=t(tab_s1[,map7$sample_type==sample_type[y1]&!is.na(map7[,y2])])
      map_y2=map7[map7$sample_type==sample_type[y1]&!is.na(map7[,y2]),y2]
      permanova_meta[x1,(y1*2-1):(y1*2)]=as.numeric(adonis(tab_s1_t~factor(map_y2),permutations=999)$aov.tab[1,5:6])
    }
    x1=x1+1
  }
  write.csv(permanova_meta,file=paste0(output, "adonis_meta_l",n,".csv",sep=""))
  
  ## Figure 1b
  fname=paste0(output, "tax_pcoa12_type_l",n,".pdf",sep="")
  pdf(fname,width=8,height=8)
  par(mfrow=c(1,1),mar=c(5,5,5,5))
  pcoa12=ordiplot(gen_pcoa,display="sites",choices=c(1,2),type="none",cex.lab=1.5,xlab=paste0("PCoA1 (",mds[[1]],"%)"), ylab=paste0("PCoA2 (",mds[[2]],"%)"), main = paste0("PERMANOVA R2=",round(ado1$aov.tab$R2[1], 3), ", P=", ado1$aov.tab$`Pr(>F)`[1]))
  col3=col11[factor(map7$sample_type)]
  points(pcoa12,"sites",col=adjustcolor(col3, alpha.f = 0.2),pch=16,cex=2)
  legend("topright",levels(factor(map7$sample_type)), cex=1.8, bty="n", col=col11[1:length(levels(factor(map7$sample_type)))], pch=20)
  dev.off()
  
  ## Figure 1b insets
  pdf(paste0(output, "tax_pcoa1234_type_box_l",n,".pdf"),width=8,height=8)
  par(mar=c(5,5,5,5),mfrow=c(1,4),bty="n")
  type_new=factor(map7$sample_type,levels=c("stool","swab","tissue"))
  boxplot(summary(gen_pcoa)$sites[,1]~type_new, lwd = 2, xlab = "", ylab = 'distance',main="",cex.lab=1,cex.main=2,cex.axis=1,col="white",border=col11,las=2)
  stripchart(summary(gen_pcoa)$sites[,1]~type_new, vertical = TRUE,  method = "jitter", add = TRUE, pch = 16, col = col11)
  boxplot(summary(gen_pcoa)$sites[,2]~type_new, lwd = 2, xlab = "", ylab = 'distance',main="",cex.lab=1,cex.main=2,cex.axis=1,col="white",border=col11,las=2)
  stripchart(summary(gen_pcoa)$sites[,2]~type_new, vertical = TRUE,  method = "jitter", add = TRUE, pch = 16, col = col11)
  boxplot(summary(gen_pcoa)$sites[,3]~type_new, lwd = 2, xlab = "", ylab = 'distance',main="",cex.lab=1,cex.main=2,cex.axis=1,col="white",border=col11,las=2)
  stripchart(summary(gen_pcoa)$sites[,3]~type_new, vertical = TRUE,  method = "jitter", add = TRUE, pch = 16, col = col11)
  boxplot(summary(gen_pcoa)$sites[,4]~type_new, lwd = 2, xlab = "", ylab = 'distance',main="",cex.lab=1,cex.main=2,cex.axis=1,col="white",border=col11,las=2)
  stripchart(summary(gen_pcoa)$sites[,4]~type_new, vertical = TRUE,  method = "jitter", add = TRUE, pch = 16, col = col11)
  dev.off()
  

  #only stool and swab
  gen_pcoa=capscale(t(tab_s1[,!map7$sample_type=="tissue"])~1,distance="bray")
  map_s=map7[!map7$sample_type=="tissue",]
  mds=list()
  mds[[1]]=round((gen_pcoa$CA$eig/sum(gen_pcoa$CA$eig))[1]*100,2)
  mds[[2]]=round((gen_pcoa$CA$eig/sum(gen_pcoa$CA$eig))[2]*100,2)
  mds[[3]]=round((gen_pcoa$CA$eig/sum(gen_pcoa$CA$eig))[3]*100,2)
  mds[[4]]=round((gen_pcoa$CA$eig/sum(gen_pcoa$CA$eig))[4]*100,2)
  
  ado1=adonis(t(tab_s1[,!map7$sample_type=="tissue"])~factor(map_s$sample_type),permutations=999)
  permanova_mat[n,3]=ado1$aov.tab[1,5]
  permanova_mat[n,4]=ado1$aov.tab[1,6]
  
  ## Figure 1c
  fname=paste0(output, "two_tax_pcoa12_type_l",n,".pdf",sep="")
  pdf(fname,width=8,height=8)
  par(mfrow=c(1,1),mar=c(5,5,5,5))
  pcoa12=ordiplot(gen_pcoa,display="sites",choices=c(1,2),type="none",cex.lab=1.5,xlab=paste0("PCoA1 (",mds[[1]],"%)"),ylab=paste0("PCoA2 (",mds[[2]],"%)"), main = paste0("PERMANOVA R2=",round(ado1$aov.tab$R2[1], 3), ", P=", ado1$aov.tab$`Pr(>F)`[1]))
  col3=col11[factor(map_s$sample_type)]
  points(pcoa12,"sites",col=adjustcolor(col3, alpha.f = 0.2),pch=16,cex=2)
  legend("topright",levels(factor(map_s$sample_type)), cex=1.8, bty="n", col=col11[1:length(levels(factor(map_s$sample_type)))], pch=20)
  dev.off()
  
  
  ## Figure 1c insets
  pdf(paste0(output, "two_tax_pcoa1234_type_box_ss_l",n,".pdf"),width=8,height=8)
  par(mar=c(5,5,5,5),mfrow=c(1,4),bty="n")
  type_new=factor(map7$sample_type[!map7$sample_type=="tissue"],levels=c("stool","swab"))
  boxplot(summary(gen_pcoa)$sites[,1]~type_new, lwd = 2, xlab = "", ylab = 'distance',main="",cex.lab=1,cex.main=2,cex.axis=1,col="white",border=col11,las=2)
  stripchart(summary(gen_pcoa)$sites[,1]~type_new, vertical = TRUE,  method = "jitter", add = TRUE, pch = 16, col = col11)
  boxplot(summary(gen_pcoa)$sites[,2]~type_new, lwd = 2, xlab = "", ylab = 'distance',main="",cex.lab=1,cex.main=2,cex.axis=1,col="white",border=col11,las=2)
  stripchart(summary(gen_pcoa)$sites[,2]~type_new, vertical = TRUE,  method = "jitter", add = TRUE, pch = 16, col = col11)
  boxplot(summary(gen_pcoa)$sites[,3]~type_new, lwd = 2, xlab = "", ylab = 'distance',main="",cex.lab=1,cex.main=2,cex.axis=1,col="white",border=col11,las=2)
  stripchart(summary(gen_pcoa)$sites[,3]~type_new, vertical = TRUE,  method = "jitter", add = TRUE, pch = 16, col = col11)
  boxplot(summary(gen_pcoa)$sites[,4]~type_new, lwd = 2, xlab = "", ylab = 'distance',main="",cex.lab=1,cex.main=2,cex.axis=1,col="white",border=col11,las=2)
  stripchart(summary(gen_pcoa)$sites[,4]~type_new, vertical = TRUE,  method = "jitter", add = TRUE, pch = 16, col = col11)
  dev.off()
  
}

colnames(permanova_mat)=c("R2","P","2R2","2P")
rownames(permanova_mat)=tax_levels
write.csv(permanova_mat,file=paste0(output, "tax_permanova.csv"))
