library(vegan)
library(ggplot2)
library(tidyr)
library(ggthemes)
library(reshape2)
library(gridExtra)
library(RGraphics)
library(pheatmap)
library(ggrepel)
library(ggpubr)

setwd("/Users/shansun/Google\ Drive/Vanderbilt/SST_test/")

col11=c("red","blue","green","orange","purple","pink","black","lightblue","lightgreen","hotpink","cyan","orchid","tan","grey","gold")
draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...))
  return(res)}

## 'Overwrite' default draw_colnames with your own version 
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))

map1=read.csv(file="metadata/study\ ID_sample\ ID_link.csv",header=T)
which(map1$visit=="") #726
map1[726,5:7]=rep(NA,3)#change missing values to NA
map1$name_v=paste(map1$study_id,map1$visit)

map2=read.csv(file="metadata/metadata_20190503.csv",header=T)
map3=read.csv(file="metadata/IHC_data_20190503.csv",header=T)
map3$name_v=paste(map3$studyid,map3$timepoints)

bug_tab=read.table(file="van_re_metaphlan2.txt",header=T,row.names=1,sep="\t")
map4=map1[match(colnames(bug_tab),map1$sample_id),]
map5=map2[match(map4$study_id,map2$study_id),]
map6=map3[match(map4$name_v,map3$name_v),]

map7=cbind(map4[,-8],map5[,-1],map6[,-c(1,2,10)])
map7$treat_visit=paste(map7$treatment,map7$visit,sep="_")
map7$treat_visit=factor(map7$treat_visit,levels=levels(factor(map7$treat_visit))[c(2,1,5,4)])


antib=read.csv(file="metadata/Antibiotic_data_20201022.csv",header=T)
antib1=data.frame(sapply(by(antib$pq_16,antib$study_id,mean),identity))
colnames(antib1)="antib"
antib1$antibiotics=antib1$antib
antib1$antibiotics[which(antib1$antib<2)]="yes"
antib1$antibiotics[which(antib1$antib==2)]="No"
antib2=antib1[match(map7$study_id,rownames(antib1)),]

days=read.csv(file="metadata/days.csv",header=T)
days1=days[match(map7$study_id,days$study_id),]

map7=cbind(map7,antib2,days1)
#write.csv(map7,file="metadata/metadata_all.csv")

sample_type=c("stool","swab","tissue")

#taxa
bug_tab=read.table(file="van_re_metaphlan2.txt",header=T,row.names=1,sep="\t")
#bug_tab=read.table(file="van_re_kraken2.txt",header=T,row.names=1,sep="\t",quote="")
map7=map7[match(colnames(bug_tab),map7$sample_id),]

tax_level=sapply(rownames(bug_tab),function(i){length(strsplit(i,"\\|")[[1]])})
tax_levels=c("kingdom","phylum","class","order","family","genus","species","strain")
permanova_mat=matrix(nrow=8,ncol=4)
for (n in 6:6){
  bug_tab2=bug_tab[tax_level==n,]
  bug_tab3=t(t(bug_tab2)/colSums(bug_tab2))*mean(colSums(bug_tab2))
  
  print(dim(bug_tab3))
  fname=paste("results/",tax_levels[n],".csv",sep="")
  write.csv(bug_tab3,fname)
  tab_s=bug_tab3
  tab_s1=log10(tab_s+1) #the average sequencing depth based on pathway data
  
  shannon=vegan::diversity(tab_s,index = "shannon", MARGIN = 2, base = exp(1))
  wilcox.test(shannon[map7$sample_type!="tissue"]~map7$sample_type[map7$sample_type!="tissue"])
  wilcox.test(shannon[map7$sample_type!="stool"]~map7$sample_type[map7$sample_type!="stool"])
  wilcox.test(shannon[map7$sample_type!="swab"]~map7$sample_type[map7$sample_type!="swab"])
  
  pdf(paste0("results/tax_alpha_box_l",n,".pdf"),width=5,height=8)
  par(mar=c(5,5,5,5),mfrow=c(1,1),bty="n")
  type_new=factor(map7$sample_type,levels=c("stool","swab","tissue"))
  boxplot(shannon~type_new, lwd = 2, ylab = 'Shannon Index',main="",cex.lab=2,cex.main=2,cex.axis=2,border=col11,las=2)
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
  write.csv(permanova_meta,file=paste0("results/adonis_meta_l",n,".csv",sep=""))
  
  fname=paste0("results/tax_pcoa12_type_l",n,".pdf",sep="")
  pdf(fname,width=8,height=8)
  par(mfrow=c(1,1),mar=c(5,5,5,5))
  pcoa12=ordiplot(gen_pcoa,display="sites",choices=c(1,2),type="none",cex.lab=1.5,xlab=paste0("PCoA1 (",mds[[1]],"%)"),ylab=paste0("PCoA2 (",mds[[2]],"%)"))
  col3=col11[factor(map7$sample_type)]
  points(pcoa12,"sites",col=adjustcolor(col3, alpha.f = 0.2),pch=16,cex=2)
  legend("topright",levels(factor(map7$sample_type)), cex=1.8, bty="n", col=col11[1:length(levels(factor(map7$sample_type)))], pch=20)
  dev.off()
  
  
  pdf(paste0("results/tax_pcoa1234_type_box_l",n,".pdf"),width=8,height=8)
  par(mar=c(5,5,5,5),mfrow=c(1,4),bty="n")
  type_new=factor(map7$sample_type,levels=c("stool","swab","tissue"))
  boxplot(summary(gen_pcoa)$sites[,1]~type_new, lwd = 2, ylab = 'distance',main="",cex.lab=1,cex.main=2,cex.axis=1,border=col11,las=2)
  stripchart(summary(gen_pcoa)$sites[,1]~type_new, vertical = TRUE,  method = "jitter", add = TRUE, pch = 16, col = col11)
  boxplot(summary(gen_pcoa)$sites[,2]~type_new, lwd = 2, ylab = 'distance',main="",cex.lab=1,cex.main=2,cex.axis=1,border=col11,las=2)
  stripchart(summary(gen_pcoa)$sites[,2]~type_new, vertical = TRUE,  method = "jitter", add = TRUE, pch = 16, col = col11)
  boxplot(summary(gen_pcoa)$sites[,3]~type_new, lwd = 2, ylab = 'distance',main="",cex.lab=1,cex.main=2,cex.axis=1,border=col11,las=2)
  stripchart(summary(gen_pcoa)$sites[,3]~type_new, vertical = TRUE,  method = "jitter", add = TRUE, pch = 16, col = col11)
  boxplot(summary(gen_pcoa)$sites[,4]~type_new, lwd = 2, ylab = 'distance',main="",cex.lab=1,cex.main=2,cex.axis=1,border=col11,las=2)
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
  
  fname=paste0("results/two_tax_pcoa12_type_l",n,".pdf",sep="")
  pdf(fname,width=8,height=8)
  par(mfrow=c(1,1),mar=c(5,5,5,5))
  pcoa12=ordiplot(gen_pcoa,display="sites",choices=c(1,2),type="none",cex.lab=1.5,xlab=paste0("PCoA1 (",mds[[1]],"%)"),ylab=paste0("PCoA2 (",mds[[2]],"%)"))
  col3=col11[factor(map_s$sample_type)]
  points(pcoa12,"sites",col=adjustcolor(col3, alpha.f = 0.2),pch=16,cex=2)
  legend("topright",levels(factor(map_s$sample_type)), cex=1.8, bty="n", col=col11[1:length(levels(factor(map_s$sample_type)))], pch=20)
  dev.off()
  
  
  pdf(paste0("results/two_tax_pcoa12_type_box_ss_l",n,".pdf"),width=8,height=8)
  par(mar=c(5,5,5,5),mfrow=c(1,4),bty="n")
  type_new=factor(map7$sample_type[!map7$sample_type=="tissue"],levels=c("stool","swab"))
  boxplot(summary(gen_pcoa)$sites[,1]~type_new, lwd = 2, ylab = 'distance',main="",cex.lab=1,cex.main=2,cex.axis=1,border=col11,las=2)
  stripchart(summary(gen_pcoa)$sites[,1]~type_new, vertical = TRUE,  method = "jitter", add = TRUE, pch = 16, col = col11)
  boxplot(summary(gen_pcoa)$sites[,2]~type_new, lwd = 2, ylab = 'distance',main="",cex.lab=1,cex.main=2,cex.axis=1,border=col11,las=2)
  stripchart(summary(gen_pcoa)$sites[,2]~type_new, vertical = TRUE,  method = "jitter", add = TRUE, pch = 16, col = col11)
  boxplot(summary(gen_pcoa)$sites[,3]~type_new, lwd = 2, ylab = 'distance',main="",cex.lab=1,cex.main=2,cex.axis=1,border=col11,las=2)
  stripchart(summary(gen_pcoa)$sites[,3]~type_new, vertical = TRUE,  method = "jitter", add = TRUE, pch = 16, col = col11)
  boxplot(summary(gen_pcoa)$sites[,4]~type_new, lwd = 2, ylab = 'distance',main="",cex.lab=1,cex.main=2,cex.axis=1,border=col11,las=2)
  stripchart(summary(gen_pcoa)$sites[,4]~type_new, vertical = TRUE,  method = "jitter", add = TRUE, pch = 16, col = col11)
  dev.off()
  
}
colnames(permanova_mat)=c("R2","P","2R2","2P")
rownames(permanova_mat)=tax_levels
write.csv(permanova_mat,file="results/tax_permanova.csv")


library(lme4)
library(lmerTest)
library(nlme)

#pairwise two types
n=6
bug_tab2=bug_tab[tax_level==n,]
map7=map7[match(colnames(bug_tab2),map7$sample_id),]

bug_tab2_1=t(t(bug_tab2)/colSums(bug_tab2))
gen_type_1=data.frame(sapply(by(t(bug_tab2_1),map7$sample_type,colMeans),identity))
gen_type_1=gen_type_1[order(gen_type_1[,2],decreasing = T),]

bug_tab3=t(t(bug_tab2)/colSums(bug_tab2))*mean(colSums(bug_tab2))

bug_tab4=log10(bug_tab3+1)
bug_tab5=bug_tab4[apply(bug_tab3!=0,1,sum)/1397>0.1,]
dim(bug_tab3) #210 1397
dim(bug_tab4)#210 1397
dim(bug_tab5)#60 1397
sample_type=c("stool","swab","tissue")
glm_tax=list()
tab_n1=bug_tab5
map_n1=map7
rownames(tab_n1)=gsub("\\|","_",rownames(tab_n1))
rownames(tab_n1)=sapply(strsplit(rownames(tab_n1),"g__"),"[[",2)

bacteriaMeta1=cbind(t(tab_n1),map_n1)
which(is.na(bacteriaMeta1$visit))#970 remove
bacteriaMeta1=cbind(t(tab_n1),map_n1)[-which(is.na(bacteriaMeta1$visit)),]#remove one sample with missing data

pVals=matrix(ncol=30,nrow=nrow(tab_n1))
for (j in 1:3){
  type1=sample_type[j]
  bacteriaMeta=bacteriaMeta1[bacteriaMeta1$sample_type!=type1,]
  bacteriaMeta$sample_type=droplevels(bacteriaMeta$sample_type)
  bacteriaMeta=bacteriaMeta[!is.na(bacteriaMeta$antibiotics)&!is.na(bacteriaMeta$visit),]
  for( i in 1:dim(tab_n1)[1]){
    print(i)
    if(sd(tab_n1[i,])==0){
      next
    }
    model <- as.formula(paste(rownames(tab_n1)[i],"~","factor(sample_type)+factor(treatment)*factor(visit)+factor(antibiotics)+age+sex+NSAIDS_use+BMI"))
    
    simpleMod <- try(gls(model,method="REML",data=bacteriaMeta))
    mixedMod <- try(lme(model,method="REML",random=~1|study_id,data=bacteriaMeta))
    if(class(simpleMod)=="try-error" | class(mixedMod)=="try-error"){
      next
    }
    
    for (m in 2:10){
      pVals[i,10*(j-1)+m-1] <- -log10(summary(mixedMod)$tTable[m,5])*sign(summary(mixedMod)$tTable[m,4])
    }
    pVals[i,10*(j-1)+10]  <- -log10(anova(simpleMod,mixedMod)$"p-value"[2])
  }
}

names1=c(rownames(summary(simpleMod)$tTable)[-1],"subject")
names1[1]="sample_type"
names1[2]="treatment"
names1[3]="visit"
names1[4]="antibiotics"
names1[6]="sex"
names1[9]="interaction"
colnames(pVals)=c(paste(names1,sample_type[1],sep="_"),paste(names1,sample_type[2],sep="_"),paste(names1,sample_type[3],sep="_"))
rownames(pVals)=sapply(strsplit(rownames(bug_tab5),"g__"),"[[",2)
pVals=data.frame(na.omit(pVals))
write.csv(pVals,file=paste0("results/l",n,"tax_lmer_P1_pair.csv"))

FDRs=10^-abs(pVals)
FDRs=apply(FDRs,1,as.character)
FDRs=apply(FDRs,1,as.numeric)
FDRs=matrix(p.adjust(FDRs,method = "fdr"),ncol=30)
colnames(FDRs)=colnames(pVals)
rownames(FDRs)=rownames(pVals)
FDRs=data.frame(na.omit(FDRs))

length(which(FDRs[,1]<0.05))#51
length(which(FDRs[,11]<0.05))#53
length(which(FDRs[,21]<0.05))#35
length(which(apply(FDRs[,c(1,11,21)],1,min)<0.05)) #56

rownames(FDRs)[which(FDRs[,1]<0.05)]
rownames(FDRs)[which(FDRs[,11]<0.05)]
rownames(FDRs)[which(FDRs[,21]<0.05)]

setdiff(rownames(FDRs)[which(FDRs[,1]<0.05)],rownames(FDRs)[which(FDRs[,11]<0.05)])#42


P_min=apply(FDRs[,c(1,11,21)],1,function(i){min(abs(i))})#all these taxa 
tab_n2=bug_tab3[apply(bug_tab3!=0,1,sum)/1397>0.1,]
rownames(tab_n2)=gsub("\\|","_",rownames(tab_n2))
rownames(tab_n2)=sapply(strsplit(rownames(tab_n2),"g__"),"[[",2)

gen_type1=data.frame(sapply(by(t(tab_n2[which(P_min<0.05),]),map7$sample_type,colMeans),identity))
gen_type1[order(gen_type1[,1]/gen_type1[,2],decreasing = T),]

length(which(P_min<0.05))
gen_type=t(scale(t(gen_type1)))
pdf("results/gen_heat.pdf",width=8,height=8)
pheatmap(gen_type,border_color = "white",fontsize_col =15)
dev.off()

pdf("results/gen_heat.pdf",width=8,height=8)
pheatmap(gen_type,border_color = "white",fontsize_col =15)
dev.off()



#one sample
pVals=matrix(ncol=27,nrow=dim(tab_n1)[1])
for (j in 1:3){
  type1=sample_type[j]
  bacteriaMeta=bacteriaMeta1[bacteriaMeta1$sample_type==type1,]
  bacteriaMeta$sample_type=droplevels(bacteriaMeta$sample_type)
  bacteriaMeta=bacteriaMeta[!is.na(bacteriaMeta$antibiotics)&!is.na(bacteriaMeta$visit),]
  
  for( i in 1:dim(tab_n1)[1]){
    print(i)
    if(sd(tab_n1[i,])==0){
      next
    }
    model <- as.formula(paste(rownames(tab_n1)[i],"~","factor(treatment)*factor(visit)+factor(antibiotics)+age+sex+NSAIDS_use+BMI"))
    
    simpleMod <- try(gls(model,method="REML",data=bacteriaMeta))
    mixedMod <- try(lme(model,method="REML",random=~1|study_id,data=bacteriaMeta))
    
    if(class(simpleMod)=="try-error"|class(mixedMod)=="try-error"){
      for (m in 2:9){
        pVals[i,9*(j-1)+m-1] =NA
      }
      pVals[i,9*(j-1)+9]=NA
    }else{
      for (m in 2:9){
        pVals[i,9*(j-1)+m-1] = -log10(summary(mixedMod)$tTable[m,5])*sign(summary(mixedMod)$tTable[m,4])
      }
      pVals[i,9*(j-1)+9]  = -log10(anova(simpleMod,mixedMod)$"p-value"[2])
    }
  }
}
pVals=data.frame(pVals)
names1=c(rownames(summary(simpleMod)$tTable)[-1],"subject")
names1[1]="treatment"
names1[2]="visit"
names1[3]="antibiotics"
names1[5]="sex"
names1[8]="interaction"
colnames(pVals)=c(paste(names1,sample_type[1],sep="_"),paste(names1,sample_type[2],sep="_"),paste(names1,sample_type[3],sep="_"))
rownames(pVals)=sapply(strsplit(rownames(bug_tab5),"g__"),"[[",2)
write.csv(pVals,file="results/onesample_l6_tax_lmer.csv")

pVals=read.csv(file="results/onesample_l6_tax_lmer.csv",row.names=1)
cor_meta=matrix(nrow=5,ncol=6)
j=0
require("ggrepel")
for (n in c(4,5,6,7,3)){
  j=j+1
  p1=list()
  fname=paste0("results/gen_tax_onesample_",names1[n],"_pair.pdf")
  pdf(fname,width=20,height=8)
  list1=c(0,1,2)
  for (m in 3:1){
    cor=cor.test(pVals[,9*list1[-m][1]+n],pVals[,9*list1[-m][2]+n],method="spearman")$estimate
    pv=as.numeric(cor.test(pVals[,9*list1[-m][1]+n],pVals[,9*list1[-m][2]+n],method="spearman")$p.value)
    if (pv==0){
      pv="< 2.2e-16"
    }else if (pv<0.001){
      pv=paste("=",formatC(pv, format = "e", digits = 2))
    }else{
      pv=paste("=",round(pv,3))
    }
    
    cor_meta[j,c(3:1)[m]]=cor
    cor_meta[j,c(3:1)[m]+3]=pv
    
    if (m==3){
      t1=paste(names1[n],":",paste(sample_type[list1[-m]+1],collapse=" vs "),"\n","rho =",round(cor,3),"P",pv)
    }else{
      t1=paste(paste(sample_type[list1[-m]+1],collapse=" vs "),"\n","rho =",round(cor,3),"P",pv)
    }
    
    x1=paste0("-log10(P) ",sample_type[list1[-m][1]+1])
    y1=paste0("-log10(P) ",sample_type[list1[-m][2]+1])
    p=ggplot(pVals, mapping=aes_string(x=colnames(pVals)[9*list1[-m][1]+n], y=colnames(pVals)[9*list1[-m][2]+n])) +geom_point(color = 'red')+ theme_classic(base_size = 20) + labs(title=t1,x =x1 , y = y1)
    tax_lab=rownames(pVals)
    #tax_lab[abs(pVals[,9*(m-1)+n])<(-log10(0.05))]=NA
    p1[[c(3:1)[m]]]=p+geom_text_repel(aes(label =tax_lab),size = 3.5)+geom_vline(xintercept=0, linetype="dotted")+geom_hline(yintercept=0, linetype="dotted")
  }
  grid.arrange(p1[[1]],p1[[2]],p1[[3]],ncol=3)
  dev.off()
}
rownames(cor_meta)=sapply(strsplit(colnames(pVals),"_"),"[[",1)[c(4,5,6,7,3)]
colnames(cor_meta)=c(paste("rho",c(1:3)),paste("P",c(1:3)))
write.csv(cor_meta,file="results/gen_tax_onesample_cor.csv")


#1397 wgs samples 
#pathway abundance
path_tab=read.table(file="van_re_humann2_unstratified.tsv",sep="\t",header=T,row.names=1, quote = "")

dim(path_tab) #484 1397
as.character(path_tab[1:5,1:5])
path_tab1=apply(path_tab,1, as.character)
path_tab1=apply(path_tab1,1, as.numeric)
as.character(path_tab1[1:5,1:5])
colnames(path_tab1)=colnames(path_tab)
rownames(path_tab1)=rownames(path_tab)
dim(path_tab1)
path_tab_ag=data.frame(sapply(by(t(path_tab1),map7$sample_type,colMeans),identity))
(path_tab_ag["UNMAPPED",]+path_tab_ag["UNINTEGRATED",])/colSums(path_tab_ag)
#stool      swab    tissue
#UNMAPPED 0.9602758 0.9625294 0.9866912
path_tab_ag["UNMAPPED",]/colSums(path_tab_ag)
#stool      swab    tissue
#UNMAPPED 0.3026347 0.2808287 0.5402413

path_tab1=path_tab1[-grep("UNMAPPED",rownames(path_tab1)),]
path_tab1=path_tab1[-grep("UNINTEGRATED",rownames(path_tab1)),]

path_tab2=t(t(path_tab1)/colSums(path_tab1))*mean(colSums(path_tab1))
colnames(path_tab2)=unlist(lapply(colnames(path_tab),function(x)strsplit(as.character(x),"_",fixed=TRUE)[[1]][1]))

path_tab4=path_tab2[apply(path_tab2!=0,1,sum)/1397>0.01,]
dim(path_tab4) #342 1397
path_tab5=log10(path_tab4+1)
map7=map7[match(colnames(path_tab4),map7$sample_id),]

richness=specnumber(t(path_tab4))

wilcox.test(richness[map7$sample_type!="tissue"]~map7$sample_type[map7$sample_type!="tissue"])
wilcox.test(richness[map7$sample_type!="stool"]~map7$sample_type[map7$sample_type!="stool"])
wilcox.test(richness[map7$sample_type!="swab"]~map7$sample_type[map7$sample_type!="swab"])

pdf(paste0("results/path_alpha_box.pdf"),width=4,height=8)
par(mar=c(10,8,3,3),mfrow=c(1,1),bty="n",mgp = c(4, 1, 0))
type_new=factor(map7$sample_type,levels=c("stool","swab","tissue"))
boxplot(richness~type_new, lwd = 2, ylab = 'Shannon Index',main="",cex.lab=2,cex.main=2,cex.axis=2,border=col11,las=2)
stripchart(richness~type_new, vertical = TRUE,  method = "jitter", add = TRUE, pch = 16, col = col11)
dev.off()

#pcoa with treament and visit plotted. 
adonis(t(path_tab5)~factor(map7$sample_type),permutations=999)
"                            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
factor(map7$sample_type)    2    14.141  7.0706  261.53 0.27284  0.001 ***
Residuals                1394    37.688  0.0270         0.72716           
Total                    1396    51.829                 1.00000               "

#gen_pcoa=metaMDS(t(path_tab5),distance="bray")
gen_pcoa=capscale(t(path_tab5)~1,distance="bray")
mds=list()
mds[[1]]=round((gen_pcoa$CA$eig/sum(gen_pcoa$CA$eig))[1]*100,2)
mds[[2]]=round((gen_pcoa$CA$eig/sum(gen_pcoa$CA$eig))[2]*100,2)
mds[[3]]=round((gen_pcoa$CA$eig/sum(gen_pcoa$CA$eig))[3]*100,2)
mds[[4]]=round((gen_pcoa$CA$eig/sum(gen_pcoa$CA$eig))[4]*100,2)
#pcoa with treament and visit plotted. 
col11=c("red","blue","green","orange","purple","pink","black","lightblue","lightgreen","hotpink","cyan","orchid","tan","grey","gold")
fname=paste("results/path_pcoa12_type.pdf",sep="")
pdf(fname,width=8,height=8)
par(mfrow=c(1,1),mar=c(5,5,5,5))
pcoa12=ordiplot(gen_pcoa,display="sites",choices=c(1,2),type="none",cex.lab=1.5,xlab=paste0("PCoA1 (",mds[[1]],"%)"),ylab=paste0("PCoA2 (",mds[[2]],"%)"))
col3=col11[factor(map7$sample_type)]
points(pcoa12,"sites",col=adjustcolor(col3, alpha.f = 0.2),pch=16,cex=2)
# for (n in 1:length(levels(factor(map7$sample_type)))){
#   ordiellipse(pcoa12, map7$sample_type, kind="se", conf=0.95, lwd=4, draw = "lines", col=col11[n],show.groups=levels(factor(map7$sample_type))[n],label=F,font=2,cex=1.3)
# }
legend("topright",levels(factor(map7$sample_type)), cex=1.8, bty="n", col=col11[1:length(levels(factor(map7$sample_type)))], pch=20)
dev.off()

pdf("results/path_pcoa1234_type_box.pdf",width=8,height=8)
par(mar=c(5,5,5,5),mfrow=c(1,4),bty="n")
boxplot(summary(gen_pcoa)$sites[,1]~map7$sample_type, lwd = 2, ylab = 'distance',main="",cex.lab=1,cex.main=2,cex.axis=1,border=col11,las=2)
stripchart(summary(gen_pcoa)$sites[,1]~map7$sample_type, vertical = TRUE,  method = "jitter", add = TRUE, pch = 16, col = col11)
boxplot(summary(gen_pcoa)$sites[,2]~map7$sample_type, lwd = 2, ylab = 'distance',main="",cex.lab=1,cex.main=2,cex.axis=1,border=col11,las=2,ylim=c(-2,2))
stripchart(summary(gen_pcoa)$sites[,2]~map7$sample_type, vertical = TRUE,  method = "jitter", add = TRUE, pch = 16, col = col11)
boxplot(summary(gen_pcoa)$sites[,1]~map7$sample_type, lwd = 2, ylab = 'distance',main="",cex.lab=1,cex.main=2,cex.axis=1,border=col11,las=2)
stripchart(summary(gen_pcoa)$sites[,1]~map7$sample_type, vertical = TRUE,  method = "jitter", add = TRUE, pch = 16, col = col11)
boxplot(summary(gen_pcoa)$sites[,2]~map7$sample_type, lwd = 2, ylab = 'distance',main="",cex.lab=1,cex.main=2,cex.axis=1,border=col11,las=2,ylim=c(-2,2))
stripchart(summary(gen_pcoa)$sites[,2]~map7$sample_type, vertical = TRUE,  method = "jitter", add = TRUE, pch = 16, col = col11)
dev.off()


#only stool and swab
gen_pcoa=capscale(t(path_tab5[,!map7$sample_type=="tissue"])~1,distance="bray")
map_s=map7[!map7$sample_type=="tissue",]
mds=list()
mds[[1]]=round((gen_pcoa$CA$eig/sum(gen_pcoa$CA$eig))[1]*100,2)
mds[[2]]=round((gen_pcoa$CA$eig/sum(gen_pcoa$CA$eig))[2]*100,2)
mds[[3]]=round((gen_pcoa$CA$eig/sum(gen_pcoa$CA$eig))[3]*100,2)
mds[[4]]=round((gen_pcoa$CA$eig/sum(gen_pcoa$CA$eig))[4]*100,2)
adonis(t(path_tab5[,!map7$sample_type=="tissue"])~factor(map7$sample_type[!map7$sample_type=="tissue"]),permutations=999)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# factor(map7$sample_type[!map7$sample_type == "tissue"])   1    0.5673 0.56729  36.644 0.03795  0.001 ***
#   Residuals                                               929   14.3821 0.01548         0.96205           
# Total                                                   930   14.9494                 1.00000           

col11=c("red","blue","green","orange","purple","pink","black","lightblue","lightgreen","hotpink","cyan","orchid","tan","grey","gold")
fname=paste("results/two_path_pcoa12_type.pdf",sep="")
pdf(fname,width=8,height=8)
par(mfrow=c(1,1),mar=c(5,5,5,5))
pcoa12=ordiplot(gen_pcoa,display="sites",choices=c(1,2),type="none",cex.lab=1.5,xlab=paste0("PCoA1 (",mds[[1]],"%)"),ylab=paste0("PCoA2 (",mds[[2]],"%)"),xlim=c(-0.6,1.6))
col3=col11[factor(map7$sample_type[!map7$sample_type=="tissue"])]
points(pcoa12,"sites",col=adjustcolor(col3, alpha.f = 0.2),pch=16,cex=2)
legend("topright",levels(factor(map7$sample_type[!map7$sample_type=="tissue"])), cex=1.8, bty="n", col=col11[1:length(levels(factor(map7$sample_type[!map7$sample_type=="tissue"])))], pch=20)
dev.off()

pdf("results/two_path_pcoa1234_type_box.pdf",width=8,height=8)
par(mar=c(5,5,5,5),mfrow=c(1,4),bty="n")
boxplot(summary(gen_pcoa)$sites[,1]~droplevels(map_s$sample_type), lwd = 2, ylab = 'distance',main="",cex.lab=1,cex.main=2,cex.axis=1,border=col11,las=2)
stripchart(summary(gen_pcoa)$sites[,1]~droplevels(map_s$sample_type), vertical = TRUE,  method = "jitter", add = TRUE, pch = 16, col = col11)
boxplot(summary(gen_pcoa)$sites[,2]~droplevels(map_s$sample_type), lwd = 2, ylab = 'distance',main="",cex.lab=1,cex.main=2,cex.axis=1,border=col11,las=2,ylim=c(-2,2))
stripchart(summary(gen_pcoa)$sites[,2]~droplevels(map_s$sample_type), vertical = TRUE,  method = "jitter", add = TRUE, pch = 16, col = col11)
boxplot(summary(gen_pcoa)$sites[,1]~droplevels(map_s$sample_type), lwd = 2, ylab = 'distance',main="",cex.lab=1,cex.main=2,cex.axis=1,border=col11,las=2)
stripchart(summary(gen_pcoa)$sites[,1]~droplevels(map_s$sample_type), vertical = TRUE,  method = "jitter", add = TRUE, pch = 16, col = col11)
boxplot(summary(gen_pcoa)$sites[,2]~droplevels(map_s$sample_type), lwd = 2, ylab = 'distance',main="",cex.lab=1,cex.main=2,cex.axis=1,border=col11,las=2,ylim=c(-2,2))
stripchart(summary(gen_pcoa)$sites[,2]~droplevels(map_s$sample_type), vertical = TRUE,  method = "jitter", add = TRUE, pch = 16, col = col11)
dev.off()


library(lme4)
library(lmerTest)
#pairwise two types
path_tab6=path_tab4[apply(path_tab4!=0,1,sum)/1397>0.1,]
dim(path_tab6)#343 1397
tab_n1=path_tab6
map_n1=map7
sample_type=c("stool","swab","tissue")
rownames(tab_n1)=paste0("a",c(1:nrow(path_tab6)))
pVals=matrix(ncol=18,nrow=nrow(tab_n1))
bacteriaMeta1=cbind(t(tab_n1),map_n1)
which(is.na(bacteriaMeta1$visit))#740
bacteriaMeta1=cbind(t(tab_n1),map_n1)[-which(is.na(bacteriaMeta1$visit)),]

pVals=matrix(ncol=30,nrow=nrow(tab_n1))
for (j in 1:3){
  type1=sample_type[j]
  bacteriaMeta=bacteriaMeta1[bacteriaMeta1$sample_type!=type1,]
  bacteriaMeta=bacteriaMeta[!is.na(bacteriaMeta$antibiotics) & !is.na(bacteriaMeta$visit),]
  for( i in 1:dim(tab_n1)[1]){
    print(i)
    if(sd(tab_n1[i,])==0){
      next
    }
    model <- as.formula(paste(rownames(tab_n1)[i],"~","factor(sample_type)+factor(treatment)*factor(visit)+factor(antibiotics)+age+sex+NSAIDS_use+BMI"))
    
    simpleMod <- try(gls(model,method="REML",data=bacteriaMeta))
    mixedMod <- try(lme(model,method="REML",random=~1|study_id,data=bacteriaMeta))
    
    if(class( mixedMod )=="try-error"){
      pVals[i,10*(j-1)+10]  <-NA
    }else{
      pVals[i,10*(j-1)+10]  <- -log10(anova(simpleMod,mixedMod)$"p-value"[2])
    }
    
    for (m in 2:10){
      pVals[i,10*(j-1)+m-1] <- -log10(summary(mixedMod)$tTable[m,5])*sign(summary(mixedMod)$tTable[m,4])
    }
    
  }
}
names1=c(rownames(summary(simpleMod)$tTable)[-1],"subject")
names1[1]="sample_type"
names1[2]="treatment"
names1[3]="visit"
names1[4]="antibiotics"
names1[6]="sex"
names1[9]="interaction"
colnames(pVals)=c(paste(names1,sample_type[1],sep="_"),paste(names1,sample_type[2],sep="_"),paste(names1,sample_type[3],sep="_"))
rownames(pVals)=rownames(path_tab6)
pVals=data.frame(na.omit(pVals))
write.csv(pVals,file="results/path_lmer_P1_pair.csv")
pVals1=read.csv(file="results/path_lmer_P1_pair.csv",header=T,row.names=1)

path_tab7=path_tab4[which(rowMeans(path_tab4)>500),]
gen_type_1=data.frame(sapply(by(t(path_tab7),map7$sample_type,colMeans),identity))
gen_type=t(scale(t(gen_type_1)))
rownames(gen_type)=rownames(path_tab7)
pdf("results/path_heat_full.pdf",width=12,height=15)
pheatmap(gen_type,border_color = "white",fontsize_col = 15)
dev.off()

FDRs=10^-abs(pVals1)
FDRs=apply(FDRs,1,as.character)
FDRs=apply(FDRs,1,as.numeric)
FDRs=matrix(p.adjust(FDRs,method = "fdr"),ncol=30)
colnames(FDRs)=colnames(pVals1)
rownames(FDRs)=rownames(pVals1)
FDRs=data.frame(na.omit(FDRs))


length(which(FDRs[,1]<0.05))#233
length(which(FDRs[,11]<0.05))#222
length(which(FDRs[,21]<0.05))#269
length(which(apply(FDRs[,c(1,11,21)],1,min)<0.05)) #318

rownames(FDRs)[which(FDRs[,1]<0.05)]
rownames(FDRs)[which(FDRs[,11]<0.05)]
rownames(FDRs)[which(FDRs[,21]<0.05)]

setdiff(rownames(FDRs)[which(FDRs[,1]<0.05)],rownames(FDRs)[which(FDRs[,11]<0.05)])#42


#one sample analysis
pVals=matrix(ncol=27,nrow=dim(tab_n1)[1])
for (j in 1:3){
  type1=sample_type[j]
  bacteriaMeta=bacteriaMeta1[bacteriaMeta1$sample_type==type1,]
  bacteriaMeta=bacteriaMeta[!is.na(bacteriaMeta$antibiotics) & !is.na(bacteriaMeta$visit),]
  
  for( i in 1:dim(tab_n1)[1]){
    print(i)
    if(sd(tab_n1[i,])==0){
      next
    }
    model <- as.formula(paste(rownames(tab_n1)[i],"~","factor(treatment)*factor(visit)+factor(antibiotics)+age+sex+NSAIDS_use+BMI"))
    
    simpleMod <- try(gls(model,method="REML",data=bacteriaMeta))
    mixedMod <- try(lme(model,method="REML",random=~1|study_id,data=bacteriaMeta))
    
    if(class(simpleMod)=="try-error"|class(mixedMod)=="try-error"){
      for (m in 2:9){
        pVals[i,9*(j-1)+m-1] =NA
      }
      pVals[i,9*(j-1)+9]=NA
    }else{
      for (m in 2:9){
        pVals[i,9*(j-1)+m-1] = -log10(summary(mixedMod)$tTable[m,5])*sign(summary(mixedMod)$tTable[m,4])
      }
      pVals[i,9*(j-1)+9]  = -log10(anova(simpleMod,mixedMod)$"p-value"[2])
    }
  }
}
pVals=data.frame(pVals)
names1=c(rownames(summary(simpleMod)$tTable)[-1],"subject")
names1[1]="treatment"
names1[2]="visit"
names1[3]="antibiotics"
names1[5]="sex"
names1[8]="interaction"
colnames(pVals)=c(paste(names1,sample_type[1],sep="_"),paste(names1,sample_type[2],sep="_"),paste(names1,sample_type[3],sep="_"))
rownames(pVals)=rownames(path_tab6)
write.csv(pVals,file=paste0("results/path_onesample_lmer.csv"))
pVals=read.csv(file=paste0("results/path_onesample_lmer.csv"),header=T,row.names=1)

cor_meta=matrix(nrow=5,ncol=6)
require("ggrepel")
j=0
for (n in c(4,5,6,7,3)){
  j=j+1
  p1=list()
  fname=paste0("results/path_tax_onesample_",names1[n],"_pair.pdf")
  pdf(fname,width=20,height=8)
  list1=c(0,1,2)
  for (m in 3:1){
    cor=cor.test(pVals[,9*list1[-m][1]+n],pVals[,9*list1[-m][2]+n],method="spearman")$estimate
    pv=as.numeric(cor.test(pVals[,9*list1[-m][1]+n],pVals[,9*list1[-m][2]+n],method="spearman")$p.value)
    if (pv==0){
      pv
    }else if (pv<0.001){
      pv=paste("=",formatC(pv, format = "e", digits = 2))
    }else{
      pv=paste("=",round(pv,3))
    }
    
    cor_meta[j,c(3:1)[m]]=cor
    cor_meta[j,c(3:1)[m]+3]=pv
    
    if (m==3){
      t1=paste(names1[n],":",paste(sample_type[list1[-m]+1],collapse=" vs "),"\n","rho =",round(cor,3),"P",pv)
    }else{
      t1=paste(paste(sample_type[list1[-m]+1],collapse=" vs "),"\n","rho =",round(cor,3),"P",pv)
    }
    
    x1=paste0("-log10(P) ",sample_type[list1[-m][1]+1])
    y1=paste0("-log10(P) ",sample_type[list1[-m][2]+1])
    p=ggplot(pVals, mapping=aes_string(x=colnames(pVals)[9*list1[-m][1]+n], y=colnames(pVals)[9*list1[-m][2]+n])) +geom_point(color = 'red')+ theme_classic(base_size = 20) + labs(title=t1,x =x1 , y = y1)
    #tax_lab=rownames(pVals)
    #tax_lab[abs(pVals[,8*list1[-m][1]+n])<(-log10(0.05))]=""
    p1[[c(3:1)[m]]]=p+geom_vline(xintercept=0, linetype="dotted")+geom_hline(yintercept=0, linetype="dotted")
  }
  grid.arrange(p1[[1]],p1[[2]],p1[[3]],ncol=3)
  dev.off()
}

rownames(cor_meta)=c("age","sex","NSAIDS_use","BMI","antibiotics")
write.csv(cor_meta,file="results/path_cor_meta_onesample.csv")

seq_dep=read.csv(file="seq_depth.csv",row.names = 1)
aggregate(seq_dep$perc~seq_dep[,4],FUN=mean)
"  seq_dep[, 4] seq_dep$perc
1           FT   0.03214532
2           ST   0.98683141
3           SW   0.52139533"
aggregate(seq_dep$perc~seq_dep[,4],FUN=sd)
" seq_dep[, 4] seq_dep$perc
1           FT   0.06565989
2           ST   0.03767827
3           SW   0.31542548"

seq_dep_m=data.frame(rbind(as.matrix(seq_dep[,c(1,3,4)]),as.matrix(seq_dep[,c(1,2,4)])))
seq_dep_m$Process=c(rep("raw",1367),rep("decontam",1367))
seq_dep_m$Process=factor(seq_dep_m$Process,levels=c("raw","decontam"))
colnames(seq_dep_m)=c("Sample","NumReads","SampleType","Process")
seq_dep_m$NumReads_log=log10(as.numeric(as.character(seq_dep_m$NumReads)))
seq_dep_m$SampleType=c("tissue","stool","swab")[factor(seq_dep_m$SampleType)]
seq_dep_m$SampleType=factor(seq_dep_m$SampleType,levels=c("stool","swab","tissue"))


pdf("results/seq_depth.pdf",width=8,height=8)
ggboxplot(seq_dep_m, x = "SampleType", y = "NumReads_log", color = "Process", palette = c("red","blue"), add = "jitter", lwd = 0.8, xlab = "SampleType",ylab = "log10(Number of reads)",
          ggtheme = theme_bw()+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
                                     panel.background = element_blank(), axis.line = element_blank(),axis.text = element_text(size=20),title=element_text(size=5),
                                     axis.title=element_text(size=20),legend.text=element_text(size=20),legend.title=element_text(size=20)))+rotate_x_text(45)
dev.off()




