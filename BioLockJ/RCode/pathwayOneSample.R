#Author: Shan Sun

#BioLockJ configuration: Ali Sorgen
#Date: 12-03-20
#Description: 

#Libraries
library(stringr)
library(gridExtra)
library(nlme)

rm(list = ls())

pipeRoot = dirname(dirname(getwd()))
moduleDir <- dirname(getwd())
input = file.path(pipeRoot,"input/")
output = file.path(moduleDir,"output/")
# pathToMeta<-paste0(pipeRoot,"/",str_subset(dir(pipeRoot), "Metadata"),"/output/")
pathToPvals<-paste0(pipeRoot,"/",str_subset(dir(pipeRoot), "PathwayHeatmap"),"/output/")

pVals1=read.csv(file=paste0(pathToPvals, "path_lmer_P1_pair.csv"),header=T,row.names=1)


#1397 wgs samples 
#pathway abundance
path_tab=read.table(file=paste0(input, "counts/van_re_humann2_unstratified.tsv"),sep="\t",header=T,row.names=1, quote = "")

dim(path_tab) #484 1397
as.character(path_tab[1:5,1:5])
path_tab1=apply(path_tab,1, as.character)
path_tab1=apply(path_tab1,1, as.numeric)
as.character(path_tab1[1:5,1:5])
colnames(path_tab1)=colnames(path_tab)
rownames(path_tab1)=rownames(path_tab)
dim(path_tab1)
# path_tab_ag=data.frame(sapply(by(t(path_tab1),map7$sample_type,colMeans),identity))
# (path_tab_ag["UNMAPPED",]+path_tab_ag["UNINTEGRATED",])/colSums(path_tab_ag)
# #stool      swab    tissue
# #UNMAPPED 0.9602758 0.9625294 0.9866912
# path_tab_ag["UNMAPPED",]/colSums(path_tab_ag)
# #stool      swab    tissue
# #UNMAPPED 0.3026347 0.2808287 0.5402413

path_tab1=path_tab1[-grep("UNMAPPED",rownames(path_tab1)),]
path_tab1=path_tab1[-grep("UNINTEGRATED",rownames(path_tab1)),]

path_tab2=t(t(path_tab1)/colSums(path_tab1))*mean(colSums(path_tab1))
colnames(path_tab2)=unlist(lapply(colnames(path_tab),function(x)strsplit(as.character(x),"_",fixed=TRUE)[[1]][1]))

path_tab4=path_tab2[apply(path_tab2!=0,1,sum)/1397>0.01,]

#pairwise two types
path_tab6=path_tab4[apply(path_tab4!=0,1,sum)/1397>0.1,]
dim(path_tab6)#343 1397
tab_n1=path_tab6
map_n1=read.csv(file=paste0(pathToPvals, "map7.csv"),header=T,row.names=1)
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

write.csv(pVals,file=paste0(output, "path_onesample_lmer.csv"))
pVals=read.csv(file=paste0(output, "path_onesample_lmer.csv"),header=T,row.names=1)

cor_meta=matrix(nrow=5,ncol=6)
require("ggrepel")
j=0
for (n in c(4,5,6,7,3)){
  j=j+1
  p1=list()
  
  ## Figure 6
  fname=paste0(output, "path_tax_onesample_",names1[n],"_pair.pdf")
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
write.csv(cor_meta,file=paste0(output, "/path_cor_meta_onesample.csv"))
