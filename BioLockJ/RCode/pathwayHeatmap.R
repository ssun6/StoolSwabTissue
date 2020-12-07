#Author: Shan Sun

#BioLockJ configuration: Ali Sorgen
#Date: 12-03-20
#Description: 

#Libraries
library(stringr)
library(pheatmap)
library(nlme)

rm(list = ls())

pipeRoot = dirname(dirname(getwd()))
moduleDir <- dirname(getwd())
input = file.path(pipeRoot,"input/")
output = file.path(moduleDir,"output/")
pathToMeta<-paste0(pipeRoot,"/",str_subset(dir(pipeRoot), "Metadata"),"/output/")


bug_tab=read.table(file=paste0(input, "/counts/van_re_metaphlan2.txt"),header=T,row.names=1,sep="\t")
tax_level=sapply(rownames(bug_tab),function(i){length(strsplit(i,"\\|")[[1]])})
map7=read.table(file=paste0(pathToMeta, "metadata_all.csv"),header=T,sep=",")
map7=map7[match(colnames(bug_tab),map7$sample_id),]

#pairwise two types
n=6
bug_tab2=bug_tab[tax_level==n,]
map7=map7[match(colnames(bug_tab2),map7$sample_id),]

#pathway abundance
path_tab=read.table(file=paste0(input, "counts/van_re_humann2_unstratified.tsv"),sep="\t",header=T,row.names=1, quote = "")

dim(path_tab) #484 1397
as.character(path_tab[1:5,1:5])
path_tab1=apply(path_tab,1, as.character)
path_tab1=apply(path_tab1,1, as.numeric)
as.character(path_tab1[1:5,1:5])
colnames(path_tab1)=colnames(path_tab)
rownames(path_tab1)=rownames(path_tab)
dim(path_tab1) #484 1397

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
write.csv(map7,file=paste0(output, "map7.csv"))

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

write.csv(pVals,file=paste0(output, "path_lmer_P1_pair.csv"))
pVals1=read.csv(file=paste0(output, "path_lmer_P1_pair.csv"),header=T,row.names=1)

path_tab7=path_tab4[which(rowMeans(path_tab4)>500),]
gen_type_1=data.frame(sapply(by(t(path_tab7),map7$sample_type,colMeans),identity))
gen_type=t(scale(t(gen_type_1)))
rownames(gen_type)=rownames(path_tab7)


## Figure 4
pdf(paste0(output, "path_heat_full.pdf"),width=12,height=15)
pheatmap(gen_type,border_color = "white",fontsize_col = 15, angle_col = 45)
dev.off()
