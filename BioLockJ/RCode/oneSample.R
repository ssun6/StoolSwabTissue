#Author: Shan Sun

#BioLockJ configuration: Ali Sorgen
#Date: 12-03-20
#Description: Generate taxa scatter plots based on various factors.

#Libraries
library(stringr)
library(nlme)
library(gridExtra)

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


#one sample
pVals=matrix(ncol=27,nrow=dim(tab_n1)[1])
for (j in 1:3){
  type1=sample_type[j]
  bacteriaMeta=bacteriaMeta1[bacteriaMeta1$sample_type==type1,]
  
  bacteriaMeta$sample_type <- factor(bacteriaMeta$sample_type)
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

write.csv(pVals,file=paste0(output, "onesample_l6_tax_lmer.csv"))

pVals=read.csv(file=paste0(output, "onesample_l6_tax_lmer.csv"),row.names=1)
cor_meta=matrix(nrow=5,ncol=6)
j=0
require("ggrepel")
for (n in c(4,5,6,7,3)){
  j=j+1
  p1=list()
  
  ## Figure 5
  fname=paste0(output, "gen_tax_onesample_",names1[n],"_pair.pdf")
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
write.csv(cor_meta,file=paste0(output, "gen_tax_onesample_cor.csv"))
