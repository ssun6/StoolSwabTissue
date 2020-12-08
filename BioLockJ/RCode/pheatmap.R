#Author: Shan Sun

#BioLockJ configuration: Ali Sorgen
#Date: 12-03-20
#Description: Generate heat map.

#Libraries
library(stringr)
library(nlme)
library(pheatmap)

rm(list = ls())

pipeRoot = dirname(dirname(getwd()))
moduleDir <- dirname(getwd())
input = file.path(pipeRoot,"input/")
output = file.path(moduleDir,"output/")
pathToMeta<-paste0(pipeRoot,"/",str_subset(dir(pipeRoot), "Metadata"),"/output/")

# draw_colnames_45 <- function (coln, gaps, ...) {
#   coord = pheatmap:::find_coordinates(length(coln), gaps)
#   x = coord$coord - 0.5 * coord$size
#   res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...))
#   return(res)}
# 
# ## 'Overwrite' default draw_colnames with your own version
# assignInNamespace(x="draw_colnames", value="draw_colnames_45",
#                   ns=asNamespace("pheatmap"))

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
write.csv(bacteriaMeta1,file=paste0(output, "bacteriaMeta1.csv"))

pVals=matrix(ncol=30,nrow=nrow(tab_n1))
for (j in 1:3){
  type1=sample_type[j]
  bacteriaMeta=bacteriaMeta1[bacteriaMeta1$sample_type!=type1,]
  
  bacteriaMeta$sample_type <- factor(bacteriaMeta$sample_type)
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

write.csv(pVals,file=paste0(output, "/l",n,"tax_lmer_P1_pair.csv"))

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

## Figure 2
pdf(paste0(output, "/gen_heat.pdf"),width=8,height=8)
pheatmap(gen_type,border_color = "white",fontsize_col =15, angle_col = 45)
dev.off()

write.csv(tab_n1,file=paste0(output, "tab_n1.csv"))
write.csv(bacteriaMeta,file=paste0(output, "bacteriaMeta.csv"))
