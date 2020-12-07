#Author: Shan Sun

#BioLockJ configuration: Ali Sorgen
#Date: 12-03-20
#Description: 

#Libraries
library(stringr)
library(vegan)

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

richness=specnumber(t(path_tab4))

wilcox.test(richness[map7$sample_type!="tissue"]~map7$sample_type[map7$sample_type!="tissue"])
wilcox.test(richness[map7$sample_type!="stool"]~map7$sample_type[map7$sample_type!="stool"])
wilcox.test(richness[map7$sample_type!="swab"]~map7$sample_type[map7$sample_type!="swab"])

col11=c("red","blue","green","orange","purple","pink","black","lightblue","lightgreen","hotpink","cyan","orchid","tan","grey","gold")

## Figure 3a
pdf(paste0(output, "path_alpha_box.pdf"),width=4,height=8)
par(mar=c(10,8,3,3),mfrow=c(1,1),bty="n",mgp = c(4, 1, 0))
type_new=factor(map7$sample_type,levels=c("stool","swab","tissue"))
boxplot(richness~type_new, lwd = 2, xlab = '', ylab = 'Number of pathways',main="",cex.lab=2,cex.main=2,cex.axis=2, col = "white", border=col11,las=2)
stripchart(richness~type_new, vertical = TRUE,  method = "jitter", add = TRUE, pch = 16, col = col11)
dev.off()

#pcoa with treament and visit plotted. 
adonis <- adonis(t(path_tab5)~factor(map7$sample_type),permutations=999)
# "                            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# factor(map7$sample_type)    2    14.141  7.0706  261.53 0.27284  0.001 ***
# Residuals                1394    37.688  0.0270         0.72716           
# Total                    1396    51.829                 1.00000               "

#gen_pcoa=metaMDS(t(path_tab5),distance="bray")
gen_pcoa=capscale(t(path_tab5)~1,distance="bray")
mds=list()
mds[[1]]=round((gen_pcoa$CA$eig/sum(gen_pcoa$CA$eig))[1]*100,2)
mds[[2]]=round((gen_pcoa$CA$eig/sum(gen_pcoa$CA$eig))[2]*100,2)
mds[[3]]=round((gen_pcoa$CA$eig/sum(gen_pcoa$CA$eig))[3]*100,2)
mds[[4]]=round((gen_pcoa$CA$eig/sum(gen_pcoa$CA$eig))[4]*100,2)
#pcoa with treament and visit plotted. 
col11=c("red","blue","green","orange","purple","pink","black","lightblue","lightgreen","hotpink","cyan","orchid","tan","grey","gold")

## Figure 3b
fname=paste(output, "path_pcoa12_type.pdf",sep="")
pdf(fname,width=8,height=8)
par(mfrow=c(1,1),mar=c(5,5,5,5))
pcoa12=ordiplot(gen_pcoa,display="sites",choices=c(1,2),type="none",cex.lab=1.5,xlab=paste0("PCoA1 (",mds[[1]],"%)"),ylab=paste0("PCoA2 (",mds[[2]],"%)"), 
                main = paste0("PERMANOVA R2=",round(adonis$aov.tab$R2[1], 3), ", P=", adonis$aov.tab$`Pr(>F)`[1]))
col3=col11[factor(map7$sample_type)]
points(pcoa12,"sites",col=adjustcolor(col3, alpha.f = 0.2),pch=16,cex=2)
# for (n in 1:length(levels(factor(map7$sample_type)))){
#   ordiellipse(pcoa12, map7$sample_type, kind="se", conf=0.95, lwd=4, draw = "lines", col=col11[n],show.groups=levels(factor(map7$sample_type))[n],label=F,font=2,cex=1.3)
# }
legend("topright",levels(factor(map7$sample_type)), cex=1.8, bty="n", col=col11[1:length(levels(factor(map7$sample_type)))], pch=20)
dev.off()

## Figure 3b insets
pdf(paste0(output, "path_pcoa1234_type_box.pdf"),width=8,height=8)
par(mar=c(5,5,5,5),mfrow=c(1,4),bty="n")
boxplot(summary(gen_pcoa)$sites[,1]~map7$sample_type, lwd = 2, xlab = '', ylab = 'distance',main="",cex.lab=1,cex.main=2,cex.axis=1, col = "white", border=col11,las=2)
stripchart(summary(gen_pcoa)$sites[,1]~map7$sample_type, vertical = TRUE,  method = "jitter", add = TRUE, pch = 16, col = col11)
boxplot(summary(gen_pcoa)$sites[,2]~map7$sample_type, lwd = 2, xlab = '', ylab = 'distance',main="",cex.lab=1,cex.main=2,cex.axis=1, col = "white", border=col11,las=2,ylim=c(-2,2))
stripchart(summary(gen_pcoa)$sites[,2]~map7$sample_type, vertical = TRUE,  method = "jitter", add = TRUE, pch = 16, col = col11)
boxplot(summary(gen_pcoa)$sites[,1]~map7$sample_type, lwd = 2, xlab = '', ylab = 'distance',main="",cex.lab=1,cex.main=2,cex.axis=1, col = "white", border=col11,las=2)
stripchart(summary(gen_pcoa)$sites[,1]~map7$sample_type, vertical = TRUE,  method = "jitter", add = TRUE, pch = 16, col = col11)
boxplot(summary(gen_pcoa)$sites[,2]~map7$sample_type, lwd = 2, xlab = '', ylab = 'distance',main="",cex.lab=1,cex.main=2,cex.axis=1, col = "white", border=col11,las=2,ylim=c(-2,2))
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
adonis <- adonis(t(path_tab5[,!map7$sample_type=="tissue"])~factor(map7$sample_type[!map7$sample_type=="tissue"]),permutations=999)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# factor(map7$sample_type[!map7$sample_type == "tissue"])   1    0.5673 0.56729  36.644 0.03795  0.001 ***
#   Residuals                                               929   14.3821 0.01548         0.96205           
# Total                                                   930   14.9494                 1.00000           

col11=c("red","blue","green","orange","purple","pink","black","lightblue","lightgreen","hotpink","cyan","orchid","tan","grey","gold")

## Figure 3c
fname=paste(output, "two_path_pcoa12_type.pdf",sep="")
pdf(fname,width=8,height=8)
par(mfrow=c(1,1),mar=c(5,5,5,5))
pcoa12=ordiplot(gen_pcoa,display="sites",choices=c(1,2),type="none",cex.lab=1.5,xlab=paste0("PCoA1 (",mds[[1]],"%)"),ylab=paste0("PCoA2 (",mds[[2]],"%)"),
                xlim=c(-0.6,1.6), main = paste0("PERMANOVA R2=",round(adonis$aov.tab$R2[1], 3), ", P=", adonis$aov.tab$`Pr(>F)`[1]))
col3=col11[factor(map7$sample_type[!map7$sample_type=="tissue"])]
points(pcoa12,"sites",col=adjustcolor(col3, alpha.f = 0.2),pch=16,cex=2)
legend("topright",levels(factor(map7$sample_type[!map7$sample_type=="tissue"])), cex=1.8, bty="n", col=col11[1:length(levels(factor(map7$sample_type[!map7$sample_type=="tissue"])))], pch=20)
dev.off()

## Figure 3c insets
pdf(paste0(output, "two_path_pcoa1234_type_box.pdf"),width=8,height=8)
par(mar=c(5,5,5,5),mfrow=c(1,4),bty="n")

boxplot(summary(gen_pcoa)$sites[,1]~map_s$sample_type, lwd = 2, xlab = '', ylab = 'distance',main="",cex.lab=1,cex.main=2,cex.axis=1, col = "white", border=col11,las=2)
stripchart(summary(gen_pcoa)$sites[,1]~map_s$sample_type, vertical = TRUE,  method = "jitter", add = TRUE, pch = 16, col = col11)

boxplot(summary(gen_pcoa)$sites[,2]~map_s$sample_type, lwd = 2, xlab = '', ylab = 'distance',main="",cex.lab=1,cex.main=2,cex.axis=1, col = "white", border=col11,las=2,ylim=c(-2,2))
stripchart(summary(gen_pcoa)$sites[,2]~map_s$sample_type, vertical = TRUE,  method = "jitter", add = TRUE, pch = 16, col = col11)

boxplot(summary(gen_pcoa)$sites[,1]~map_s$sample_type, lwd = 2, xlab = '', ylab = 'distance',main="",cex.lab=1,cex.main=2,cex.axis=1, col = "white", border=col11,las=2)
stripchart(summary(gen_pcoa)$sites[,1]~map_s$sample_type, vertical = TRUE,  method = "jitter", add = TRUE, pch = 16, col = col11)

boxplot(summary(gen_pcoa)$sites[,2]~map_s$sample_type, lwd = 2, xlab = '', ylab = 'distance',main="",cex.lab=1,cex.main=2,cex.axis=1, col = "white", border=col11,las=2,ylim=c(-2,2))
stripchart(summary(gen_pcoa)$sites[,2]~map_s$sample_type, vertical = TRUE,  method = "jitter", add = TRUE, pch = 16, col = col11)
dev.off()


