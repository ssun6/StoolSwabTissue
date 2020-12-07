#Author: Shan Sun

#BioLockJ configuration: Ali Sorgen
#Date: 12-03-20
#Description: 

#Libraries
library(ggplot2)
library(ggpubr)

rm(list = ls())

pipeRoot = dirname(dirname(getwd()))
moduleDir <- dirname(getwd())
input = file.path(pipeRoot,"input/")
output = file.path(moduleDir,"output/")

seq_dep=read.csv(file=paste0(input, "seq_depth.csv"),row.names = 1)
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


pdf(paste0(output, "seq_depth.pdf"),width=8,height=8)
ggboxplot(seq_dep_m, x = "SampleType", y = "NumReads_log", color = "Process", palette = c("red","blue"), add = "jitter", lwd = 0.8, xlab = "SampleType",ylab = "log10(Number of reads)",
          ggtheme = theme_bw()+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
                                     panel.background = element_blank(), axis.line = element_blank(),axis.text = element_text(size=20),title=element_text(size=5),
                                     axis.title=element_text(size=20),legend.text=element_text(size=20),legend.title=element_text(size=20)))+rotate_x_text(45)
dev.off()




