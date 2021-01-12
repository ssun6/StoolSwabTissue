#Author: Shan Sun

#BioLockJ configuration: Ali Sorgen
#Date: 12-02-20
#Description: Combine metadata files.

#Libraries

rm(list = ls())

pipeRoot = dirname(dirname(getwd()))
moduleDir <- dirname(getwd())
input = file.path(pipeRoot,"input/")
output = file.path(moduleDir,"output/")

col11=c("red","blue","green","orange","purple","pink","black","lightblue","lightgreen","hotpink","cyan","orchid","tan","grey","gold")

draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...))
  return(res)}

## 'Overwrite' default draw_colnames with your own version 
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))

message(paste0(input, "metadata/study_ID_sample_ID_link.csv"))
map1=read.csv(file=paste0(input, "metadata/study_ID_sample_ID_link.csv"),header=T)
which(map1$visit=="") #726
map1[726,5:7]=rep(NA,3)#change missing values to NA
map1$name_v=paste(map1$study_id,map1$visit)

map2=read.csv(file=paste0(input, "metadata/metadata_20190503.csv"),header=T)
map3=read.csv(file=paste0(input, "metadata/IHC_data_20190503.csv"),header=T)
map3$name_v=paste(map3$studyid,map3$timepoints)

bug_tab=read.table(file=paste0(input, "counts/van_re_metaphlan2.txt"),header=T,row.names=1,sep="\t")
map4=map1[match(colnames(bug_tab),map1$sample_id),]
map5=map2[match(map4$study_id,map2$study_id),]
map6=map3[match(map4$name_v,map3$name_v),]

map7=cbind(map4[,-8],map5[,-1],map6[,-c(1,2,10)])
map7$treat_visit=paste(map7$treatment,map7$visit,sep="_")
map7$treat_visit=factor(map7$treat_visit,levels=levels(factor(map7$treat_visit))[c(2,1,5,4)])


antib=read.csv(file=paste0(input, "metadata/Antibiotic_data_20201022.csv"),header=T)
antib1=data.frame(sapply(by(antib$pq_16,antib$study_id,mean),identity))
colnames(antib1)="antib"
antib1$antibiotics=antib1$antib
antib1$antibiotics[which(antib1$antib<2)]="yes"
antib1$antibiotics[which(antib1$antib==2)]="No"
antib2=antib1[match(map7$study_id,rownames(antib1)),]

days=read.csv(file=paste0(input, "metadata/days.csv"),header=T)
days1=days[match(map7$study_id,days$study_id),]

map7=cbind(map7,antib2,days1)
write.csv(map7,file=paste0(output, "metadata_all.csv"))