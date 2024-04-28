######### telomere calculate
library(data.table)
library(dplyr)
args <- commandArgs(trailingOnly = TRUE)
hmtelo<-fread(args[2])
hmtelo$V2<-hmtelo$V2+1
chrpc<-fread(args[3])
if(args[1]=="REF"){
chr<-chrpc[,c(1,2)]
chr<- distinct(chr)
}else{
chr<-chrpc[,c(3,4)]
chr<- distinct(chr) 
}


colnames(chr)=c("chr","chr_len")
colnames(hmtelo)[1]="chr"

hmall<-merge(hmtelo,chr,by="chr")
hmall$endtelo<-hmall$chr_len-hmall$V2
hmend<-hmall[hmall$V3<100000 |hmall$endtelo<100000,]
hmend$id='zero'
hmend[hmend$V3<100000,]$id="first"
hmend[hmend$endtelo<100000,]$id="end"


result_hmtelo<-hmend %>% group_by(chr,id)%>%
  summarise(loc=ifelse(
    id=="first",max(V3),min(V2)
  ),loc2=ifelse(
    id=="first",1,chr_len
  ),
  ,telolen=ifelse(
    id=="first",max(V3),max(endtelo)
  ))
result_hmtelo<- distinct(result_hmtelo)


result_hmtelo <- within(result_hmtelo, {
  temp <- loc
  loc[loc2 == 1] <- loc2[loc2 == 1]
  loc2[loc2 == 1] <- temp[loc2 == 1]
  rm(temp)
})

colnames(result_hmtelo)<-c("chr","id","ref_start","ref_end","telolen")
write.table(result_hmtelo, args[4], quote = FALSE, sep = "\t", row.names = FALSE)
