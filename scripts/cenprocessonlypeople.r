######### telomere calculate
library(data.table)
library(dplyr)
chrpc<-fread("D:/MS/SDR/centromere/hg19chrpc.txt")
if(args[1]=="REF"){
  chr<-chrpc[,c(1,2)]
  chr<- distinct(chr)
}else{
  chr<-chrpc[,c(3,4)]
  chr<- distinct(chr) 
}

colnames(chr)=c("chr","chr_len")

######### centromere calculate
pancen<-fread("D:/MS/SDR/centromere/hg19centr.txt", fill=TRUE)

colnames(pancen)[5]='chr'
panall<-merge(pancen,chr,by="chr")
panall<-panall[,c(1,6,7,10,11,17)]
panall<-panall[panall$V10=="ALR/Alpha",]
for (i in 1: dim(panall)[1]){
  panall$len[i]<-panall[i+1,2]-panall[i,3]
}

panall$name<-0
panall[panall$len<500000,]$name<-1 ##moify 1M
print(panall)
end<-as.data.frame(matrix(nrow=length(unique(panall$chr)),ncol=4))
a<-1
for (m in unique(panall$chr)){
  #现在更改为找HSATII的区域和satellite的区域  
  xxx=panall[panall$chr==m,]
  sa_cen<-xxx[xxx$V11=="Satellite/centr"]
  print(sa_cen)
  if(nrow(sa_cen)!=0){
    x1=rle(sa_cen$name)
    #loc1<-which(x1$lengths==max(x1$lengths[which(x1$values==1)]))
    loc1<-which(x1$values==1)
    alllen<-0
    al<-0
    for(k in loc1){
      test1<-sa_cen[c((sum(x1$lengths[1:k-1])+1):sum(x1$lengths[1:k])),]
      if(nrow(test1)>=2){
        if(test1[nrow(test1),]$V7-test1[1,]$V6 >alllen){
          alllen<-test1[nrow(test1),]$V7-test1[1,]$V6 
          al<-k
        }
      }
    }
    test1<-sa_cen[c((sum(x1$lengths[1:al-1])+1):sum(x1$lengths[1:al])),]
    print(test1)
    end$V1[a]<-m
    end$V2[a]<-min(test1$V6)
    end$V3[a]<-max(test1$V7)
    end$V4[a]<-test1$chr_len[1]
    a<-a+1
  }
}
colnames(end)<-c("chr","query_start","query_end","chr_len")
write.table(end, "D:/MS/SDR/centromere/hg19cenend.txt", quote = FALSE, sep = "\t", row.names = FALSE)
