# 
library(dbscan)
library(dplyr)
library(data.table)
library(IRanges)
args <- commandArgs(trailingOnly = TRUE)
print(args)
alignintersect<-function(align1,align2,align3){
  del.list<-c()
  cat("-------开始去除着丝粒和端粒区域的比对---------\n")
  cat("去除前的：",dim(align1)[1],"\n")
  ir1ref <- IRanges(start = align1$V8, end = align1$V9)
  ir1que  <- IRanges(start = align1$V3, end = align1$V4)
  ir2 <- IRanges(start =align2$ref_start , end =align2$ref_end)
  overlapsref<-findOverlaps(ir1ref,ir2) ##我们的比对序列的参考基因组序列和人的端粒|着丝粒的交集
  if(length(overlapsref@from)!=0){  
    del.list<-append(del.list,overlapsref@from)
  }
  if(!missing(align3)){
    ir3 <- IRanges(start = align3$query_start, end = align3$query_end)
    overlapsque<-findOverlaps(ir1que,ir3)
    if(length(overlapsque@from)!=0){
      del.list<-append(del.list,overlapsque@from)
    }
  }
  if(length(del.list)!=0)
  {align1<-align1[-del.list,]}
  cat("去除后的：",dim(align1)[1],"\n")
  return(align1)
}
#-------------------function2 把pos改为saffire文件函数 -----------------------------------------------------
flit2saffire<-function(pos,name){
  cat("-------开始转saffire--------\n")
  col_names <- colnames(pos)
  col_names[c(1,3, 4, 8, 9,6,5)] <- c("query_chr","query_start", "query_end", "ref_start", "ref_end","ref_chr","orient")
  colnames(pos) <- col_names
  chrpc<-fread(args[1])
  chrpc<-distinct(chrpc)
  colnames(chrpc)<-c("ref_chr","ref_len","query_chr","query_len")
  saffire<-merge(pos,chrpc,by=c("ref_chr","query_chr"))
  saffire<-saffire[,c("ref_chr","ref_start","ref_end","ref_len","orient","query_chr","query_start","query_end","query_len")]
  saffire_name<-c("#reference_name","reference_start","reference_end","reference_length","strand",  "query_name","query_start","query_end","query_length","perID_by_matches","perID_by_events", "perID_by_all","matches","mismatches","deletion_events", "insertion_events", "deletions","insertions")
  saffire[,c("perID_by_matches","perID_by_events", "perID_by_all","matches","mismatches","deletion_events", "insertion_events", "deletions","insertions")]<-0
  colnames(saffire)<-saffire_name
  saffire$matches<-100
  saffire$mismatches<-100
  write.table(saffire, file = paste(args[2],name,sep = ""), sep = "\t", row.names = FALSE,quote = FALSE)
}
#-------------------function4 split 并聚类 ----------------------



## 分为三种情况：
## 什么都没有
## 只有参考基因组的端粒和着丝粒
## 两者都有
## cta=TRUE 都计算 FALSE都不计算
## cts=TRUE 计算参考基因组的
## ctn=TRUE 计算参考基因组的
cts=TRUE
if(cts){  ##只计算参考基因组的
  result_reftelo<-fread(args[5])
  result_quetelo<-""
  result_refcentr<-fread(args[6])
  result_quecentr<-""
}


pos<-fread(args[3],fill=TRUE )
cat("去除前的：",dim(pos)[1],"\n")
flit2saffire(pos,"init.saffire") ## 看看过滤了以后的比对结果
chrnames <- unique(pos$V6)
numeric_part <- as.numeric(gsub("\\D", "", chrnames))
sorted_chrnames <- chrnames[order(numeric_part)] #染色体排序
## 染色体
## -----------对文件进行过滤：删除距离大于某某G的比对,相当于不考虑translocation的情况
##把着丝粒和端粒区域的比对删掉
for(chrid in sorted_chrnames){
  pos.chr<-pos[pos$V6==chrid,]
  print(chrid)
  align1<-pos.chr  ##我们的比对

  if(cts){
    align2<-result_refcentr[result_refcentr$chr==chrid,]          ##人着丝粒
    intersect<-which(pos.chr$V8>=align2$ref_start & pos.chr$V9<=align2$ref_end)
    if(length(intersect)!=0){
      print("have")
      pos.chr<-pos.chr[-intersect,]
    }
    align2<-result_reftelo[result_reftelo$chr==chrid,]          ##人端粒
    for(j in dim(align2)[1]){
      intersect<-which(pos.chr$V8>=align2[j,]$ref_start & pos.chr$V9<=align2[j,]$ref_end)
      if(length(intersect)!=0){
      pos.chr<-pos.chr[-intersect,]
    }
    }
    assign(chrid,pos.chr)
  }  
}


pos<-do.call(rbind,mget(chrnames))
pos$sourse<-1:dim(pos)[1]
rm(list=chrnames) ##删除变量
cat("去除着丝粒后的：",dim(pos)[1],"\n")
flit2saffire(pos,"filt_process.saffire") ## 看看过滤了以后的比对结果
all<-data.frame()
for(i in unique(pos$V1)){
  for(j in unique(pos$V6)){
     data<-pos[(pos$V1==i) & (pos$V6==j),]
    if(nrow(data)==0){
      next
    }else{
x<-c(i,j,as.numeric(sum(data[data$V5 == "-", "V4"] - data[data$V5 == "-", "V3"])),as.numeric(sum(data[data$V5 == "+", "V4"] - data[data$V5 == "+", "V3"])))
all<-rbind(all,x)
    }
  }
}

# 将字符串列转换为数字列
colnames(all)<-c("V1","V2","V3","V4")
all$V5 <- as.numeric(all$V3)
all$V6 <- as.numeric(all$V4)

print("a")
initdata<-all
data<-all
## 筛选较大的错误组装----------------------------------------------
data$len<-data$V5+data$V6
grouped_vectors <- tapply(data$len, data$V1, FUN = function(x) x)
all=list()
for(m in names(grouped_vectors)){
  if(length(grouped_vectors[[m]])!=1){
    all[[m]]<-grouped_vectors[[m]]
  }
}
all1=list()
for(m in names(all)){
  if(length(which(all[[m]]>max(all[[m]])/5))!=1){
    all1[[m]]<-all[[m]]
  }
}
all2=list()
for(m in names(all1)){
  if(max(all1[[m]])>1000000){
    all2[[m]]<-all1[[m]]
  }
}

saf<-fread(paste(args[2],"filt_process.saffire",sep = ""))
allista<-data.frame(matrix(ncol = 5, nrow = 0))
for(j in names(all2)){
  a<-data[data$V1==j,]
  allist<-data.frame(matrix(ncol = 5, nrow = 0))
  for(m in a$V2){
    x<-saf[saf$query_name==j & saf$`#reference_name`==m,]
    count_minus <- sum(x$strand == "-")
    count_plus <- sum(x$strand == "+")
    if(count_minus>count_plus){
      colnames(allist)<-c("contig","chr","start","end","strand")
      allist<-rbind(allist,c(j,m,min(x$query_start),max(x$query_end),"-"))
    }
    if(count_plus>=count_minus){
      colnames(allist)<-c("contig","chr","start","end","strand")
      allist<-rbind(allist,c(j,m,min(x$query_start),max(x$query_end),"+"))
    }
    
  }
  breakpos<-c()
  allist<-allist[order(allist$start),]
  apenddat<-a[,c(1,2,7)]
  colnames(apenddat)<-c("contig","chr","len")
  allist<-merge(allist,apenddat)
  #allist$len<-as.numeric(allist$end)-as.numeric(allist$start)
  allist<-allist[allist$len>=1000000,]
  if(nrow(allist)>=2){
    for(line in 2:dim(allist)[1]){
      if(allist[line,]$start>allist[line-1,]$end){
        if(line==2){
          allist[1,]$start=0
        }
        if(line==dim(allist)[1]){
          allist[line,]$end=unique(saf[saf$query_name==unique(allist$contig),]$query_length)
        }
        allist[line,]$start<-allist[line-1,]$end
      }
    }
    allista<-rbind(allista,allist)
  }

  
}

write.csv(allista, paste(dirname(args[4]),"/add.csv",sep=""), quote = FALSE, row.names = FALSE)


## 重新定义
all<-data.frame()
for(i in unique(pos$V1)){
     data<-pos[(pos$V1==i),]
    if(nrow(data)==0){
      next
    }else{
x<-c(i,as.numeric(sum(data[data$V5 == "-", "V4"] - data[data$V5 == "-", "V3"])),as.numeric(sum(data[data$V5 == "+", "V4"] - data[data$V5 == "+", "V3"])))
all<-rbind(all,x)
    }
}
print(all)
# 将字符串列转换为数字列
colnames(all)<-c("V1","V2","V3")
all$V5 <- as.numeric(all$V2)
all$V6 <- as.numeric(all$V3)
print(all)


write.table(pos, file = args[4], sep = "\t", row.names = FALSE, col.names = FALSE,quote=FALSE)
write.table(setdiff(all[all$V5 > all$V6,]$V1,allista$contig), file = paste(dirname(args[4]),"/minus.txt",sep=""), sep = "\t", row.names = FALSE, col.names = FALSE,quote=FALSE)
write.table(setdiff(all[all$V5 <= all$V6,]$V1,allista$contig), file = paste(dirname(args[4]),"/plus.txt",sep=""), sep = "\t", row.names = FALSE, col.names = FALSE,quote=FALSE)

print("complete!!!!")
