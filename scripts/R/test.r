source("/home/jmhan/SDR/scripts/R/SDRlib.r")
library(ggplot2)
library(dbscan)
library(dplyr)
library(data.table)
library(IRanges)
library(plotly)
#args <- commandArgs(trailingOnly = TRUE)
pos <- read.delim("D:/MS/SDR/SDRall/macaque/h1syntenic_blocks.tsv", header = FALSE, col.names = c("ref_chr", "ref_start", "ref_end", "ref_pos", "query_chr", "query_start", "query_end", "query_pos", "orient"))
chrpc<-fread("D:/MS/SDR/SDRall/macaque/p_c_chrlen1.txt")
chrpc<-distinct(chrpc)
colnames(chrpc)<-c("ref_chr","ref_len","query_chr","query_len")
## 染色体
chrnames <- unique(pos$ref_chr)
numeric_part <- as.numeric(gsub("\\D", "", chrnames))
sorted_chrnames <- chrnames[order(numeric_part)] #染色体排序


## -----------对文件进行过滤：删除距离大于某某G的比对,相当于不考虑translocation的情况
# pos$ref_start<-pos$ref_start+1 ##发现起始位点是从0开始的，所以加1
# pos$query_start<-pos$query_start+1

# store$len1<-store$ref_end-store$ref_start
# store$len2<-store$query_end-store$query_start
## 用来存储一些小的translocation


cluster0paras<-700000  ##如果比对比较松散，就需要设置的大一点，尽量设置的大一点，太小会更分散
cluster1paras<-700000
transpara=10 # 如果比对较为松散，即比对较多，需要提高该参数从而最终找到大范围的translocation
for(chrid in sorted_chrnames){
  ## store用来存取所有的结构变异NM_SDR
  ## duplication用来存储所有的Duplication
  ## inversion用来存储所有的inversion
  store<- data.frame(ref_chr=0,ref_start = 0, ref_end = 0,query_chr=0,query_start = 0, query_end = 0,anno='ss')
  storesmall<- data.frame(ref_chr=0,ref_start = 0, ref_end = 0,query_chr=0,query_start = 0, query_end = 0,anno='ss')
  duplication<- data.frame(ref_chr=0,ref_start = 0, ref_end = 0,query_chr=0,query_start = 0, query_end = 0,orient=0,cluster=0)
  inversion<-data.frame(ref_chr=0,ref_start = 0, ref_end = 0,query_chr=0,query_start = 0, query_end = 0)
  smalltrans<- data.frame(ref_chr=0,ref_start = 0, ref_end = 0,query_chr=0,query_start = 0, query_end = 0,orient=0,cluster=0)
  pos.chr<-pos[pos$ref_chr==chrid,]
  pos.chr<-pos.chr[order(pos.chr$ref_start),]
  

  ## 1.第一步，大聚类参数，将大聚类结果sytenic中的translocation提取出来并删掉
  endcluster0<-split_region(pos.chr,0,cluster0paras)
  aftertans<-smalltransf(endcluster0,transpara)
  #dotplot_cluster(endcluster0)
  ###!!!!这里需要加一步提取内部的SDR
  for(i in unique((aftertans$tranbefore)$cluster)){
    a<-aftertans$tranbefore[(aftertans$tranbefore)$cluster==i,]
    a$cluster<-as.character(1:dim(a)[1])
    if(sum(a$orient == "-") > 0.8 * nrow(a)){
      reverse_end<-reverse.region(a,chrid,3)
    }
    if(sum(a$orient == "+") > 0.8 * nrow(a)){
      reverse_end<-reverse.region(a,chrid,2)
    }
    middle<-reverse_end$reverse
    middle<-intersect.unit(middle,unique(middle$query_chr)) #将among中复杂片段作为一个大SDR
    storesmall<-rbind(storesmall,middle)
  }
  if(length(aftertans$tran)!=0){
    smalltrans<-rbind(smalltrans,aftertans$tran)
  }
  
  endcluster0<-aftertans$endcluster0
  #\dotplot_cluster(endcluster1)
  ## 2.第二步，inversion提取（-）
  endcluster1<-split_region(endcluster0,1,cluster1paras)
  inver<-inversion.extract(endcluster1,chrid)  ### 先把inversion提取出来
  ##先把很乱的区域存储为SDR
  list<-which(rle(rle(endcluster1$cluster)$lengths)$values==1 &rle(rle(endcluster1$cluster)$lengths)$lengths >5)
  for(k in list){
    chaoval<-rle(rle(endcluster1$cluster)$lengths)$values
    chaolen<-rle(rle(endcluster1$cluster)$lengths)$lengths
    start<-sum(chaoval[1:k-1])+1
    endcluster1[start:(start+chaolen[k]-1),]$cluster<-"1000"
    print(endcluster1[start:(start+chaolen[k]-1),])
    chaos<-endcluster1[start:(start+chaolen[k]-1),]
    chaos$query_start<-abs(chaos$query_start)
    chaos$query_end<-abs(chaos$query_end)
    chaos<-cluster(chaos)
    chaos$anno<-"BIG"
    store<-rbind(store,chaos[,colnames(store)])
  }
  
  reverse_end<-reverse.region(endcluster1,chrid,0)
  duplic<-reverse_end$dup
  if(!is.character(duplic)){
    duplication<-rbind(duplication,duplic[,colnames(duplication)])
  }
  minimap<-reverse_end$minimaploc  ### !!!!!!!!不知道怎么求
  
  
    for (value in unique(reverse_end$reverse$query_chr)) {
      datarev<-reverse_end$reverse[reverse_end$reverse$query_chr==value,]
      if(nrow(datarev)!=0){
      # 找到query_chr等于当前唯一值的行的索引
      rows_to_remove <- which(datarev$query_chr == value)
      # 如果有多行，删除第一行和最后一行
      if (length(rows_to_remove) > 1) {
        startend<-datarev[c(rows_to_remove[1], tail(rows_to_remove, 1)), ]
        store<-rbind(store,datarev[-c(rows_to_remove[1], tail(rows_to_remove, 1)), ]) ##为了防止有一条query对应多个ref最终导致前后计算的端粒过长，在去除store中含有的比对时误删
        storesmall<-rbind(storesmall,startend)
        }
    }
  }

  ##处理一下大cluster之间的store
  vectore<-store[(store$ref_start-store$ref_end)==2,]
  if(nrow(vectore)!=0){
    vectore$ref_start<-vectore$ref_start-1
    vectore$ref_end<-vectore$ref_end+1
    
  }
  vectore<-store[(store$query_start-store$query_end)==2,]
  if(nrow(vectore)!=0){
    vectore$query_start<-vectore$query_start-1
    vectore$query_end<-vectore$query_end+1
  }
  if(nrow(store[store$ref_start>store$ref_end,])!=0){
    store[store$ref_start>store$ref_end,]$ref_end=store[store$ref_start>store$ref_end,]$ref_start #如果大片段有重叠，直接将其重叠删掉
  }
  if(length(unique(store$query_chr)[-1])!=0){
    for(chr_child in unique(store$query_chr)[-1]){
      new_row<-intersect.unit(store,chr_child) #将among中复杂片段作为一个大SDR
      assign(chr_child,new_row)
    }
    store<-docall(store)
    rm(list=unique(store$query_chr)) ##删除变量
  }
  
  if(dim(inver)[1]!=0){
    inver<-inver[,-1]
    inver<-inver[,-7]
    inversion<-rbind(inversion,inver)
  }
  
  cat("第一次大聚类后（Big cluster among）获得的的SDR数目是",dim(store)[1],"\n")  ##但是这个计算了每条染色体的最开始，如果是多个染色体的比对需要改一下
  cat("获得的的inversion数目是",dim(inversion)[1]-1,"\n") 
  cat("获得的的translocation数目是",dim(smalltrans)[1]-1,"\n") 
  
  ## 现在需要二轮的聚类啦
  
  ## 如果比对中有找到的杂乱区域，就把这些比对去掉
  endcluster1<-reverse_end$endcluster1
  count_minus <- table(endcluster1$cluster[endcluster1$orient == '-'])
  count_total <- table(endcluster1$cluster[endcluster1$cluster%in% as.numeric(names(count_minus))])
  ## 在这里先找一下正负cluster之间的重复，相当于内部重复
  
  
  
  SDRminudlist<-as.numeric(names(count_minus[(count_minus / count_total) >0.6])) ##cluster
  ## 先对大inversion计算一下 SDR
  for(k in SDRminudlist){
    region<-endcluster1[endcluster1$cluster==k,]
    if(dim(region)[1]==1){
      next
    }
    #endcluster2<-split_region(region,2)
    region$cluster<-1:dim(region)[1]
    #dotplot_cluster(endcluster2)
    #endcluster2<-inte.minud(endcluster2,2)## 对其中的inversion进行整合
    claster2<-(repeat.integrate(region,2)) ##将inversion返回去了
    endcluster2<-claster2$afterdup
    if(!is.null(claster2$repeat.region)){
      dup<-distinct(claster2$repeat.region)
      duplication<-rbind(duplication,dup)
    }
    ## 从一堆-的中寻找为+的各自进行聚类
    endcluster2before<-endcluster2
    orientid="+"
    endcluster2<-smallcluster(endcluster2,orientid) ##得到把负链中的相邻正链进行整合后计算
    reverse_end<-reverse.region(endcluster2,chrid,3)
    middle<-reverse_end$reverse
    if(nrow(middle)!=0){
      middle[middle$ref_start>middle$ref_end,]$ref_end=middle[middle$ref_start>middle$ref_end,]$ref_start #如果大片段有重叠，直接将其重叠删掉
      middle<-intersect.unit(middle,unique(middle$query_chr)) #将among中复杂片段作为一个大SDR
      middle$orient<-"s"
      middle$cluster<-"1000000"
      storesmall<-rbind(storesmall,middle[,colnames(storesmall)]) ## 存到smalltrans里，防止存在store中会将synetenic中的区域删掉
    }
    storesmall<-insertsmall(endcluster2before,storesmall,orientid)
    }
    

  
  
  for(chr_child in unique(endcluster1$query_chr)){
    new_row<-alignintersect(endcluster1[endcluster1$query_chr==chr_child,],store[store$query_chr==chr_child,],store[store$query_chr==chr_child,]) 
    ## 这里重新写一下，如果交集大于70%则去掉
    assign(chr_child,new_row)
  }
  endcluster1<-docall(endcluster1)
  rm(list=unique(endcluster1$query_chr)) ##删除变量
  
  for(k in setdiff(unique(endcluster1$cluster),SDRminudlist)){
    region<-endcluster1[endcluster1$cluster==k,]
    region$cluster<-1:dim(region)[1]
    if(dim(region)[1]==1){
      next
    }
    #endcluster2<-split_region(region,2)
    claster2<-(repeat.integrate(region,1))
    endcluster2<-claster2$afterdup
    if(!is.null(claster2$repeat.region)){
      dup<-distinct(claster2$repeat.region)
      duplication<-rbind(duplication,dup)
    }
    endcluster2before<-endcluster2
    orientid="-"
    endcluster2<-smallcluster(endcluster2,orientid) ##得到把负链中的相邻正链进行整合后计算
    reverse_end<-reverse.region(endcluster2,chrid,2)
    middle<-reverse_end$reverse
    middle[middle$ref_start>middle$ref_end,]$ref_end=middle[middle$ref_start>middle$ref_end,]$ref_start #如果大片段有重叠，直接将其重叠删掉
    middle<-intersect.unit(middle,unique(middle$query_chr)) #将among中复杂片段作为一个大SDR
    store<-rbind(store,middle)
    storesmall<-insertsmall(endcluster2before,storesmall,orientid)

  }
  cat("第二次小聚类后一起获得的的SDR数目是",dim(store)[1],"\n")
  cat("第二次小聚类后一起获得的的INV数目是",dim(inversion)[1],"\n")
  cat("第二次小聚类后一起获得的的DUP数目是",dim(duplication)[1],"\n")
  cat("第二次小聚类后一起获得的TRAN数目是",dim(smalltrans)[1],"\n")
  inversion$anno<-"INV"
  duplication$anno<-"DUP"
  smalltrans$anno<-"TRANS"
  storesmall$anno<-"SDR_NM"
  duplication<-duplication[,colnames(store)]
  smalltrans<-smalltrans[-1,colnames(store)]
  storesmall<-storesmall[-1,colnames(store)]
  all<-rbind(store,inversion[-1,],duplication[-1,],smalltrans,storesmall)
  all<-all[order(all$ref_start),]
  for (chr_child in unique(all$query_chr)){
  assign(chr_child,endfilter(all[all$query_chr==chr_child,],chrid,chr_child))
  
  }
  data<-docall(all)
  write.table(data, paste("D:/MS/saffire测试/R/result/",chrid,"SDRend.tsv",sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)
  write.table(minimap, paste("D:/MS/saffire测试/R/minimap/",chrid,"minimap.tsv",sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)
  }

write.csv(result_data[result_data$SV=="DUP",], paste("D:/MS/saffire测试/R/result/SDRend.csv",sep = ""), quote = FALSE, row.names = FALSE)