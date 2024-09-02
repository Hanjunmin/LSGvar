library(dbscan)
library(dplyr)
library(data.table)
library(IRanges)
args <- commandArgs(trailingOnly = TRUE)
source(args[1])
# ## 端粒
# result_hmtelo<-fread("D:/MS/SDR/telomere/hm_teloend.tsv")
# result_chantelo<-fread("D:/MS/SDR/telomere/chan_teloend.tsv")
#
# ## 着丝粒
# result_hmcentr<-fread("D:/MS/SDR/centromere/hm_centroend.tsv")
# result_chancentr<-fread("D:/MS/SDR/centromere/chan_centroend.tsv")
#transpara=10
#transpara=5
pos <- read.delim(args[2], header = FALSE, col.names = c("ref_chr", "ref_start", "ref_end", "ref_pos", "query_chr", "query_start", "query_end", "query_pos", "orient"))
chrpc<-fread(args[3])
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


INVball<-data.frame()
cluster0paras<-700000  ##如果比对比较松散，就需要设置的大一点，尽量设置的大一点，太小会更分散
cluster1paras<-700000
#transpara=10 # 如果比对较为松散，即比对较多，需要提高该参数从而最终找到大范围的translocation
#sorted_chrnames<- sorted_chrnames[14:23]
for(chrid in sorted_chrnames){
  ## store用来存取所有的结构变异NM_SDR
  ## duplication用来存储所有的Duplication
  ## inversion用来存储所有的inversion
  store<- data.frame(ref_chr=0,ref_start = 0, ref_end = 0,query_chr=0,query_start = 0, query_end = 0,anno='ss')
  storesmall<- data.frame(ref_chr=0,ref_start = 0, ref_end = 0,query_chr=0,query_start = 0, query_end = 0,anno='ss',orient="+")
  duplication<- data.frame(ref_chr=0,ref_start = 0, ref_end = 0,query_chr=0,query_start = 0, query_end = 0,orient=0,cluster="0")
  inversion<-data.frame(ref_chr=0,ref_start = 0, ref_end = 0,query_chr=0,query_start = 0, query_end = 0)
  smalltrans<- data.frame(ref_chr=0,ref_start = 0, ref_end = 0,query_chr=0,query_start = 0, query_end = 0,orient=0,cluster="0")
  pos.chr<-pos[pos$ref_chr==chrid,]
  pos.chr<-pos.chr[order(pos.chr$ref_start),]
  highdupall<-data.frame(ref_chr=0,ref_start = 0, ref_end = 0)
  storehighdup_sdr<- data.frame(ref_chr=0,ref_start = 0, ref_end = 0,query_chr=0,query_start = 0, query_end = 0,anno='ss')
  
  ## 想计算一个高度重复区域计算
  ## 
  list<-c()
  dupset<-highdupregion(pos.chr)
  if(length(dupset)!=1 ||length(dupset[[1]]>=5)){
    for (dupsetname in names(dupset)){
      print(length(dupset[[dupsetname]]))
      if(length(dupset[[dupsetname]])>=5){
        vectore<-append(dupset[[dupsetname]],as.numeric(dupsetname))
        higndup<-pos.chr[vectore,]
        list<-append(list,vectore)
        higndup$cluster=1
        duplication<-rbind(duplication,higndup[,colnames(duplication)])
        highdupall<-rbind(highdupall,c(chrid,min(higndup$ref_start),max(higndup$ref_end)))
        ## 防止误删，找到重复的边界看看周围有没有insertion和deletion
      
        dup_boundary<-pos.chr[min(vectore),]
        dup_boundary$anno="SDR_NM"
        if(min(vectore)!=1){
          dup_boundary$ref_end<-dup_boundary$ref_start
          dup_boundary$ref_start<-pos.chr[min(vectore)-1,]$ref_end
          if(pos.chr[min(vectore)-1,"orient"]=="-" & pos.chr[min(vectore),"orient"]=="-"){
            dup_boundary$query_start<-dup_boundary$query_end
            lis<-pos.chr[pos.chr$query_start>=dup_boundary$query_end,]$query_start
            dup_boundary$query_end<-lis[which.min(abs(lis-dup_boundary$query_end))]
            storehighdup_sdr<-rbind(storehighdup_sdr,dup_boundary[,colnames(storehighdup_sdr)])
            
            
          }else{
            dup_boundary$query_end<-dup_boundary$query_start
            if(nrow(pos.chr[pos.chr$query_end<=dup_boundary$query_start,])==0){
              next
            }
            lis<-pos.chr[pos.chr$query_end<=dup_boundary$query_start,]$query_end
            dup_boundary$query_start<-lis[which.min(abs(lis-dup_boundary$query_start))]
            storehighdup_sdr<-rbind(storehighdup_sdr,dup_boundary[,colnames(storehighdup_sdr)])
            ##
            dup_boundary$query_start<-pos.chr[min(vectore)-1,]$query_end
            lis<-pos.chr[pos.chr$query_start>=dup_boundary$query_start,]$query_start
            dup_boundary$query_end<-lis[which.min(abs(lis-dup_boundary$query_start))]
            storehighdup_sdr<-rbind(storehighdup_sdr,dup_boundary[,colnames(storehighdup_sdr)])
            
          }
          
        }
          
          dup_boundary<-pos.chr[max(vectore),]
          dup_boundary$anno="SDR_NM"
        if(max(vectore)!=dim(pos.chr)[1]){
          dup_boundary$ref_start<-dup_boundary$ref_end
          dup_boundary$ref_end<-pos.chr[max(vectore)+1,]$ref_start
          
          if(pos.chr[max(vectore)+1,"orient"]=="-" & pos.chr[max(vectore),"orient"]=="-"){
            dup_boundary$query_end<-dup_boundary$query_start
            lis<-pos.chr[pos.chr$query_end<=dup_boundary$query_start,]$query_end
            dup_boundary$query_end<-lis[which.min(abs(lis-dup_boundary$query_start))]
            storehighdup_sdr<-rbind(storehighdup_sdr,dup_boundary[,colnames(storehighdup_sdr)])
            

          }else{
            dup_boundary$query_start<-dup_boundary$query_end
            lis<-pos.chr[pos.chr$query_start>=dup_boundary$query_end,]$query_start
            dup_boundary$query_end<-lis[which.min(abs(lis-dup_boundary$query_end))]
            storehighdup_sdr<-rbind(storehighdup_sdr,dup_boundary[,colnames(storehighdup_sdr)])
            
            ##
            dup_boundary$query_end<-pos.chr[max(vectore)+1,]$query_start
            lis<-pos.chr[pos.chr$query_end<=dup_boundary$query_end,]$query_end
            dup_boundary$query_start<-lis[which.min(abs(lis-dup_boundary$query_end))]
            storehighdup_sdr<-rbind(storehighdup_sdr,dup_boundary[,colnames(storehighdup_sdr)])
            
          }
          
          
        }
        
        }
    }
    if(length(list)!=0){
      pos.chr<-pos.chr[-list,]
    }
  }
  storehighdup_sdr<-distinct(storehighdup_sdr)
 
  ## 1.第一步，大聚类参数，将大聚类结果sytenic中的translocation提取出来并删掉
  endcluster0<-split_region(pos.chr,0,cluster0paras)
  # Dupex<-clusterall(endcluster0)
  # Dupex<-Dupex[order(Dupex$ref_start),]
  # if(nrow(Dupex)!=1){
  #   for(i in 2:dim(Dupex)[1]){
  #     if(Dupex[i,]$ref_start<Dupex[i-1,]$ref_end & Dupex[i,]$ref_start>Dupex[i-1,]$ref_start){
  #       a<-Dupex[i,]$ref_start
  #       b<-min(Dupex[c(i,i-1),]$ref_end)
  #       if(nrow(endcluster0[endcluster0$ref_start>=a & endcluster0$ref_end<=b,])!=0)
  #       {
  #         c<-clusterall(endcluster0[endcluster0$ref_start>=a & endcluster0$ref_end<=b,])
  #         duplication<-rbind(duplication,c)
  #       }
  #       else{
  #         if((b-a)>10000){
  #           c<-Dupex[i,]
  #           c$ref_start<-a
  #           c$ref_end<-b
  #           c$query_end<-Dupex[i,]$query_start
  #           cc<-Dupex[which((c$query_end-Dupex$query_end)>0),]
  #           c$query_start<-cc$query_end[which.min(c$query_end - cc$query_end)]
  #           duplication<-rbind(duplication,c)
  #         }
  #       }
  #       
  #     }
  #   }
  # }
  # 
  aftertans<-smalltransf(endcluster0)
  ## 这里加一步，如果translocation在dup中，就把dup对应的序列删掉
  # if(length(aftertans$tran)!=0){
  #   data1<-aftertans$tran
  #   for(j in 2:dim(duplication)[1]){
  #     data2<-data1[data1$ref_start>=duplication[j,]$ref_start & data1$ref_end<=duplication[j,]$ref_end & data1$query_start>=duplication[j,]$query_start &data1$query_end<=duplication[j,]$query_end,]
  #     if(nrow(data2)!=0){
  #       duplication<-duplication[-j,]
  #     }
  #   }
  # }
  # 
  
  
  
  
  #dotplot_cluster(endcluster0)
  ###!!!!这里需要加一步提取内部的SDR
  for(i in unique((aftertans$tranbefore)$cluster)){
    a<-aftertans$tranbefore[(aftertans$tranbefore)$cluster==i,]
    a$cluster<-as.character(1:dim(a)[1])
    if(sum(a$orient == "-") > 0.5* nrow(a)){
      reverse_end<-reverse.region(a,chrid,3,"init")
    }
    if(sum(a$orient == "+") > 0.5 * nrow(a) |sum(a$orient == "-") == 0.5 * nrow(a)){
      reverse_end<-reverse.region(a,chrid,2,"init")
    }
    if(exists("reverse_end")){
      middle<-reverse_end$reverse
      if(nrow(middle)!=0){
        middle<-intersect.unit(middle,unique(middle$query_chr)) #将among中复杂片段作为一个大SDR
      }
      if(sum(a$orient == "+") > 0.5 * nrow(a) & nrow(middle)!=0){
        middle$orient<-"+"
      }
      if(sum(a$orient == "-") > 0.5 * nrow(a) & nrow(middle)!=0){
        middle$orient<-"-"
      }
      if(sum(a$orient == "-") == 0.5 * nrow(a) & nrow(middle)!=0){
        middle$orient<-"no"
      }
      
      storesmall<-rbind(storesmall,middle)
      minusmiddle<-reverse_end$reverse[reverse_end$reverse$ref_start<=reverse_end$reverse$ref_end & reverse_end$reverse$query_start<=reverse_end$reverse$query_end,]
      if(nrow(minusmiddle)!=0){
        minusmiddle$orient<-"no"
      }
      
      storesmall<-rbind(storesmall,minusmiddle)
      storesmall<-distinct(storesmall)
    }
  }
  if(length(aftertans$tran)!=0){
    smalltrans<-rbind(smalltrans,aftertans$tran)
  }
  if(!isEmpty(aftertans$tran)){
    inversion<-rbind(inversion,aftertans$tranbefore[aftertans$tranbefore$orient=="-",colnames(inversion)])
    transdup<-transduplication_extract(endcluster0,aftertans$tran)
    if(nrow(transdup)!=0){
      duplication<-rbind(duplication,transdup[,colnames(duplication)])
    }
    
 }
  
  
  
  endcluster0<-aftertans$endcluster0
  #\dotplot_cluster(endcluster1)
  ## 2.第二步，inversion提取（-）

  endcluster1<-split_region(endcluster0,1,cluster1paras)
  
  ## 构建一个函数重新调整minus的聚类
  endcluster1<-minus.next(endcluster1)
  
  i=0
  for(miusclu in unique(endcluster1[endcluster1$orient=='-',]$cluster)){
    print(miusclu)
    list<-which(endcluster1$cluster==miusclu)
    if(length(list)>=2){
      for(line in 2:length(list)){
        if(endcluster1[list[line],]$ref_start-endcluster1[list[line-1],]$ref_end >50000){
          endcluster1[list[line],]$cluster<-paste(endcluster1[list[line-1],]$cluster,i,sep="")
          i<-i+1
        }else{
          endcluster1[list[line],]$cluster<-endcluster1[list[line-1],]$cluster
        }
      }
    }
   
    endcluster1[endcluster1$cluster==miusclu,]
  }
  
  

  ## 如果两个inversion之间距离大于某一个值，就把它分开为两类
  
  
  
  inver<-inversion.extract(endcluster1,chrid)  ### 先把inversion提取出来
  if(nrow(inver)!=0){
    colnames(inver)[8]="orient"
  }
  
  if(nrow(duplication)!=1){
    duplication$cluster<-as.character(duplication$cluster)
    inver<-rbind(inver,duplication[duplication$orient=="-",colnames(inver)])
  }
  ##先把很乱的区域存储为SDR
  list<-which(rle(rle(endcluster1$cluster)$lengths)$values==1 &rle(rle(endcluster1$cluster)$lengths)$lengths >5)
  chaoval<-rle(rle(endcluster1$cluster)$lengths)$values
  chaolen<-rle(rle(endcluster1$cluster)$lengths)$lengths
  cross.region<-cross.calcu(endcluster1)
  initclsuer<-1000
  for(k in list){
    
    start<-sum(chaoval[1:k-1]*chaolen[1:k-1])+1
    endcluster1[intersect(which(!endcluster1$cluster %in% cross.region),start:(start+chaolen[k]-1)),]$cluster<-as.character(initclsuer)
    print(endcluster1[intersect(which(!endcluster1$cluster %in% cross.region),start:(start+chaolen[k]-1)),])
    chaos<-endcluster1[intersect(which(!endcluster1$cluster %in% cross.region),start:(start+chaolen[k]-1)),]
    chaos$query_start<-abs(chaos$query_start)
    chaos$query_end<-abs(chaos$query_end)
    if(nrow(chaos)!=0){
      chaos<-clusterall(chaos)
    }
    chaos$anno<-"COMPLEX"
    store<-rbind(store,chaos[,colnames(store)])
    initclsuer<-initclsuer+1
  }
  
  reverse_end<-reverse.region(endcluster1,chrid,0,"init")
  duplic<-reverse_end$dup
  if(!is.character(duplic)){
    if(nrow(duplic)!=0){
      duplication<-rbind(duplication,duplic[,colnames(duplication)])
    }
    
  }
  minimap<-reverse_end$minimaploc  ### !!!!!!!!不知道怎么求
  
  
  if ("delsytenic" %in% names(reverse_end)){
    for(cluster in unique(reverse_end$delsytenic$cluster)){
      datalist<-reverse_end$delsytenic[reverse_end$delsytenic$cluster==cluster,]
      datalist$cluster<-1:dim(datalist)[1]
      reverse_end1<-reverse.region(datalist,unique(datalist$ref_chr),2,"init")
      middle<-reverse_end1$reverse
      if(nrow(middle)!=0){
        middle<-intersect.unit(middle,unique(middle$query_chr)) #将among中复杂片段作为一个大SDR
        middle$orient<-"no"
        storesmall<-rbind(storesmall,middle)
      }
    }
  }
  
  
  for (value in unique(reverse_end$reverse$query_chr)) {
    datarev<-reverse_end$reverse[reverse_end$reverse$query_chr==value,]
    if(nrow(datarev)!=0){
      # 找到query_chr等于当前唯一值的行的索引
      rows_to_remove <- which(datarev$query_chr == value)
      # 如果有多行，删除第一行和最后一行
      if (length(rows_to_remove) > 1) {
        startend<-datarev[c(rows_to_remove[1], tail(rows_to_remove, 1)), ]
        store<-rbind(store,datarev[-c(rows_to_remove[1], tail(rows_to_remove, 1)), ]) ##为了防止有一条query对应多个ref最终导致前后计算的端粒过长，在去除store中含有的比对时误删
        startend$orient<-"+"
        storesmall<-rbind(storesmall,startend)
      }
    }
  }
  
  ##处理一下大cluster之间的store
  vectore<-which(store$ref_start-store$ref_end==2)
  if(length(vectore)!=0){
    for(l in vectore){
      store[l,]$ref_start<-store[l,]$ref_start-1
      store[l,]$ref_end<-store[l,]$ref_end+1
    }
  }
  vectore<-which((store$query_start-store$query_end)==2)
  if(length(vectore)!=0){
    for(l in vectore){
      store[l,]$query_start<-store[l,]$query_start-1
      store[l,]$query_end<-store[l,]$query_end+1
    }
  }
  
  if(nrow(store[store$ref_start>store$ref_end,])!=0){
    store[store$ref_start>store$ref_end,]$ref_end=store[store$ref_start>store$ref_end,]$ref_start #如果大片段有重叠，直接将其重叠删掉
  }
  initstore<-store[store$ref_chr!=0,]
  if(length(unique(store$query_chr)[-1])!=0){
    for(chr_child in unique(store$query_chr)[-1]){
      new_row<-intersect.unit(store,chr_child) #将among中复杂片段作为一个大SDR
      if(nrow(new_row)!=1){
        new_row<-complex(new_row)
      }
      percen<-max(new_row[new_row$anno=="COMPLEX",]$ref_end-new_row[new_row$anno=="COMPLEX",]$ref_start)
      if(percen>100000000){ ##说明全基因组都有跨染色体变异
        new_row<-intersect.unit(store,chr_child) #将among中复杂片段作为一个大SDR
      }
      assign(chr_child,new_row)
    }
    store<-docall(store)
    rm(list=unique(store$query_chr)) ##删除变量
  }
  store<-rbind(store,anti_join(initstore,store))
  
 if(dim(inver)[1]!=0){
    inver<-inver[,colnames(inversion)]
    inversion<-rbind(inversion,inver)
    inversion<-distinct(inversion)
  }
  
  cat("第一次大聚类后（Big cluster among）获得的的SDR数目是",dim(store)[1],"\n")  ##但是这个计算了每条染色体的最开始，如果是多个染色体的比对需要改一下
  cat("获得的的inversion数目是",dim(inversion)[1]-1,"\n") 
  cat("获得的的translocation数目是",dim(smalltrans)[1]-1,"\n") 

  endcluster1<-reverse_end$endcluster1
  #endcluster1<-endcluster1[endcluster1$cluster!="1000",]
  # statistic<-rle(endcluster1[endcluster1$ref_start>=store[store$anno=="COMPLEX",]$ref_start&
  #               endcluster1$ref_end<=store[store$anno=="COMPLEX",]$ref_end&
  #               endcluster1$query_start>=store[store$anno=="COMPLEX",]$query_start&
  #               endcluster1$query_end<=store[store$anno=="COMPLEX",]$query_end,]$cluster)
  ## 下面为如果complex区域中大多数都是比较分散的，就把比对中的含有COMPLEX的区域去掉
  # if (sum(statistic$values=="1" |statistic$values=="2")>0.8*length(statistic)){
  #   for(chr_child in unique(endcluster1$query_chr)){
  #     new_row<-alignintersect(endcluster1[endcluster1$query_chr==chr_child,],store[store$query_chr==chr_child,],store[store$query_chr==chr_child,]) 
  #     ## 这里重新写一下，如果交集大于70%则去掉
  #     assign(chr_child,new_row)
  #   }
  #   endcluster1<-docall(endcluster1)
  #   rm(list=unique(endcluster1$query_chr)) ##删除变量
  # }
  
  ## 现在需要二轮的聚类啦
  
  ## 如果比对中有找到的杂乱区域，就把这些比对去掉
  ## 在计算大SDR的时候把INV小的cluster聚到一起了，在这里想更精确一点
  
  m=0
  count_minus <- table(endcluster1$cluster[endcluster1$orient == '-'])
  count_total <- table(endcluster1$cluster[endcluster1$cluster%in% names(count_minus)])
  ##cluster
  for(clusid in names(count_minus[(count_minus / count_total) >0.6]) ){
    print(clusid)
    miuscluster<-split_region(endcluster1[endcluster1$cluster==clusid,],1,100000)
    minusduplic<-duplication_extract(endcluster1,miuscluster)  
    if(!isEmpty(minusduplic$dupli)){
      duplication<-rbind(duplication,minusduplic$dupli[,colnames(duplication)])
    }
    
    ##计算一下这些miuscluster中的duplication，因为之后可能直接分开成不同的cluster了，造成最终的结果有差异
    reverse_end<-reverse.region(miuscluster,chrid,3,"mius")
    middle<-reverse_end$reverse
    if(nrow(middle)!=0){
      #middle[middle$ref_start>middle$ref_end,]$ref_end=middle[middle$ref_start>middle$ref_end,]$ref_start #如果大片段有重叠，直接将其重叠删掉
      local<-intersect.unit(middle,unique(middle$query_chr)) #将among中复杂片段作为一个大SDR
      changed_rows <- anti_join(local,middle)
      store<-rbind(store,inner_join(local,middle))
      if(nrow(changed_rows)!=0){
        changed_rows$anno<-'COMPLEX'
        store<-rbind(store,changed_rows[,colnames(store)])
        store<-rbind(store,middle)
        store<-distinct(store)
      }}
    
    if(length(unique(miuscluster$cluster))!=1){
      miuscluster<-reverse_end$endcluster1
      miuscluster$cluster<-paste(miuscluster$cluster, "8",as.character(m), sep = "")
      endcluster1<-endcluster1[endcluster1$cluster!=clusid,]
      endcluster1<-rbind(endcluster1,miuscluster)
    }
    
    m=m+1
  }
  
  
  count_minus <- table(endcluster1$cluster[endcluster1$orient == '-'])
  count_total <- table(endcluster1$cluster[endcluster1$cluster %in% names(count_minus)])
  ## 在这里先找一下正负cluster之间的重复，相当于内部重复
  
  
  
  SDRminudlist<-names(count_minus[(count_minus / count_total) >0.6]) ##cluster
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
    #claster2<-(repeat.integrate(region,2)) ##将inversion返回去了
    # endcluster2<-claster2$afterdup
    # if(!is.null(claster2$repeat.region)){
    #   dup<-distinct(claster2$repeat.region)
    #   duplication<-rbind(duplication,dup)
    # }
    ## 从一堆-的中寻找为+的各自进行聚类
    if(dim(region)[1]==1){
      next
    }
    
    endcluster2<-duplication_extract(region,region)
    if(nrow(endcluster2$dupli)!=0){
      duplication<-rbind(duplication,endcluster2$dupli[,colnames(duplication)])
    }

    endcluster2before<-endcluster2$pos_end
    orientid="+"
    endcluster2<-smallcluster(endcluster2$pos_end,orientid) ##得到把负链中的相邻正链进行整合后计算
    reverse_end<-reverse.region(endcluster2,chrid,3,"init")
    middle<-reverse_end$reverse
    if(nrow(middle)!=0){
      #middle[middle$ref_start>middle$ref_end,]$ref_end=middle[middle$ref_start>middle$ref_end,]$ref_start #如果大片段有重叠，直接将其重叠删掉
      local<-intersect.unit(middle,unique(middle$query_chr)) #将among中复杂片段作为一个大SDR
      changed_rows <- anti_join(local,middle)
      if(nrow(changed_rows)!=0){
        changed_rows$anno<-'COMPLEX'
        store<-rbind(store,changed_rows[,colnames(store)])
        store<-rbind(store,middle)
        store<-distinct(store)
      }
      
      if(length(which(middle$query_start>middle$query_end))!=0){
        middle<-middle[-(which(middle$query_start>middle$query_end)),]
        
      }
      if(nrow(middle)!=0){
        middle$orient<-"-"
        storesmall<-rbind(storesmall,middle[,colnames(storesmall)]) ## 存到smalltrans里，防止存在store中会将synetenic中的区域删掉
      }
      
    }
    storesmall<-insertsmall(endcluster2before,storesmall,orientid)
  }
  
  
  
  
  for(k in setdiff(unique(endcluster1$cluster),SDRminudlist)){
    region<-endcluster1[endcluster1$cluster==k,]
    region$cluster<-1:dim(region)[1]
    if(dim(region)[1]==1){
      next
    }
    
    endcluster2<-duplication_extract(region,region)
    if(nrow(endcluster2$dupli)!=0){
      duplication<-rbind(duplication,endcluster2$dupli[,colnames(duplication)])
    }

    #claster2<-(repeat.integrate(region,1))
    #endcluster2<-claster2$afterdup
    # if(!is.null(claster2$repeat.region)){
    #   dup<-distinct(claster2$repeat.region)
    #   duplication<-rbind(duplication,dup)
    # }
    endcluster2before<-endcluster2$pos_end
    orientid="-"
    endcluster2<-smallcluster(endcluster2$pos_end,orientid) ##得到把负链中的相邻正链进行整合后计算
    reverse_end<-reverse.region(endcluster2,chrid,2,"init")
    middle<-reverse_end$reverse
    local<-intersect.unit(middle,unique(middle$query_chr)) #将among中复杂片段作为一个大SDR
    changed_rows <- anti_join(local,middle)
    if(nrow(changed_rows)!=0){
      print(k)
      changed_rows$anno<-'COMPLEX'
      store<-rbind(store,changed_rows[,colnames(store)])
      
    }
    # middle[middle$ref_start>middle$ref_end,]$ref_end=middle[middle$ref_start>middle$ref_end,]$ref_start #如果大片段有重叠，直接将其重叠删掉
    # middle<-intersect.unit(middle,unique(middle$query_chr)) #将among中复杂片段作为一个大SDR
    if(length(which(middle$query_start>middle$query_end))!=0){
      middle<-middle[-(which(middle$query_start>middle$query_end)),]
    }
    if(nrow(middle)!=0){
      middle$orient<-"+"
      storesmall<-rbind(storesmall,middle[,colnames(storesmall)])
    }
    
    storesmall<-insertsmall(endcluster2before,storesmall,orientid)
    
  }
  cat("第二次小聚类后一起获得的的SDR数目是",dim(store)[1],"\n")
  cat("第二次小聚类后一起获得的的INV数目是",dim(inversion)[1],"\n")
  cat("第二次小聚类后一起获得的的DUP数目是",dim(duplication)[1],"\n")
  cat("第二次小聚类后一起获得的TRAN数目是",dim(smalltrans)[1],"\n")
  inversion$anno<-"INV"
  inversion$orient<-"-"
  duplication<-distinct(duplication)
  duplication$anno<-"DUP"
  duplication$orient<-"no"
  smalltrans$anno<-"TRANS"
  smalltrans$orient<-"no"
  storesmall$anno<-"SDR_NM"
  store$orient<-"no"
  storehighdup_sdr$orient<-"no"
  duplication<-duplication[,colnames(store)]
  smalltrans<-smalltrans[,colnames(store)]
  storesmall<-storesmall[-1,colnames(store)]
  store<-store[store$ref_chr!=0,]
  storehighdup_sdr<-storehighdup_sdr[storehighdup_sdr$ref_chr!=0,]
  all<-rbind(store,inversion[-1,],distinct(duplication[-1,]),smalltrans[-1,],storesmall,storehighdup_sdr)
  all<-all[order(all$ref_start),]
  for (chr_child in unique(all$query_chr[all$query_chr!=0])){
    assign(chr_child,endfilter(all[all$query_chr==chr_child,],chrid,chr_child))
    
  }
  data<-docall(all)
  if(length(which(data$ref_start>data$ref_end))!=0){
    data<-data[-which(data$ref_start>data$ref_end),]
  }
  data<-distinct(data)
  if(length(which(data$reflen==0 & data$querylen==0))!=0){
    data<-data[-which(data$reflen==0 & data$querylen==0),]
  }
  COMPLEX<-data[data$anno=="COMPLEX",]
  if(nrow(data[data$anno=="COMPLEX",])!=0){
    data<-data[-which(data$anno=="COMPLEX"),]
    newcomplex<-complexinte(COMPLEX)
    newcomplex$anno<-"COMPLEX"
    newcomplex$orient<-"no"
    newcomplex$reflen<-newcomplex$ref_end-newcomplex$ref_start
    newcomplex$querylen<-newcomplex$query_end-newcomplex$query_start
    data<-rbind(data,newcomplex[,colnames(data)])
  }
  ## 加一步计算一下breakpoints
  INVb<-inversion[inversion$ref_chr!=0,]
  if(nrow(INVb)!=0){
    INVb$ref_p_1s<-0
    INVb$ref_p_1e<-0
    INVb$que_p_1s<-0
    INVb$que_p_1e<-0
    INVb$ref_p_2s<-0
    INVb$ref_p_2e<-0
    INVb$que_p_2s<-0
    INVb$que_p_2e<-0
    for(i in 1:dim(INVb)[1]){
      invpos<-max(which(pos.chr$ref_start==INVb[i,]$ref_start & pos.chr$orient=='-'))
      if(invpos!=1){
        INVb[i,]$ref_p_1s<-pos.chr[invpos-1,]$ref_end
        INVb[i,]$ref_p_1e<-pos.chr[invpos,]$ref_start
        INVb[i,]$que_p_1s<-pos.chr[invpos-1,]$query_end
        INVb[i,]$que_p_1e<-pos.chr[invpos,]$query_start
      }
     
      invpos<-max(which(pos.chr$ref_end==INVb[i,]$ref_end & pos.chr$orient=='-'))
      if(invpos!=dim(pos.chr)[1]){
        INVb[i,]$ref_p_2s<-pos.chr[invpos,]$ref_end
        INVb[i,]$ref_p_2e<-pos.chr[invpos+1,]$ref_start
        INVb[i,]$que_p_2s<-pos.chr[invpos,]$query_end
        INVb[i,]$que_p_2e<-pos.chr[invpos+1,]$query_start
      }
      
    }
    INVball<-rbind(INVball,INVb)
  }
data<-data[data$ref_start<=data$ref_end & data$query_start<=data$query_end,]
if(nrow(highdupall[highdupall$ref_chr!=0,])!=0){
  hidup<-highdupall[highdupall$ref_chr!=0,]
  hidup$query_chr=""
  hidup$query_start=""
  hidup$query_end=""
  hidup$anno="high-dup"
  hidup$orient=""
  hidup$reflen=""
  hidup$querylen=""
  data<-rbind(data,hidup)
}
if(nrow(data[data$ref_start==0 & data$ref_end==0,])!=0){
  data<-data[!(data$ref_start==0 & data$ref_end==0),]
}
if(nrow(data[data$query_start==0 & data$query_end==0,])!=0){
  data<-data[!(data$query_start==0 & data$query_end==0),]
}
write.table(data, paste(args[4],chrid,"end.tsv",sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)
#write.table(minimap, paste("D:/MS/saffire测试/R/minimap/",chrid,"minimap.tsv",sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)
}
colnames(INVball)<-c("ref_chr","ref_start","ref_end","query_chr","query_start",
"query_end"  , "anno"     ,   "orient"   ,   "ref_p_1s"  ,  "ref_p_1e"  ,
"que_p_1s"  ,  "que_p_1e"  ,  "ref_p_2s"   , "ref_p_2e"  ,  "que_p_2s"  ,
"que_p_2e")
write.table(INVball, paste(args[4],"INVresult.tsv",sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)
#write.csv(result_data[result_data$SV=="DUP",], paste("D:/MS/saffire测试/R/result/SDRend.csv",sep = ""), quote = FALSE, row.names = FALSE)
