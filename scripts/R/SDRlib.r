## 定义两个小函数
exchange<-function(data){  ##这个函数用于交换开始和结束不是小于关系的比对
  data <- data[data$ref_start != 0, ]
  data[data$ref_start > data$ref_end, c('ref_start', 'ref_end')] <-  data[data$ref_start > data$ref_end, c('ref_end','ref_start')]
  data[data$query_start > data$query_end, c('query_start', 'query_end')] <-  data[data$query_start > data$query_end, c('query_end','query_start')] 
return(data)
}
## 这个函数用于将有交集的SDR结果整合为大片段SDR

intersect.unit<-function(data,chr_child){
  rm(inte)
  data<-data[data$query_chr==chr_child,]
  data<-exchange(data)
  ir1ref <- IRanges(start = data$ref_start, end = data$ref_end)
  ir1que  <- IRanges(start = data$query_start, end = data$query_end)
  overlapsref<-findOverlaps(ir1ref,ir1ref) 
  overlapsque<-findOverlaps(ir1que,ir1que) 
  if(length(which(overlapsref@from!=overlapsref@to))!=0){
    inte<-as.data.frame(overlapsref[which(overlapsref@from!=overlapsref@to)])
    }
  if(length(which(overlapsque@from!=overlapsque@to))!=0){
  if(exists("inte")){
    inte<-rbind(inte,as.data.frame(overlapsque[which(overlapsque@from!=overlapsque@to)]))
  }
  else{
    inte<-as.data.frame(overlapsque[which(overlapsque@from!=overlapsque@to)])
  }
  }
  if(exists("inte")){
    df_sorted <- t(apply(inte, 1, function(x) sort(x)))
    unique_df <- unique(df_sorted)
    inte <- as.data.frame(unique_df)
    for(j in 1:dim(inte)[1]){
      a<-c()
      a<-append(a,inte[j,1])
      a<-append(a,inte[j,2])
      for(k in 1:dim(inte)[1]){
        if(inte[k,1] %in% ((min(a)-1):(max(a)+1))){
          a<-unique(append(a,inte[k,1]))
        }
        if(inte[k,2] %in% ((min(a)-1):(max(a)+1))){
          a<-unique(append(a,inte[k,2]))
        }
      }
      if(j==1){
        x<-a
        data[min(x),]$ref_start<-min(data[x,]$ref_start)
        data[min(x),]$ref_end<-max(data[x,]$ref_end)
        data[min(x),]$query_start<-min(data[x,]$query_start)
        data[min(x),]$query_end<-max(data[x,]$query_end)
        data[x,]<-data[min(x),]
        print(x)
      }
      else{
        if(length(intersect(a,x))==0){
          x<-a
          data[min(x),]$ref_start<-min(data[x,]$ref_start)
          data[min(x),]$ref_end<-max(data[x,]$ref_end)
          data[min(x),]$query_start<-min(data[x,]$query_start)
          data[min(x),]$query_end<-max(data[x,]$query_end)
          data[x,]<-data[min(x),]
        }
      }
    }
    data<-distinct(data)
  }
  data<-data[order(data$ref_start),]
  return(data)
}


##----------------------function1 定义一个函数：去除端粒和着丝粒区域的比对（有交集就去除）
# align1<-endcluster1[endcluster1$query_chr==chr_child,]
# align2<-store[store$query_chr==chr_child,]
# align3<-store[store$query_chr==chr_child,]
alignintersect<-function(align1,align2,align3){
del.list<-c()
cat("-------开始去除着丝粒和端粒区域的比对---------\n")

ir1ref <- IRanges(start = align1$ref_start, end = align1$ref_end)
ir1que  <- IRanges(start = align1$query_start, end = align1$query_end)
ir2 <- IRanges(start =align2$ref_start , end =align2$ref_end)
overlapsref<-findOverlaps(ir1ref,ir2) ##我们的比对序列的参考基因组序列和人的端粒|着丝粒的交集
if(length(overlapsref@from)!=0){  
  xx<-cbind(overlapsref@from,
            align1[overlapsref@from,]$ref_chr,
            align1[overlapsref@from,]$ref_start,
            align1[overlapsref@from,]$ref_end,
            overlapsref@to,
            align2[overlapsref@to,]$ref_chr,
            align2[overlapsref@to,]$ref_start,
            align2[overlapsref@to,]$ref_end)
  xx<-as.data.frame(xx)
  for(i in 1:dim(xx)[1]){
    if(xx$V3[i]>=xx$V7[i] & xx$V4[i]<=xx$V8[i]){  ##说明store中的结果把比对的结果覆盖了
      del.list<-append(del.list,xx$V1[i])
    }
  }
  
}

ir3 <- IRanges(start = align3$query_start, end = align3$query_end)
overlapsque<-findOverlaps(ir1que,ir3)
if(length(overlapsque@from)!=0){  
  xx<-cbind(overlapsque@from,
            align1[overlapsque@from,]$query_chr,
            align1[overlapsque@from,]$query_start,
            align1[overlapsque@from,]$query_end,
            overlapsque@to,
            align3[overlapsque@to,]$query_chr,
            align3[overlapsque@to,]$query_start,
            align3[overlapsque@to,]$query_end)
  xx<-as.data.frame(xx)
  for(i in 1:dim(xx)[1]){
    if(xx$V3[i]>=xx$V7[i] & xx$V4[i]<=xx$V8[i]){  ##说明store中的结果把比对的结果覆盖了
      del.list<-append(del.list,xx$V1[i])
    }
  }
  
}
 ##我们的比对序列的参考基因组序列和猩猩的端粒|着丝粒的交集
if(length(del.list)!=0)
{align1<-align1[-as.numeric(del.list),]}
return(align1)
}

reverse_xy<-function(data){
  selected_rows <- data$orient == '-'
  new_query_start <- data[selected_rows, ]$query_end
  new_query_end <- data[selected_rows, ]$query_start
  data[selected_rows, ]$query_start <- new_query_start
  data[selected_rows, ]$query_end <- new_query_end
  return(data)
}
#-------------------function4 split 并聚类 ----------------------
split_region<-function(pos.chr.region,cluster.id,clusterparas){
  pos.chr.region<-reverse_xy(pos.chr.region)  #负链反向
  df <- pos.chr.region %>%
    rowwise() %>%
    mutate(length = sqrt((ref_end - ref_start)^2 + (query_end- query_start)^2),
           segments = ceiling(length / 5000))
  df_segments <- df[rep(seq_len(nrow(df)), df$segments), ]
  df_segments <- df_segments %>%
    group_by(ref_start, ref_end, query_start, query_end) %>%
    mutate(segment_id = row_number(),
           segment_mid_x = ref_start + (ref_end - ref_start) * (segment_id - 0.5) / segments,
           segment_mid_y = query_start + (query_end - query_start) * (segment_id - 0.5) / segments)
  
  if(cluster.id==2){
    dbs.para<-10000
  }
  if(cluster.id==0){
    dbs.para<-clusterparas
  
  }
  if(cluster.id==1){
    df_segments[df_segments$orient=='-',]$segment_mid_y<--df_segments[df_segments$orient=='-',]$segment_mid_y
    df_segments[df_segments$orient=='-',]$query_start<--df_segments[df_segments$orient=='-',]$query_start
    df_segments[df_segments$orient=='-',]$query_end<--df_segments[df_segments$orient=='-',]$query_end
    pos.chr.region[pos.chr.region$orient=='-',]$query_start<--pos.chr.region[pos.chr.region$orient=='-',]$query_start
    pos.chr.region[pos.chr.region$orient=='-',]$query_end<--pos.chr.region[pos.chr.region$orient=='-',]$query_end
    dbs.para<-clusterparas}
  dbscan_result <- dbscan(df_segments[,c("segment_mid_x", "segment_mid_y")], eps = dbs.para, minPts =1)
  df_segments$cluster <- dbscan_result$cluster
  df_segments <- df_segments %>%
    group_by(ref_chr,ref_start, ref_end, ref_pos,query_chr,query_start, query_end,query_pos,orient) %>%
    summarise(cluster = names(which.max(table(cluster))))
  df_segments$query_start<-abs(df_segments$query_start)
  df_segments$query_end<-abs(df_segments$query_end)
  df_segments<-reverse_xy(df_segments)  #负链反向
  return(df_segments)
}

#-------------------function5 画图 ----------------------
# cluster.id标识是大聚类还是小聚类，大聚类1考虑大片段的替换，不考虑inversion,小聚类2考虑小片段，考虑inversion
dotplot_cluster<-function(plotpos,region){
  plotpos<-reverse_xy(plotpos) 
if(!missing(region)){plotpos<-plotpos[plotpos$ref_end<region[2] & plotpos$ref_start>region[1],]}

  # 创建图形并绘制线段，根据 cluster 列进行着色
  # plot <- ggplot(segments, aes(x = NULL, y = NULL)) +
  #   geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2, color = as.factor(cluster)), size = 1, arrow = arrow(length = unit(0.1, "cm"))) +
  #   # 添加标题和标签
  #   ggtitle("Segments") +
  #   xlab("X-axis") +
  #   ylab("Y-axis")
  # 
  # # 显示图形
  # print(plot)
  plot <- ggplot(plotpos, aes(x = NULL, y = NULL)) +
    geom_segment(aes(x = ref_start, y = query_start, xend = ref_end, yend = query_end,color = as.factor(cluster)), size = 1, arrow = arrow(length = unit(0.3, "cm"))) +
    ggtitle("Segments") +
    xlab("X-axis") +
    ylab("Y-axis")
  plotly_obj <- ggplotly(plot, tooltip = c("ref_start", "query_start", "ref_end", "query_end"))
  print(plotly_obj)
}


### function 将小cluster对应的负链连续出现的位置进行整合，如果负链想整合的区域中有小片段的align就删掉
inte.minud<-function(endcluster1,id){
  if(id==1){
    if(length( which(endcluster1$orient=="-"))!=0){
      endcluster1<-endcluster1[-(which(endcluster1$orient=="-")),]
    }
  }
  if(id==2){
    if(length( which(endcluster1$orient=="+"))!=0){
      endcluster1<-endcluster1[-(which(endcluster1$orient=="+")),]
    }
  }
  
  
  # len<-rle(endcluster1[['orient']])$lengths
  # mask<-rle(endcluster1[['orient']])$values
  # ## 整合负链
  # num<-which(mask=='+' & len<3)[which(mask=='+' & len<3)!=1 & which(mask=='+' & len<3)!=length(len)]
  # for(i in num){
  #   len[i-1]<-len[i-1]+len[i]
  # }
  # minus_rows <- which(mask=='-' & len>=1)
  # minus_rows <-minus_rows[minus_rows!=1 & minus_rows!=length(len)]
  # for(i in minus_rows){
  #   i = minus_rows
  #   from=sum(len[1:i-1])+1
  #   loc<- seq(from = from, by = 1, length.out = len[i])
  #   endcluster1[loc,]$ref_start<-min(endcluster1[loc,]$ref_start)
  #   endcluster1[loc,]$ref_end<-max(endcluster1[loc,]$ref_end)
  #   endcluster1[loc,]$query_start<-min(endcluster1[loc,]$query_start)
  #   endcluster1[loc,]$query_end<-max(endcluster1[loc,]$query_end)
  #   endcluster1[loc,]$cluster<-0
  #   endcluster1[loc,]$orient<-'-'
  # }
  # endcluster1<-distinct(endcluster1)
  return(endcluster1)
}



docall<-function(data){
  if(length(unique(data$query_chr)[unique(data$query_chr)!='0'])!=1){
    data<-do.call(rbind,mget(unique(data$query_chr)[unique(data$query_chr)!='0'], envir = .GlobalEnv))
  }
  else{
    data<-get( unique(data$query_chr)[unique(data$query_chr)!='0'])
  }
  return(data)
}



##----function 专门提取cluster2对应的inversion
inversion.extract<-function(endcluster1,chrid){
  endcluster1$query_start<-abs(endcluster1$query_start)
  endcluster1$query_end<-abs(endcluster1$query_end)
  inversion<-endcluster1[endcluster1$orient=="-",] %>% group_by(cluster)%>%  #聚cluster
    summarise(ref_chr=ref_chr,ref_start=min(ref_start),
              ref_end=max(ref_end),
              query_chr=query_chr,
              query_starttem=min(pmin(query_end,query_start)),
              query_endtem=max(pmax(query_start,query_end)),
              (names(which.max(table(orient))))) ###正链负链？？重新算一下
  colnames(inversion)[6]<-"query_start"
  colnames(inversion)[7]<-"query_end"
  inversion<-distinct(inversion) ##在进行处理之前先把inversion找到，已经找了全部的了
  return(inversion)
  
}

## 这里加一个函数，如果前后都是-并且是反转的，就把这两个类聚在一起
clusterbigminus<-function(endcluster1,cluster.id){
  rm("xx")
  pos_end<-endcluster1 %>% group_by(cluster)%>%  #聚cluster
    summarise(ref_chr=ref_chr,ref_start=min(ref_start),ref_end=max(ref_end),query_chr=query_chr,query_start=min(query_start),query_end=max(query_end),(names(which.max(table(orient))))) ###正链负链？？重新算一下
  # pos_end<-pos_end[,-1]
  # pos_end<-pos_end[,-7]
  pos_end<-distinct(pos_end)
  pos_end<-pos_end[order(pos_end$ref_start),]
  if(cluster.id==0){
    newpos_end<-pos_end
    for(chr_child in unique(pos_end$query_chr)){
      minus<-which(pos_end$`(names(which.max(table(orient))))`=="-" & pos_end$query_chr==chr_child)
      if(length(minus)!=1 &length(minus)!=0){
        for(i in 2:length(minus)){
          if((minus[i]==minus[i-1]+1) & pos_end[minus[i],]$query_end<=pos_end[minus[i-1],]$query_start ){ ## 说明这两个-相邻且为反转
            endcluster1[endcluster1$cluster==newpos_end[minus[i],]$cluster,]$cluster<-newpos_end[minus[i-1],]$cluster
            newpos_end[minus[i],]$cluster<-newpos_end[minus[i-1],]$cluster
            colnames(newpos_end)[8]<-"orient"
            if(exists("xx")){
              xx<-rbind(xx,newpos_end %>% group_by(cluster)%>%  #聚cluster
                          summarise(ref_chr=ref_chr,ref_start=min(ref_start),ref_end=max(ref_end),query_chr=query_chr,query_start=min(query_start),query_end=max(query_end),(names(which.max(table(orient))))))
            }else{
              xx<-newpos_end %>% group_by(cluster)%>%  #聚cluster
                summarise(ref_chr=ref_chr,ref_start=min(ref_start),ref_end=max(ref_end),query_chr=query_chr,query_start=min(query_start),query_end=max(query_end),(names(which.max(table(orient))))) ###
              xx<-distinct(xx)
            }
            }
        }
        
      }
    }
  }
  if(exists("xx")){
    return(list(endcluster1=endcluster1,pos_end=xx))
  }else{
    return(list(endcluster1=endcluster1,pos_end=pos_end))
  }
}



#-------------------function6 获取大cluster前中后的区域SDR以及小cluster的中间区域! ----------------------
reverse.region<-function(endcluster1,chrid,cluster.id){
  endcluster1$query_start<-abs(endcluster1$query_start)
  endcluster1$query_end<-abs(endcluster1$query_end)
  ## cluster.id筛选的为大聚类还是小聚类
  if(cluster.id==0){
    ## 如果大聚类的cluster中包含了一些小聚类我们把它划为大的
    for(chr_child in unique(endcluster1$query_chr)){
      x=table(endcluster1[endcluster1$query_chr==chr_child,]$cluster)
      thereshod<-quantile(sort(as.vector(x)),0.8)/5
      if((thereshod)<5){ ## 说明比对数较少eg，1  10  11   2   3   4   5   6   7   8   9 
                                        # 390   1   2   1   2 838   1   1   2   1   4 
        thereshod=5
      }
      for(k in names(x)[x<=thereshod]){
        print(k)
        if(k==1){
          next
        }
        for(j in names(x)[x>thereshod]){
          a<-which(endcluster1$cluster==k & endcluster1$query_chr==chr_child )
          b<-which(endcluster1$cluster==j & endcluster1$query_chr==chr_child)
          if(min(b)<min(a) & max(b)>max(a)){
            if(unique(endcluster1[endcluster1$cluster==k & endcluster1$query_chr==chr_child,]$orient)=='-'){
              endcluster1[endcluster1$cluster==k & endcluster1$query_chr==chr_child,]$ref_pos<-(min(endcluster1[endcluster1$cluster==k,]$ref_start)+max(endcluster1[endcluster1$cluster==k,]$ref_end))/2
              endcluster1[endcluster1$cluster==k & endcluster1$query_chr==chr_child,]$ref_start<-min(endcluster1[endcluster1$cluster==k,]$ref_start)
              endcluster1[endcluster1$cluster==k & endcluster1$query_chr==chr_child,]$ref_end<-max(endcluster1[endcluster1$cluster==k,]$ref_end)
              minval=min(pmin(endcluster1[endcluster1$cluster==k & endcluster1$query_chr==chr_child,]$query_end,endcluster1[endcluster1$cluster==k & endcluster1$query_chr==chr_child,]$query_start))
              maxval=max(pmax(endcluster1[endcluster1$cluster==k & endcluster1$query_chr==chr_child,]$query_end,endcluster1[endcluster1$cluster==k & endcluster1$query_chr==chr_child,]$query_start))
              endcluster1[endcluster1$cluster==k & endcluster1$query_chr==chr_child,]$query_start<-minval
              endcluster1[endcluster1$cluster==k & endcluster1$query_chr==chr_child,]$query_end<-maxval
              endcluster1[endcluster1$cluster==k & endcluster1$query_chr==chr_child,]$query_pos<-(minval+maxval)/2
            }
            endcluster1[(which(endcluster1$cluster==k & endcluster1$query_chr==chr_child)),]$cluster<-j
            endcluster1<-distinct(endcluster1)
            break
          }
        }
      }
    }
    
   
    ## 处理一下边缘
    rm("dupli")
    for (iteration in 1:3) { ##多循环几次，防止出现删除之后后面还有重复的情况
      del.list<-c()
    for(chr_child in unique(endcluster1$query_chr)){
      for(k in unique(endcluster1[endcluster1$query_chr==chr_child,]$cluster)[-1]){
        datalist<-which(endcluster1$cluster==k  & endcluster1$query_chr==chr_child)
        if(endcluster1[min(datalist),]$ref_start<endcluster1[min(datalist)-1,]$ref_end){
          intervalue<-abs(endcluster1[min(datalist),]$ref_start-endcluster1[min(datalist)-1,]$ref_end)
          allval<-endcluster1[min(datalist),]$ref_end-endcluster1[min(datalist),]$ref_start
          allval2<-endcluster1[min(datalist)-1,]$ref_end-endcluster1[min(datalist)-1,]$ref_start
          if(intervalue/allval>0.7){
            del.list<-append(del.list,min(datalist))
          }
          else if(intervalue/allval2>0.7){
            del.list<-append(del.list,min(datalist)-1)
          }
          else{
            endcluster1[min(datalist),]$ref_start<-endcluster1[min(datalist)-1,]$ref_end+1
          }
        }
      }
      
    }
      if(length(del.list)!=0){
        if(exists("dupli")){
          dupli<-rbind(dupli,endcluster1[del.list,])
        }else{dupli<-endcluster1[del.list,]}
        endcluster1<-endcluster1[-del.list,]}
    }
    if(!exists("dupli")){dupli<-"no"}
    
    # for(k in unique(endcluster1[endcluster1$orient=='-',]$cluster)){
    #   if(k==1){
    #     next
    #   }
    #   a<-endcluster1[min(which(endcluster1$cluster==k))-1,]$cluster
    #   print('start')
    #   print(k)
    #   print(a)
    #   b<-endcluster1[max(which(endcluster1$cluster==k))+1,]$cluster
    #   print(b)
    #   if((!is.na(b)& length(b)!=0)  &  (!is.na(a) & length(a)!=0 ) ){
    #     if(a==b){
    #       
    #     }
    #   }
    #   
    #   #if(nrow(endcluster1[endcluster1$cluster==k,])!=0){endcluster1[min(which(endcluster1$cluster==k,)):max(which(endcluster1$cluster==k,)),"cluster"]=k}
    # }
    # for(k in unique(endcluster1[endcluster1$orient=='+',]$cluster)){
    #   if(k==1){
    #     next
    #   }
    #   a<-endcluster1[min(which(endcluster1$cluster==k))-1,]$cluster
    #   print('start')
    #   print(k)
    #   print(a)
    #   b<-endcluster1[max(which(endcluster1$cluster==k))+1,]$cluster
    #   print(b)
    #   if((!is.na(b)& length(b)!=0)  &  (!is.na(a) & length(a)!=0 ) ){
    #     if(a==b){
    #       endcluster1[(which(endcluster1$cluster==k)),]$cluster<-a
    #     }
    #   }
    #   #if(nrow(endcluster1[endcluster1$cluster==k,])!=0){endcluster1[min(which(endcluster1$cluster==k,)):max(which(endcluster1$cluster==k,)),"cluster"]=k}
    # }
    
    ## 处理一下边缘的inversion：1用
    if(cluster.id==1){
      invdata<-endcluster1[endcluster1$orient=="-",]
      for(k in unique(invdata$cluster)){
        if(k==1){
          next
        }
        temdata<-invdata[invdata$cluster==k,]
        mark<-which(endcluster1$cluster==k & endcluster1$orient=="-")
        if(endcluster1[min(mark),]$ref_start<endcluster1[min(mark)-1,]$ref_end | endcluster1[min(mark),]$query_start<endcluster1[min(mark)-1,]$query_end){
          endcluster1[mark,]$cluster<-endcluster1[min(mark)-1,]$cluster
        }
        if(endcluster1[max(mark),]$ref_end>endcluster1[max(mark)+1,]$ref_start | endcluster1[min(mark),]$query_end>endcluster1[min(mark)+1,]$query_start){
          endcluster1[mark,]$cluster<-endcluster1[max(mark)+1,]$cluster
        }
      }
    }

    #dotplot_cluster(endcluster1)
    }

  testminus<-clusterbigminus(endcluster1,cluster.id)
  endcluster1<-testminus$endcluster1
  pos_end<-testminus$pos_end
  pos_end<-pos_end[order(pos_end$ref_start),]
  reall<-merge(pos_end,chrpc,by=c("ref_chr","query_chr"))
  reall<-reall[order(reall$ref_start),]
  for(chr_child in unique(endcluster1$query_chr)){
  ## 第一行
  chrchch<-reall[reall$query_chr==chr_child,]
  new_row <- data.frame(ref_chr=chrid,ref_start = 1, ref_end = as.numeric(chrchch$ref_start[1])-1,query_chr=chr_child,query_start = 1, query_end = as.numeric(min(as.numeric(chrchch$query_start))-1)) ##改为初始queryend值为querystart最小值减1
  if(cluster.id==0){
    if(dim(chrchch)[1]!=1){
      for (i in 1:dim(chrchch)[1]){
        ## 序列前后相同的情况
        ## ref还是原来的方法，query寻找最近的query
          nowvalue=chrchch$query_end[i]
          ## 如果当前的queryend是最大值，则把它作为端粒前的部分
          if(nowvalue==max(chrchch$query_end)){
            if(i==dim(chrchch)[1]){ ##如果是位于最后一个的queryend
              xx=data.frame(ref_chr=chrchch$ref_chr[i],ref_start = as.numeric(chrchch$ref_end[i])+1, ref_end = chrchch$ref_len[i],query_chr=chrchch$query_chr[i],query_start = as.numeric(chrchch$query_end[i])+1,  query_end = chrchch$query_len[i])
            }
            else{
              xx=data.frame(ref_chr=chrchch$ref_chr[i],ref_start = as.numeric(chrchch$ref_end[i])+1, ref_end = as.numeric(chrchch$ref_start[i+1])-1,query_chr=chrchch$query_chr[i],query_start = as.numeric(chrchch$query_end[i])+1,  query_end = chrchch$query_len[i])
            }
          }else{
            query_endlist<-chrchch$query_start[which(chrchch$query_start>=nowvalue)]
            queendval<-query_endlist[which.min(abs(query_endlist -nowvalue))]
            if(i==dim(chrchch)[1]){
              xx=data.frame(ref_chr=chrchch$ref_chr[i],ref_start = as.numeric(chrchch$ref_end[i])+1, ref_end = chrchch$ref_len[i],query_chr=chrchch$query_chr[i],query_start = as.numeric(chrchch$query_end[i])+1, query_end = queendval-1)
            }else{
              xx=data.frame(ref_chr=chrchch$ref_chr[i],ref_start = as.numeric(chrchch$ref_end[i])+1, ref_end = as.numeric(chrchch$ref_start[i+1])-1,query_chr=chrchch$query_chr[i],query_start = as.numeric(chrchch$query_end[i])+1, query_end = queendval-1)
            }
          }
          colnames(xx)<-colnames(new_row)
          new_row <-rbind(new_row,xx)
        }
        
    }
    else{
      for (i in 1:dim(chrchch)[1]){
        xx=data.frame(ref_chr=chrchch$ref_chr[i],ref_start = as.numeric(chrchch$ref_end[i])+1, ref_end = chrchch$ref_len[i],query_chr=chrchch$query_chr[i],query_start = as.numeric(chrchch$query_end[i])+1,  query_end = chrchch$query_len[i])
        colnames(xx)<-colnames(new_row)
        new_row <-rbind(new_row,xx)
        }
    }
    
    }
  if(cluster.id==2){
    if(dim(chrchch)[1]!=1){
      for (i in 2:dim(chrchch)[1]-1){
        ## 序列前后相同的情况
        
        ## ref还是原来的方法，query寻找最近的query
        if(as.numeric(chrchch$ref_end[i])==chrchch$ref_start[i+1]){
          xx=data.frame(ref_chr=chrchch$ref_chr[i],ref_start = as.numeric(chrchch$ref_end[i]), ref_end =as.numeric(chrchch$ref_start[i+1]),query_chr=chrchch$query_chr[i],query_start = as.numeric(chrchch$query_end[i])+1, query_end = as.numeric(chrchch$query_start[i+1])-1)
        }else if(as.numeric(chrchch$query_end[i])==chrchch$query_start[i+1]){
          xx=data.frame(ref_chr=chrchch$ref_chr[i],ref_start = as.numeric(chrchch$ref_end[i])+1, ref_end = as.numeric(chrchch$ref_start[i+1])-1,query_chr=chrchch$query_chr[i],query_start = as.numeric(chrchch$query_end[i]), query_end = as.numeric(chrchch$query_start[i+1]))
        }else{
            xx=data.frame(ref_chr=chrchch$ref_chr[i],ref_start = as.numeric(chrchch$ref_end[i])+1, ref_end = as.numeric(chrchch$ref_start[i+1])-1,query_chr=chrchch$query_chr[i],query_start = as.numeric(chrchch$query_end[i])+1, query_end =as.numeric(chrchch$query_start[i+1])-1)
        }
        colnames(xx)<-colnames(new_row)
        new_row <-rbind(new_row,xx)
      }
    }
  }
  
  
  
  if(cluster.id==3){
    ## 3代表是大inversion
    if(dim(chrchch)[1]!=1){
      for (i in 2:dim(chrchch)[1]-1){
        ## 序列前后相同的情况
        if(as.numeric(chrchch$ref_end[i])==chrchch$ref_start[i+1]){
          xx=data.frame(ref_chr=chrchch$ref_chr[i],ref_start = as.numeric(chrchch$ref_end[i]), ref_end =as.numeric(chrchch$ref_start[i+1]),query_chr=chrchch$query_chr[i],query_start = as.numeric(chrchch$query_end[i+1])+1, query_end = as.numeric(chrchch$query_start[i])-1)
        }else if(as.numeric(chrchch$query_start[i])==chrchch$query_end[i+1]){
          xx=data.frame(ref_chr=chrchch$ref_chr[i],ref_start = as.numeric(chrchch$ref_end[i])+1, ref_end = as.numeric(chrchch$ref_start[i+1])-1,query_chr=chrchch$query_chr[i],query_start = as.numeric(chrchch$query_end[i+1]), query_end = as.numeric(chrchch$query_start[i]))
        }else{
          xx=data.frame(ref_chr=chrchch$ref_chr[i],ref_start = as.numeric(chrchch$ref_end[i])+1, ref_end = as.numeric(chrchch$ref_start[i+1])-1,query_chr=chrchch$query_chr[i],query_start = as.numeric(chrchch$query_end[i+1])+1, query_end = as.numeric(chrchch$query_start[i])-1)
        }
        colnames(xx)<-colnames(new_row)
        new_row <-rbind(new_row,xx)
        
      }
      
    }
    
  }
  

  assign(chr_child,new_row)
  }
  allchr.reverse<-do.call(rbind,mget(unique(endcluster1$query_chr)))
  pointer<-which(allchr.reverse$ref_start==1)[which(allchr.reverse$ref_start==1)!=1]
  allchr.reverse[pointer,]$ref_start<-allchr.reverse[pointer-1,]$ref_start
  allchr.reverse[pointer-1,]$ref_end<-allchr.reverse[pointer,]$ref_end
  rm(list=unique(endcluster1$query_chr)) ##删除变量
  

  if(cluster.id==2 |cluster.id==3){
  allchr.reverse<-allchr.reverse[-1,]
  }
  if(nrow(allchr.reverse)!=0){
    allchr.reverse$anno<-"SDR_NM"
  }
  if(cluster.id==0){
    return(list(reverse=allchr.reverse,endcluster1=endcluster1,minimaploc=reall,dup=dupli))
  }else{
    return(list(reverse=allchr.reverse,endcluster1=endcluster1))
  }
  
}




#-------------------function7 将重复区域进行整合 ----------------------
repeat.integrate<-function(data,repeatid){
  repeat_region<-c()
  data$query_start<-abs(data$query_start)
  data$query_end<-abs(data$query_end)
  data<-exchange(data)
  ir1ref <- IRanges(start = data$ref_start, end = data$ref_end)
  overlapsref<-findOverlaps(ir1ref,ir1ref) 
   
  if(length(which(overlapsref@from!=overlapsref@to))!=0){
    inte<-as.data.frame(overlapsref[which(overlapsref@from!=overlapsref@to)])
    if(exists("inte")){
      df_sorted <- t(apply(inte, 1, function(x) sort(x)))
      unique_df <- unique(df_sorted)
      inte <- as.data.frame(unique_df)
      for(j in 1:dim(inte)[1]){
        a<-c()
        a<-append(a,inte[j,1])
        a<-append(a,inte[j,2])
        for(k in 1:dim(inte)[1]){
          if(inte[k,1] %in% ((min(a)):(max(a)))){
            a<-unique(append(a,inte[k,1]))
          }
          if(inte[k,2] %in% ((min(a)):(max(a)))){
            a<-unique(append(a,inte[k,2]))
          }
        }
        x<-a
        if(abs(max(data[x,]$ref_start)-min(data[x,]$ref_end))>10000){
          repeat_region<-append(repeat_region,min(x))
          data[min(x),]$ref_start<-min(data[x,]$ref_start)
          data[min(x),]$ref_end<-max(data[x,]$ref_end)
          data[min(x),]$query_start<-min(data[x,]$query_start)
          data[min(x),]$query_end<-max(data[x,]$query_end)
          data[x,]<-data[min(x),]
        }
        if(abs(max(data[x,]$ref_start)-min(data[x,]$ref_end))<10000){
          data[min(x)+1,]$ref_start<- data[min(x),]$ref_end
        }
      }
    }
  }
  data<-distinct(data)
  ir1que  <- IRanges(start = data$query_start, end = data$query_end)
  overlapsque<-findOverlaps(ir1que,ir1que)
  if(length(which(overlapsque@from!=overlapsque@to))!=0){
    inte<-as.data.frame(overlapsque[which(overlapsque@from!=overlapsque@to)])
    if(exists("inte")){
      df_sorted <- t(apply(inte, 1, function(x) sort(x)))
      unique_df <- unique(df_sorted)
      inte <- as.data.frame(unique_df)
      for(j in 1:dim(inte)[1]){
        a<-c()
        a<-append(a,inte[j,1])
        a<-append(a,inte[j,2])
        for(k in 1:dim(inte)[1]){
          if(inte[k,1] %in% ((min(a)):(max(a)))){
            a<-unique(append(a,inte[k,1]))
          }
          if(inte[k,2] %in% ((min(a)):(max(a)))){
            a<-unique(append(a,inte[k,2]))
          }
        }
        x<-a
        if(abs(max(data[x,]$query_start)-min(data[x,]$query_end))>10000){
          repeat_region<-append(repeat_region,min(x))
          data[min(x),]$ref_start<-min(data[x,]$ref_start)
          data[min(x),]$ref_end<-max(data[x,]$ref_end)
          data[min(x),]$query_start<-min(data[x,]$query_start)
          data[min(x),]$query_end<-max(data[x,]$query_end)
          data[x,]<-data[min(x),]
        }
        if(abs(max(data[x,]$query_start)-min(data[x,]$query_end))<10000){
          if(repeatid==2){
            data[min(x)+1,]$query_end<- data[min(x),]$query_start
          }
          if(repeatid==1){
            data[min(x)+1,]$query_start<- data[min(x),]$query_end
          }
         
        }
      }
    }
  }
  
##处理一下query跨染色体的情况:
if(repeatid==1){
for(i in 2:dim(data)[1]){
  list<-which(data$query_start[i]<data[1:(i-1),]$query_end)
  if(i>dim(data)[1]){
    break
  }
  if(length(list)!=0){
    print(i)
    data[i,]$ref_start<-min(data[c(list,i),]$ref_start)
    data[i,]$ref_end<-max(data[c(list,i),]$ref_end)
    data[i,]$query_start<-min(data[c(list,i),]$query_start)
    data[i,]$query_end<-max(data[c(list,i),]$query_end)
    data<-data[-list,]
    ##没有改位置哦
  }
}}
dup<-data[unique(repeat_region),]
data<-data[order(data$ref_start),]
data<-distinct(data)
dup<-dup[c("ref_chr","ref_start","ref_end","query_chr","query_start","query_end","orient","cluster")]
  return(list(afterdup=data,repeat.region=dup))

}
#-------------------function7 整理最终的结果（去除大于10k的片段，去除端粒和着丝粒区域） ----------------------
endfilter<-function(all,chrid,chr_child){
  data<-all
  data$reflen<-data$ref_end-data$ref_start
  data$querylen<-data$query_end-data$query_start
  data<-data[data$reflen>10000 | data$querylen>10000,]
  if(length(which(data$reflen==0))!=0){
    data[data$reflen==0,]$anno<-"INS"
  }
  if(length(which(data$querylen==0))!=0){
    data[data$querylen==0,]$anno<-"DEL"
  }

  # ##去掉端粒
  # teloquery<-result_chantelo[result_chantelo$chr==chr_child,]
  # m<-teloquery[teloquery$id=="end",]
  # n<-teloquery[teloquery$id=="first",]
  # if(nrow(m)!=0){ ##说明是端粒末端
  #   data[nrow(data),6]<-m$query_start
  # }
  # if(nrow(n)!=0){
  #   data[1,5]<-n$query_end
  # }
  # teloref<-result_hmtelo[result_hmtelo$chr==chrid,]
  # m<-teloref[teloref$id=="end",]
  # n<-teloref[teloref$id=="first",]
  # if(nrow(m)!=0){ ##说明是端粒末端
  #   data[nrow(data),3]<-m$ref_start
  # }
  # if(nrow(n)!=0){
  #   data[1,2]<-n$ref_end
  # }
  # ####去掉着丝粒
  # centroref<-result_hmcentr[result_hmcentr$chr==chrid,]
  # centroquery<-result_chancentr[result_chancentr$chr==chr_child,]
  # #data<-data[!((data$ref_start %in% (centroref$ref_start:centroref$ref_end)) & (data$ref_end %in% (centroref$ref_start:centroref$ref_end))),]
  # #data<-data[!((data$query_start %in% (centroquery$query_start:centroquery$query_end)) & (data$query_end %in% (centroquery$query_start:centroquery$query_end))),]
  # ir1ref <- IRanges(start = data$ref_start, end = data$ref_end)
  # ir2 <- IRanges(start =centroref$ref_start , end =centroref$ref_end)
  # overlapsref<-findOverlaps(ir1ref,ir2) ##我们的比对序列的参考基因组
  # ir1que  <- IRanges(start = data$query_start, end = data$query_end)
  # ir3 <- IRanges(start =centroquery$query_start , end =centroquery$query_end)
  # overlapsquery<-findOverlaps(ir1que,ir3) ##我们的比对序列的参考基因组
  # ref<-overlapsref@from
  # que<-overlapsquery@from
  # print(overlapsquery)
  # data$note<-"a"
  # ## 分为： 同一个 SDR的ref和query都与着丝粒有交集
  # if(length(intersect(ref,que))!=0){
  #   for(l in intersect(ref,que)){
  #     ## 然后继续分，先看参考基因组
  #     ## 1.如果位于着丝粒区域内
  #     if(data[l,]$ref_start> centroref$ref_start & data[l,]$ref_end<centroref$ref_end)
  #     {data$note[l]<-'del'}
  #     ## 2.如果位于着丝粒区域外
  #     if(data[l,]$ref_start< centroref$ref_start & data[l,]$ref_end>centroref$ref_end)
  #     {data$note[l]<-'del'
  #     data<-rbind(data,data[l,])
  #     data<-rbind(data,data[l,])
  #     data[dim(data)[1]-1,]$ref_end<-centroref$ref_start-1
  #     data[dim(data)[1],]$ref_start<-centroref$ref_end+1
  #     if(data[l,]$query_start< centroquery$query_start & data[l,]$query_end>centroquery$query_end){
  #       ## ref和query都位于着丝粒区域外，包括了着丝粒区域
  #       data[dim(data)[1]-1,]$query_end<-centroquery$query_start-1
  #       data[dim(data)[1]-1,]$note<-"a"
  #       data[dim(data)[1],]$query_start<-centroquery$query_end+1
  #       data[dim(data)[1],]$note<-"a"
  #     }
  #     }
  #   }
  # }
  # ## 分为： ref与着丝粒有交集
  # setdiff(ref,que)
  # ## 分为： query与着丝粒有交集
  # setdiff(que,ref)
  # 
  # if(length(which(data$note=="del"))!=0){
  #   data<-data[-which(data$note=="del"),]
  # }
  # 
  # data<-data[order(data$ref_start),]
  # data<-data[,-10]
  # 

  return(data)
}

insertsmall<-function(endcluster2before,storesmall,orientid){
  orientlist<-rle(endcluster2before[["orient"]] == orientid)
  if(orientid=="+"){ #说明的大inversion
    reverseid=2
  }
  else{
    reverseid=3
  }
  for (i in which(orientlist$values)){
    if(i!=1){
      start<-(sum(orientlist$lengths[1:(i-1)])+1)
      end<-(sum(orientlist$lengths[1:(i-1)]))+orientlist$lengths[i]
      insertpos=endcluster2before[start:end,]
      reverse_end<-reverse.region(insertpos,chrid,reverseid)
      if(nrow(reverse_end$reverse)!=0){
        storesmall<-rbind(storesmall,reverse_end$reverse[,colnames(storesmall)])
        print(storesmall)
      }
     
    }
    else{
      insertpos=data[1:orientlist$lengths[i],]
      reverse_end<-reverse.region(insertpos,chrid,reverseid)
      if(nrow(reverse_end$reverse)!=0){
        storesmall<-rbind(storesmall,reverse_end$reverse[,colnames(storesmall)])
        print(storesmall)
      }
      
    }
  }
  return(storesmall)
}


smalltransf<-function(endcluster0,transpara){
  x=table(endcluster0$cluster)
  dellist<-c()
  for(k in names(x)[x<=transpara]){
    if(k==1){
      next
    }
    for(j in names(x)[x>transpara]){
      a<-which(endcluster0$cluster==k)
      b<-which(endcluster0$cluster==j)
      if(min(b)<min(a) & max(b)>max(a)){
        dellist<-append(dellist,k)
      }
    }
  }
  if(length(dellist)!=0){
    tran<-endcluster0[endcluster0$cluster %in% dellist,]
    tranbefore<-tran
    tran<-tran %>% group_by(cluster)%>%  #聚cluster
      summarise(ref_chr=ref_chr,ref_start=min(ref_start),
                ref_end=max(ref_end),
                query_chr=query_chr,
                query_starttem=min(pmin(query_end,query_start)),
                query_endtem=max(pmax(query_start,query_end)),
                (names(which.max(table(orient))))) ###正链负链？？重新算一下
    colnames(tran)[6]<-"query_start"
    colnames(tran)[7]<-"query_end"
    endcluster0<-endcluster0[!endcluster0$cluster %in% dellist,]
    colnames(tran)[8]<-"orient"
    return(list(tran=distinct(tran),endcluster0=endcluster0,tranbefore=tranbefore))
  }
  else{
    return(list(tran=NULL,endcluster0=endcluster0))
  }
  
}
cluster<-function(data){
  data<-data %>% group_by(cluster)%>%  #聚cluster
    summarise(ref_chr=ref_chr,ref_start=min(ref_start),
              ref_end=max(ref_end),
              query_chr=query_chr,
              query_starttem=min(pmin(query_end,query_start)),
              query_endtem=max(pmax(query_start,query_end)),
              (names(which.max(table(orient))))) ###正链负链？？重新算一下
  colnames(data)[6]<-"query_start"
  colnames(data)[7]<-"query_end"
  colnames(data)[8]<-"orient"
  data<-distinct(data)
  return(data)
}
## 这个函数用于把cluster为-中相邻的+进行合并
smallcluster<-function(data,orientid){
  orientlist<-rle(endcluster2[["orient"]] == orientid)
  for (i in which(orientlist$values)){
    if(i!=1){
      start<-(sum(orientlist$lengths[1:(i-1)])+1)
      end<-(sum(orientlist$lengths[1:(i-1)]))+orientlist$lengths[i]
    }
    else{
      start<-1
      end<-orientlist$lengths[i]
    }
    data[start,]$ref_start<-min(data[start:end,]$ref_start)
    data[start,]$ref_end<-max(data[start:end,]$ref_end)
    data[start,]$query_start<-min(data[start:end,]$query_start)
    data[start,]$query_end<-max(data[start:end,]$query_end)
    data[start,]$ref_pos<-(data[start,]$ref_start+data[start,]$ref_end)/2
    data[start,]$query_pos<-(data[start,]$query_start+data[start,]$query_end)/2
    data[start:end,]<-data[start,]
  }
  return(distinct(data))
}
