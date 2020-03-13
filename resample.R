##之前因为需要自己写了一个resample的函数。
#因为传统实现resample的方法好像没有做迭代，只会重抽一次。这就导致了每次重抽会有一些差别。于是我加入了迭代。
#懒得写成独立的函数了，就这样放出来，可以看到我每一步的想法。
#思路是对于每个样本，先将每一个OTU和其对应的序列数相乘，从这个结果中进行重抽，并加入迭代。最后把迭代结果取平均并取整，即为该样本最终结果。

#虽然用了几种方法提高速度：并行；提前建好最后的数据框；利用foreach；每次循环清空内存。
#但是本身方法比较笨，算得特别慢，加入迭代之后就更慢了。不推荐平时使用。但是需要迭代的时候可以试试。

rm(list=ls());gc()

##设置并行，提高速度
library(doParallel) 
core <- makeCluster(6)  ##设置集群数。默认为3
registerDoParallel(core)

###读入OTU
otu = read.table(file="OTU.txt",sep="\t",header=T,row.names=1)

ptm = proc.time()

#提前建立空数据框，提高速度
total = as.data.frame(array(NA ,dim=c(length(rownames(otu)),length(colnames(otu))))) 

#使用foreach方法，提高速度
library(foreach)
total = foreach(i=1:length(colnames(otu)),.combine = "cbind") %dopar% {
  read.num = sum(as.numeric(otu[,i])) #序列数
  count = c()
  for (j in 1:length(rownames(otu))){ #计算reads数
    times_read = as.numeric(re[j,i])
    reppp = rep(rownames(re)[j],times_read)
    count = c(count,reppp) }

  #设置resample数，这里为10000，可根据要求自己修改
  k = 10000
  total.summary = c()
  for (m in 1:1000){                #迭代1000次，无放回抽样并统计物种，可自己修改
    sample_read = sample(count,k,replace=FALSE)
    summary = as.data.frame(table(sample_read))  
    diff = setdiff(rownames(re),summary[,1])  ##补集
    diff_otu = cbind(sample_read = as.data.frame(diff),Freq = rep(0,length(diff))) ##添加0
    each_otu = as.data.frame(t (cbind(t(summary),t(diff_otu))))
    each_otu_sort = each_otu[order(each_otu[,1]),]
    total.summary  = as.data.frame(cbind(total.summary,as.vector(each_otu_sort[,2])))
  }  
  colnames(total.summary)=1:1000
  total.summary = as.data.frame(t(total.summary)) ##注意转置
  f<-function(x){mean(as.numeric(as.vector(x)))}   
  mean_read = apply(total.summary,2,f)      #对1000次迭代的结果取平均后再取整作为最后结果
  round(mean_read)
  #gc()
}
rownames(total) = rownames(otu)[ order(rownames(otu)) ]
colnames(total) = colnames(otu)

proc.time() - ptm

total

# 关闭集群
stopCluster(cl)

write.table(total,file="resample_10000.txt",col.names = NA,sep="\t")
