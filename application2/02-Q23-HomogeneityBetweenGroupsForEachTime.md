

Homogeneity between groups for each time - Q23
========================================================

For *Woman Q23*, we may test homogeneity of groups for each time  (dates: 7/93, 6/94, 11/94, 2/95, 9/95), i.e., we compute $BG_t, t=1, \ldots, 5$. Then, we can test if the homogeneity of groups is the same for all times using  $\mathbf{T}_n=\mathbf{C}\mathbf{B}_n = (T_{n1}, T_{n2}, T_{n3}, T_{n4})$, with $\mathbf{B}_n = (BG_1, BG_2, BG_3, BG_4, BG_5)$ and 
\[ \mathbf{C} = \left ( \begin{array}{ccccc}
-1 & 1 & 0 & 0 & 0\\
-1 & 0 & 1 & 0 & 0 \\
-1 & 0 & 0 & 1 & 0 \\
-1 & 0 & 0 & 0 & 1  
\end{array} 
\right ) . 
\]


If one wants to test the group effect taking out the time effect, the structure of times should be preserved. So, for each $b=1,\ldots,B$: 

 1. sample with replacement (w.r.) from the $m_{i\cdot g}$ sequences (i.e., the total number of sequences for individuals $i$ at time $t$) ignoring the group structure. For each time, the first $m_{i1t}$ sequences will be from group 1, the second $m_{i2t}$ from group 2 and so on.
2. compute the vector $\mathbf{B}_{nb} = (BG_{1_b}, BG_{2_b}, \ldots , BG_{Tb})$ for each $b=1, \ldots B$.
3. standardize each $BG_{t_b}$ ($t=1, \ldots , T$), i.e., subtract the mean and standard deviation of the $B$ values of $BG_{t_b}'s$. 
4. compute the appropriate $\mathbf{T}_{nb}=\mathbf{C}\mathbf{B}_{nb}$.

The code for this procedure is:




```r
if (!file.exists("./data/DataBootQ23.rda"))
  {
  load("./data/AllSeqNew.rda")
  load("./data/DataInfoNew.rda")
  
  DataInfoNewQ23 <- DataInfoNew[DataInfoNew$Patient=="Q23",]
  AllSeqNewQ23 <- AllSeqNew[DataInfoNew$Patient=="Q23",]
  
  save(DataInfoNewQ23,file="./data/DataInfoNewQ23.rda")
  save(AllSeqNewQ23,file="./data/AllSeqNewQ23.rda")


  DataBoot <- DataInfoNewQ23[,c("Group","Patient","Visit")]
  tmp <- unique(DataBoot[,c("Patient","Group","Visit")])
  B=2000
  IndexBoot <- matrix(NA,ncol=B,nrow=dim(DataInfoNewQ23)[1])
  set.seed(2014)
  for (b in 1:B)
      {
          tmpIndexBoot <- c()
          for (i in 1:dim(tmp)[1])
              {
                  index <- which(DataInfoNewQ23$Visit==tmp[i,"Visit"] &
                                 DataInfoNewQ23$Patient==tmp[i,"Patient"])
                  if (length(index)>1)
                      {          
                          tmpIndexBoot <- sample(index,replace=TRUE)
                      }
                  else
                      {
                          tmpIndexBoot <- index
                      }
                  IndexBoot[index,b] <- tmpIndexBoot
              }
      }
  DataBoot$IndexOriginal <- 1:dim(DataBoot)[1]
  DataBoot <- cbind(DataBoot,IndexBoot)
  save(DataBoot,file="./data/DataBootQ23.rda")
  }
```

After generating all boostrap indexes, we need to compute the statistics for each bootstrap sample.


```r
if (!file.exists("./data/RbootQ23.rda")) {
  load("./data/AllSeqNewQ23.rda")
  load("./data/DataInfoNewQ23.rda")
  load("./data/DataBootQ23.rda")
  library(ape)
  
  
  IndexOrig=DataBoot[,4]
  SequenceMatrix <- AllSeqNewQ23[IndexOrig,]
  
  AllDistances <- as.matrix(dist.gene(SequenceMatrix,method="percentage",pairwise.deletion=TRUE))
  
  source("./function/distances.R")
  
  Rboot <- list()
  B <- dim(DataBoot)[2]-3
  # checking if the group/visit label are in the same order
  DataBoot$Group==DataInfoNewQ23$Group # should be true
  DataBoot$Visit==DataInfoNewQ23$Visit # should be true
  
  for (b in 1:B) {
    Rboot[[b]] <- distances(AllDistances=AllDistances,Group=DataBoot$Group,Time=DataBoot$Visit,Index=DataBoot[,b+3])
    
  }
  
  save(Rboot,file="./data/RbootQ23.rda")
}
```

We now calculate $\mathbf{B}_n=(BG_1,BG_2,BG_3,BG_4,BG_5)$ for the original sample and the bootstrap samples.

$BG_t=\sum_{g<g'}\frac{m_{.gt}m_{.g't}}{N(N-1)}(2\bar{R}_{gt;g't}-\bar{R}_{gt;gt}-\bar{R}_{g't;g't})$, $t=1,\ldots, 5$ for testing the homogeneity of groups (Plasma, PBMC and Cervical) at each time.


```r
load("./data/RbootQ23.rda")
load("./data/DataBootQ23.rda")

tmp <- as.data.frame(table(DataBoot$Group))
minSeq <- which(tmp$Freq==min(tmp$Freq))
n <- tmp$Freq[minSeq]

if(! file.exists("./data/BGt23.rda"))
  {
  groups <- unique(DataBoot$Group)
  time <- unique(DataBoot$Visit)

  g <- expand.grid(groups,time) # always use groups and time in this order (keep same order as used in distance.R
  colnames(g) <- c("groups","time")

  BGt <- matrix(NA,ncol=length(time),nrow=length(Rboot))

  library(combinat)

  for (b in 1:length(Rboot))
      {
          N <- dim(DataBoot)[1] # total number of sequences
          Migt <- Rboot[[b]]$Migt # m.gt
          R <- Rboot[[b]]$R
          GroupTime <- Rboot[[b]]$GroupTime

          for (j in 1:length(time))
              {
                  d <- c()
                  tmp <- as.character(g[g[,2]==time[j],1])
                  aa <- combn(tmp,2) # all pairwise combination of the groups

                  d <- apply(aa,2,function (x)
                                        {
                                            tmp1 <- which(Migt[,2]==time[j] & Migt[,1]==x[1])
                                            tmp2 <- which(Migt[,2]==time[j] & Migt[,1]==x[2])
                                            
                                            group1time <- paste(x[1],time[j],sep="|")
                                            group1timegroup1time <- paste(group1time,group1time,sep="_")

                                            group2time <- paste(x[2],time[j],sep="|")
                                            group2timegroup2time <- paste(group2time,group2time,sep="_")

                                            group1timegroup2timeA <- paste(group1time,group2time,sep="_")
                                            group1timegroupp2timeB <- paste(group2time,group1time,sep="_")
                                            
                                            tmp3 <- which(GroupTime %in% c(group1timegroup2timeA,group1timegroupp2timeB),arr.ind=TRUE)
                                            tmp4 <- which(GroupTime %in% group1timegroup1time,arr.ind=TRUE)
                                            tmp5 <- which(GroupTime %in% group2timegroup2time,arr.ind=TRUE)
                                            
                                            rr <- Migt[tmp1,3]*Migt[tmp2,3]*(2*R[tmp3]-R[tmp4]-R[tmp5])
                                            return(rr)
                                        }
                             )
                 BGt[b,j] <- sum(d)/(N*(N-1))
              }
      }
  colnames(BGt) <- time

  save(BGt,file="./data/BGtQ23.rda")
}
```




Empirical distribution of $BG_t$: 


![](02-Q23-HomogeneityBetweenGroupsForEachTime_files/figure-html/unnamed-chunk-5-1.png) ![](02-Q23-HomogeneityBetweenGroupsForEachTime_files/figure-html/unnamed-chunk-5-2.png) ![](02-Q23-HomogeneityBetweenGroupsForEachTime_files/figure-html/unnamed-chunk-5-3.png) 



Empirical distribution of $BG_t$: 

![](02-Q23-HomogeneityBetweenGroupsForEachTime_files/figure-html/unnamed-chunk-7-1.png) 






We then calculate $\mathbf{T}_n=\mathbf{C}\mathbf{B}_n$. First, we standardize each $BG_t$.


```r
Tn1 <- (BGt[,2]-mean(BGt[,2]))/sd(BGt[,2]) - (BGt[,1]-mean(BGt[,1]))/sd(BGt[,1])
Tn2 <- (BGt[,3]-mean(BGt[,3]))/sd(BGt[,3]) - (BGt[,1]-mean(BGt[,1]))/sd(BGt[,1])
Tn3 <- (BGt[,4]-mean(BGt[,4]))/sd(BGt[,4]) - (BGt[,1]-mean(BGt[,1]))/sd(BGt[,1])
Tn4 <- (BGt[,5]-mean(BGt[,5]))/sd(BGt[,5]) - (BGt[,1]-mean(BGt[,1]))/sd(BGt[,1])
```


Empirical distribution of $T_{n1}$-$T_{n4}$.


![](02-Q23-HomogeneityBetweenGroupsForEachTime_files/figure-html/unnamed-chunk-10-1.png) 



$T_{n1}$ p-value is 0.916.

$T_{n2}$ p-value is 0.938.

$T_{n3}$ p-value is 0.554.

$T_{n4}$ p-value is 0.756.



![](02-Q23-HomogeneityBetweenGroupsForEachTime_files/figure-html/unnamed-chunk-12-1.png) 








Next, we calculate $\cal{L}_n$:



```r
Tn <- cbind(Tn1,Tn2,Tn3,Tn4)
source("./function/Ln.R")
Ln1 <- Ln(Tn,n)  
pvLn <- round(length(which(Ln1[-1]>Ln1[1]))/(length(Ln1)-1),3)
```


The observed $\cal{L}_n$ is 1203.77.

Empirical distribution of $\cal{L}_n$:


![](02-Q23-HomogeneityBetweenGroupsForEachTime_files/figure-html/unnamed-chunk-16-1.png) 





