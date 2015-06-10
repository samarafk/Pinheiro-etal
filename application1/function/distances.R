### AllDistances: distance matrix for all pairs of sequences in the original dataset. The element [i,j] represents the distance between sequence i and sequence j. 

### Index indicates the rows from AllDistances that should be considered, default being 1:(Total number of sequences)

### Group: represents the group of the sequence

### Time: represents the time of the sequence


### This function returns a list with:
###    - R : symmetric matrix with elements \bar{R}_{gt}{g't'} 
###    - NumPair : symmetric matrix with elements being the total pairwise comparisons involved in the calculation of \bar{R}_{gt}{g't'}
###    - Migt : total number of sequences from group g time t


distances <- function(AllDistances=AllDistances,Group=DataBoot$Group,Time=DataBoot$Visit,Index=1:dim(AllDistances)[1])
    {
        require(dplyr)
        require(bitops)

        groups <- unique(Group)
        time <- unique(Time)
        g <- expand.grid(groups,time)

        AllDistances1 <- AllDistances[Index,Index]

        AllDistancesData <- as.data.frame( which( row(AllDistances1) <= col(AllDistances1) , arr.ind=TRUE) )
        
        AllDistancesData$ij <- paste(AllDistancesData$row,AllDistancesData$col,sep="-")

        AllDistancesData$GroupI <- Group[AllDistancesData[,1]]
        AllDistancesData$TimeI <- Time[AllDistancesData[,1]]
        AllDistancesData$GroupJ <- Group[AllDistancesData[,2]]
        AllDistancesData$TimeJ <- Time[AllDistancesData[,2]]

        AllDistancesData$Dist <- AllDistances1[which( row(AllDistances1) <= col(AllDistances1) , arr.ind=TRUE)]
        
        AllDistancesDataFinal <- AllDistancesData[-which(AllDistancesData$row==AllDistancesData$col),]

        AllDistancesData2 <- AllDistancesDataFinal %>%
                mutate(GroupTime1 = paste0(as.hexmode(as.integer(as.factor(GroupI)) +
                                          bitShiftL(TimeI, 2)),
                             as.hexmode(as.integer(as.factor(GroupJ)) +
                                          bitShiftL(TimeJ, 2)))) %>%
                  rowwise() %>%
                mutate(GroupTime2 = do.call(paste0, as.list(sort(unlist(strsplit(GroupTime1, split = NULL))))))

        MyGroups <- unique(as.factor(AllDistancesDataFinal$GroupI))

        AllDistancesAgg <- AllDistancesData2 %>%
           group_by(GroupTime2) %>%
           summarize(mean = mean(Dist), num = n()) %>%
           mutate(GTnum = strtoi(paste0("0x", GroupTime2))) %>%
           mutate(Group1 = MyGroups[bitwAnd(GTnum, 0x3)]) %>%
           mutate(Time1 = bitwAnd(bitShiftR(GTnum, 2), 0x3)) %>%
           mutate(Group2 = MyGroups[bitwAnd(bitShiftR(GTnum, 4), 0x3)]) %>%
           mutate(Time2 = bitShiftR(GTnum, 6))

      AllDistancesMat <- matrix(nrow = dim(g)[1], ncol = dim(g)[1])
      AllDistancesMat[lower.tri(AllDistancesMat, diag = TRUE)] <- AllDistancesAgg$mean
      AllDistancesMat <- t(AllDistancesMat)

      NumPair <- matrix(nrow = dim(g)[1], ncol = dim(g)[1])
      NumPair[lower.tri(NumPair, diag = TRUE)] <- AllDistancesAgg$num
      NumPair <- t(NumPair)

      Migt <- aggregate(DataInfoNew$Accession,by=list(DataInfoNew$Group,DataInfoNew$Visit),length)$x  #m.gt: total number of sequences from group g time t 
            
      return(list(R=AllDistancesMat,NumPair=NumPair,Migt=Migt))
    }

