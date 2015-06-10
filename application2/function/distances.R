### AllDistances: distance matrix for all pairs of sequences in the original dataset. The element [i,j] represents the distance between sequence i and sequence j.

### Index indicates the rows from AllDistances that should be considered, default being 1:(Total number of sequences)

### Group: represents the group of the sequence

### Time: represents the time of the sequence


### This function returns a list with:
###    - R : symmetric matrix with elements \bar{R}_{gt}{g't'}
###    - NumPair : symmetric matrix with elements being the total pairwise comparisons involved in the calculation of \bar{R}_{gt}{g't'}
###    - Migt : total number of sequences from group g time t


distances <- function(AllDistances = AllDistances,
                      Group = DataBoot$Group,
                      Time = DataBoot$Visit,
                      Index = 1:dim(AllDistances)[1]) {

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

  paste_ <- function(...) {
    paste(..., sep = "_")
  }

  AllDistancesData2 <- AllDistancesDataFinal %>%
    mutate(GroupTime1 = paste(paste(GroupI, TimeI, sep = "|"),
                              paste(GroupJ, TimeJ, sep = "|"), sep = "_")) %>%
    rowwise() %>%
    mutate(GroupTime2 = do.call(paste_,
                                as.list(sort(unlist(strsplit(GroupTime1, "_"))))))

  AllDistancesAgg <- AllDistancesData2 %>%
    group_by(GroupTime2) %>%
    summarize(mean = mean(Dist), num = n())

  #m.gt: total number of sequences from group g time t
  Migt <- aggregate(Index, by = list(Group, Time), length)

  AllDistancesMat <- matrix(nrow = dim(g)[1], ncol = dim(g)[1])
  AllDistancesMat[lower.tri(AllDistancesMat, diag = TRUE)] <- AllDistancesAgg$mean
  AllDistancesMat <- t(AllDistancesMat)

  NumPair <- matrix(nrow = dim(g)[1], ncol = dim(g)[1])
  NumPair[lower.tri(NumPair, diag = TRUE)] <- AllDistancesAgg$num
  NumPair <- t(NumPair)

  GroupTimeMat <- matrix(nrow = dim(g)[1], ncol = dim(g)[1])
  GroupTimeMat[lower.tri(GroupTimeMat, diag = TRUE)] <- AllDistancesAgg$GroupTime2
  GroupTimeMat <- t(GroupTimeMat)

  return(list(R = AllDistancesMat,
              NumPair = NumPair,
              Migt = Migt,
              GroupTime = GroupTimeMat))
}

