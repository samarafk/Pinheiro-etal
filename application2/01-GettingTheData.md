
Getting the Dataset
========================================================

We need the list of all accession numbers from [J. Virol. October 1998 vol. 72 no. 10 8240-8251](http://jvi.asm.org/content/72/10/8240.short). Then, for each accession number we need the corresponding sequence and annotation.


```r
if (!file.exists("./data/Acc.csv") | !file.exists("./data/Seq.rda") | !file.exists("./data/InfoSeq.rda"))
  {
  library("seqinr")
  choosebank("genbank")
  tmp <- query("allsequences", "AU=Poss AND Y=1998") #707
  t <- sapply(tmp$req, getName, as.string = TRUE)
  write.csv(t,file="./data/Acc.csv",row.names = FALSE)
  Seq <- getSequence(tmp$req)
  InfoSeq <- getAnnot(tmp$req,"features")
  closebank()
  save(Seq,file="./data/Seq.rda")
  save(InfoSeq,file="./data/InfoSeq.rda")
  }
```




There is a total of 707 sequences.

We created a dataframe, `DataInfo`, with the following columns: Patient, Isolate, Clone, Region, Accession, SampleNumber, Location, Date.




We will work only with sequences from Region V3.


```r
if (!file.exists("./data/V3.rda"))

  {
  load("./data/DataInfo.rda")
  library(Biostrings)
  library(DECIPHER)
  tmp <- which(DataInfo[,"Region"]=="V3")
  V3seq <- sapply(tmp,function(x) paste0(Seq[[x]],collapse=""))

  #BrowseSequences(DNAStringSet(V3seq))
  alignedV3 <- AlignSeqs(DNAStringSet(V3seq),verbose = FALSE)
  #BrowseSequences(alignedV3)

  # each row is a sequence
  SeqV3 <- t(sapply(1:length(alignedV3),function(x) unlist(strsplit(as.character(alignedV3[[x]]),split=character(0)))))

  SeqV3[SeqV3=="-"] <- NA

  DataInfoV3 <- DataInfo[DataInfo[,"Region"]=="V3",]

  save(SeqV3,DataInfoV3,file="./data/V3.rda")
}
```

We consider only dates with sequences available for all locations (Plasma, PBMC and Cervical).


```r
load("./data/V3.rda")
bb <- as.data.frame(xtabs(~Location+Date, data=DataInfoV3,drop.unused.levels = TRUE))
cc <- unique(bb[which(bb$Freq==0),"Date"])

DataInfoV3Final <- DataInfoV3[-which(DataInfoV3$Date %in% cc),]
SeqV3Final <- SeqV3[-which(DataInfoV3$Date %in% cc),]

N <- dim(SeqV3Final)[1]

DataInfoNew <- DataInfoV3Final
AllSeqNew <- SeqV3Final

DataInfoNew$Group <- DataInfoNew$Location
DataInfoNew$Visit <- DataInfoNew$Date

save(DataInfoNew,file="./data/DataInfoNew.rda")
save(AllSeqNew,file="./data/AllSeqNew.rda")
```


The dataset used in the analysis contains 249 sequences.


```r
xtabs(~Date+Location+Patient, data=DataInfoV3Final,drop.unused.levels = TRUE)
```

```
## , , Patient = Q23
## 
##                  Location
## Date              Cervical PBMC Plasma
##   August 1993"           0    0      0
##   February 1994"         0    0      0
##   February 1995"         6    9      7
##   July 1993"             9   10     12
##   July 1995"             0    0      0
##   June 1994"             6   11     11
##   November 1994"         9   11     12
##   October 1994"          0    0      0
##   September 1995"        8   11      9
## 
## , , Patient = Q45
## 
##                  Location
## Date              Cervical PBMC Plasma
##   August 1993"           5    3      8
##   February 1994"         0    0      0
##   February 1995"         0    0      0
##   July 1993"             0    0      0
##   July 1995"             0    0      0
##   June 1994"             0    0      0
##   November 1994"         0    0      0
##   October 1994"          3    9     10
##   September 1995"        0    0      0
## 
## , , Patient = Q47
## 
##                  Location
## Date              Cervical PBMC Plasma
##   August 1993"           0    0      0
##   February 1994"         8   11      7
##   February 1995"         0    0      0
##   July 1993"             0    0      0
##   July 1995"             5   10      5
##   June 1994"             0    0      0
##   November 1994"        10   10      4
##   October 1994"          0    0      0
##   September 1995"        0    0      0
```


The following table presents the number of sequences per patient, date, location and region.


|Patient |Date            |Location | Freq|
|:-------|:---------------|:--------|----:|
|Q23     |February 1995"  |Cervical |    6|
|Q23     |July 1993"      |Cervical |    9|
|Q23     |June 1994"      |Cervical |    6|
|Q23     |November 1994"  |Cervical |    9|
|Q23     |September 1995" |Cervical |    8|
|Q23     |February 1995"  |PBMC     |    9|
|Q23     |July 1993"      |PBMC     |   10|
|Q23     |June 1994"      |PBMC     |   11|
|Q23     |November 1994"  |PBMC     |   11|
|Q23     |September 1995" |PBMC     |   11|
|Q23     |February 1995"  |Plasma   |    7|
|Q23     |July 1993"      |Plasma   |   12|
|Q23     |June 1994"      |Plasma   |   11|
|Q23     |November 1994"  |Plasma   |   12|
|Q23     |September 1995" |Plasma   |    9|
|Q45     |August 1993"    |Cervical |    5|
|Q45     |October 1994"   |Cervical |    3|
|Q45     |August 1993"    |PBMC     |    3|
|Q45     |October 1994"   |PBMC     |    9|
|Q45     |August 1993"    |Plasma   |    8|
|Q45     |October 1994"   |Plasma   |   10|
|Q47     |February 1994"  |Cervical |    8|
|Q47     |July 1995"      |Cervical |    5|
|Q47     |November 1994"  |Cervical |   10|
|Q47     |February 1994"  |PBMC     |   11|
|Q47     |July 1995"      |PBMC     |   10|
|Q47     |November 1994"  |PBMC     |   10|
|Q47     |February 1994"  |Plasma   |    7|
|Q47     |July 1995"      |Plasma   |    5|
|Q47     |November 1994"  |Plasma   |    4|


The next table presents the number of sequences per patient, location and size (number of nucleotides in the sequence):


|Patient |Location | TamanhoSeq| Freq|
|:-------|:--------|----------:|----:|
|Q23     |Cervical |        105|   38|
|Q23     |PBMC     |        105|   52|
|Q23     |Plasma   |        105|   51|
|Q45     |Cervical |        102|    8|
|Q45     |PBMC     |        102|   12|
|Q45     |Plasma   |        102|   18|
|Q47     |Cervical |        105|   23|
|Q47     |PBMC     |        105|   31|
|Q47     |Plasma   |        105|   16|
