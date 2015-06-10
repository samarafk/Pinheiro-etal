

Getting the Dataset
========================================================

Code to reproduce the analyses for the first example (application 1) in the manuscript is available in this repository.

The dataset contains multiple-time points HIV sequences from the CAPRISA 002 Study Team [(PLoS Pathog. 2008 Mar 21;4(3):e1000033)](http://www.ncbi.nlm.nih.gov/pubmed?term=18369479).

To download the data, go to http://www.hiv.lanl.gov/components/sequence/HIV/search/search.html and fill in the PUBMEDID (18369479) in __Publication Information__ menu. In the __Output__ menu, select to list 400 records per page. Then click on __Search__. The next screen lists all the sequences. Click on __Select All__ and then on __Download Sequences__.

In __Download Options__, select __Align__, __Format__ Table, __Gap__ none, __Sequence type__ nucleotides.

In __Compose sequence labels__, select: 1 - Accession, 2 - Patient Code, 3 - Sampling Year, 4 - Subtype, 5 - Name, 6 - HXB2/MAC239 start, 7 - HXB2/MAC239 stop, 8 - Sequence Length. Click __Ok__ e then __Download__.

The file is available for download in this [link](http://www.hiv.lanl.gov/cgi-bin/common_code/download.cgi?/tmp/download/f57248d8/hiv-db.table). We provide a copy of this file in the subfolder __/data__.

Also, by clicking on each Patients's id, we have more information, such as the progression classification. The patients are classified as __Slow Progressors__ (SP), __Rapid Progressors__ (RP) or __Normal Progressors__ (P).


Reading the dataset:


```r
data <- read.delim("./data/hiv-db.table",as.is=TRUE,header=FALSE)

InfoSeq <- data[,1]
Seq <- data[,2]

save(Seq,file="./data/Seq.rda")
```

Creating a data set with all the information from the sequences:


```r
Accession <- sapply(InfoSeq,function(x) strsplit(x, "\\.")[[1]][1])
names(Accession) <- NULL

Patient <- sapply(InfoSeq,function(x) strsplit(x, "\\.")[[1]][2])
names(Patient) <- NULL

Date <- sapply(InfoSeq,function(x) strsplit(x, "\\.")[[1]][3])
names(Date) <- NULL

Isolate <- sapply(InfoSeq,function(x) strsplit(x, "\\.")[[1]][5])
names(Isolate) <- NULL

Start <- as.numeric(sapply(InfoSeq,function(x) strsplit(x, "\\.")[[1]][6]))
names(Start) <- NULL

End <- as.numeric(sapply(InfoSeq,function(x) strsplit(x, "\\.")[[1]][7]))
names(End) <- NULL

Length <- as.numeric(sapply(InfoSeq,function(x) strsplit(x, "\\.")[[1]][8]))
names(Length) <- NULL

DataInfo <- data.frame(Accession=Accession,Patient=Patient,Date=Date,
                       Isolate=Isolate,Start=Start,End=End,Length=Length,stringsAsFactors=FALSE)
```

The individuals are classified as __Slow Progressors__ (SP), __Rapid Progressors__ (RP) or __Normal Progressors__ (P). The following code includes the group information into the dataset:


```r
DataInfo$Group=NA
DataInfo$Group[DataInfo$Patient%in% c("CAP45","CAP61","CAP228")] <- "SP"
DataInfo$Group[DataInfo$Patient %in% c("CAP85","CAP88","CAP225","CAP255","CAP30","CAP244","CAP257")] <- "P"
DataInfo$Group[DataInfo$Patient%in% c("CAP8","CAP69","CAP174","CAP256","CAP258","CAP65")] <- "RP"
DataInfo$Group[DataInfo$Patient %in% c("CAP248","CAP200","CAP84","CAP206","CAP210")] <- NA
save(DataInfo,file="./data/DataInfo.rda")

# http://uctscholar.uct.ac.za/PDF/92521_Ntale_R.pdf
# http://jvi.asm.org/content/83/1/470.full.pdf
# table 3
```

Checking alingment:

```r
library(DECIPHER)
BrowseSequences(DNAStringSet(Seq))
```

Saving the sequences in matrix form (each column is a site, each row is a sequence) and checking inconsistencies:


```r
Seqq <- list()

for (i in 1:length(Seq))
    {
        Seqq[[i]] <- strsplit(Seq[i],character(0))
    }

tmp <- sapply(Seqq,function(x) length(x[[1]]))

AllSeq <- matrix(NA,nrow=length(Seq),ncol=max(tmp))
for (i in 1:length(Seq))
  {
      tmp1 <- length(Seqq[[i]][[1]])
      AllSeq[i,1:tmp1] <- Seqq[[i]][[1]]
  }

AllSeq[AllSeq=="-"] <- NA
AllSeq[AllSeq=="K"] <- NA
AllSeq[AllSeq=="S"] <- NA
AllSeq[AllSeq=="t"] <- "T"
AllSeq[AllSeq=="a"] <- "A"

save(AllSeq,file="./data/AllSeq.rda")
```


We exclude patients without group information. We keep only patients that have information for at least 2 visits. We will consider only the first and second visits.


```r
# Exclude patients without group information
index <- !is.na(DataInfo$Group)
DataInfo <- DataInfo[index,]
AllSeq <- AllSeq[index,]

# Keep only patients that have info for at least 2 time points

tmp <- aggregate(DataInfo$Accession, by=list(DataInfo$Patient,DataInfo$Date),function(x) length(x))
tmp1 <- as.data.frame(table(tmp[,1]))
pat <- tmp1[tmp1[,2]==1,1]

DataInfoNew <- DataInfo[-which(DataInfo$Patient %in% pat),]
AllSeqNew <- AllSeq[-which(DataInfo$Patient %in% pat),]

# Exclude info from third visit, if available

tmp <- aggregate(DataInfoNew$Accession, by=list(DataInfoNew$Patient,DataInfoNew$Date),function(x) length(x))
tmp1 <- as.data.frame(table(tmp[,1]))

pat <- tmp1[tmp1[,2]==3,1]

patdate <- sapply(pat,function(x) sort(as.numeric(as.character(unique(DataInfoNew[DataInfoNew$Patient==x,"Date"]))),decreasing=TRUE)[1])

tmp2 <- data.frame(as.character(pat),as.character(patdate))

tmp3 <- list()
for (i in 1:dim(tmp2)[1])
    {
        tmp3[[i]] <- which(DataInfoNew$Patient==tmp2[i,1] & DataInfoNew$Date==tmp2[i,2])
    }

tmp4 <- unlist(tmp3)

DataInfoNew <- DataInfoNew[-tmp4,]
AllSeqNew <- AllSeqNew[-tmp4,]

### Defining visits 1 and 2
DataInfoNew$Visit <- NA
pat <- unique(DataInfoNew$Patient)

for (i in 1:length(pat))
    {
        Datatmp <- DataInfoNew[DataInfoNew$Patient==pat[i],]
        tmp <- sort(unique(Datatmp$Date))
        DataInfoNew$Visit[DataInfoNew$Patient==pat[i] & DataInfoNew$Date==tmp[1]] <- 1
        DataInfoNew$Visit[DataInfoNew$Patient==pat[i] & DataInfoNew$Date==tmp[2]] <- 2
    }

index <- order(DataInfoNew$Group,DataInfoNew$Patient,DataInfoNew$Visit)

DataInfoNew <- DataInfoNew[index,]
AllSeqNew <- AllSeqNew[index,]

save(DataInfoNew,file="./data/DataInfoNew.rda")
save(AllSeqNew,file="./data/AllSeqNew.rda")
```

There is a total of 224 sequences from 13 patients (3 SP, 4 RP and 6 P).

Table 4 of the manuscript: Distribution of sequences per individual:

|Group |Date |Patient | Freq|
|:-----|:----|:-------|----:|
|P     |2005 |CAP225  |   18|
|P     |2006 |CAP225  |   32|
|P     |2005 |CAP255  |    3|
|P     |2006 |CAP255  |    1|
|P     |2005 |CAP257  |    2|
|P     |2006 |CAP257  |    1|
|P     |2004 |CAP30   |    2|
|P     |2005 |CAP30   |    2|
|P     |2005 |CAP85   |   31|
|P     |2006 |CAP85   |   13|
|P     |2005 |CAP88   |   23|
|P     |2006 |CAP88   |   25|
|RP    |2005 |CAP174  |    2|
|RP    |2006 |CAP174  |    1|
|RP    |2005 |CAP256  |    2|
|RP    |2006 |CAP256  |    1|
|RP    |2005 |CAP65   |    3|
|RP    |2006 |CAP65   |    1|
|RP    |2005 |CAP69   |    1|
|RP    |2006 |CAP69   |    2|
|SP    |2005 |CAP228  |   20|
|SP    |2006 |CAP228  |    1|
|SP    |2005 |CAP45   |    2|
|SP    |2006 |CAP45   |    1|
|SP    |2004 |CAP61   |   20|
|SP    |2005 |CAP61   |   14|

