# install.packages("seqinr")
library("seqinr", lib.loc="~/R/win-library/3.2")
# preparation of data -----------------------------------------------------
train_plus  <-  read.delim("train_plus.data.txt", header = FALSE, sep = ">")
train_plus_seq  <-  as.vector(train_plus[[2]])

train_minus  <-  read.delim("train_minus.data.txt", header = FALSE, sep = ">")
train_minus_seq <- as.vector(train_minus[[2]])

test_plus <-  read.delim("test_plus.data.txt", header = FALSE, sep = ">")
test_plus_seq <- as.vector(test_plus[[2]])

test_minus  <-  read.delim("test_minus.data.txt", header = FALSE, sep = ">")
test_minus_seq <- as.vector(test_minus[[2]])


# seqAnalysis -------------------------------------------------------------
seqAnalyais <- function(sequence){
  x <- count(s2c(sequence), word = 3, start = 0, by = 3, freq = FALSE, alphabet = s2c("ACGT"))
  # startsite
  k <- 0
  if (x["ATG"] > 0) {
    while (substr(sequence, 3*k+1, 3*k+3)!="ATG" && 3*k+3<= nchar(sequence)) {
      k <- k+1
    }
    startsite <- 3*k+1
  } else {
    startsite <- -1;
  }
  
  m <- 0
  if (x["TAA"] >0 || x["TAG"]>0 ||x["TGA"]>0) {
    stopCodon <- 1;
    while (substr(sequence, 3*m+1, 3*m+3)!="TAA" && substr(sequence, 3*m+1, 3*m+3)!="TAG" 
           && substr(sequence, 3*m+1, 3*m+3)!="TGA" && 3*m+3<=nchar(sequence)) {
      m <- m+1
    }
    endsite <- 3*m+3
  } else{
    stopCodon <- 0
    endsite <- nchar(sequence)
  }
  # anaseq <- substr(sequence, startsite, endsite)
  result <- c(startsite, endsite, stopCodon)
  return(result)
}


# training minus model ----------------------------------------------------
mModel <- function(train_minus_seq){
  transp_minus <- matrix(0, 6, 6)
  
  for (i in 1:length(train_minus_seq)) {
    dinuCount <- count(s2c(train_minus_seq[i]),2, alphabet = s2c("ACGT"))
    
    transp_minus[2,2] <- transp_minus[2,2] + dinuCount["AA"]
    transp_minus[2,3] <- transp_minus[2,3] + dinuCount["AC"]
    transp_minus[2,4] <- transp_minus[2,4] + dinuCount["AG"]
    transp_minus[2,5] <- transp_minus[2,5] + dinuCount["AT"]
    
    transp_minus[3,2] <- transp_minus[3,2] + dinuCount["CA"]
    transp_minus[3,3] <- transp_minus[3,3] + dinuCount["CC"]
    transp_minus[3,4] <- transp_minus[3,4] + dinuCount["CG"]
    transp_minus[3,5] <- transp_minus[3,5] + dinuCount["CT"]
    
    transp_minus[4,2] <- transp_minus[4,2] + dinuCount["GA"]
    transp_minus[4,3] <- transp_minus[4,3] + dinuCount["GC"]
    transp_minus[4,4] <- transp_minus[4,4] + dinuCount["GG"]
    transp_minus[4,5] <- transp_minus[4,5] + dinuCount["GT"]
    
    transp_minus[5,2] <- transp_minus[5,2] + dinuCount["TA"]
    transp_minus[5,3] <- transp_minus[5,3] + dinuCount["TC"]
    transp_minus[5,4] <- transp_minus[5,4] + dinuCount["TG"]
    transp_minus[5,5] <- transp_minus[5,5] + dinuCount["TT"]
    
    transp_minus[1,s2n(substr(train_minus_seq[i],1,1))+2] <- 
      transp_minus[1,s2n(substr(train_minus_seq[1],1,1))+2]+1
    transp_minus[s2n(substr(train_minus_seq[i],nchar(train_minus_seq[i]),nchar(train_minus_seq[i])))+2,6] <- 
      transp_minus[s2n(substr(train_minus_seq[i],nchar(train_minus_seq[i]),nchar(train_minus_seq[i])))+2,6]+1
  }
  
  for (p in 1:6) {
    if (sum(transp_minus[p, ])>0) {
      transp_minus[p, ]=transp_minus[p, ]/sum(transp_minus[p, ])
    } else{
      transp_minus[p, ]=0;
    }
  }
  return(transp_minus)
}


# training plus model -----------------------------------------------------
pModel <- function(train_plus_seq){
  transp_plus <- matrix(0,21,21)
  
  for (i in 1:length(train_plus_seq)) {
    result <- seqAnalyais(train_plus_seq[i])
    
    if (result[1]>0 && result[2]>result[1] && result[3]==1) {
      codingRegion <- substr(train_plus_seq[i],result[1]+3, result[2]-3)
      y <- count(s2c(codingRegion),word = 3, start = 0, by = 3, freq = FALSE, alphabet = s2c("ACGT"))
    } else {
      codingRegion <- substr(train_plus_seq[i],result[1]+3, result[2])
      y <- count(s2c(codingRegion),word = 3, start = 0, by = 3, freq = FALSE,alphabet = s2c("ACGT"))
    }
    
    if (result[1]==1) {
      transp_plus[1,2] <- transp_plus[1,2]+1
      transp_plus[2,3] <- transp_plus[2,3]+1
      transp_plus[3,4] <- transp_plus[3,4]+1
      transp_plus[4,s2n(substr(train_plus_seq[i],4,4))+5]<- transp_plus[4,s2n(substr(train_plus_seq[i],4,4))+5]+1
    }
    
    n <- 0
    if (result[3]==0) {
      while (3*n+4<=nchar(codingRegion)) {
        transp_plus[s2n(substr(codingRegion,3*n+3,3*n+3))+15,s2n(substr(codingRegion,3*n+4,3*n+4))+5] <- transp_plus[s2n(substr(codingRegion,3*n+3,3*n+3))+15,s2n(substr(codingRegion,3*n+4,3*n+4))+5]+1
        n <- n+1
      }
      transp_plus[s2n(substr(codingRegion, nchar(codingRegion), nchar(codingRegion)))+15,21] <- transp_plus[s2n(substr(codingRegion, nchar(codingRegion), nchar(codingRegion)))+15,21]+1
      
    } else if (result[3]==1) {
      while (3*n+4<=nchar(codingRegion)){
        transp_plus[s2n(substr(codingRegion,3*n+3,3*n+3))+15,s2n(substr(codingRegion,3*n+4,3*n+4))+5] <- transp_plus[s2n(substr(codingRegion,3*n+3,3*n+3))+15,s2n(substr(codingRegion,3*n+4,3*n+4))+5]+1
        n <- n+1
      }
      transp_plus[s2n(substr(codingRegion, nchar(codingRegion), nchar(codingRegion)))+15,8] <- transp_plus[s2n(substr(codingRegion, nchar(codingRegion), nchar(codingRegion)))+15,8]+1
      
      if (substr(train_plus_seq[i],result[2]-2, result[2])=="TAA") {
        transp_plus[8,13] <- transp_plus[8,13]+1
        transp_plus[13,19] <- transp_plus[13,19]+1
        transp_plus[19,21] <- transp_plus[19,21]+1
      } else if (substr(train_plus_seq[i],result[2]-2, result[2])=="TAG") {
        transp_plus[8,13] <- transp_plus[8,13]+1
        transp_plus[13,20] <- transp_plus[13,20]+1
        transp_plus[20,21] <- transp_plus[20,21]+1
      } else if (substr(train_plus_seq[i],result[2]-2, result[2])=="TGA") {
        transp_plus[8,14] <- transp_plus[8,14]+1
        transp_plus[14,19] <- transp_plus[14,19]+1
        transp_plus[19,21] <- transp_plus[19,21]+1
      }
    }
    
    transp_plus[5,9] <- y["AAA"]+y["AAC"]+y["AAG"]+y["AAT"]+transp_plus[5,9]
    transp_plus[5,10] <- y["ACA"]+y["ACC"]+y["ACG"]+y["ACT"]+transp_plus[5,10]
    transp_plus[5,11] <- y["AGA"]+y["AGC"]+y["AGG"]+y["AGT"]+transp_plus[5,11]
    transp_plus[5,12] <- y["ATA"]+y["ATC"]+y["ATG"]+y["ATT"]+transp_plus[5,12]
    
    transp_plus[6,9] <- y["CAA"]+y["CAC"]+y["CAG"]+y["CAT"]+transp_plus[6,9]
    transp_plus[6,10] <- y["CCA"]+y["CCC"]+y["CCG"]+y["CCT"]+transp_plus[6,10]
    transp_plus[6,11] <- y["CGA"]+y["CGC"]+y["CGG"]+y["CGT"]+transp_plus[6,11]
    transp_plus[6,12] <- y["CTA"]+y["CTC"]+y["CTG"]+y["CTT"]+transp_plus[6,12]
    
    transp_plus[7,9] <- y["GAA"]+y["GAC"]+y["GAG"]+y["GAT"]+transp_plus[7,9]
    transp_plus[7,10] <- y["GCA"]+y["GCC"]+y["GCG"]+y["GCT"]+transp_plus[7,10]
    transp_plus[7,11] <- y["GGA"]+y["GGC"]+y["GGG"]+y["GGT"]+transp_plus[7,11]
    transp_plus[7,12] <- y["GTA"]+y["GTC"]+y["GTG"]+y["GTT"]+transp_plus[7,12]
    
    transp_plus[8,10] <- y["TCA"]+y["TCC"]+y["TCG"]+y["TCT"]+transp_plus[8,10]
    transp_plus[8,12] <- y["TTA"]+y["TTC"]+y["TTG"]+y["TTT"]+transp_plus[8,12]
    transp_plus[8,13] <- y["TAA"]+y["TAC"]+y["TAG"]+y["TAT"]+transp_plus[8,13]
    transp_plus[8,14] <- y["TGA"]+y["TGC"]+y["TGG"]+y["TGT"]+transp_plus[8,14]
    
    transp_plus[9,15] <- y["AAA"]+y["CAA"]+y["GAA"]+transp_plus[9,15]
    transp_plus[9,16] <- y["AAC"]+y["CAC"]+y["GAC"]+transp_plus[9,16]
    transp_plus[9,17] <- y["AAG"]+y["CAG"]+y["GAG"]+transp_plus[9,17]
    transp_plus[9,18] <- y["AAT"]+y["CAT"]+y["GAT"]+transp_plus[9,18]
    
    transp_plus[10,15] <- y["ACA"]+y["CCA"]+y["GCA"]+y["TCA"]+transp_plus[10,15]
    transp_plus[10,16] <- y["ACC"]+y["CCC"]+y["GCC"]+y["TCC"]+transp_plus[10,16]
    transp_plus[10,17] <- y["ACG"]+y["CCG"]+y["GCG"]+y["TCG"]+transp_plus[10,17]
    transp_plus[10,18] <- y["ACT"]+y["CCT"]+y["GCT"]+y["TCT"]+transp_plus[10,18]
    
    transp_plus[11,15] <- y["AGA"]+y["CGA"]+y["GGA"]+transp_plus[11,15]
    transp_plus[11,16] <- y["AGC"]+y["CGC"]+y["GGC"]+transp_plus[11,16]
    transp_plus[11,17] <- y["AGG"]+y["CGG"]+y["GGG"]+transp_plus[11,17]
    transp_plus[11,18] <- y["AGT"]+y["CGT"]+y["GGT"]+transp_plus[11,18]
    
    transp_plus[12,15] <- y["ATA"]+y["CTA"]+y["GTA"]+y["TTA"]+transp_plus[12,15]
    transp_plus[12,16] <- y["ATC"]+y["CTC"]+y["GTC"]+y["TTC"]+transp_plus[12,16]
    transp_plus[12,17] <- y["ATG"]+y["CTG"]+y["GTG"]+y["TTG"]+transp_plus[12,17]
    transp_plus[12,18] <- y["ATT"]+y["CTT"]+y["GTT"]+y["TTT"]+transp_plus[12,18]
    
    transp_plus[13,16] <- y["TAC"]+transp_plus[13,16]
    transp_plus[13,18] <- y["TAT"]+transp_plus[13,18]
    transp_plus[13,19] <- y["TAA"]+transp_plus[13,19]
    transp_plus[13,10] <- y["TAG"]+transp_plus[13,20]
    
    transp_plus[14,16] <- y["TGC"]+transp_plus[14,16]
    transp_plus[14,17] <- y["TGG"]+transp_plus[14,17]
    transp_plus[14,18] <- y["TGT"]+transp_plus[14,18]
    transp_plus[14,19] <- y["TGA"]+transp_plus[14,19]
    
  }
  if (transp_plus[13,19] == 0) {
    transp_plus[13,19] <- transp_plus[13,19]+1
  }
  if (transp_plus[14,19] == 0) {
    transp_plus[14,19] <- transp_plus[14,19]+1
  }
  if (transp_plus[13,20] == 0) {
    transp_plus[13,20] <- transp_plus[13,20]+1
  }
  if (transp_plus[19,21] == 0) {
    transp_plus[19,21] <- transp_plus[19,21]+1
  }
  if (transp_plus[20,21] == 0) {
    transp_plus[20,21] <- transp_plus[20,21]+1
  }
  
  for (p in 1:21) {
    if (sum(transp_plus[p, ])>0) {
      transp_plus[p,] <- transp_plus[p,]/sum(transp_plus[p,])
    }else {
      transp_plus[p,] <- 0
    }
  }
  return(transp_plus)
}



# calculation of log-odd ratio for a given sequence -----------------------
logORatio <- function(sequence, transp_plus, transp_minus){
  result <- seqAnalyais(sequence)
  logtrp_p <- log(transp_plus)
  logtrp_m <- log(transp_minus)
  
  # P(seq|mModel)
  n <- 1
  pseq_m <- 0
  pseq_m <- pseq_m + logtrp_m[1,s2n(substr(sequence, 1,1))+2]
  while (n+1<=nchar(sequence)) {
    pseq_m <- pseq_m + logtrp_m[s2n(substr(sequence, n, n))+2, s2n(substr(sequence, n+1, n+1))+2]
    n <- n+1
  }
  pseq_m <- pseq_m + logtrp_m[s2n(substr(sequence, nchar(sequence), nchar(sequence)))+2,6]

  
  # P(seq|pModel)
  pseq_p <- 0

  if (result[1]!=1) {
    pseq_p <- 0
  }else {
    pseq_p <- pseq_p+logtrp_p[1,2]+logtrp_p[2,3]+logtrp_p[3,4]
    pseq_p <- pseq_p+logtrp_p[4,s2n(substr(sequence, 4, 4))+5]

    m <- 1
    while (3*m+4<=result[2]) {
      pseq_p <- pseq_p+logtrp_p[s2n(substr(sequence, 3*m+3, 3*m+3))+15, s2n(substr(sequence,3*m+4,3*m+4))+5]
      m <- m+1
    }
    
    if (result[3]==0) {
      codingRegion <- substr(sequence,result[1]+3, result[2])
      pseq_p <- pseq_p+logtrp_p[s2n(substr(sequence,result[2],result[2]))+15,21]
      
    } else {
      if (substr(sequence,result[2]-1,result[2]-1)=="A") {
        if (substr(sequence,result[2],result[2])=="A") {
          pseq_p <- pseq_p+logtrp_p[8,13]+logtrp_p[13,19]+logtrp_p[19,21]
        } else if (substr(sequence,result[2],result[2])=="G") {
          pseq_p <- pseq_p+logtrp_p[8,13]+logtrp_p[13,20]+logtrp_p[20,21]
        }
      } else if (substr(sequence,result[2]-1,result[2]-1)=="G") {
        pseq_p <- pseq_p+logtrp_p[8,14]+logtrp_p[14,19]+logtrp_p[19,21]
      }
      
      if (result[2]-result[1]>6) {
        codingRegion <- substr(sequence,result[1]+3, result[2]-3)
      } else {
        codingRegion <- ""
      }
      
    }
    
    if (codingRegion!="") {
      y <- count(s2c(codingRegion),word = 3, start = 0, by = 3, freq = FALSE,alphabet = s2c("ACGT"))      
      pseq_p <- pseq_p+logtrp_p[5,9]*(y["AAA"]+y["AAC"]+y["AAG"]+y["AAT"])
      pseq_p <- pseq_p+logtrp_p[5,10]*(y["ACA"]+y["ACC"]+y["ACG"]+y["ACT"])
      pseq_p <- pseq_p+logtrp_p[5,11]*(y["AGA"]+y["AGC"]+y["AGG"]+y["AGT"])
      pseq_p <- pseq_p+logtrp_p[5,12]*(y["ATA"]+y["ATC"]+y["ATG"]+y["ATT"])
      
      pseq_p <- pseq_p+logtrp_p[6,9]*(y["CAA"]+y["CAC"]+y["CAG"]+y["CAT"])
      pseq_p <- pseq_p+logtrp_p[6,10]*(y["CCA"]+y["CCC"]+y["CCG"]+y["CCT"])
      pseq_p <- pseq_p+logtrp_p[6,11]*(y["CGA"]+y["CGC"]+y["CGG"]+y["CGT"])
      pseq_p <- pseq_p+logtrp_p[6,12]*(y["CTA"]+y["CTC"]+y["CTG"]+y["CTT"])
      
      pseq_p <- pseq_p+logtrp_p[7,9]*(y["GAA"]+y["GAC"]+y["GAG"]+y["GAT"])
      pseq_p <- pseq_p+logtrp_p[7,10]*(y["GCA"]+y["GCC"]+y["GCG"]+y["GCT"])
      pseq_p <- pseq_p+logtrp_p[7,11]*(y["GGA"]+y["GGC"]+y["GGG"]+y["GGT"])
      pseq_p <- pseq_p+logtrp_p[7,12]*(y["GTA"]+y["GTC"]+y["GTG"]+y["GTT"])
      
      pseq_p <- pseq_p+logtrp_p[8,10]*(y["TCA"]+y["TCC"]+y["TCG"]+y["TCT"])
      pseq_p <- pseq_p+logtrp_p[8,12]*(y["TTA"]+y["TTC"]+y["TTG"]+y["TTT"])
      pseq_p <- pseq_p+logtrp_p[8,13]*(y["TAC"]+y["TAT"])
      pseq_p <- pseq_p+logtrp_p[8,14]*(y["TGC"]+y["TGG"]+y["TGT"])
      
      pseq_p <- pseq_p+logtrp_p[9,15]*(y["AAA"]+y["CAA"]+y["GAA"]);
      pseq_p <- pseq_p+logtrp_p[9,16]*(y["AAC"]+y["CAC"]+y["GAC"]);
      pseq_p <- pseq_p+logtrp_p[9,17]*(y["AAG"]+y["CAG"]+y["GAG"]);
      pseq_p <- pseq_p+logtrp_p[9,18]*(y["AAT"]+y["CAT"]+y["GAT"]);
      
      pseq_p <- pseq_p+logtrp_p[10,15]*(y["ACA"]+y["CCA"]+y["GCA"]+y["TCA"]);
      pseq_p <- pseq_p+logtrp_p[10,16]*(y["ACC"]+y["CCC"]+y["GCC"]+y["TCC"]);
      pseq_p <- pseq_p+logtrp_p[10,17]*(y["ACG"]+y["CCG"]+y["GCG"]+y["TCG"]);
      pseq_p <- pseq_p+logtrp_p[10,18]*(y["ACT"]+y["CCT"]+y["GCT"]+y["TCT"]);
      
      pseq_p <- pseq_p+logtrp_p[11,15]*(y["AGA"]+y["CGA"]+y["GGA"]);
      pseq_p <- pseq_p+logtrp_p[11,16]*(y["AGC"]+y["CGC"]+y["GGC"]);
      pseq_p <- pseq_p+logtrp_p[11,17]*(y["AGG"]+y["CGG"]+y["GGG"]);
      pseq_p <- pseq_p+logtrp_p[11,18]*(y["AGT"]+y["CGT"]+y["GGT"]);
      
      pseq_p <- pseq_p+logtrp_p[12,15]*(y["ATA"]+y["CTA"]+y["GTA"]+y["TTA"]);
      pseq_p <- pseq_p+logtrp_p[12,16]*(y["ATC"]+y["CTC"]+y["GTC"]+y["TTC"]);
      pseq_p <- pseq_p+logtrp_p[12,17]*(y["ATG"]+y["CTG"]+y["GTG"]+y["TTG"]);
      pseq_p <- pseq_p+logtrp_p[12,18]*(y["ATT"]+y["CTT"]+y["GTT"]+y["TTT"]);
      
      pseq_p <- pseq_p+logtrp_p[13,16] * y["TAC"];
      pseq_p <- pseq_p+logtrp_p[13,18] * y["TAT"];
      
      pseq_p <- pseq_p+logtrp_p[14,16] * y["TGC"];
      pseq_p <- pseq_p+logtrp_p[14,17] * y["TGG"];
      pseq_p <- pseq_p+logtrp_p[14,18] * y["TGT"];
    } 
  }
  logr <- (pseq_p-pseq_m)/nchar(sequence)
  return(logr)
}

transp_minus <- mModel(train_minus_seq)
transp_plus <- pModel(train_plus_seq)

# log-odd ratios for all training sequences -------------------------------
logr_trainp <- matrix(0, length(train_plus_seq),1)
logr_trainm <- matrix(0, length(train_minus_seq), 1)

i <- 1
while (i<=length(train_plus_seq)) {
  logr_trainp[i] <- logORatio(sequence = train_plus_seq[i], transp_plus, transp_minus)
  i <- i+1
}

i <- 1
while (i<=length(train_minus_seq)) {
  logr_trainm[i] <- logORatio(sequence = train_minus_seq[i], transp_plus, transp_minus)
  i <- i+1
}

# plot the histogram of logr_trainp and logr_trainm -----------------------
library(ggplot2)
dat <- data.frame(classtype = factor(rep(c("train_plus","train_minus"),c(length(train_plus_seq),length(train_minus_seq)))), Length_normalized_log_odd_ratio_trainingSeq = c(logr_trainp, logr_trainm))
ggplot(dat, aes(x=Length_normalized_log_odd_ratio_trainingSeq, fill=classtype))+geom_histogram(binwidth = 0.06, position = "dodge")


# log-odd ratios for testing sequences ------------------------------------
logr_testp <- matrix(0,length(test_plus_seq),1)
logr_testm <- matrix(0,length(test_minus_seq),1)

i <- 1
while (i<=length(test_plus_seq)) {
  logr_testp[i] <- logORatio(test_plus_seq[i], transp_plus, transp_minus)
  i <- i+1
}

i <- 1
while (i<=length(test_minus_seq)) {
  logr_testm[i] <- logORatio(test_minus_seq[i], transp_plus, transp_minus)
  i <- i+1
}


# plot the histogram of logr_testp and logr_testm -------------------------
dat <- data.frame(classtype = factor(rep(c("test_plus","test_minus"),c(length(test_plus_seq),length(test_minus_seq)))),Length_normalized_log_odd_ratio_testingSeq = c(logr_testp, logr_testm))
ggplot(dat, aes(x=Length_normalized_log_odd_ratio_testingSeq, fill=classtype))+geom_histogram(binwidth = 0.06, position = "dodge")


# determination of threshold and evaluation of model with training --------
mEvalue <- matrix(0,5,200)
for (loop in 1:200) {
  TP <- 0
  FP <- 0
  FN <- 0
  TN <- 0
  threshold <- -0.20+0.005*loop
  
  i <- 1
  while (i<=length(train_plus_seq)) {
    if (logr_trainp[i]< threshold) {
      TP <- TP+1
    } else {
      FN <- FN+1
    }
    i <- i+1
  }
  i <- 1
  while (i<=length(train_minus_seq)) {
    if (logr_trainm[i]< threshold) {
      FP <- FP+1
    } else {
      TN <- TN+1
    }
    i <- i+1
  }
  mEvalue[1,loop] <- threshold
  mEvalue[2,loop] <- (TP+TN)/(TP+TN+FP+FN)
  mEvalue[3,loop] <- TP/(TP+FN)
  mEvalue[4,loop] <- TN/(TN+FP)
  mEvalue[5,loop] <- 2*TP/(2*TP+FN+FP)
}
thresholdIndex <- which.max(mEvalue[5,])
threshold <- mEvalue[1,thresholdIndex]
print(mEvalue[,thresholdIndex])

# Model evaluation using testing sequences --------------------------------
testEvaluation <- matrix(5,1)

TP <- 0
FP <- 0
FN <- 0
TN <- 0

i <- 1
while (i<=length(test_plus_seq)) {
  if (logr_testp[i]< threshold) {
    TP <- TP+1
  } else {
    FN <- FN+1
  }
  i <- i+1
}
i <- 1
while (i<=length(test_minus_seq)) {
  if (logr_testm[i]< threshold) {
    FP <- FP+1
  } else {
    TN <- TN+1
  }
  i <- i+1
}

testEvaluation[1] <- threshold
testEvaluation[2] <- (TP+TN)/(TP+TN+FP+FN)
testEvaluation[3] <- TP/(TP+FN)
testEvaluation[4] <- TN/(TN+FP)
testEvaluation[5] <- 2*TP/(2*TP+FN+FP)
print(testEvaluation)