# sequence analysis
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
    while (substr(sequence, 3*m+1, 3*m+3)!="TAA" && substr(sequence, 3*m+1, 3*m+3)!="TAG" && substr(sequence, 3*m+1, 3*m+3)!="TGA" && 3*m+3<=nchar(sequence)) {
      m <- m+1
    }
    endsite <- 3*m+3
  } else{
    stopCodon <- 0
    endsite <- nchar(sequence)
  }
  anaseq <- substr(sequence, startsite, endsite)
  result <- c(startsite, endsite, stopCodon, anaseq)
  return(result)
}
