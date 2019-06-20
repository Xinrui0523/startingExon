# training plus model -----------------------------------------------------
pModel <- function(train_plus_seq){
  transp_plus <- matrix(0,21,21)
  
  for (i in 1:length(train_plus_seq)) {
    result <- seqAnalyais(train_plus_seq[i])
    
    if (result[1]>0 && result[2]>result[1] && result[3]==1) {
      codingRegion <- substr(train_plus_seq[i],result[1]+3, result[2]-3)
      trinCount_coding <- count(s2c(codingRegion),word = 3, start = 0, by = 3, freq = FALSE, alphabet = s2c("ACGT"))
    } else {
      codingRegion <- substr(train_plus_seq[i],result[1]+3, result[2])
      trinCount_coding <- count(s2c(codingRegion),word = 3, start = 0, by = 3, freq = FALSE,alphabet = s2c("ACGT"))
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
    
    transp_plus[5,9] <- trinCount_coding["AAA"]+trinCount_coding["AAC"]+trinCount_coding["AAG"]+trinCount_coding["AAT"]+transp_plus[5,9]
    transp_plus[5,10] <- trinCount_coding["ACA"]+trinCount_coding["ACC"]+trinCount_coding["ACG"]+trinCount_coding["ACT"]+transp_plus[5,10]
    transp_plus[5,11] <- trinCount_coding["AGA"]+trinCount_coding["AGC"]+trinCount_coding["AGG"]+trinCount_coding["AGT"]+transp_plus[5,11]
    transp_plus[5,12] <- trinCount_coding["ATA"]+trinCount_coding["ATC"]+trinCount_coding["ATG"]+trinCount_coding["ATT"]+transp_plus[5,12]
    
    transp_plus[6,9] <- trinCount_coding["CAA"]+trinCount_coding["CAC"]+trinCount_coding["CAG"]+trinCount_coding["CAT"]+transp_plus[6,9]
    transp_plus[6,10] <- trinCount_coding["CCA"]+trinCount_coding["CCC"]+trinCount_coding["CCG"]+trinCount_coding["CCT"]+transp_plus[6,10]
    transp_plus[6,11] <- trinCount_coding["CGA"]+trinCount_coding["CGC"]+trinCount_coding["CGG"]+trinCount_coding["CGT"]+transp_plus[6,11]
    transp_plus[6,12] <- trinCount_coding["CTA"]+trinCount_coding["CTC"]+trinCount_coding["CTG"]+trinCount_coding["CTT"]+transp_plus[6,12]
    
    transp_plus[7,9] <- trinCount_coding["GAA"]+trinCount_coding["GAC"]+trinCount_coding["GAG"]+trinCount_coding["GAT"]+transp_plus[7,9]
    transp_plus[7,10] <- trinCount_coding["GCA"]+trinCount_coding["GCC"]+trinCount_coding["GCG"]+trinCount_coding["GCT"]+transp_plus[7,10]
    transp_plus[7,11] <- trinCount_coding["GGA"]+trinCount_coding["GGC"]+trinCount_coding["GGG"]+trinCount_coding["GGT"]+transp_plus[7,11]
    transp_plus[7,12] <- trinCount_coding["GTA"]+trinCount_coding["GTC"]+trinCount_coding["GTG"]+trinCount_coding["GTT"]+transp_plus[7,12]
    
    transp_plus[8,10] <- trinCount_coding["TCA"]+trinCount_coding["TCC"]+trinCount_coding["TCG"]+trinCount_coding["TCT"]+transp_plus[8,10]
    transp_plus[8,12] <- trinCount_coding["TTA"]+trinCount_coding["TTC"]+trinCount_coding["TTG"]+trinCount_coding["TTT"]+transp_plus[8,12]
    transp_plus[8,13] <- trinCount_coding["TAA"]+trinCount_coding["TAC"]+trinCount_coding["TAG"]+trinCount_coding["TAT"]+transp_plus[8,13]
    transp_plus[8,14] <- trinCount_coding["TGA"]+trinCount_coding["TGC"]+trinCount_coding["TGG"]+trinCount_coding["TGT"]+transp_plus[8,14]
    
    transp_plus[9,15] <- trinCount_coding["AAA"]+trinCount_coding["CAA"]+trinCount_coding["GAA"]+transp_plus[9,15]
    transp_plus[9,16] <- trinCount_coding["AAC"]+trinCount_coding["CAC"]+trinCount_coding["GAC"]+transp_plus[9,16]
    transp_plus[9,17] <- trinCount_coding["AAG"]+trinCount_coding["CAG"]+trinCount_coding["GAG"]+transp_plus[9,17]
    transp_plus[9,18] <- trinCount_coding["AAT"]+trinCount_coding["CAT"]+trinCount_coding["GAT"]+transp_plus[9,18]
    
    transp_plus[10,15] <- trinCount_coding["ACA"]+trinCount_coding["CCA"]+trinCount_coding["GCA"]+trinCount_coding["TCA"]+transp_plus[10,15]
    transp_plus[10,16] <- trinCount_coding["ACC"]+trinCount_coding["CCC"]+trinCount_coding["GCC"]+trinCount_coding["TCC"]+transp_plus[10,16]
    transp_plus[10,17] <- trinCount_coding["ACG"]+trinCount_coding["CCG"]+trinCount_coding["GCG"]+trinCount_coding["TCG"]+transp_plus[10,17]
    transp_plus[10,18] <- trinCount_coding["ACT"]+trinCount_coding["CCT"]+trinCount_coding["GCT"]+trinCount_coding["TCT"]+transp_plus[10,18]
    
    transp_plus[11,15] <- trinCount_coding["AGA"]+trinCount_coding["CGA"]+trinCount_coding["GGA"]+transp_plus[11,15]
    transp_plus[11,16] <- trinCount_coding["AGC"]+trinCount_coding["CGC"]+trinCount_coding["GGC"]+transp_plus[11,16]
    transp_plus[11,17] <- trinCount_coding["AGG"]+trinCount_coding["CGG"]+trinCount_coding["GGG"]+transp_plus[11,17]
    transp_plus[11,18] <- trinCount_coding["AGT"]+trinCount_coding["CGT"]+trinCount_coding["GGT"]+transp_plus[11,18]
    
    transp_plus[12,15] <- trinCount_coding["ATA"]+trinCount_coding["CTA"]+trinCount_coding["GTA"]+trinCount_coding["TTA"]+transp_plus[12,15]
    transp_plus[12,16] <- trinCount_coding["ATC"]+trinCount_coding["CTC"]+trinCount_coding["GTC"]+trinCount_coding["TTC"]+transp_plus[12,16]
    transp_plus[12,17] <- trinCount_coding["ATG"]+trinCount_coding["CTG"]+trinCount_coding["GTG"]+trinCount_coding["TTG"]+transp_plus[12,17]
    transp_plus[12,18] <- trinCount_coding["ATT"]+trinCount_coding["CTT"]+trinCount_coding["GTT"]+trinCount_coding["TTT"]+transp_plus[12,18]
    
    transp_plus[13,16] <- trinCount_coding["TAC"]+transp_plus[13,16]
    transp_plus[13,18] <- trinCount_coding["TAT"]+transp_plus[13,18]
    transp_plus[13,19] <- trinCount_coding["TAA"]+transp_plus[13,19]
    transp_plus[13,10] <- trinCount_coding["TAG"]+transp_plus[13,20]
    
    transp_plus[14,16] <- trinCount_coding["TGC"]+transp_plus[14,16]
    transp_plus[14,17] <- trinCount_coding["TGG"]+transp_plus[14,17]
    transp_plus[14,18] <- trinCount_coding["TGT"]+transp_plus[14,18]
    transp_plus[14,19] <- trinCount_coding["TGA"]+transp_plus[14,19]
    
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