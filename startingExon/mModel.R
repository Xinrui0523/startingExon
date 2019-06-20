# training minus model
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
    
    transp_minus[1,s2n(substr(train_minus_seq[i],1,1))+2] <- transp_minus[1,s2n(substr(train_minus_seq[1],1,1))+2]+1
    transp_minus[s2n(substr(train_minus_seq[i],nchar(train_minus_seq[i]),nchar(train_minus_seq[i])))+2,6] <- transp_minus[s2n(substr(train_minus_seq[i],nchar(train_minus_seq[i]),nchar(train_minus_seq[i])))+2,6]+1
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