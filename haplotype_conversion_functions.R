### Stacks haplotype conversion - functions ###
# Mark Ravinet - October 2013

# file_convert - converts from original TSV to a haplotype matrix
file_convert <- function(infile){
  
  haplo <- t(infile[, -1:-2])
  colnames(haplo) <- infile$Catalog.ID

  return(haplo)
}

# Calculate nucleotide diversity using Nei & Li's (1979) method
nuc.div <- function(x, seq_length) {
  
  haplotype_freq <- table(x)/length(x)
  
  if(length(unique(x)) <= 1){
    nd <- 0
  } else{  
    
    haplotype_comb <- combn(unique(x), 2)
    # Calculate pairwise nucleotide mismatches
    pairwise_mm <- apply(haplotype_comb, 2, function(z) {
      mapply(function(x, y) sum(x!=y), strsplit(z[1], ""), strsplit(z[2], ""))
    })
    # Calculate nucleotide diversity from frequency and mismatch data
    multi_freq <- apply(haplotype_comb, 2, function(x) haplotype_freq[x[1]] * haplotype_freq[x[2]])
    nd <- 2*sum(multi_freq*pairwise_mm)/seq_length
  }
  
  return(nd)
}

# Basic method
# # Based on formula in Hartl and Clark - poss computationally slower
# L <- 92
# pair_comp <- combn(seq, 2)
# pair_mm <- apply(pair_comp, 2, function(z) {
#   mapply(function(x, y) sum(x!=y), strsplit(z[1], ""), strsplit(z[2], ""))
# })
# 
# nuc_mm <- sum(pair_mm)/2346
# 
# ((length(seq)^2)-length(seq))
# 
# ncol(pair_comp)
# nuc_div <- nuc_mm/L 
