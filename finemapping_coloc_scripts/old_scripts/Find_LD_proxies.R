# Add if statement here for palindromic snps
# Need to only filter out those pallindromic with intermediate allele frequencies

# Obtain snps in LD
find_prox = function(pal_snps, LD){
  
  # Transpose the subset of LD_coloc and filter
  fil_LD <- transpose(LD[rownames(LD) == pal_snps, .SD, .SDcols = colnames(LD) != pal_snps])
  fil_LD <- as.data.table(fil_LD)
  setnames(fil_LD, colnames(LD)[colnames(LD) != pal_snps])
  fil_LD <- fil_LD[, SNP := .I][rowSums(.SD > 0.7) > 0, .(SNP)]
  
  # Filter d2 for SNPs in fil_LD and perform transformations
  d2_filtered <- d2[SNP %in% fil_LD$SNP]
  d2_filtered[, c("effect_allele", "other_allele", "MAF", "beta", "N") := .(EA, NEA, EAF, Beta, n)]
  d2_filtered[, trait := "T2D"]
  d2_filtered[, palindromic := fifelse(effect_allele %in% c("A", "T") & other_allele %in% c("A", "T") | 
                                         effect_allele %in% c("C", "G") & other_allele %in% c("C", "G"), 1, 0)]
  d2_filtered <- d2_filtered[palindromic == 0]
  
  # Join fil_LD with d2_filtered to get SNPs
  result <- d2_filtered[fil_LD, on = "SNP"]
  result[, hit2 := SNP]
  
  # Extract the row with the maximum value in the 15th column (assumed here based on your original code)
  # Ensure you adjust the column index if necessary
  extr_prx <- result[which.max(unlist(result[[15]])), .SD, .SDcols = -"palindromic"]
  extr_prx[, original_snp := pal]
  
  return(extr_prx)
  
}

#find_prox(df2$hit2)
# sapply(df2$hit2, find_prox)
# ld_prox = sapply(fil_LD, max, na.rm=T)
# data.frame(palidromic)
# sapply(colnames(fil_LD), function(x) rownames(fil_LD)[max_ld[names(max_ld) == x]])