rm(list=ls())
library(data.table)

# Load gencode data
annot_file = "/gpfs/commons/groups/lappalainen_lab/woliveros/231005_OneK1K/data/Gencode/gencode.v33lift37.GRCh38.genes.gtf"
# Annotation file
annot <- read.table(annot_file, header = F, sep = "\t", stringsAsFactors = F)
## Keep only genes from chr1-22
annot <- annot[annot$V1 %in% c(paste0("chr", 1:22)), ]
annot <- annot[annot$V3 %in% "gene", ]
annot$gene_id <- sapply(annot$V9, function(x) {
  unlist(strsplit(sub("gene_id ", "", unlist(strsplit(x, ";"))[1]), "[.]"))[1]
})
annot$gene_name <- sapply(annot$V9, function(x) {
  unlist(strsplit(sub(".*gene_name ", "", unlist(strsplit(x, ";"))[4]), "[.]"))[1]
})
## Add start (TSS -1) and end (TSS)
## Note: if strand == +, then (start - 1, start)
## Note: if strand == -, then (end -1, end)
annot$start <- ifelse(annot$V7 %in% "+", annot$V4 - 1, annot$V5 - 1)
annot$end <- ifelse(annot$V7 %in% "+", annot$V4, annot$V5)
annot$chr_number <- as.numeric(sub("chr", "", annot$V1))
annot <- annot[order(annot$chr_number, annot$start),]

# Load GTEx whole blood tpm
gtex_tpm = fread("data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz")
gtex_tpm = gtex_tpm[`Whole Blood` > 0,]
gtex_tpm$ensembl_id = str_split_fixed(gtex_tpm$Name, "\\.",2)[,1]


# Filter gencode to only expressed genes
annot.filtered = annot[annot$gene_id %in% gtex_tpm$ensembl_id,]

fwrite(annot.filtered[,c("V1","start","end","gene_id")], "data/gencode.v33lift37.annotation.bed", quote = F, row.names = F, col.names = F, sep = "\t")
# load GWAS variants with gene distances

system("bash scripts/06.calculate_gwas_gene_distance.sh")

all_gwas = data.frame()

for (num in seq(30000,30300,10)){
    
  all_chr = data.frame()
  
  for (i in c(1:22)){
    
    if(file.exists(paste0("../results/UKBB_SuSiE_finemap/",num,"/",num,"_chr",i,"_gene_distance.bed"))){
      
      gwas = fread(paste0("../results/UKBB_SuSiE_finemap/",num,"/",num,"_chr",i,"_gene_distance.bed")) 
      
      if(nrow(gwas) > 0){
        
        gwas = gwas %>%
          arrange(abs(V12)) %>% distinct(V5,.keep_all = T)
        
        all_chr = bind_rows(all_chr, gwas)
        
      } else{
        next()
      }
    } else{
      next()
    }
  }
  
  all_gwas = bind_rows(all_gwas,all_chr)
  
}
