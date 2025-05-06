libraries <- c("coloc",
               "susieR",
               "data.table",
               "tidyverse",
               "glue",
               "purrr",
               "jsonlite",
               "ggrepel",
               "httr")

invisible(suppressMessages(lapply(libraries, require, character.only = TRUE)))

setDTthreads(8)

args = commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Supply necessary files")
} else {
  gwas_id = args[1]
}

wd = "/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/eQTL_catalogue/"

setwd(wd)

# eQTL studies to include
bulk_study_ids = c("QTS000001", "QTS000002","QTS000012", "QTS000013","QTS000021","QTS000024","QTS000029","QTS000026",
                   "QTS000019","QTS000006","QTS000003")

sc_study_ids = c("QTS000036","QTS000037","QTS000040","QTS000041")

# Change parameters
max_pulled_rows = 1000 #All datasets will be pulled if this parameter is bigger than the actual number of datasets

URL = glue("https://www.ebi.ac.uk/eqtl/api/v2/datasets/?size={max_pulled_rows}")

# Make a request
r <- GET(URL, accept_json())
# Check status
status_code(r)
# Extract content
cont <- content(r, "text", encoding = "UTF-8")
# Convert content to dataframe
datasets <- fromJSON(cont)
datasets = datasets[datasets$quant_method == "ge" & datasets$study_id %in% bulk_study_ids,]
rm_tissues = c("fat","skin")
datasets = datasets[!datasets$sample_group %in% rm_tissues,]

dict = c(
  '30000' = 'White blood cell (leukocyte) count',
  '30010' = 'Red blood cell (erythrocyte) count',
  '30020' = 'Haemoglobin concentration',
  '30030' = 'Haematocrit percentage',
  '30040' = 'Mean corpuscular volume',
  '30050' = 'Mean corpuscular haemoglobin',
  '30060' = 'Mean corpuscular haemoglobin concentration',
  '30070' = 'Red blood cell (erythrocyte) distribution width',
  '30080' = 'Platelet count',
  '30090' = 'Platelet crit',
  '30100' = 'Mean platelet (thrombocyte) volume',
  '30110' = 'Platelet distribution width',
  '30120' = 'Lymphocyte count',
  '30130' = 'Monocyte count',
  '30140' = 'Neutrophil count',
  '30150' = 'Eosinophil count',
  '30160' = 'Basophil count',
  '30180' = 'Lymphocyte percentage',
  '30190' = 'Monocyte percentage',
  '30200' = 'Neutrophil percentage',
  '30210' = 'Eosinophil percentage',
  '30220' = 'Basophil percentage',
  '30240' = 'Reticulocyte percentage',
  '30250' = 'Reticulocyte count',
  '30260' = 'Mean reticulocyte volume',
  '30270' = 'Mean sphered cell volume',
  '30280' = 'Immature reticulocyte fraction',
  '30290' = 'High light scatter reticulocyte percentage',
  '30300' = 'High light scatter reticulocyte count'
)

gwas_id = "30070"
study_id = datasets$study_id[15]
dataset_id = datasets$dataset_id[15]

run_coloc = function(gwas_id,study_id,dataset_id){
  
  translated_gwas_name = dict[gwas_id]
  
  ## Import lbf variable values from eQTL catalogue into R
  eqtl_file = paste0(study_id,"/",dataset_id,"/",dataset_id,".lbf_variable.txt.gz")
  eqtl_name = datasets$sample_group[datasets$dataset_id == dataset_id]
  eqtl_base = fread(eqtl_file)
  
  coloc_res = data.frame()
  
  for (i in 1:22){
    
    chr = paste0("chr",i)
    
    gwas_file = paste0("UKBB_susie_finemap_hg38/",gwas_id,"/",gwas_id,".",chr,"_finemap_results.hg38.txt")
    
    cat("Running coloc between GWAS",gwas_id,"&",study_id,dataset_id,chr,"\n")
    
    ukbb_lbf = fread(gwas_file) %>%
      select(-variant.x,-chr.x,-chr.y) %>% mutate(variant_hg38 = paste(chr,hg38_pos,a1,a0,sep="_"))
    
    #Filter eqtl lbf for variants in gwas lbf region
    eqtl = eqtl_base[variant %in% ukbb_lbf$variant_hg38,]
    
    genes = unique(eqtl$molecular_trait_id)
    
    for (gene in genes){
      
      eqtl2 = eqtl %>% filter(molecular_trait_id == gene)
      
      ukbb_lbf2 = ukbb_lbf[variant_hg38 %in% eqtl2$variant,]
      
      if(nrow(ukbb_lbf2) < 500){
        cat("Number of overlapping variants too low",nrow(ukbb_lbf2), "\n")
        next
      }
      
      ## Convert LBF data frames into matrixes suitable for the coloc.bf_bf() method
      eqtl_mat = as.matrix(dplyr::select(eqtl2, lbf_variable1:lbf_variable10))
      row.names(eqtl_mat) = eqtl2$variant
      eqtl_mat = t(eqtl_mat)
      
      ukbb_mat = as.matrix(dplyr::select(ukbb_lbf2, lbf_1:lbf_10))
      row.names(ukbb_mat) = ukbb_lbf2$variant_hg38
      ukbb_mat = t(ukbb_mat)
      row.names(ukbb_mat) = row.names(eqtl_mat)
      
      coloc = coloc::coloc.bf_bf(ukbb_mat,eqtl_mat)
      coloc = dplyr::as_tibble(coloc$summary) %>% dplyr::filter(PP.H4.abf > 0.5) %>% 
        rename(gwas_hit = hit1, eqtl_hit = hit2)
      
      if(nrow(coloc) == 0){
        cat(gene,"No significant colocalisation \n")
        next
      }
      
      cat(gene,"Significant colocalisation! \n")
      
      temp_df = ukbb_lbf %>% select(gwas_hit_hg19 = variant, gwas_hit = variant_hg38, gwas_pval = pval)
      temp_df2 = ukbb_lbf %>% select(eqtl_hit_hg19 = variant, eqtl_hit = variant_hg38)
      
      coloc = coloc %>% mutate(eqtl_name = eqtl_name, gwas = translated_gwas_name, study_id = study_id,
                               dataset_id = dataset_id, molecular_id = gene) %>%
        left_join(temp_df, "gwas_hit") %>% left_join(temp_df2,"eqtl_hit")
      
      coloc_res = bind_rows(coloc_res,coloc)
    }
    
  }
  
  return(coloc_res)
}

options(warn=1)
### In case of an error return NA ### 
run_coloc = purrr::possibly(run_coloc, otherwise = NA, quiet = F)
testie = lapply(1:nrow(datasets), function(x) run_coloc(gwas_id, datasets$study_id[x], datasets$dataset_id[x]))
t2 = do.call(rbind, testie) %>% filter(gwas_pval < 1e-05)
# Write results
res_dir = paste0("coloc_results/")
system(paste0("bash -c 'mkdir -p ", res_dir, "'"))

fwrite(t2, paste0(res_dir,gwas_id,"_coloc_results.txt"), sep = ",", row.names = F, quote=F)