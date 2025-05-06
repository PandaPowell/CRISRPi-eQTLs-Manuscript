library("httr")
library("glue")
library("ggrepel")
library("jsonlite")
library("data.table")

wd = "/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/eQTL_catalogue/"

setwd(wd)

# eQTL studies to include
bulk_study_ids = c("QTS000001", "QTS000002","QTS000012", "QTS000013","QTS000021","QTS000024","QTS000029","QTS000026",
                   "QTS000019","QTS000006","QTS000003")

sc_study_ids = c("QTS000036","QTS000037","QTS000040","QTS000041")

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

fwrite(datasets,"eQTL_catalogue_datasets.txt", row.names=F, quote =F)

# Download lbf variables from eQTL catalogue
options(timeout = 30000)
# download files if not already present

remote = "http://ftp.ebi.ac.uk/pub/databases/spot/eQTL/susie/"

for (i in 1:nrow(datasets)) {
  
  study_id = datasets$study_id[i]
  dataset_id = datasets$dataset_id[i]
  
  full_remote = paste0(remote,study_id,"/",dataset_id,"/")
  
  filename = paste0(dataset_id,".lbf_variable.txt.gz")
  
  data_dir = paste0(study_id,"/",dataset_id)
  system(paste0("mkdir -p ", data_dir))
  
  if (!file.exists(paste0(data_dir, "/", filename))) {
    
    cat(paste0("Downloading ", filename, "\n"))
    
    source <- paste0(full_remote, filename)
    
    dest <- paste0(data_dir, "/", filename)
    system(paste0("wget -t 3 -T 300 -P ", data_dir, " ", source))
    
    
  }
}
