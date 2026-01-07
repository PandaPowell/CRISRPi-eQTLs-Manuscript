library(sceptre)
library(sceptredata)
library(fst)

setwd("Post_coloc_analysis/SCEPTRE/")

resample_results = read.fst("resampling_results.fst") %>% filter(p_value<=0.1,quality_rank_grna != "control") %>%
  filter(rejected == "FALSE", outlier_gene =="FALSE", quality_rank_grna == "top_two") %>%
  mutate(distance = abs(TSS - target_site.start)) %>% unique()

resample_results = read.fst("resampling_results.fst") %>% filter(rejected == T)

# load the data, creating a sceptre_object
directories <- paste0(
  system.file("extdata", package = "sceptredata"),
  "/highmoi_example/gem_group_", 1:2
)
data(grna_target_data_frame_highmoi)
sceptre_object <- import_data_from_cellranger(
  directories = directories,
  moi = "high",
  grna_target_data_frame = grna_target_data_frame_highmoi
)

# construct the grna-gene pairs to analyze
# returns the set of response-target pairs located on the same chromosome within distance_threshold
positive_control_pairs <- construct_positive_control_pairs(sceptre_object)
discovery_pairs <- construct_cis_pairs(
  sceptre_object, 
  positive_control_pairs = positive_control_pairs
)

# A left-tailed test is the most appropriate choice for a CRISPRi screen of enhancers as it tests a decrease
side <- "left"

sceptre_object <- set_analysis_parameters(
  sceptre_object = sceptre_object,
  discovery_pairs = discovery_pairs,
  positive_control_pairs = positive_control_pairs,
  side = side
)

sceptre_object <- run_discovery_analysis(sceptre_object, parallel = F)

result <- get_result(
  sceptre_object = sceptre_object,
  analysis = "run_discovery_analysis"
)

sig_results = result[result$p_value<0.1,]

run_calibration_check()
# apply the pipeline functions to the sceptre_object in order
sceptre_object <- sceptre_object |> # |> is R's base pipe, similar to %>%
  set_analysis_parameters(discovery_pairs, positive_control_pairs) |>
  run_calibration_check() |>
  run_power_check() |>
  run_discovery_analysis()

# write the results to disk
write_outputs_to_directory(sceptre_object, "~/sceptre_outputs")