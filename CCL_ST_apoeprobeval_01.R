########################################################################
# Name: CCL_ST_apoe_probe_val_01.R
# Project: 
# Purpose: verifying that Xenium apoe SNP probes work 
# scripts:
#     - 01. merging samples together and cleaning (nCount > 5)
# Input Files: Xenium slide output folders 
# Output Files: individual subset objs, and merged cleaned obj (20260915_merged_clean_apoeprobeval_01.qs2)
# Date created: 9/9/26
# Last updated: 9/15/26
# Author: Chloe Lucido
########################################################################


# Load Libraries ----
BiocManager::install("spacexr")
library(spacexr)
library(Seurat)
library(SeuratDisk)
library(future)
library(ggplot2)
library(arrow)
library(hdf5r)
library(presto)
library(glmGamPoi)
library(readr)
library(dplyr)
library(data.table)
library(tidyverse)
library(readxl)
library(patchwork)
library(sceasy)
library(reticulate)
library(SPLIT)
library(qs2) 
library(RColorBrewer)
library(Polychrome)
library(purrr)
options(future.globals.maxSize = 100 *1024^3)


# Vignette: https://satijalab.org/seurat/articles/seurat5_spatial_vignette_2 

################# PATHS and IMPORTANT VECTORS/VARIABLES ###########################
# slide paths ----
## run 2 slide 0021991  
path_991 <- '/Volumes/JOHNSON-L/07. Xenium Spatial Transcriptomics/2 Aging_Metab_Xenium/Run2_20250606__202944__20250606_AgingXMetabolism_2/Slide1_output-XETG00118__0021991__Region_1__20250606__202953'

## run 5 slide 0069002 
path_002 <- '/Volumes/JOHNSON-L/07. Xenium Spatial Transcriptomics/2 Aging_Metab_Xenium/Run5_20260306_210910_20260306_AgingXMetabolism_5/slide2_output-XETG00118__0069002__Region_1__20260306__210910'

## run 6 slide 0069130
path_130 <- '/Volumes/JOHNSON-L/07. Xenium Spatial Transcriptomics/2 Aging_Metab_Xenium/Run6_20260408__203018__20260408_AgingXMet_6_redo/output-XETG00118__0069130__Region_1__20260408__203028'

## run 7 slide 0069080
path_080 <- '/Volumes/JOHNSON-L/07. Xenium Spatial Transcriptomics/2 Aging_Metab_Xenium/Run7_20260417__210457__20260417_AgeXMetabolism_7/output-XETG00118__0069080__Region_1__20260417__210507'

## full brain XE annotations ----
XE_ann_list <- c(
  "/Users/cclu223/Desktop/Xenium_panel_verification/APOE_probe_verification/csv/XE_annotations/cell_stats/run2_slide0021991_36wk_F1_cells_stats.csv", # 36wk_F1
  "/Users/cclu223/Desktop/Xenium_panel_verification/APOE_probe_verification/csv/XE_annotations/cell_stats/run5_slide0069002_E2FAD_F3_cells_stats.csv", # E2FAD_F3
  "/Users/cclu223/Desktop/Xenium_panel_verification/APOE_probe_verification/csv/XE_annotations/cell_stats/run6_slide0069130_E3FAD_F2_cells_stats.csv", # E3FAD_F2
  "/Users/cclu223/Desktop/Xenium_panel_verification/APOE_probe_verification/csv/XE_annotations/cell_stats/run6_slide0069130_E4FAD_F3_cells_stats.csv", # E4FAD_F3
  "/Users/cclu223/Desktop/Xenium_panel_verification/APOE_probe_verification/csv/XE_annotations/cell_stats/run7_slide0069080_4s2_F1_cells_stats.csv", # 4s2_F1
  "/Users/cclu223/Desktop/Xenium_panel_verification/APOE_probe_verification/csv/XE_annotations/cell_stats/run7_slide0069080_4s2_F2_cells_stats.csv", # 4s2_F2
  "/Users/cclu223/Desktop/Xenium_panel_verification/APOE_probe_verification/csv/XE_annotations/cell_stats/run7_slide0069080_4s2M_F1_cells_stats.csv", # 4s2M_F1
  "/Users/cclu223/Desktop/Xenium_panel_verification/APOE_probe_verification/csv/XE_annotations/cell_stats/run7_slide0069080_4s2M_F2_cells_stats.csv" # 4s2M_F2
)

# slide vectors ----
slides <- c("0021991", "0069002", "0069130", "0069080")
cellID_suffixes <- c("slide991", "slide002", "slide130", "slide080")

# prefix map ----
prefix_map <- c(
  "36wk_F1" = "slide991_", 
  "E2FAD_F3" = "slide002_", 
  "E3FAD_F2" = "slide130_", 
  "E4FAD_F3" = "slide130_", 
  "4s2_F1" = "slide080_", 
  "4s2_F2" = "slide080_", 
  "4s2M_F1" = "slide080_", 
  "4s2M_F2" = "slide080_" 
) 

# metadata, dir paths, and date ----
metadata_xlsx <- "/Users/cclu223/Desktop/Xenium_panel_verification/APOE_probe_verification/csv/Xenium_apoeprobe_validation_metadata.xlsx"
output_path <- "/Users/cclu223/Desktop/Xenium_panel_verification/APOE_probe_verification/obj_files/"
date <- "20260915_"
#####################################################################

# 01. Loading in Xenium data and attaching cell ID suffixes ----

## putting paths into vector 
paths <- c(path_991, path_002, path_130, path_080)

## load xenium data
xenium_list <- map(paths, function(p) {LoadXenium(p, fov = "fov", segmentations = "cell")})

## adding slide number suffix to cell names 
### NOTE: .x is the element of xenium_list, .y is the element of cellID_suffixes list
xenium_list <- map2(xenium_list, cellID_suffixes, ~ {
  .x <- RenameCells(.x, add.cell.id = .y)
  .x 
})

### SANITY CHECK: make sure there are no duplicated cell IDs
all_ids <- unlist(map(xenium_list, colnames))

any(duplicated(all_ids))
# output: FALSE

# SANITY CHECK: map through each object and print the dims
xenium_list <- map(xenium_list, function(x) {
  print(dim(x))
  x
})

# output:
# [[1]] 480 413608
# [[2]] 480 426681
# [[3]] 480 465484
# [[4]] 480 415342

# 02. attaching sample level information from XE ann ----

## 02a. creating list of all read XE ann files ----
annotation_list <- map(XE_ann_list, function(x) { 
  
  ### reading each csv file, skipping first 2 metadata lines
  read.csv(x, skip = 2, blank.lines.skip = TRUE)
})

## 02b. changing annotation_list names to corresponding sample id ----
names(annotation_list) <- names(prefix_map)

## 02c. iterating through annotation_list and adding the proper suffix to the cell ids ----

annotation_list <- imap(annotation_list, function(df, name) {
  
  df <- annotation_list[[name]] # extracting each dataframe from annotation_list
  
  ### appending suffix to cell ids
  df <- df %>% 
    mutate(Cell.ID = paste0(prefix_map[[name]], Cell.ID), 
           #### adding sample id column
           Sample.ID = name)  
  
  return(df)
})

## 02d. merge together XE annotations belonging to same slide ----
### slide 0021991
ann_991 <- annotation_list[["36wk_F1"]]

unique(ann_991$Sample.ID) # SANITY CHECK
# output: [1] "36wk_F1"

dim(ann_991) # SANITY CHECK: rows (cells) X cols
# output:[1] 68063     5

### slide 0069002
ann_002 <- annotation_list[["E2FAD_F3"]]

unique(ann_002$Sample.ID) # SANITY CHECK
# output: [1] "E2FAD_F3"

dim(ann_002) #  SANITY CHECK
# output: [1] 73967     5

### slide 0069130
ann_130 <- bind_rows(
  annotation_list[["E3FAD_F2"]], 
  annotation_list[["E4FAD_F3"]]
)

unique(ann_130$Sample.ID) # SANITY CHECK
# output:[1] "E3FAD_F2" "E4FAD_F3"

dim(ann_130) # SANITY CHECK
# output: [1] 129168      5

### slide 0069080
ann_080 <- bind_rows(
  annotation_list[["4s2_F1"]], 
  annotation_list[["4s2_F2"]], 
  annotation_list[["4s2M_F1"]], 
  annotation_list[["4s2M_F2"]]
)

unique(ann_080$Sample.ID) # SANITY CHECK
# output: [1] "4s2_F1"  "4s2_F2"  "4s2M_F1" "4s2M_F2"

dim(ann_080) # SANITY CHECK
# output: [1] 285030      5

slide_ann_list <- list(ann_991, ann_002, ann_130, ann_080) # can't use c() to make list of dataframes because it will concatenate the dfs together instead of keeping them separate

## 02e. set rownames to cellID to match seurat obj ----
slide_ann_list <- map(slide_ann_list, function(x) {
  rownames(x) <- x$Cell.ID
  x
})

# 03. add metadata to each obj in xenium_list ----

## 03a. reorder rows to match corresponding xenium obj ----
### NOTE: .x is the element of slide_ann_list, .y is the element of xenium_list
slide_ann_list <- map2(slide_ann_list, xenium_list, ~{
   .x <- .x[Cells(.y), , drop = T] 
   print(dim(.x))
   .x
})

# output:
# [[1]] 413608      5
# [[2]] 426681      5
# [[3]] 465484      5
# [[4]] 415342      5

## 03b. add to metadata ----
### NOTE: .x is the element of slide_ann_list, .y is the element of xenium_list
xenium_list <- map2(slide_ann_list, xenium_list, ~{
  .y <- AddMetaData(.y, .x$Sample.ID, col.name = "sample_ID") 
  .y
})

# 04. subset to brains of interest -----
## 04a. remove cells with NA in sample_ID col ----
xenium_list <- map(xenium_list, function(x) {
  cat("Before subsetting: ", dim(x), "\n") # SANITY CHECK
  
  x <- subset(x, cells = colnames(x)[!is.na(x$sample_ID)])
  cat("After subsetting: ", dim(x), "\n") # SANITY CHECK
  
  print(unique(x$sample_ID)) # SANITY CHECK
  
  x
})

## 04b. save each subset file ----
xenium_list <- map2(xenium_list, slides, ~{
  qs_save(.x, paste0(output_path, date, "slide", .y, ".qs2"))
  .x
})



# 05. merge subset objs ----
merged.obj <- merge(xenium_list[[1]], y = xenium_list[-1])

# 06. attach metadata from excel sheet ----
## 06a. read in xlsx ----
meta_xlsx <- data.frame(read_excel(metadata_xlsx))

## 06b. extract metadata from seurat obj ----
meta_seurat <- merged.obj@meta.data

## 06c. convert rownames to columns to preserve cell IDs ----
meta_seurat <- meta_seurat %>%
  rownames_to_column(var = "cell_id")

## 06d. joining metadata from xlsx to metadata in merged seurat obj ----
meta_seurat <- meta_seurat %>%
  left_join(meta_xlsx, by = "sample_ID"
            #, keep = T
  )

# SANITY CHECK
table(meta_seurat$sample_ID, meta_seurat$slide_ID)
# output: 
#           0021991 0069002 0069080 0069130
#36wk_F1    68063       0       0       0
#4s2_F1         0       0   72914       0
#4s2_F2         0       0   71978       0
#4s2M_F1        0       0   73828       0
#4s2M_F2        0       0   66310       0
#E2FAD_F3       0   73967       0       0
#E3FAD_F2       0       0       0   62889
#E4FAD_F3       0       0       0   66279

## 06e. convert column back to rownames to preserve cell IDs ----
meta_seurat <- column_to_rownames(meta_seurat, var = "cell_id")

## 06f. add metadata back to merged obj ----
merged.obj@meta.data <- meta_seurat

dim(merged.obj) # SANITY CHECK
# output: 480 556228

# 07. clean obj (nCount>5) ----
merged.obj <- subset(merged.obj, subset = nCount_Xenium > 5)

dim(merged.obj) # SANITY CHECK
# output: 480 554889

# 08. save merged clean obj ----
qs_save(merged.obj, paste0(output_path, date, "merged_clean_apoeprobeval_01.qs2"))

