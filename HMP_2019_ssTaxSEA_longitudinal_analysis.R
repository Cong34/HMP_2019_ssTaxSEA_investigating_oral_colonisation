################################################################################
# Metagenomic Analysis of IBD using curatedMetagenomicData
#
# This script demonstrates how to:
#   1. Explore available datasets in curatedMetagenomicData (CMD)
#   2. Extract and filter count data from the HMP_2019_ibdmdb dataset
#   3. Run a linear model for differential abundance using LinDA (MicrobiomeStat)
#   4. Perform taxon set enrichment analysis using TaxSEA
#
# Dependencies: curatedMetagenomicData, TaxSEA, MicrobiomeStat
################################################################################

library(curatedMetagenomicData)
library(TaxSEA)
library(MicrobiomeStat)
library(grid)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(dplyr)
library(hrbrthemes)

################################################################################
# 1. Explore available data in curatedMetagenomicData
################################################################################
# # Preview:
# head(sampleMetadata[, 1:10])
#
# # List available studies:
# table(sampleMetadata$study_name) # USing this code to find HMP_2012
#
# Example_data <- sampleMetadata[sampleMetadata$study_name == "HMP_2012",]
#
# # List disease/condition labels:
# table(sampleMetadata$study_condition)
#
# Example_data <- sampleMetadata[sampleMetadata$study_condition == "schizophrenia",]



################################################################################
# 2. Load and pre-process the HMP_2019_ibdmdb dataset:
################################################################################

# Download the count data for HMP_2019 ibd database.
HMP_2019_ibdmdb_cmd_object <- curatedMetagenomicData(
  pattern = "HMP_2019_ibdmdb.relative_abundance",
  counts   = TRUE, # counts = TRUE returns raw counts rather than relative abundances
  dryrun   = FALSE # dryrun = FALSE performs the actual download
)

### curatedMetagenomicData(pattern = "HMP_2019_ibdmdb") # Use this command to check the database.

# Extract the SummarizedExperiment Object from the returned list
# and pull out the count matrix (taxa x samples)
HMP_2019_counts.df <- assay(
  HMP_2019_ibdmdb_cmd_object$
    `2021-10-14.HMP_2019_ibdmdb.relative_abundance`)

### HMP_2019_ibdmdb_cmd_object$#press_tab # Use this command to find the object.

# QC: retain only taxa present at > 1000 counts in at least 4 samples:
HMP_2019_counts.df <- HMP_2019_counts.df[
  apply(HMP_2019_counts.df > 1000, 1, sum) > 3,
]

# Subset the global sample metadata to only samples present in our count matrix
HMP_2019_md.df <- sampleMetadata[
  sampleMetadata$sample_id %in% colnames(HMP_2019_counts.df),
]
rownames(HMP_2019_md.df) <- HMP_2019_md.df$sample_id

# Synchronize the count matrix and the metadata to the same set of sample IDs
s2k <- intersect(rownames(HMP_2019_md.df), colnames(HMP_2019_counts.df))
HMP_2019_counts.df <- HMP_2019_counts.df[, s2k]
HMP_2019_md.df     <- HMP_2019_md.df[s2k, ]

# Recode disease subtype: samples with no subtype annotation are controls
HMP_2019_md.df$disease_subtype[is.na(HMP_2019_md.df$disease_subtype)] <- "Control"

# Set factor levels so Control is the reference group in downstream models
HMP_2019_md.df$disease_subtype <- factor(
  HMP_2019_md.df$disease_subtype,
  levels = c("Control", "UC", "CD")
)

# The longitudinal analysis require the multiple visit data to be retained, 

# Simplify row names from full taxonomic strings to species-level names
# CMD rownames are pipe-delimited (k__|p__|...|s__Genus_species)
# We extract the 7th element (species level) and strip the "s__" prefix
vector_simplified <- sapply(
  rownames(HMP_2019_counts.df), function(y) strsplit(y, split = "\\|")[[1]][7]
)
vector_simplified <- gsub("s__", "", vector_simplified)
names(vector_simplified) <- NULL
rownames(HMP_2019_counts.df) <- vector_simplified

# Update matrix: 
HMP_2019_counts.df <- HMP_2019_counts.df[
  , rownames(HMP_2019_md.df)
]

################################################################################
# 3. Sorting for longitudinal analysis: 
################################################################################
### Check to see if there are any full 0 columns: 
zero_samples <- colSums(HMP_2019_counts.df) == 0
sum(zero_samples)

# Remove the 0 columns from count and update the mt matrix:
# Removing the zero columns:
HMP_2019_counts.df <- HMP_2019_counts.df[rowSums(is.na(HMP_2019_counts.df)) == 0, ]
HMP_2019_counts.df <- HMP_2019_counts.df[, colSums(HMP_2019_counts.df) > 0]

# Check for similar samples between the 2 matrices: 
samples <- intersect(colnames(HMP_2019_counts.df), HMP_2019_md.df$sample_id)

# update the 2 matrix for similar samples: 
HMP_2019_counts.df   <- HMP_2019_counts.df[, samples]
HMP_2019_md.df       <- HMP_2019_md.df[HMP_2019_md.df$sample_id %in% samples, ]

# Check if the dimension is accurate: 
dim(HMP_2019_counts.df)      # rows = features, cols = samples
dim(HMP_2019_md.df)         # rows = samples

### As there are a lot of patients data, a few main metadata df will be created 
### First, to make the last day only metadata: 
last_day_md <- HMP_2019_md.df[
  order(HMP_2019_md.df$visit_number, decreasing = TRUE),
]
last_day_md <- last_day_md[
  !duplicated(last_day_md$subject_id), 
]
last_day_md <- last_day_md[!(last_day_md$visit_number == 1),] # Remove patients that only come in once. 

### Making a function that can Isolate data set by specified time point: 
### First curate the patient list for only patients with more than (or equal to) 20 visits
last_day_md_only_high_visit <- last_day_md[last_day_md$visit_number >=20,] # Last_day_md is the md of the last visit of all patients
List_of_longitudinal_subjects <- last_day_md_only_high_visit$subject_id
List_of_longitudinal_subjects <- unique(List_of_longitudinal_subjects)

# Curate the metadata so only these patients are left: 
HMP_2019_md_longitudinal.df <- HMP_2019_md.df[HMP_2019_md.df$subject_id %in% List_of_longitudinal_subjects,]

# Function that create a metadata df base on the specified visit number: 
make_md_longitudinal <- function(time_point) {
  timepoint_list <- c(time_point, time_point-2, time_point-1,  time_point+1, time_point+2, time_point+3)
  df <- HMP_2019_md_longitudinal.df[HMP_2019_md_longitudinal.df$visit_number %in% timepoint_list,]
  df <- df[order(abs(df$visit_number - time_point)), ]
  df <- df[ !duplicated(df$subject_id),]
  return(df)
}

################################################################################
# 4. Differential abundance analysis with LinDA (MicrobiomeStat)
################################################################################

# LinDA fits a linear model to log-ratio transformed counts.
run_linda_study_condition_function <- function(mt_df, count_df) {
  # Running Linda: 
  linda_res <- linda(
    count_df,
    mt_df,
    formula = '~study_condition'
  )
  
  linda_res <- linda_res$output$study_conditionIBD
  linda_rank        <- linda_res$log2FoldChange
  names(linda_rank) <- rownames(linda_res)
  
  return(linda_rank)
}


################################################################################
# 5. Run TaxSEA and ssTaxSEA: 
################################################################################
# Function that generate ssTaxSEA output given an input metadata: 
runssTaxSEA_Oral_taxon_only <- function(md) {
  
  count_df <- HMP_2019_counts.df[, rownames(md)]
  
  rank_md <- run_linda_study_condition_function(md, count_df)
  
  TaxSEA_results <- TaxSEA(taxon_ranks = rank_md)
  
  write_important_taxonset <- function( df , x) {
    # df is the input TaxSEA results
    # x is the total of important taxa extract needed.
    
    # df <- subset(df, PValue < 0.05) # make sure data is statistically significant.
    
    # Important list must include significant results,
    high_list <- head(df[order(-df$median_rank_of_set_members), "taxonSetName"], x)    # meaning highest median results
    low_list  <- head(df[order( df$median_rank_of_set_members), "taxonSetName"], x)    # and also lowest median results.
    important_list <- c(high_list, low_list)
    
    df_important <- subset(df, taxonSetName %in% important_list)
    
    named_list <- setNames(
      lapply(df_important$TaxonSet, function(taxa_string) {
        strsplit(taxa_string, ", ")[[1]]
      }),
      df_important$taxonSetName
    )
    return(named_list) #The output is a namelist of important taxonsets and all of the bacteria taxa
  }
  
  metabolites_res <- write_important_taxonset(TaxSEA_results$Metabolite_producers, 10)
  disease_res     <- write_important_taxonset(TaxSEA_results$Health_associations, 20)
  bacdive_res     <- write_important_taxonset(TaxSEA_results$BacDive_bacterial_physiology, 10)
  
  Important_taxonset <- c(metabolites_res, bacdive_res, disease_res)
  
  # RUn ssTaxSEA with custom DB
  TaxSEA_results <- ssTaxSEA(counts = count_df, custom_db = Important_taxonset)
  TaxSEA_results <- TaxSEA_results$scores
  
  TaxSEA_results <- TaxSEA_results[, colnames(TaxSEA_results)
                                       %in% c("mBodyMap_Oral"
                                              # ,"BacDive_facultative anaerobe"
                                              # ,"MiMeDB_producers_of_Serotonin"
                                              # ,"MiMeDB_producers_of_Dopamine"
                                              ),
                                       drop=FALSE]
  # Curate from the metadata: 
  Relevant_samples = subset(md, select=c(disease_subtype, sample_id, subject_id))
  rownames(Relevant_samples) = Relevant_samples$sample_id
  
  TaxSEA_results_merged <- merge(Relevant_samples, TaxSEA_results, by = "row.names", all = TRUE)
  TaxSEA_results_merged <- TaxSEA_results_merged[, !colnames(TaxSEA_results_merged) %in% c("Row.names", "sample_id")]

  return(TaxSEA_results_merged)
}



################################################################################
# 6. Curating data for the time series: 
################################################################################
md_1 <- make_md_longitudinal(1)
md_5 <- make_md_longitudinal(5)
md_10 <- make_md_longitudinal(10)
md_15 <- make_md_longitudinal(15)
md_20 <- make_md_longitudinal(20)
md_last <- last_day_md_only_high_visit

# count_df <- HMP_2019_counts.df[, rownames(md_1)]
# rank_md <- run_linda_study_condition_function(md_1, count_df)
# TaxSEA_results <- TaxSEA(taxon_ranks = rank_md)
# TaxSEA_results$Health_associations


ssTaxSEA_output_visit_1 <- runssTaxSEA_Oral_taxon_only(md_1)
ssTaxSEA_output_visit_5 <- runssTaxSEA_Oral_taxon_only(md_5)
ssTaxSEA_output_visit_10 <- runssTaxSEA_Oral_taxon_only(md_10)
ssTaxSEA_output_visit_15 <- runssTaxSEA_Oral_taxon_only(md_15)
ssTaxSEA_output_visit_20 <- runssTaxSEA_Oral_taxon_only(md_20)
ssTaxSEA_output_visit_last <- runssTaxSEA_Oral_taxon_only(md_last)

TaxSEA_results_merged <- Reduce(function(x, y) merge(x, y, by = c("subject_id", "disease_subtype"), all = TRUE),
                                list(ssTaxSEA_output_visit_1,
                                     ssTaxSEA_output_visit_5,
                                     ssTaxSEA_output_visit_10,
                                     ssTaxSEA_output_visit_15,
                                     ssTaxSEA_output_visit_20,
                                     ssTaxSEA_output_visit_last)
)
colnames(TaxSEA_results_merged) <- c("subject_id ", "disease_subtype ", "Timepoint_1", "Timepoint_2", "Timepoint_3", "Timepoint_4", "Timepoint_5", "Timepoint_6")
rownames(TaxSEA_results_merged) <- TaxSEA_results_merged$`subject_id `
TaxSEA_results_merged <- TaxSEA_results_merged[order(TaxSEA_results_merged$`disease_subtype `),]

TaxSEA_results_dropped <- TaxSEA_results_merged[, -c(1,2)]

# ha = HeatmapAnnotation(study_condition = TaxSEA_results_merged$`disease_subtype `,
#                        col = list(study_condition = setNames(
#                          c("#bcc8ff","#aa6f73", "#8e6d3d" ),
#                          c("Control", "UC", "CD")  
#                        )),
#                        annotation_label = "Disease Subtype",
#                        # annotation_name_side = "left",
#                        annotation_name_gp =gpar(fontsize = 8, fontface="bold"),
#                        show_legend = TRUE,
#                        simple_anno_size = unit(0.4, "cm"),
#                        border=FALSE
# )
# 
# 
# 
# 
# Heatmap(t(TaxSEA_results_dropped),
#         name="Enrichment Score",
#         # show_column_names = FALSE,
#         cluster_rows = FALSE,
#         # cluster_columns = FALSE,
#         row_names_gp = gpar(fontsize = 8),
#         top_annotation = ha,
# )


################################################################################
# 7. Making the time series: 
################################################################################

### Using the heatmap above, a list of potentially meaningful patients data were chosen: 
### Now, we will show case the different groups: 
### Increasing group: "P6018", "H4015", "C3011", "C3037", "H4016", "E5013", "C3009", "H4030", "H4027", "H4038"
### Decreasing group: "H4042", "P6016", "P6010", "H4032", "E5004", "C3013", "M2077", "M2027", "H4008", "M2085", "M2047",
###                   "H4014", "H4023", "H4040", "M2068", "H4024"
### Control group: "M2075", "H4045", "M2061"
### Disease signature group: "M2034", "E5001", "M2072", "P6024", "H4006", "H4004", "H4018", "P6012", "M2069"



list_of_patients_increasing = c("P6018", "H4015", "C3011", "C3037", "H4016", "E5013", "C3009", "H4030", "H4027", "H4038")
list_of_patients_decreasing = c("H4042", "P6016", "P6010", "H4032", "E5004", "C3013", "M2077", "M2027", "H4008", "M2085", "M2047",
                                                  "H4014", "H4023", "H4040", "M2068", "H4024")
list_of_disease_n_control =  c("M2034", "E5001", "M2072", "P6024", "H4006", "H4004", "H4018", "P6012", "M2069",
                               "M2075", "H4045", "M2061")

Timeseries_function <- function(list_of_patients) {
  df <- TaxSEA_results_merged[TaxSEA_results_merged$`subject_id ` %in% list_of_patients,]
  rownames(df) <- df$`subject_id `
  df_dropped <- df[,-c(1,2)]
  
  df_long <- df %>%
    pivot_longer(
      cols = starts_with("Timepoint"),
      names_to = "Timepoint",
      values_to = "Value"
    ) %>%
    mutate(
      Timepoint = factor(
        gsub("Timepoint_", "T", Timepoint),
        levels = c("T1", "T2", "T3", "T4", "T5", "T6"),
        ordered = TRUE
      )
    )
  
  ggplot(df_long, aes(x = Timepoint, y = Value,
                      group = `subject_id `, color = `disease_subtype `)) +
    geom_line(alpha = 0.3) +
    geom_point(alpha = 0.3) +
    scale_color_brewer(palette = "Set2") +
    stat_summary( aes(group = `disease_subtype `, color = `disease_subtype `),
                  fun = mean, geom = "line", linewidth = 1.5) + 
    stat_summary( aes(group = `disease_subtype `, color = `disease_subtype `),
                  fun = mean, geom = "point") +
    theme_minimal() +
    labs(title = "Patient oral colonisation Over Time",
         x = "Time Point", y = "Enrichment Score", color = "Disease")
  }

Timeseries_function(list_of_patients = c(list_of_patients_increasing, list_of_patients_decreasing, list_of_disease_n_control))
