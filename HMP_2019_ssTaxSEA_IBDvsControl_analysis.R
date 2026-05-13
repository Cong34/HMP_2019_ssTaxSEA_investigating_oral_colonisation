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
# Set factor level so it is ordered into Oral-samples and Gut-Samples
# # Delete all the non-oral and non-gut samples:
HMP_2019_md.df$disease_subtype <- factor(
  HMP_2019_md.df$disease_subtype,
  levels = c("Control", "UC", "CD")
)

# For subjects with multiple visits, retain only the latest visit
# (highest visit_number) to avoid pseudoreplication
HMP_2019_md.df <- HMP_2019_md.df[
  order(HMP_2019_md.df$visit_number, decreasing = TRUE),
]
HMP_2019_md.df <- HMP_2019_md.df[
  !duplicated(HMP_2019_md.df$subject_id),
]

# Update count matrix:
HMP_2019_counts.df <- HMP_2019_counts.df[
  , rownames(HMP_2019_md.df)
]

# Simplify row names from full taxonomic strings to species-level names
# CMD rownames are pipe-delimited (k__|p__|...|s__Genus_species)
# We extract the 7th element (species level) and strip the "s__" prefix
vector_simplified <- sapply(
  rownames(HMP_2019_counts.df), function(y) strsplit(y, split = "\\|")[[1]][7]
)
vector_simplified <- gsub("s__", "", vector_simplified)
names(vector_simplified) <- NULL
rownames(HMP_2019_counts.df) <- vector_simplified

################################################################################
# 3. Differential abundance analysis with LinDA (MicrobiomeStat)
################################################################################

# LinDA fits a linear model to log-ratio transformed counts.
# Here we test for an effect of study_condition (IBD vs. Control)
HMP_2019_res <- linda(
  HMP_2019_counts.df,
  HMP_2019_md.df,
  formula = '~study_condition'
)

HMP_2019_lindaresult = HMP_2019_res$output$study_conditionIBD
HMP_2019_ranks <- HMP_2019_res$output$study_conditionIBD$log2FoldChange
names(HMP_2019_ranks) <- rownames(
  HMP_2019_res$output$study_conditionIBD
)


TaxSEA_results_IBD <- TaxSEA(taxon_ranks = HMP_2019_ranks)

write_important_list <- function( df , x) {
  # df is the input TaxSEA results
  # x is the total of important taxa extract needed.
  
  df_sig <- subset(df, PValue < 0.05) # make sure data is statistically significant.
  
  # Important list must include significant results,
  # meaning highest median results
  # and also lowest median results.
  high_list <- head(df_sig[order(-df_sig$median_rank_of_set_members), "taxonSetName"], x)
  low_list  <- head(df_sig[order( df_sig$median_rank_of_set_members), "taxonSetName"], x)
  important_list <- c(high_list, low_list)
  
  df_important <- subset(df_sig, taxonSetName %in% important_list)
  
  named_list <- setNames(
    lapply(df_important$TaxonSet, function(taxa_string) {
      strsplit(taxa_string, ", ")[[1]]
    }),
    df_important$taxonSetName
  )
  return(named_list)
}

# Combining the most important Taxon sets in both
metabolite_res  <- TaxSEA_results_IBD$Metabolite_producers
bacdive_res     <- TaxSEA_results_IBD$BacDive_bacterial_physiology
disease_res     <- TaxSEA_results_IBD$Health_associations





metabolite_res  <- write_important_list(metabolite_res , 4)
bacdive_res     <- write_important_list(bacdive_res    , 4)

Important_taxonset <- c(metabolite_res, bacdive_res)

ssTaxSEA_results_IBD <- ssTaxSEA(counts = HMP_2019_counts.df, custom_db = Important_taxonset, min_set_size = 3)

Relevant_samples = subset(HMP_2019_md.df, select=c(study_condition, sample_id))
ssTaxSEA_results_IBD_data <- ssTaxSEA_results_IBD$score
ssTaxSEA_results_IBD_data <- ssTaxSEA_results_IBD_data[Relevant_samples$sample_id,]


ha = HeatmapAnnotation(study_condition = Relevant_samples$study_condition,
                         col = list(study_condition = setNames(
                           c("#aa6f73", "#bcc8ff"),
                           c("IBD", "control")  
                         )),
                         annotation_label = "Disease Subtype",
                         # annotation_name_side = "left",
                         annotation_name_gp =gpar(fontsize = 8, fontface="bold"),
                         show_legend = TRUE,
                         simple_anno_size = unit(0.4, "cm"),
                         border=FALSE
  )
  
Heatmap(t(ssTaxSEA_results_IBD_data),
          name="Enrichment Score",
          show_column_names = FALSE,
          row_names_gp = gpar(fontsize = 8),
          top_annotation = ha,
  )




