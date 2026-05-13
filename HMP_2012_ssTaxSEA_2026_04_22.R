# A script to run and execute TaxSEA basic code and practice.

### 1. Explore the available data in curatedMetagenomicData:
library(curatedMetagenomicData)
library(TaxSEA)
library(MicrobiomeStat)
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

### 2. Load and pre-process the HMP_2019_ibdmdb dataset:
# counts = TRUE returns raw counts rather than relative abundances
# dryrun = FALSE performs the actual download

HMP_2012_cmd_object <- curatedMetagenomicData(
  pattern = "HMP_2012.relative_abundance",
  counts   = TRUE,
  dryrun   = FALSE
)

# Extract the SummarizedExperiment Object from the returned list.
HMP_2012_counts.df <- assay(
  HMP_2012_cmd_object$
    `2021-03-31.HMP_2012.relative_abundance`)

# QC: retain only taxa present at > 1000 counts in at least 4 samples:
HMP_2012_counts.df <- HMP_2012_counts.df[
  apply(HMP_2012_counts.df > 1000, 1, sum) > 3,
]

# Subset the global sample metadata to only samples present in our count matrix
HMP_2012_md.df <- sampleMetadata[
  sampleMetadata$sample_id %in% colnames(HMP_2012_counts.df),
]
rownames(HMP_2012_md.df) <- HMP_2012_md.df$sample_id

# Synchronize the count matrix and the metadata to the same set of sample IDs
s2k <- intersect(rownames(HMP_2012_md.df), colnames(HMP_2012_counts.df))
HMP_2012_counts.df <- HMP_2012_counts.df[, s2k]
HMP_2012_md.df     <- HMP_2012_md.df[s2k, ]

# Set factor level so it is ordered into Oral-samples and Gut-Samples
# # Delete all the non-oral and non-gut samples:
# HMP_2012_md.df <- HMP_2012_md.df[
#   HMP_2012_md.df$body_site %in% c("oralcavity", "stool"),
HMP_2012_md.df$body_site <- factor(
  HMP_2012_md.df$body_site,
  levels = c("stool", "oralcavity")
)

# Remove NA rows: 
HMP_2012_md.df <- HMP_2012_md.df[!is.na(HMP_2012_md.df$body_site), ]

# For subjects with multiple visits, retain only the latest visit
# (highest visit_number) to avoid pseudoreplication
HMP_2012_md.df <- HMP_2012_md.df[
  order(HMP_2012_md.df$visit_number, decreasing = TRUE),
]
HMP_2012_md.df <- HMP_2012_md.df[
  !duplicated(HMP_2012_md.df$subject_id),
]

# Update count matrix:
HMP_2012_counts.df <- HMP_2012_counts.df[
  , rownames(HMP_2012_md.df)
]

vector_simplified <- sapply(
  rownames(HMP_2012_counts.df), function(y) strsplit(y, split = "\\|")[[1]][7]
)
vector_simplified <- gsub("s__", "", vector_simplified)
names(vector_simplified) <- NULL
rownames(HMP_2012_counts.df) <- vector_simplified

# saveRDS(HMP_2012_counts.df , "count.rds")
# saveRDS(HMP_2012_md.df , "metadata.rds")


### Linda Differential abundance analysis:
HMP_2012_res <- linda(
  HMP_2012_counts.df,
  HMP_2012_md.df,
  formula = '~body_site'
)

HMP_2012_lindaresult = HMP_2012_res$output$body_siteoralcavity
HMP_2012_ranks <- HMP_2012_res$output$body_siteoralcavity$log2FoldChange
names(HMP_2012_ranks) <- rownames(
  HMP_2012_res$output$body_siteoralcavity
)

TaxSEA_results_HMP_2012 <- TaxSEA(taxon_ranks = HMP_2012_ranks)

metabolites.df <- TaxSEA_results_HMP_2012$Metabolite_producers
disease.df <- TaxSEA_results_HMP_2012$Health_associations
Bacdive.df = TaxSEA_results_HMP_2012$BacDive_bacterial_physiology


#select important taxon sets. 
#Perhaps 4 per subtype
#have a clear scientific rationale for why you selected them
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

ssTaxSEA_metabolites = write_important_list(metabolites.df, 5)
ssTaxSEA_disease = write_important_list(disease.df, 5)
ssTaxSEA_Bacdive = write_important_list(Bacdive.df, 5)

List_of_bacteria_oral <- readRDS("~/List_of_bacteria.rds")

ssTaxSEA_Bacdive_oral <- lapply(ssTaxSEA_Bacdive, function(x) x[x %in% List_of_bacteria_oral])
ssTaxSEA_Bacdive_oral <- Filter(function(x) length(x) > 0, ssTaxSEA_Bacdive_oral)
ssTaxSEA_metabolites_oral <- lapply(ssTaxSEA_metabolites, function(x) x[x %in% List_of_bacteria_oral])
ssTaxSEA_metabolites_oral <- Filter(function(x) length(x) > 9, ssTaxSEA_metabolites_oral)

df_sig <- subset(disease.df, PValue < 0.05) # make sure data is statistically significant.

high_list <- head(df_sig[order(-df_sig$median_rank_of_set_members), "taxonSetName"], 5)
df_important <- subset(df_sig, taxonSetName %in% high_list)
df_important <- setNames(
  lapply(df_important$TaxonSet, function(taxa_string) {
    strsplit(taxa_string, ", ")[[1]]
  }),
  df_important$taxonSetName
)

List_of_bacteria <- c(df_important$mBodyMap_Oral)
saveRDS(unique(List_of_bacteria), "List_of_bacteria.rds")
readRDS("~/List_of_bacteria.rds")

ssTaxSEA_metabolites_1 <- Filter(function(x) any(x %in% List_of_bacteria), ssTaxSEA_metabolites)
ssTaxSEA_metabolites_1 <- lapply(ssTaxSEA_metabolites_1, function(x) x[x %in% List_of_bacteria])

ssTaxSEA_metabolites$MiMeDB_producers_of_Tyramine
ssTaxSEA_metabolites_1
# Run ssTaxSEA to with the curated database above. 
HMP_2012_ssTaxSEA_oral_res_metabolites = ssTaxSEA(counts = HMP_2012_counts.df, custom_db = ssTaxSEA_metabolites_oral,  min_set_size = 3)
HMP_2012_ssTaxSEA_oral_res_Bacdive = ssTaxSEA(counts = HMP_2012_counts.df, custom_db = ssTaxSEA_Bacdive,  min_set_size = 3)
# VIP_taxa <- apply(HMP_2012_ssTaxSEA_oral_res$pvalues, 1, function(x) {x < 0.05})
# HMP_2012_ssTaxSEA_oral_res_metabolites$pvalues
# HMP_2012_ssTaxSEA_oral_res_Bacdive$pvalues


Relevant_samples = subset(HMP_2012_md.df, select=c(body_site, sample_id))
Relevant_samples = Relevant_samples[order(Relevant_samples$body_site),]

ssTaxSEA_set_size <- function(df) { 
  df = df$scores
  df = df[Relevant_samples$sample_id,]
  df =t(df)
  }

HMP_2012_ssTaxSEA_oral_res_metabolites = ssTaxSEA_set_size(HMP_2012_ssTaxSEA_oral_res_metabolites)

HMP_2012_ssTaxSEA_oral_res_Bacdive = ssTaxSEA_set_size(HMP_2012_ssTaxSEA_oral_res_Bacdive)

#complexheatmap: 
library(grid)
library(ComplexHeatmap)
library(circlize)

ha = HeatmapAnnotation(body_site = Relevant_samples$body_site,
                       # gp = gpar(col = "black"),
                       col = list(body_site = c("oralcavity" = "#6faaa6", "stool" = "#aa6f73")),
                       annotation_label = "Body site",
                       annotation_name_side = "left",
                       annotation_name_gp =gpar(fontsize = 8, fontface="bold"),
                       show_legend = TRUE,
                       # legend_gp = gpar(fontsize = 8),
                       simple_anno_size = unit(0.4, "cm"),
                       border=FALSE
                       )
Heatmap(HMP_2012_ssTaxSEA_oral_res_Bacdive, 
        name="Enrichment Score",
        show_column_names = FALSE,
        row_names_gp = gpar(fontsize = 8),
        # bottom_annotation = NULL,
        top_annotation = ha, 
        column_km = 2,
)

