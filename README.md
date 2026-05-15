# Oral colonisation in IBD patients
> HMP-2019 microbiome analysis using ssTaxSEA
## Overview

This project applies **ssTaxSEA** (subject-specific TaxSEA) to characterise oral microbiota in inflammatory bowel disease (IBD) patients using the HMP-2019 dataset. Where TaxSEA identifies important taxon sets from differential abundance analysis, ssTaxSEA extends this to pinpoint the specific subjects harbouring those taxa — enabling personalised microbial profiling.

## Goals

- Identify taxon sets and individual taxa associated with IBD.
- Characterise oral colonisation patterns in IBD patients over time.
- Assess the stability of oral colonisation through longitudinal time-series analysis.

## Dataset

| Dataset | Purpose |
|---|---|
| **HMP-2019** | Primary — longitudinal oral microbiota samples from IBD patients |
| **HMP-2012** | Secondary — oral vs gut microbiome differential abundance |

## Tools & technologies

| Tool | Role |
|---|---|
| TaxSEA / ssTaxSEA | Taxon set enrichment and individual-level analysis |
| LinDA | Differential abundance analysis |
| R | Data processing and visualisation |

## Methods

### Step 1 — Oral vs gut microbiome (HMP-2012)

Using the HMP-2012 dataset, we compared microbiome composition across body sites, focusing on differences between the oral and gut compartments. This produced a heatmap showing enrichment of different taxon sets in oral vs gut samples.

### Step 2 — IBD-enriched taxon sets (HMP-2019)

TaxSEA was applied to last-day samples from all IBD patients to identify enriched taxon sets. Cross-referencing with the HMP-2012 oral enrichment results, four taxon sets emerged with oral bacterial relevance:

- Producers of dopamine
- Producers of serotonin
- Facultative anaerobes
- Nitrate utilisers

> The `mBodyMap_Oral` taxon set was also applied to directly identify oral species detected in the gut of IBD patients.

### Step 3 — Longitudinal stability

A time-series analysis tracks how oral colonisation changes (or persists) across sampling time points, providing a measure of colonisation stability over the study period.
