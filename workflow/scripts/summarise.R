library(tidyverse)
library(jsonlite)

args = commandArgs(trailingOnly=TRUE)

## samples

# read sample names
samples_dat <- read.csv(args[1], header = T) %>%
  select(ID)

## fastp

# get fastp json paths
fastp_files <- list.files("results/fastp/", pattern = ".json", full.names = T)

# get read counts from json files
fastp_list <- lapply(fastp_files, function(x) {
  fastp_json <- fromJSON(x, flatten = T)
  sample_name <- gsub("results/fastp//|_fastp.json", "", x)
  reads_before <- fastp_json$summary$before_filtering$total_reads
  reads_after  <- fastp_json$summary$after_filtering$total_reads
  return(data.frame(ID = sample_name, 
                    reads_raw = reads_before, 
                    reads_qc = reads_after))
})

# rbind list
fastp_dat <- do.call("rbind", fastp_list)

## seqkit

# get seqkit file paths
seqkit_files <- list.files("results/seqkit/", pattern = ".txt", full.names = T)

# read as list
seqkit_list <- lapply(seqkit_files, function(x) read.table(x, header = T))

# rbind list
seqkit_dat <- do.call("rbind", seqkit_list)

# remove commas
seqkit_dat[] <- lapply(seqkit_dat, function(x) {
  gsub(",","", x)
})

# format seqkit data
seqkit_dat <- seqkit_dat %>%
  mutate(file = gsub("results/seqkit/|\\.fasta", "", file)) %>%
  select(file, num_seqs, sum_len, min_len, avg_len, max_len) %>%
  rename("ID" = "file", 
         "N seqs" = "num_seqs", 
         "Sum length" = "sum_len", 
         "Min. length" = "min_len", 
         "Avg. length" = "avg_len",
         "Max. length" = "max_len")

## blobtools

# get blobtools paths
blobtools_files <- list.files("results/blobtools/", recursive = T, pattern = "table.tsv", full.names = T)

# read as list
blobtools_list <- lapply(blobtools_files, function(x) {
  read.table(x, header = T, sep = "\t", 
             col.names = c("Index","Contig","GC","Length","Coverage",
                           "Superkingdom","Kingdom","Phylum","Class",
                           "Order","Family","Species"))
  })

# rbind list
blobtools_dat <- do.call("rbind", blobtools_list)

# add sample names
blobtools_dat <- cbind(ID = gsub("_circular$|_contig\\d*$","", blobtools_dat$Contig), 
      blobtools_dat)

## annotations

# get annotations file paths
annotations_files <- list.files("results/annotations/", recursive = T, pattern = "result.bed", full.names = T)

# read as list
annotations_list <- lapply(annotations_files, function(x) {
  read.table(x, header = F, sep = "\t", 
             col.names =c("identifiers", "start", "stop", "name", "score", "orientation") )
})

# rbind list
annotations_dat <- do.call("rbind", annotations_list)

# add sample names
annotations_dat <- cbind(ID = gsub("_circular$|_contig\\d*$","", annotations_dat$identifiers), 
                       annotations_dat)

# remove OH, OL and trn* annotations
annotations_dat <- annotations_dat %>% 
  filter(!str_detect(name, "OH|OL|trn")) %>% 
  arrange(ID, name)

# summarise annotations per contig
annotations_dat <- annotations_dat %>%
  group_by(identifiers) %>%
  summarise(genes = n(), 
            cox1 = sum(name == "cox1"), 
            genes_list = paste0(name, collapse = ",")) %>%
  rename("Contig" = "identifiers", 
         "N. genes" = "genes", 
         "Cox1" = "cox1", 
         "Genes list" = "genes_list")

## join 

samples_out <- left_join(samples_dat, fastp_dat, by = "ID") %>%
  rename("Reads raw" = "reads_raw", 
         "Reads QC" = "reads_qc") %>%
  left_join(., seqkit_dat, by = "ID")

samples_out [ is.na(samples_out) ] <- "NA"

contigs_out <- left_join(samples_dat, blobtools_dat, by = "ID") %>%
  left_join(., annotations_dat, by = "Contig") 

contigs_out [ is.na(contigs_out) ] <- "NA"

## write outputs

dir.create("results/summary/", showWarnings = FALSE)
write.table(samples_out, "results/summary/summary_samples_mqc.txt", sep = "\t", col.names = T, row.names = F, quote = F)
write.table(contigs_out, "results/summary/summary_contigs_mqc.txt", sep = "\t", col.names = T, row.names = F, quote = F)

