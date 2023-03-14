#!/bin/Rscript
args = commandArgs(trailingOnly=TRUE)

# get output dir and output_path from arguments
output_dir <- args[1]
output_type <- args[2] # mitos or barrnap
output_path <- args[3]

library(tidyverse)

# get list of samples
sample_list <- list.files(path = paste0(output_dir, "/blobtools/"))

dat_list <- lapply(sample_list, function(sample) {
  
  # test 
  # sample <- "Bath_M299691"
  
  # check blobtools file exists
  if(file.exists(paste0(output_dir, "/blobtools/", sample, "/table.tsv"))) {
    
    # set column names 
    column_names <- c("index","identifiers","gc","length","cov","superkingdom","kingdom","phylum","class","order","family","species")
    
    # read blobtools output 
    blobtools <- read.table(paste0(output_dir, "/blobtools/", sample, "/table.tsv"), header = T, sep = "\t", col.names = column_names)
    
    # add column for sample name
    blobtools <- cbind(sample = sample, blobtools)
    
    # create lineage as comma separated string
    blobtools <- mutate(blobtools, lineage = paste(superkingdom,
                                                   kingdom,
                                                   phylum,
                                                   class,
                                                   order,
                                                   family,
                                                   species, sep = ",")) 
    
    # select columns 
    blobtools <-  select(blobtools, 
                         index,
                         identifiers,
                         gc,
                         length,
                         cov, 
                         lineage)

    # one sequence 
    if(nrow(blobtools) == 1) {

      if(output_type == "mitos") {

        df <- read.table(paste0(output_dir, "/annotations/", sample, "/result.bed"),
                       col.names = c("identifiers", "start", "stop", "name", "score", "orientation"))

        df <- filter(df, !str_detect(name, "OH|OL|trn"))

        df <- arrange(df, name)

        gene_names <- paste0(df$name, collapse = ",")

      } else {

        if(output_type == "barrnap"){

          df <- read.table(paste0(output_dir, "/annotations/", sample, "/result.gff"), 
                  sep = "\t", col.names = c("sequence", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"))

          df <- arrange(df, attribute)

          df <- mutate(df, gene = gsub("Name=", "", str_split_i(attribute, ";", 1)))

          gene_names <- paste0(df$gene, collapse = ",")

        }

      }
        
      blobtools$annotations <- gene_names
    
    # more than one sequence
    } else {
    
      list_bed <- sapply(blobtools$index, function(x) {
        
        if(output_type == "mitos") {

          df <- read.table(paste0(output_dir, "/annotations/", sample, "/", x, "/result.bed"), 
                         col.names = c("identifiers", "start", "stop", "name", "score", "orientation"))

          df <- filter(df, !str_detect(name, "OH|OL|trn"))

          df <- arrange(df, name)

          gene_names <- paste0(df$name, collapse = ",")
          
        } else {

          if(output_type == "barrnap"){

            df <- read.table(paste0(output_dir, "/annotations/", sample, "/result.gff"),
                           sep = "\t", col.names = c("sequence", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"))

            df <- arrange(df, attribute)

            df <- mutate(df, gene = gsub("Name=", "", str_split_i(attribute, ";", 1)))

            gene_names <- paste0(df$gene, collapse = ",")

          }

        }
        
        gene_names
        
      })
      
      df2 <- data.frame(index = blobtools$index, annotations = list_bed)
      
      blobtools <- left_join(blobtools, df2, by = "index")
    
      } 
  
    # remove sample name from column names
    colnames(blobtools) <- gsub(paste0(sample, "_"), "", colnames(blobtools))
    
    # return blobtools output
    return(blobtools)
    
  }
  
})

# rbind to dataframe
dat <- do.call("rbind", dat_list)

# write output
write.table(dat, output_path, sep = "\t", quote = F, col.names = T, row.names = F)


