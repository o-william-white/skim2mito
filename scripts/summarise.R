#!/bin/Rscript
args = commandArgs(trailingOnly=TRUE)

# get output dir and output_path from arguments
output_dir <- args[1]
output_path <- args[2]

library(tidyverse)


# get list of samples
sample_list <- list.files(path = paste0(output_dir, "/blobtools/"))


dat_list <- lapply(sample_list, function(sample) {
  
  # test 
  # sample <- "Bath_M299691"
  
  # check blobtools file exists
  if(file.exists(paste0(output_dir, "/blobtools/", sample, "/table.tsv"))) {
    
    # read blobtools output 
    blobtools <- read.table(paste0(output_dir, "/blobtools/", sample, "/table.tsv"), header = T, sep = "\t")
    
    # add column for sample name
    blobtools <- cbind(sample = sample, blobtools)
    
    # create lineage as comma separated string
    blobtools <- mutate(blobtools, lineage = paste(bestsumorder_superkingdom,
                                                   bestsumorder_kingdom,
                                                   bestsumorder_phylum,
                                                   bestsumorder_class,
                                                   bestsumorder_order,
                                                   bestsumorder_family,
                                                   bestsumorder_species, sep = ",")) 
    
    # select columns 
    blobtools <-  select(blobtools, 
                         index,
                         identifiers,
                         gc,
                         length,
                         ends_with("cov"), 
                         lineage)
    
    # one sequence 
    if(nrow(blobtools) == 1) {
      
      df <- read.table(paste0(output_dir, "/mitos/", sample, "/result.bed"), 
                       col.names = c("identifiers", "start", "stop", "name", "score", "orientation"))
      
      df <- filter(df, !str_detect(name, "OH|OL|trn"))
      
      df <- arrange(df, name)
      
      gene_names <- paste0(df$name, collapse = ",")
      
      blobtools$annotations <- gene_names
    
    # more than one sequence
    } else {
    
      list_bed <- sapply(blobtools$index, function(x) {
        
        df <- read.table(paste0(output_dir, "/mitos/", sample, "/", x, "/result.bed"), 
                         col.names = c("identifiers", "start", "stop", "name", "score", "orientation"))
        
        df <- filter(df, !str_detect(name, "OH|OL|trn"))
        
        df <- arrange(df, name)
        
        gene_names <- paste0(df$name, collapse = ",")
        
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


