library(ggtree)

args = commandArgs(trailingOnly=TRUE)

# run as 
# plot_tree input.tree output.png height width

# read tree
tre <- read.tree(args[1])

# create ggtree object
p <- ggtree(tre) 

# write png
png(args[2], res=600, height = as.integer(args[3]), width = as.integer(args[4]), units = "cm")
p + 
  geom_tiplab(align = TRUE) + 
  geom_nodelab(aes(x=branch), vjust=-0.75) + 
  xlim(c(0,2.5*max(p$data$x)))
dev.off()

