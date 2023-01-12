library(ggtree)

args = commandArgs(trailingOnly=TRUE)

# run as 
# plot_tree input.tree output.png

# read tree
tre <- read.tree(args[1])

# create ggtree object
p <- ggtree(tre) 

# write png
png(args[2], res=600, height = 12, width = 12, units = "in")
p + 
  geom_tiplab() + 
  geom_nodelab(aes(x=branch), vjust=-0.75) + 
  xlim(c(0,1.5*max(p$data$x)))
dev.off()

