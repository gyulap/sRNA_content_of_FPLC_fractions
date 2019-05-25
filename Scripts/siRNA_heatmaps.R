setwd("./sRNA-seq")

sRNA_class = commandArgs(trailingOnly = T)
length = as.integer(substr(sRNA_class, 1, 2))

library(pheatmap)
library(RColorBrewer)

# Reading the count table

df = "Top_5000_sequences_miRBase_tasiRNA_TAIR10_annotated.txt"
a = read.table(file = df, sep = "\t", header = T, row.names = 1, check.names = F, fill = T, quote = "")

# Calculating the means of the replicates

a$Leaf_HMW = rowMeans(a[,16:17])
a$Leaf_LMW = rowMeans(a[,18:19])
a$Leaf_free = rowMeans(a[,20:21])
a$Flower_HMW = rowMeans(a[,8:9])
a$Flower_LMW = rowMeans(a[,10:11])
a$Flower_free = rowMeans(a[,12:13])

# Calculating the z-scores for the two tissues separately

siRNAs = function(length){
  b = log2(a[a$miRBase == "No hit" & !grepl("miRNA", a$TAIR10) & a$Length == length, 22:27] + 0.01)
  b$leaf_mean = rowMeans(b[,1:3])
  b$flower_mean = rowMeans(b[,4:6])
  b$leaf_sd = apply(b[,1:3], 1, sd)
  b$flower_sd = apply(b[,4:6], 1, sd)
  c = (b[,1:3] - b$leaf_mean) / b$leaf_sd
  d = (b[,4:6] - b$flower_mean) / b$flower_sd
  e = cbind(c, d)
  min = min(e[!is.na(e)])
  max = max(e[!is.na(e)])

  # Converting missing data, i.e. tissue-specific miRNAs to a value that is not anywhere else in the table

  e[is.na(e)] = -2

  # Making nicer column names

  colnames(e) = gsub("_", " ", colnames(e))

  return(e)
}

heatmap = function(length){
  png(filename = paste0(sRNA_class, ".png"), units = "cm" , width = 6, height = 23, res = 600, type = "cairo-png")
  siRNA = siRNAs(length)
  heat = pheatmap(siRNA,
         color = c("grey50", colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)),
         breaks = c(-2, seq(min - 0.001, max + 0.001, length.out = 100)),
         scale = "none",
         show_rownames = F,
         cellwidth = 18,
         treeheight_row = 0,
         treeheight_col = 0,
         cluster_cols = F,
         gaps_col = 3
         )
  print(heat)
  dev.off()

  table = a[rownames(siRNA[heat$tree_row$order,]), ]
  rownames(table) = gsub("T", "U", rownames(table))
  write.table(table, file = paste0(sRNA_class, ".txt"), sep = "\t", quote = F, row.names = T, col.names = T)
}

heatmap(length)

rm(list=ls())
