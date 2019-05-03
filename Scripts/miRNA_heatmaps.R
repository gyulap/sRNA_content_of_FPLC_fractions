setwd("./sRNA-seq")

sRNA_class = commandArgs(trailingOnly = T)

library(pheatmap)
library(RColorBrewer)

# Reading the count table

a = read.table(file="miRBase_mature_sequences_norm_count_table.txt", sep = "\t", header = T, row.names = 1, check.names = F)

# Calculating the means of the replicates

a$Flower_input = rowMeans(a[,2:3])
a$Flower_HMW = rowMeans(a[,4:5])
a$Flower_LMW = rowMeans(a[,6:7])
a$Flower_unbound = rowMeans(a[,8:9])
a$Leaf_input = rowMeans(a[,10:11])
a$Leaf_HMW = rowMeans(a[,12:13])
a$Leaf_LMW = rowMeans(a[,14:15])
a$Leaf_unbound = rowMeans(a[,16:17])
a$mean = rowMeans(a[,c(4:9, 12:17)])

if (sRNA_class == "miRNA"){
  b = log2(a[rowMeans(a[, c(23:25, 19:21)]) > 1, c(23:25, 19:21)] + 1)
  cluster_rows = T,
  cldist = "manhattan"
  cutree_rows = 6
  filename = "Fig3a_miRNAs_min1RPM"
  height = 21
} else if (sRNA_class == "5p-U") {
  b = log2(a[substr(a$Sequence, 1, 1) == "T" & rowMeans(a[, c(23:25, 19:21)]) > 0, c(23:25, 19:21)] + 1)
  cluster_rows = T,
  cldist = "euclidean"
  cutree_rows = NA
  filename = "FigS6a_miRNAs_5p-U"
  height = 32
} else if (sRNA_class == "abundance") {
  cluster_rows = F,
  cldist = NA
  a = a[order(a$mean, decreasing = T),]
  b = log2(a[rowMeans(a[, c(23:25, 19:21)]) > 1, c(23:25, 19:21)] + 1)
  cutree_rows = NA
  filename = "FigS6b_miRNAs_by_abundance"
  height = 21
}

# Calculating the z-scores for the two tissues separately

b$leaf_mean = rowMeans(b[,1:3])
b$flower_mean = rowMeans(b[,4:6])
b$leaf_sd = apply(b[,1:3], 1, sd)
b$flower_sd = apply(b[,4:6], 1, sd)
c = (b[,1:3]-b$leaf_mean)/b$leaf_sd
d = (b[,4:6]-b$flower_mean)/b$flower_sd
e = cbind(c, d)
min = min(e[!is.na(e)])
max = max(e[!is.na(e)])

# Converting missing data, i.e. tissue-specific miRNAs to a value that is not anywhere else in the table

e[is.na(e)] = -2

# Making nicer column names

colnames(e) = gsub("_", " ", colnames(e))

# Drawing miRNA heatmap

heat = pheatmap(e,
         color = c("grey50", colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)),
         breaks = c(-2, seq(min - 0.001, max + 0.001, length.out = 100)),
         cellheight = 6,
         cellwidth = 18,
         treeheight_row = 0,
         treeheight_col = 0,
         gaps_col = 3,
         cutree_rows = cutree_rows,
         scale = "none",
         cluster_cols = F,
         cluster_rows = cluster_rows,
         clustering_distance_rows = cldist,
         clustering_method = "complete",
         show_colnames = F,
         fontsize_row = 6
)
png(filename = paste0(filename, ".png"), units = "cm" , width = 9, height = height, res = 600, type = "cairo-png")
print(heat)
dev.off()

table = a[rownames(e[heat$tree_row$order,]), c(1, 10:17, 2:9)]
table$Sequence = gsub("T", "U", table$Sequence)
write.table(table, paste0(filename, ".txt"), sep="\t", quote=F, row.names=T, col.names=T)
rm(list=ls())
