setwd("~/Dropbox/Arabidopsis_FPLC")

library(pheatmap)
library(RColorBrewer)

timestamp = format(Sys.time(), "%Y.%m.%d_%H.%M.%S")

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
#a = a[order(a$mean, decreasing = T),]
# Calculating the z-scores for the two tissues separately
b = log2(a[rowMeans(a[, c(23:25, 19:21)]) > 1, c(23:25, 19:21)] + 1)
#b = log2(a[substr(a$Sequence, 1, 1) == "T" & rowMeans(a[, c(23:25, 19:21)]) > 0, c(23:25, 19:21)] + 0.01)
#b = log2(a[rowMeans(a[, c(23:25, 19:21)]) > 0, c(23:25, 19:21)] + 0.01)
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
#tiff(filename = paste0("miRNA_min1_sep", timestamp, ".tiff"), units = "cm" , width = 9, height = 23, compression = "lzw", res = 600, type = "cairo")
#ord = pheatmap(e, clustering_distance_rows = "manhattan", clustering_method = "complete")
#f = e[ord$tree_row$order[c(1:7,53:23,8:22,54:94)],]
#dend = as.dendrogram(hclust(dist(e, method = "manhattan"), method = "complete"))
#hc = as.hclust(dend)

heat = pheatmap(e,
         color = c("grey50", colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)),
         breaks = c(-2, seq(min - 0.001, max + 0.001, length.out = 100)),
         cellheight = 6,
         cellwidth = 18,
         treeheight_row = 0,
         treeheight_col = 0,
         gaps_col = 3,
         #cutree_rows = 6,
         #gaps_row = c(3, 7, 53, 74),
         scale = "none",
         cluster_cols = F,
         cluster_rows = F,
         #clustering_distance_rows = "manhattan",
         #clustering_method = "complete",
         show_colnames = F,
         fontsize_row = 6
)
png(filename = paste0("miRNA_min1_", timestamp, ".png"), units = "cm" , width = 9, height = 21, res = 600, type = "cairo-png")
#png(filename = paste0("U_miRNA_", timestamp, ".png"), units = "cm" , width = 9, height = 32, res = 600, type = "cairo-png")
print(heat)
dev.off()

table = a[rownames(e[heat$tree_row$order,]), c(1, 10:17, 2:9)]
#table = a[rownames(e), c(1, 10:17, 2:9)]
table$Sequence = gsub("T", "U", table$Sequence)
write.table(table, paste0("miRNA_heatmap_", timestamp, ".txt"), sep="\t", quote=F, row.names=T, col.names=T)
rm(list=ls())