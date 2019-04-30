setwd("~/Dropbox/Arabidopsis_FPLC")

library(ggfortify)

timestamp = format(Sys.time(), "%Y.%m.%d_%H.%M.%S")

a = read.table(file="Top_5000_sequences2.txt", sep = "\t", header = T, row.names = 1, check.names = F)
a = as.data.frame(t(a))
a = as.data.frame(a[, apply(a, 2, var) != 0])
#a$Replicates = as.factor(rep(c("Leaf HMW", "Leaf LMW", "Leaf free", "Flower HMW", "Flower LMW", "Flower free"), 1, each=2))
a$Fractions = as.factor(rep(c("HMW", "LMW", "unbound"), 2, each=2))
pca = prcomp(log2(a[,1:ncol(a)-1] + 0.01), scale. = F, center = T)
p = autoplot(pca, data = a, colour = "Fractions", label = F, scale = 0, x = 1, y = 2)
p = p + geom_text(vjust = -1, size = 3, label = rownames(a), show.legend = F, aes(colour = factor(Fractions)))
p = p + xlim(-220, 220) + ylim(-220, 220)
png(filename = paste0("PCA_plot_", timestamp, ".png"), units = "cm" , width = 16, height = 12.3, res = 600, type = "cairo-png")
plot(p)
dev.off()