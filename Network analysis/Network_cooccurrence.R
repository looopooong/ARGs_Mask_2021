
library(Hmisc)
phylum <- read.delim('otu2.txt', row.name = 1, check.names = FALSE)
ARGs <- read.delim('ARGs_final.txt', row.name = 1, check.names = FALSE)
phylum_ARGs_corr <- rcorr(as.matrix(phylum), as.matrix(ARGs), type = 'spearman')
r <- phylum_ARGs_corr$r
p <- phylum_ARGs_corr$P
r <- r[colnames(phylum),colnames(ARGs)]
p <- p[colnames(phylum),colnames(ARGs)]
r[abs(r) < 0.8] <- 0
p <- p.adjust(p, method = 'BH') 
p[p>=0.001] <- -1
p[p<0.001 & p>=0] <- 1
p[p==-1] <- 0
z <- r * p
z1 <- phylum_ARGs_corr$r
z1[z1 != 0] <- 0
z1[rownames(z),colnames(z)] <- z
z1[colnames(z),rownames(z)] <- z

write.table(data.frame(z1, check.names = FALSE), 'phylum_ARGs_corr.matrix.txt', col.names = NA, sep = '\t', quote = FALSE)
library(igraph)

g <- graph.adjacency(z1, weighted = TRUE, mode = 'undirected')
g
g <- delete.vertices(g, names(degree(g)[degree(g) == 0]))
E(g)$correlation <- E(g)$weight
E(g)$weight <- abs(E(g)$weight)
plot(g)
write.graph(g, 'network.gml', format = 'gml')
write.graph(g, 'network.graphml', format = 'graphml')
