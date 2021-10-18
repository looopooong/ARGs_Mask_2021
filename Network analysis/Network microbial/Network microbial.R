
library(Hmisc)

genus <- read.delim('genus_table.txt', row.name = 1, check.names = FALSE)


genus <- genus[which(rowSums(genus) >= 0.005), ] 

genus1 <- genus
genus1[genus1>0] <- 1
genus <- genus[which(rowSums(genus1) >= 5), ]    


genus_corr <- rcorr(t(genus), type = 'spearman')

r <- genus_corr$r
r[abs(r) < 0.7] <- 0
p <- genus_corr$P
p <- p.adjust(p, method = 'BH')    
p[p>=0.05] <- -1
p[p<0.05 & p>=0] <- 1
p[p==-1] <- 0
z <- r * p
diag(z) <- 0    
head(z)[1:6,1:6]
write.table(data.frame(z, check.names = FALSE), 'genus_corr.matrix.txt', col.names = NA, sep = '\t', quote = FALSE)
library(igraph)
g <- graph.adjacency(z, weighted = TRUE, mode = 'undirected')
g
g <- simplify(g)
g <- delete.vertices(g, names(degree(g)[degree(g) == 0]))
E(g)$correlation <- E(g)$weight
E(g)$weight <- abs(E(g)$weight)
tax <- read.delim('genus_taxonomy.txt', row.name = 1, check.names = FALSE, stringsAsFactors = FALSE)
tax <- tax[as.character(V(g)$name), ]

V(g)$kingdom <- tax$kingdom
V(g)$phylum <- tax$phylum
V(g)$class <- tax$class
V(g)$order <- tax$order
V(g)$family <- tax$family
V(g)$genus <- tax$genus
g
plot(g)
adj_matrix <- as.matrix(get.adjacency(g, attr = 'correlation'))
write.table(data.frame(adj_matrix, check.names = FALSE), 'network.adj_matrix.txt', col.names = NA, sep = '\t', quote = FALSE)
edge <- data.frame(as_edgelist(g))

edge_list <- data.frame(
  source = edge[[1]],
  target = edge[[2]],
  weight = E(g)$weight,
  correlation = E(g)$correlation
)
head(edge_list)

write.table(edge_list, 'network.edge_list.txt', sep = '\t', row.names = FALSE, quote = FALSE)
node_list <- data.frame(
  label = names(V(g)),
  kingdom = V(g)$kingdom,
  phylum = V(g)$phylum,
  class = V(g)$class,
  order = V(g)$order,
  family = V(g)$family,
  genus = V(g)$genus
)
head(node_list)

write.table(node_list, 'network.node_list.txt', sep = '\t', row.names = FALSE, quote = FALSE)
write.graph(g, 'network.graphml', format = 'graphml')
write.graph(g, 'network.gml', format = 'gml')
