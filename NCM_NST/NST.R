library(NST)

data(tda)

comm <- read.delim('allsample_wet.txt', row.names = 1, sep = '\t', check.names = FALSE)


group <- read.delim('group_wet.txt', row.names = 1, sep = '\t', check.names = FALSE)




tnst <- tNST(comm = comm, group = group, dist.method = 'jaccard', null.model = 'PF', 
             rand = 1000, nworker = 1)

nst_group <- tnst$index.pair.grp
nst_group

write.table(nst_group, 'nst_group_fungal.txt', sep = '\t', row.names = FALSE, quote = FALSE)

library(ggpubr)

ggboxplot(data = nst_group, x = 'group', y = 'MST.ij.ruzicka', color = 'group') +
  stat_compare_means(method = 'wilcox.test', comparisons = list(c('a', 'b'))) +
  labs(y = 'Modified Stochasticity Ratio (MST)')

write.table(ggboxplot, 'data.txt', sep = '\t', row.names = FALSE, quote = FALSE)

