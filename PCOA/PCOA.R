##PCOA
otu <- read.delim('otu_table.txt', row.names = 1, stringsAsFactors = FALSE)
otu <- data.frame(t(otu))
group <- read.delim('group.txt', stringsAsFactors = FALSE)

library(vegan)

adonis_drug <- adonis(otu~drug, group, distance = 'bray', permutations = 999)
adonis_drug

bray_dis <- vegdist(otu, method = 'bray')
pcoa <- cmdscale(bray_dis, k = (nrow(otu) - 1), eig = TRUE)
pcoa_exp <- pcoa$eig/sum(pcoa$eig)
pcoa1 <- paste('PCoA axis1 :', round(100*pcoa_exp[1], 2), '%')
pcoa2 <- paste('PCoA axis2 :', round(100*pcoa_exp[2], 2), '%')
site <- data.frame(pcoa$point)[1:2]
site$sample <- rownames(site)
site <- merge(site, group, by = 'sample')
names(site)[2:3] <- c('pcoa1', 'pcoa2')
library(ggplot2)
library(ggrepel)

p <- ggplot(data = site) +
  geom_point(aes(x = pcoa1, y = pcoa2, color = drug), size = 2) + 
  geom_text_repel(aes(x = pcoa1, y = pcoa2, label = sample, color = drug), size = 2.5, 
                  box.padding = unit(0.3, 'lines'), show.legend = FALSE) +  
  scale_color_manual(limits = c('Before', 'After'), values = c('#D27FB2', '#764697')) +  
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(color = 'black'), legend.key = element_blank()) +
  labs(x = pcoa1, y = pcoa2, color = '')

p

group_average <- aggregate(cbind(pcoa1, pcoa2)~drug, data = site, FUN = mean)

p1 <- p +
  stat_ellipse(aes(x = pcoa1, y = pcoa2, color = drug), level = 0.95, linetype = 2, show.legend = FALSE) +  
  geom_point(data = group_average, aes(x = pcoa1, y = pcoa2, color = drug), size = 5, show.legend = FALSE)  

p1 <- p1 +
  annotate('text', label = 'PERMANOVA', x = 0.18, y = 0.15, size = 3) +
  annotate('text', label = sprintf('italic(P) == %.3f', adonis_drug$aov.tab[1,6]), x = 0.18, y = 0.13, size = 3, parse = TRUE)

p1