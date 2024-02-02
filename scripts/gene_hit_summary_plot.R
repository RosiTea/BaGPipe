require(ggplot2)
require(ggrepel)
library(ggrepel)

gene_hits = read.table("gene_hits.txt", stringsAsFactors=FALSE, header=TRUE, sep="\t")

pdf("Gene_Hit_Plot.pdf", width = 12, height = 9) 

ggplot(gene_hits, aes(x=avg_beta, y=maxp, colour=avg_maf, size=hits, label=gene)) +
   geom_point(alpha=0.5) +
   geom_text_repel(aes(size=60), show.legend = FALSE, colour='black') +
   scale_size("Number of unitigs/k-mers", range=c(1,10)) +
   scale_colour_gradient('Average MAF') +
   theme_bw(base_size=14) +
   ggtitle("Penicillin resistance") +
   xlab("Average effect size") +
   ylab("Maximum -log10(p-value)")

dev.off()