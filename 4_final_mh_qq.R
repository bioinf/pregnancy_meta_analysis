if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")s
if (!requireNamespace("RColorBrewer", quietly = TRUE)) install.packages("RColorBrewer")
if (!requireNamespace("CMplot", quietly = TRUE)) install.packages("CMplot")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")


library("ggplot2")
library("RColorBrewer")
library("CMplot")
library("dplyr")

prepare_snp <- function(FILE) {
  table <- read.table(FILE, header = T, sep = "\t")
  
  table$chr[table$chr == "X"] <- 23
  table$chr[table$chr == "Y"] <- 24
  
  snp <- table$rsid
  chr <- as.numeric(table$chr)
  pos <- as.numeric(table$pos)
  pval <- as.numeric(table$pval)
  final_table <- data.frame(snp, chr, pos, pval)
  return(final_table)
  
}

draw_mhplot_qq <- function(t, name, format, color, SNPs, genes, max_pval = 8) {
  CMplot(t,
         plot.type = "m", col = c("grey40", "grey70"), # chr colors
         highlight = SNPs,
         highlight.text = genes,
         highlight.col = c("#fb8072"),
         highlight.cex = 1, highlight.pch = c(16),
         LOG10 = TRUE, ylim = c(0, max_pval), # limits of log pval
         threshold = c(0.05 / nrow(t), 1e-4), # cut-offs of pval
         threshold.lty = c(1, 2), threshold.lwd = c(1, 1),
         threshold.col = c("black", "grey"), # threshold colors
         amplify = TRUE, chr.den.col = NULL,
         signal.col = c("#fb8072", "#b3de69"), # colors of significant
         signal.cex = c(1.5, 1.5), signal.pch = c(19, 19),
         file = format, # file format
         memo = name, # file postfix
         dpi = 30, file.output = TRUE, verbose = TRUE,
         width = 14, height = 5
  )
  CMplot(t,
         plot.type = "q", col = color, box = FALSE, file = format, memo = name, dpi = 30,
         conf.int = TRUE, conf.int.col = NULL, threshold.col = "red", threshold.lty = 2,
         file.output = TRUE, verbose = TRUE,
         width = 5, height = 3.5
  )
}

SNPS = list(c('rs13306561', 'rs35954793', 'rs10882398', 'rs167479', 'rs259983'),
            c('rs2208589'),
            c('rs58835482,rs796221113'),
            c('rs780094', 'rs9275373', 'rs10659211', 'rs10830963'))
GENES = list(c('MTHFR,CLCN6', 'FGF5', 'PLCE1', 'RGL3', 'ZNF831'),
             c('PREX1'),
             c('GDF15'),
             c('GCKR', 'MTCO3P1,HLA-*', 'TCF7L2', 'MTNR1B'))
FILES <- c('./data/f_special/maf_fg_I9_HYPTENSPREG_hg19lifted.tsv_.tsv',
           './data/f_special/maf_fg_O15_GESTAT_HYPERT_hg19lifted.tsv_.tsv',
           './data/f_special/maf_fg_O15_EXCESS_VOMIT_PREG_hg19lifted.tsv_.tsv',
           './data/f_special/maf_fg_GEST_DIABETES_hg19lifted.tsv_.tsv')
pheno_names = c('_FG_I9_HYPTENSPREG',
                '_FG_O15_GESTAT_HYPERT',
                '_FG_O15_EXCESS_VOMIT_PREG',
                '_FG_GEST_DIABETES')
colors <- c( "#fb8072", "#fdb462", "#a1cf53", "#80b1d3", "#bebada","#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5", "#ffed6f")
colors <- c("#ffed6f", "#bbf0af","#bbb6e3", "#d090d1")
for (i in 3:3){

  pheno_name <- pheno_names[i]
  color <- colors[i]
  
  t <- prepare_snp(FILES[i])
  
  max_pval <- max(-log10(t$pval))
  max_pval <- ceiling(max_pval)
  
  snps <- SNPS[[i]]
  genes <- GENES[[i]]
  
  draw_mhplot_qq(t,
                 pheno_name,
                 "pdf",
                 color,
                 snps,
                 genes,
                 max_pval)
}



SNPS = list(c('rs35954793', 'rs10882398', 'rs167479'),
            c('rs36090025', 'rs10830963'),
            c('rs2963457'))
GENES = list(c('FGF5', 'PLCE1', 'RGL3'),
             c('TCF7L2', 'MTNR1B'),
             c('EBF1'))
FILES <- c('./data/meta_special/extended_I9_HYPTENSPREG1.TBL',
           './data/meta_special/extended_GEST_DIABETES1.TBL',
           './data/meta_special/extended_O15_PRETERM1.TBL')
pheno_names = c('_MET_I9_HYPTENSPREG',
                '_MET_GEST_DIABETES',
                '_MET_O15_PRETERM')
colors <- c( "#fb8072", "#fdb462", "#a1cf53", "#80b1d3", "#bebada","#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5", "#ffed6f")
colors <- c("#ffed6f",  "#d090d1", "#699df0")
for (i in 1:3){
  
  pheno_name <- pheno_names[i]
  color <- colors[i]
  
  t <- prepare_snp(FILES[i])
  
  max_pval <- max(-log10(t$pval))
  max_pval <- ceiling(max_pval)
  
  snps <- SNPS[[i]]
  genes <- GENES[[i]]
  
  draw_mhplot_qq(t,
                 pheno_name,
                 "pdf",
                 color,
                 snps,
                 genes,
                 max_pval)
}






