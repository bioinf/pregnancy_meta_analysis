library(ggplot2)
library(tidyverse)


# ==========================================================
# Prepare meta data
# ==========================================================

d <- read.csv("out_data/meta_feature.csv", stringsAsFactor=F)

levels <- c("HP", "GD", "PTB")

ds <- d
ds$cat_new <- d$priority
ds$p2_new <- d$p2
ds$ymin <- ds$rg-1.96*ds$se
ds$ymin <- sapply(ds$ymin, function(x) ifelse(x < -1, -1, x))
ds$ymax <- ds$rg+1.96*ds$se
ds$ymax <- sapply(ds$ymax, function(x) ifelse(x > 1, 1, x))
ds$Trait <- ordered(ds$p1,levels=levels)
ds$cat_pregnancy_trait <- ordered(ds$cat_new,levels=c("2", "1"))
tt <- aggregate(ds$rg,list(ds$p2_new),mean)
tt <- tt[order(tt$x),]
tt <- tt[c(1,2,5,3,4,6),]

ds$p2_pregnancy_trait <- ordered(ds$p2_new,levels=tt$Group.1)
ds$gsign <- factor(ds$significant)
dsO <- ds[rev(order(ds$cat_pregnancy_trait,ds$p2_pregnancy_trait,ds$Trait)),]

x_max <- ceiling(length(ds$p2) * 4 / 3)
x <- 1:x_max
dsO$pos <- x[x%%4!=0]

y_min <- floor(min(dsO$ymin)*10)/10
y_max <- ceiling(max(dsO$ymax)*10)/10




# ==========================================================
# Draw forrest plot for meta data
# ==========================================================

pdf("img/meta_gen_cor.pdf",width=11,height=5)

p1 <- ggplot(dsO, aes(x=pos, y=rg, ymin=ymin, ymax=ymax)) +
  # geom_hline(yintercept=0.05/nrow(d), linetype="dashed") + # se vuoi una linea verticale nel plot, se no togli
  geom_point(aes(color=Trait), size=4) +
  geom_errorbar(width = .9, aes(color=Trait)) + # questo usa ymin e ymax
  theme_bw() +
  coord_flip() + # mette tutto orizzontale
  scale_x_continuous('',breaks=seq(1.5,x_max + 0.5, 4), labels=unique(dsO$p2_new)) +
  scale_y_continuous('Genetic correlation', 
                     # limits=c(-0.5,1), 
                     breaks=c(seq(y_min,y_max,0.2)), 
                     labels=c(as.character(round(seq(y_min,y_max,0.2),1)))) +
  geom_text(aes(label=gsub('e-0*', ' %*% 10^-', prettyNum(p, digits=2))), 
            hjust=-.25, 
            vjust=-0.4, 
            size=3.6, 
            parse=TRUE) + 
  geom_text(aes(label=gsign), 
            col="black", 
            vjust=+0.8, 
            size=6) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.title=element_text(size=12), 
        legend.text=element_text(size=12))  + 
  scale_colour_manual(values=c("#f0d935",  "#d090d1", "#699df0"),  #"#86d9b8","#de9090", "#a07fe3"
  # scale_colour_manual(values=c("#86d9b8","#de9090", "#a07fe3"),  #
                      # guide="none" # можно добавить легенду, убрав эту строку
  ) + 
  expand_limits(y = c(-0.5, 1))

p1

dev.off()

# ==========================================================
# Draw heatmap for meta data
# ==========================================================


dat <- dsO
dat$lci <- dat$ymin
dat$uci <- dat$ymax
dat$signif <- case_when(dat$p > 0.05 ~ "NS", 
                        dat$p < 0.05 & dat$significant == "" ~ "Nominal", 
                        dat$significant != "" ~ "Significant")
dat$signif <- fct_relevel(dat$signif, "NS", "Nominal", "Significant")
dat$exposure.name <- dat$p1
dat$outcome.name <- dat$p2
dat$stars <- dat$significant

dat$exposure.name <- factor(dat$exposure.name, levels = levels)


pdf("img/meta_gen_cor_heatmap.pdf",width=5,height=12)

p2 <-ggplot(dat) +
  geom_raster(aes(x = exposure.name, y = outcome.name, fill = rg)) +
  geom_text(data = dat, size = 5, aes(label = stars, x = exposure.name, y = outcome.name)) +
  scale_fill_gradient2(low="steelblue", high="firebrick", mid = "white", na.value = "grey75", name = "rg", limits = c(-1,1)) +
  geom_vline(xintercept=seq(0.5, 40.5, 1),color="white") +
  geom_hline(yintercept=seq(0.5, 11.5, 1),color="white") +
  coord_equal() +
  theme_classic() +
  theme(legend.position = 'right', 
        legend.key.height = unit(1, "line"),
        axis.text.x = element_text(angle = 35, hjust = 0),
        legend.text = element_text(hjust = 1.5), 
        text = element_text(size=8),
        title = element_text(size=8),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank()) +
  scale_x_discrete(position = "top")

p2

dev.off()


#===============================================================================


# ==========================================================
# Prepare FG data
# ==========================================================


d_fg <- read.csv("out_data_fg/feature_supp.csv", stringsAsFactor=F)
FG_levels <- c("HP", 'GH', "EV", "GD")

ds_fg <- d_fg
# ds_fg$cat_new <- d$priority
ds_fg$p2_new <- d_fg$p2
ds_fg$ymin <- ds_fg$rg-1.96*ds_fg$se
ds_fg$ymin <- sapply(ds_fg$ymin, function(x) ifelse(x < -1, -1, x))
ds_fg$ymax <- ds_fg$rg+1.96*ds_fg$se
ds_fg$ymax <- sapply(ds_fg$ymax, function(x) ifelse(x > 1, 1, x))
ds_fg$rg <- sapply(ds_fg$rg, function(x) ifelse(x > 1, 1, x))
ds_fg$Trait <- ordered(ds_fg$p1,levels=FG_levels)
# ds_fg$cat_pregnancy_trait <- ordered(ds_fg$cat_new,levels=c("2", "1"))
tt_fg <- aggregate(ds_fg$rg,list(ds_fg$p2_new),mean)
tt_fg <- tt_fg[order(tt_fg$x),]
# tt <- tt[c(1,2,5,3,4,6),]

ds_fg$p2_pregnancy_trait <- ordered(ds_fg$p2_new,levels=tt_fg$Group.1)
ds_fg$gsign <- factor(ds_fg$significant)
# dsO <- ds[rev(order(ds_fg$cat_pregnancy_trait,ds_fg$p2_pregnancy_trait,ds_fg$Trait)),]
dsO_fg <- ds_fg[rev(order(ds_fg$p2_pregnancy_trait,ds_fg$Trait)),]

x_max_fg <- ceiling(length(ds_fg$p2) * 5 / 4)
x_fg <- 1:x_max_fg
dsO_fg$pos <- x_fg[x_fg%%5!=0]

y_min_fg <- floor(min(dsO_fg$ymin)*10)/10
y_max_fg <- ceiling(max(dsO_fg$ymax)*10)/10



# ==========================================================
# Draw forrest plot for FG data
# ==========================================================


pdf("img/gen_cor_FG_supp.pdf",width=30,height=25)

p1_fg <- ggplot(dsO_fg, aes(x=pos, y=rg, ymin=ymin, ymax=ymax)) +
  # geom_hline(yintercept=0.05/nrow(d), linetype="dashed") + # se vuoi una linea verticale nel plot, se no togli
  geom_point(aes(color=Trait), size=4) +
  geom_errorbar(width = .9, aes(color=Trait)) + # questo usa ymin e ymax
  theme_bw() +
  coord_flip() + # mette tutto orizzontale
  scale_x_continuous('',breaks=seq(1.5,x_max_fg +0.5, 5), labels=unique(dsO_fg$p2_new)) +
  scale_y_continuous('Genetic correlation', 
                     # limits=c(-0.5,1), 
                     breaks=c(seq(y_min_fg,y_max_fg,0.2)), 
                     labels=c(as.character(round(seq(y_min_fg,y_max_fg,0.2),1)))) +
  geom_text(aes(label=gsub('e-0*', ' %*% 10^-', prettyNum(p, digits=2))), 
            hjust=-.25, 
            vjust=-0.4, 
            size=3.6, 
            parse=TRUE) + 
  geom_text(aes(label=gsign), 
            col="black", 
            vjust=+0.8, 
            size=9) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.title=element_text(size=15), 
        legend.text=element_text(size=15))  + 
  scale_colour_manual(values=c("#f0d935", "#9aed87","#bbb6e3", "#d090d1"),  #"#86d9b8","#de9090", "#a07fe3"
                      # scale_colour_manual(values=c("#86d9b8","#de9090", "#a07fe3"),  #
                      # guide="none" # можно добавить легенду, убрав эту строку
  ) + 
  expand_limits(y = c(-0.5, 1))

p1_fg

dev.off()

# ==========================================================
# Draw heatmap for FG data
# ==========================================================


dat_fg <- dsO_fg
dat_fg$lci <- dat_fg$ymin
dat_fg$uci <- dat_fg$ymax
dat_fg$signif <- case_when(dat_fg$p > 0.05 ~ "NS", 
                           dat_fg$p < 0.05 & dat_fg$significant == "" ~ "Nominal", 
                           dat_fg$significant != "" ~ "Significant")
dat_fg$signif <- fct_relevel(dat_fg$signif, "NS", "Nominal", "Significant")
dat_fg$exposure.name <- dat_fg$p1
dat_fg$outcome.name <- dat_fg$p2
dat_fg$stars <- dat_fg$significant

dat_fg$exposure.name <- factor(dat_fg$exposure.name, levels = FG_levels)

pdf("img/heatmap_gen_cor_FG_supp.pdf",width=8,height=14)

p2_fg <- ggplot(dat_fg) +
  geom_raster(aes(x = exposure.name, y = outcome.name, fill = rg)) +
  geom_text(data = dat_fg, size = 10, aes(label = stars, x = exposure.name, y = outcome.name)) +
  scale_fill_gradient2(low="steelblue", high="firebrick", mid = "white", na.value = "grey75", name = "rg", limits = c(-1,1)) +
  geom_vline(xintercept=seq(0.5, 40.5, 1),color="white") +
  geom_hline(yintercept=seq(0.5, 11.5, 1),color="white") +
  coord_equal() +
  theme_classic() +
  theme(legend.position = 'right', 
        legend.key.height = unit(1, "line"),
        axis.text.x = element_text(angle = 35, hjust = 0),
        legend.text = element_text(hjust = 1.5), 
        text = element_text(size=15),
        title = element_text(size=15),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank()) +
  scale_x_discrete(position = "top")

p2_fg

dev.off()


