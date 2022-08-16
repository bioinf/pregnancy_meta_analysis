library(ggplot2)

d <- read.csv("out_data/feature.csv", stringsAsFactor=F)

ds <- d
ds$cat_new <- d$priority
ds$p2_new <- d$p2
ds$ymin <- ds$rg-1.96*ds$se
ds$ymin <- sapply(ds$ymin, function(x) ifelse(x < -1, -1, x))
ds$ymax <- ds$rg+1.96*ds$se
ds$ymax <- sapply(ds$ymax, function(x) ifelse(x > 1, 1, x))
ds$Trait <- ordered(ds$p1,levels=c("HP", "GD", "PTB"))
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


pdf("img/final_gen_cor.pdf",width=11,height=5)

p1 <- ggplot(dsO, aes(x=pos, y=rg, ymin=ymin, ymax=ymax)) +
  geom_point(aes(color=Trait), size=4) +
  geom_errorbar(width = .9, aes(color=Trait)) + 
  theme_bw() +
  coord_flip() + 
  scale_x_continuous('',breaks=seq(1.5,x_max + 0.5, 4), labels=unique(dsO$p2_new)) +
  scale_y_continuous('Genetic correlation', 
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
  scale_colour_manual(values=c("#f0d935",  "#d090d1", "#699df0"), 
  # scale_colour_manual(values=c("#86d9b8","#de9090", "#a07fe3"),
  ) + 
  expand_limits(y = c(-0.5, 1))

p1

dev.off()

