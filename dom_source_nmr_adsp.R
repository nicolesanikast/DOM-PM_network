library(ggplot2);library(ggrepel);library(reshape2)

d.data  <- read.csv("~/Documents/Projects/version_control/DOM-PM-network/newdata/DOM_propNMR_integration.csv")
d.color  <- data.frame(type = levels(d.data$TYPE), col = topo.colors(length(levels(d.data$TYPE)), alpha = 1))
d.data <- d.data[d.data$INCLUDE == "yes", ]
d.data$TYPE <- as.character(d.data$TYPE)
d.data$C.H <- d.data$C/d.data$H
d.data$C.O <- d.data$C/d.data$O
d.data$col  <- sapply(1:nrow(d.data), function(x) d.color$col[d.color$type == d.data$TYPE[x]])
d.data.red <- d.data[ ,c("TYPE", "SOURCE", "col",
                         "Aliphatic.0..110ppm.", "Aroamtic.110..165.", "Carbonyl.165..220.", "aromatic.aliphatic", "aromatic.carbonyl")]
d.data.red <- d.data.red[complete.cases(d.data.red), ]
d.data.red <- d.data.red[!duplicated(d.data.red$SOURCE),]
d.pca  <- princomp(d.data.red[ ,c("Aliphatic.0..110ppm.", "Aroamtic.110..165.", "Carbonyl.165..220."
                                  )], cor =T)
summary(d.pca, loadings = T)


d.pca.summary <- data.frame(PC1 = jitter(d.pca$scores[, 1]), PC2 = d.pca$score[, 2], type = d.data.red$TYPE, 
                            color = d.data.red$col, source = d.data.red$SOURCE)

#analysis: how does the distance between same type humic substances is compared to different type humic substances?
inter.type.dist <- c()
intra.type.dist <- c()
dist.matrix <- data.frame(as.matrix(dist(d.pca.summary[ ,1:2], upper = T, diag = T)))

for(i in 1:(nrow(dist.matrix)-1)){
  current.type <- as.character(d.pca.summary$type[i])
  for(j in (i+1):nrow(dist.matrix)){#going over only the upper half of the matrix without the diagonal
    print(as.character(d.pca.summary$type[j]))
    if(as.character(d.pca.summary$type[j]) == current.type){
      intra.type.dist <- c(intra.type.dist, dist.matrix[i, j])
    }else{inter.type.dist <- c(inter.type.dist, dist.matrix[i, j])}
  }
}
#a histogram that summerizes the differences between the inter and intra distances:
hist(inter.type.dist, xlim = range(c(intra.type.dist, inter.type.dist)), 
     freq = F, ylim = c(0,1), col = "#F39C1250", border = "grey", xlab = "euclidean distance", main = "")
hist(intra.type.dist, add = T, xlim = range(c(intra.type.dist, inter.type.dist)), freq = F, col = "#45993C70", border = "grey")


#with lables:
ggplot(d.pca.summary) +
  geom_vline(xintercept = 0, size = 0.1) +
  geom_hline(yintercept = 0, size = 0.1) +
  geom_point(aes(PC1, PC2, color = type), size = 3)+#, color = 'grey') +
  geom_label_repel(
    aes(PC1, PC2, fill = type, label = source),
    fontface = 'bold', color = 'white',
    label.size = 1.5,
    size = 2,
    box.padding = unit(0.25, "lines"),
    point.padding = unit(0.5, "lines")
  ) +
  scale_fill_discrete(guide = guide_legend(ncol =  1, title = "DOM source")) +
  #geom_text(aes(d.pca.summary[, 1], d.pca.summary[, 2], label = ""), show.legend  = F) +
  theme_classic(base_size = 7) +
  theme(axis.line.x = element_line(color="black", size = 0.1),
        axis.line.y = element_line(color="black", size = 0.1))
ggsave("FigureS1_lables.pdf", width = 10, height = 7.5)

#without labels:

ggplot(d.pca.summary) +
  geom_vline(xintercept = 0, size = 0.1) +
  geom_hline(yintercept = 0, size = 0.1) +
  geom_point(aes(PC1, PC2, color = type), size = 3)+#, color = 'grey') +
  #geom_label_repel(
  # aes(PC1, PC2, fill = type, label = source),
  #  fontface = 'bold', color = 'white',
  # label.size = 1.5,
  #  size = 2,
  #box.padding = unit(0.25, "lines"),
  # point.padding = unit(0.5, "lines")
  #) +
  scale_fill_discrete(guide = guide_legend(ncol =  1, title = "DOM source")) +
  #geom_text(aes(d.pca.summary[, 1], d.pca.summary[, 2], label = ""), show.legend  = F) +
  theme_classic(base_size = 7) +
  theme(axis.line.x = element_line(color="black", size = 0.1),
        axis.line.y = element_line(color="black", size = 0.1))
ggsave("FigureS1_nolables.pdf", width = 10, height = 7.5)



#Figure iterpertation of the pca results: 
#It depicts the representative percentages of the different carbon content for the materials that fall within each quarter.

#all the materials that fall in the upper left region, exclusing those fall directly on the horizontal and vertical lines:
upper.left <- which(d.pca.summary$PC1 <0 & d.pca.summary$PC2 >0)
#all the materials that fall in the lower left region, exclusing those fall directly on the horizontal and vertical lines:
lower.left <- which(d.pca.summary$PC1 <0 & d.pca.summary$PC2 <0)
#all the materials that fall in the lower right region, exclusing those fall directly on the horizontal and vertical lines:
lower.right <- which(d.pca.summary$PC1 >0 & d.pca.summary$PC2 <0)
#all the materials that fall in the upper left region, exclusing those fall directly on the horizontal and vertical lines:
upper.right <- which(d.pca.summary$PC1 >0 & d.pca.summary$PC2 >0)
#list of all these regions:
regions <- list(upper.left, lower.left, lower.right, upper.right)
#obtain the mean and std of the carbonyl distributiono of each region
hist.data <- data.frame()
for(i in 1:length(regions)){
  t <- d.data.red[regions[[i]], c("Aliphatic.0..110ppm.", "Aroamtic.110..165.", "Carbonyl.165..220.")]
  t <- melt(t)
  t$panel <- i
  hist.data <- rbind(hist.data, t)
}
hist.data$panel <- factor(hist.data$panel, levels = c(1,4,2,3))#reorder the panels

#density plots:
ggplot(data = hist.data, aes(x = value, fill = variable)) +  
  geom_density(alpha = 0.7, ) + 
  facet_wrap(~panel) +
  #facet_wrap(~panel, scales='free_y') + #to have each panel have it's own y scale
  theme_classic(base_size = 7) + 
  theme(axis.line.x = element_line(color="black", size = 0.1),
        #axis.line.y = element_line(color="black", size = 0.1),
        #axis.ticks = element_blank(),
        #axis.text = element_blank(),
        axis.title = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) + 
  scale_fill_discrete(guide = guide_legend(ncol =  1, title = "carbon type"))
ggsave("FigureS1_b_option1.pdf", width = 10, height = 7.5)


#bar plot with error bars:
mean.agg.data <- aggregate(hist.data$value, by = list(hist.data$variable, hist.data$panel), FUN = "mean")
sd.agg.data <- aggregate(hist.data$value, by = list(hist.data$variable, hist.data$panel), FUN = "sd")
row.names.barplot <-mean.agg.data$Group.1
barplot.data <- data.frame(mean = mean.agg.data$x, sd = sd.agg.data$x)
barplot.data$variable <- mean.agg.data$Group.1
barplot.data$panel <- mean.agg.data$Group.2


