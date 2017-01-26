#the following file analyzes the data from Deonarine et at. 2011 in order to inspect the similarity between the DOM analyzed in this paper:
data <- read.csv("../../NOM_model/deonarine_2011.csv", header = T, na.strings = "NA")
data.red <- data[-8, c(4:9, 11:17)]
col <- 
plot(data[, c("mw", "aromatic")], col = data[, "color"])
identify(data[, c("mw", "SUVA")])
pca <- princomp(data.red, cor = T)
summary(pca, loadings = T)
plot(pca$scores[ ,1], pca$scores[ ,2])

data_aromatic <- data[ , c("category", "aromatic")]
data_aromatic <- data_aromatic[complete.cases(data_aromatic), ]
data_aromatic <- data_aromatic[order(data_aromatic$category), ]
colors <- sample(colors(distinct = T), size = length(unique(data_aromatic$category)), replace = F)
d.col <- data.frame(col = colors, category = as.character(unique(data_aromatic$category)))
d.col$col <- as.character(d.col$col)
d.col$category <- as.character(d.col$category)
for(i in 1:nrow(data_aromatic)){data_aromatic$col[i] <- d.col$col[d.col$category == data_aromatic$category[i]]}
plot(data_aromatic$aromatic~data_aromatic$category, pch = 16, cex = 1, las = 2, xlab = "")
plot(data_aromatic$aromatic, pch = 16, cex = 1, las = 2, xlab = "", col = data_aromatic$col)

pca <- princomp(data_aromatic[complete.cases(data_aromatic), -c(1, 2, 5)], cor = T)
summary(pca)
plot(pca$scores[ ,1:2])
