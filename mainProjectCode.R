load.required.packages()

# Read in the input data --------------------------------------------------
##Uncomment the lines below if the database was changed and you'd like to recompile it (the file "newdata.csv 
#corresponds to the xlsx supporting information file of the publications), else the existing database stored 
#as an R object in the working directory will be loaded as "g.data.nom.type" variable

  nom.data.column <- "NOM.type.detailed"#DON'T comment this out! this variable is needed throughout the script
  # all.data <- read.csv("../newdata/newdata.csv", header =TRUE, strip.white = TRUE)#read in the entire database
  # #columns for the bipartite network
  # d.data <- all.data[c("ENP", "NOM.by.location", "Experiment.overall", "NOM.type", "ENP.type", "reference", "year",
  #  "Full.paper.title", "DOI", "NOM.type.extended", nom.data.column)]#extract only relevant columns
  # d.data <- d.data[d.data$ENP != "", ]#eliminate empty lines
  # d.data$DOI <- tolower(d.data$DOI)#convert DOI to lower case
  # d.data.nom.type <- d.data[d.data[, nom.data.column] != "", ]#double check there are no empty lines...
  # save(d.data.nom.type, file = "dataBase.Rdata")#output the dataframe into an R object. This will be stored in the working directory

load("dataBase.Rdata")#load the R object with the database from the working directory, the data is stored in the variable "d.data.nom.type"

# DOM and PM categorization -----------------------------------------------
#summary of DOM types in the database:

group1.dom <- c("humic substance", "humic acid", "fulvic acid", "total DOC",
                "hydrophilic acid", "hydrophilic neutral acid", "transphilic acid", "leachate", "hydrophobic neutral")
water.samples <- c("natural water sample", "industrial wastewater", "STP effluent", "STP influent")
group2.dom <- unique(as.character(d.data.nom.type$NOM.type))
group2.dom <- group2.dom[!group2.dom %in% group1.dom]                                    
group2.dom <- group2.dom[!group2.dom %in% water.samples]                                    
#group1 DOM:
n.group1DOM <- sum(d.data.nom.type$NOM.type%in%group1.dom)#the number of experiments that employ group-1 DOM
n.group1DOM
#write.csv(file = "~/Downloads/group1DOM.csv",unique(d.data.nom.type$NOM.by.location[d.data.nom.type$NOM.type%in%group1.dom]))
#d.data.nom.type$Full.paper.title[d.data.nom.type$NOM.by.location == "forest floor leachate"]
length(unique(as.character(d.data.nom.type$NOM.type.detailed[d.data.nom.type$NOM.type%in%group1.dom])))#number of group1-DOM types
n.group1DOM/nrow(d.data.nom.type)#percentage of group-1 DOM
#water samples:
n.waterSampleDOM <- sum(d.data.nom.type$NOM.type%in%water.samples)#the number of water samples
n.waterSampleDOM
n.waterSampleDOM/nrow(d.data.nom.type)
length(unique(as.character(d.data.nom.type$NOM.type.detailed[d.data.nom.type$NOM.type%in%water.samples])))
#group-2 DOM:
n.group2DOM <- sum(d.data.nom.type$NOM.type%in%group2.dom)#the number of experiments that employ group-2 DOM
n.group2DOM
length(unique(as.character(d.data.nom.type$NOM.type.detailed[d.data.nom.type$NOM.type%in%group2.dom])))#number of group2-DOM types
n.group2DOM/nrow(d.data.nom.type)#percentage of group-2 DOM


#PCA analysis of humic substances:
  # d.data.pca  <- read.csv("~/Documents/Projects/version_control/DOM-PM-network/newdata/DOM_propNMR_integration.csv")#read in the csv file
  # d.color  <- data.frame(type = levels(d.data.pca$TYPE), col = topo.colors(length(levels(d.data.pca$TYPE)), alpha = 1))#create color scheme
  # d.data.pca <- d.data.pca[d.data.pca$INCLUDE == "yes", ]#include the conssitent data
  # d.data.pca$TYPE <- as.character(d.data.pca$TYPE)
  # d.data.pca$col  <- sapply(1:nrow(d.data.pca), function(x) d.color$col[d.color$type == d.data.pca$TYPE[x]])
  # d.pca.red <- d.data.pca[ ,c("TYPE", "SOURCE", "col",
  #                          "Aliphatic.0..110ppm.", "Aroamtic.110..165.", "Carbonyl.165..220.")]
  # d.pca.red <- d.pca.red[complete.cases(d.pca.red), ]#keep only datapoints that have values in the above categories
  # d.pca.red <- d.pca.red[!duplicated(d.pca.red$SOURCE),]#remove duplicates of exact same substances that were reported in different studies
  # save(d.pca.red, file = "domProperties.Rdata")#save the file to the working directory
load("domProperties.Rdata")#read in DOM NMR data file for PCA analysis
nrow(d.pca.red)#number of humic substances in this analysis

d.pca.output  <- princomp(d.pca.red[ ,c("Aliphatic.0..110ppm.", "Aroamtic.110..165.", "Carbonyl.165..220."
)], cor =T)#perform PCA on the scaled data

# # #perform robust PCA, using a covariace matrix insensitive to outliers, it showed almost identical results
# # library(MASS)
# d.pca.output  <- princomp(d.pca.red[ ,c("Aliphatic.0..110ppm.", "Aroamtic.110..165.", "Carbonyl.165..220."
# )], cor =T, covmat = cov.rob(d.pca.red[ ,c("Aliphatic.0..110ppm.", "Aroamtic.110..165.", "Carbonyl.165..220."
# )]))#perform PCA on the scaled data

summary(d.pca.output, loadings = T)#PCA output including the variable loadings
d.pca.summary <- data.frame(PC1 = jitter(d.pca.output$scores[, 1]), PC2 = d.pca.output$score[, 2], type = d.pca.red$TYPE, 
                            color = d.pca.red$col, source = d.pca.red$SOURCE)#summary data frame

# Figure 1 - carbon content distribution in PC dimensions ---------------
#Figure iterpertation of the pca results: 
#It depicts the representative percentages of the different carbon content for the materials that fall within each quarter
#of the PC1-PC2 space. It excludes the material that lay exactly on the horizontal and vertical lines

#all the materials that fall in the upper left region, exclusing those fall directly on the horizontal and vertical lines:
upper.left <- which(d.pca.summary$PC1 <0 & d.pca.summary$PC2 >0)
length(upper.left)#number of matrials in this quadrant
#what is the most prevelnt material in this quadrant?
top.source.upper.left <- percentage.in.quadrant(d.pca.summary, upper.left)
#all the materials that fall in the lower left region, exclusing those fall directly on the horizontal and vertical lines:
lower.left <- which(d.pca.summary$PC1 <0 & d.pca.summary$PC2 <0)
length(lower.left)#number of matrials in this quadrant
top.source.lower.left <- percentage.in.quadrant(d.pca.summary, lower.left)

#all the materials that fall in the lower right region, exclusing those fall directly on the horizontal and vertical lines:
lower.right <- which(d.pca.summary$PC1 >0 & d.pca.summary$PC2 <0)
length(lower.right)#number of matrials in this quadrant
top.source.lower.right <- percentage.in.quadrant(d.pca.summary, lower.right)

#all the materials that fall in the upper left region, exclusing those fall directly on the horizontal and vertical lines:
upper.right <- which(d.pca.summary$PC1 >0 & d.pca.summary$PC2 >0)
length(upper.right)#number of matrials in this quadrant
top.source.upper.right <- percentage.in.quadrant(d.pca.summary, upper.right)

#list of all these regions:
regions <- list(upper.left, lower.left, lower.right, upper.right)
#obtain the mean and std of the carbonyl distributiono of each region
hist.data <- data.frame()
for(i in 1:length(regions)){
  t <- d.pca.red[regions[[i]], c("Aliphatic.0..110ppm.", "Aroamtic.110..165.", "Carbonyl.165..220.")]
  t <- melt(t)
  t$panel <- i
  hist.data <- rbind(hist.data, t)
}
hist.data$panel <- factor(hist.data$panel, levels = c(1,4,2,3))#reorder the panels

#density plots:
ggplot(data = hist.data, aes(x = value, fill = variable)) +  
  geom_density(alpha = 0.7) + 
  facet_wrap(~panel) +
  #facet_wrap(~panel, scales='free_y') + #to have each panel have it's own y scale
  theme_classic(base_size = 7) + 
  theme(axis.line.x = element_line(color="black", size = 0.1),
        #axis.line.y = element_line(color="black", size = 0.1),
        #axis.ticks = element_blank(),
        #axis.text = element_blank(),
        #axis.title = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.background = element_blank(),
        legend.text = element_text(size  = 12),
        legend.title = element_text(size = 14),
        legend.position = "top",
        #legend.direction = "vertical",
        strip.text.x = element_blank()) + 
  annotate("text", label = c("PC1 < 0, PC2 > 0","PC1 > 0, PC2 > 0", "PC1 < 0, PC2 < 0", "PC1 > 0, PC2 < 0"), x =50, y = 0.15, size = 5, fontface =2) + 
  annotate("text", label = c(top.source.upper.left, top.source.upper.right, top.source.lower.left, top.source.lower.right), x =50, y = 0.10, size = 4.5) +
  annotate("text", label = c(paste(as.character(length(upper.left)), "samples", sep =  " "),
                             paste(as.character(length(upper.right)), "samples", sep =  " "), 
                             paste(as.character(length(lower.left)), "samples", sep =  " "),
                             paste(as.character(length(lower.right)), "samples", sep =  " ")), x =50, y = 0.07, size = 3.5) +
  annotate("text", label = c("a","b", "c", "d"), x =0, y = 0.15, size = 5, fontface =2) + 
  scale_fill_discrete(guide = guide_legend(ncol =  1, title = "carbon type"), labels = c("aliphatic", "aromatic", "carbonyl")) +
  labs(x = "% Carbon", y = "Probability")
 
ggsave("Figure1.pdf", width = 7.5, height = 7.5)

# Organizing the experimental data in a network ---------------------------
#nom.data.column  <- "NOM.type.detailed"#what NOM column should be used for the analysis?
g.data.nom.type <- adj.to.graph(d.data.nom.type[c("ENP", nom.data.column)], d.data.nom.type$ENP)#Converting the dataset into a graph instance
V(g.data.nom.type)$degree.lables.code <- label.node.degree(g.data.nom.type)#label nodes of highest degree, DOM with letters and PM with numbers

#save the network structure to an external file
#write.graph(graph = g.data.nom.type, format = "pajek", file = "gDataNOMtype.net")
write.graph(graph = g.data.nom.type, format = "graphml", file = "gDataNOMtype.graphml")#this file format is easilt read by gephi
# Network topology --------------------------------------------------------
#Basic properties (more properties such as degree assortativity and diameter are calculated in the code section: Supporting Information comparison to random network)
length(d.data.nom.type$DOI[!d.data.nom.type$DOI == ""])#numpber of publications
sum(E(g.data.nom.type)$weight)#number of experiments (sum of links' weights)
n.pm <- sum(V(g.data.nom.type)$type); n.pm#number of PM types (number of PM nodes)
n.dom <- sum(V(g.data.nom.type)$type == FALSE); n.dom#number of DOM types (number of DOM nodes)
n.comb.stud <- ecount(g.data.nom.type); n.comb.stud#number of unique DOM-PM combinations studied
density  <- n.comb.stud/(n.pm*n.dom); signif(density, digits = 2)#network density
theor.density <- sum(E(g.data.nom.type)$weight)/(n.pm*n.dom);signif(theor.density, digits = 2)#theoretical/potential density
n.pm*n.dom
max(E(g.data.nom.type)$weight)#the max number of times a given DOM-PM combination was studied
min(E(g.data.nom.type)$weight)#the min number of times a given combination was studied
E(g.data.nom.type)[which.max(E(g.data.nom.type)$weight)]#DOM-PM combinations studied often
E(g.data.nom.type)[E(g.data.nom.type)$weight == max(E(g.data.nom.type)$weight)]#DOM-PM combinations studied often
E(g.data.nom.type)[E(g.data.nom.type)$weight == 1]#DOM-PM combinations studied less often

# Figure 2 - network overview plus hists of degree and weight distributions ------------------
#par(xpd = NA)
m  <- rbind(
  c(0.15, 0.7, 0.3, 1),#for the network,
  c(0, 0.15, 0.1, 0.4),#for the PM legend
  c(0, 0.15, 0.4, 0.7),#for the DOM legend
  c(0.7, 0.8, 0.4, 0.7),#for the degree distribution
  c(0.85, 0.95, 0.4, 0.7))#for the link weight distribution
pdf("Figure2.pdf", width = 14, height = 7)#the figure will be saved under this name in the working directory
split.screen(m)
screen(1)
par(mar = rep(0.1, 4))
plot(g.data.nom.type, vertex.frame.color = V(g.data.nom.type)$frame.color, vertex.size = log(degree(g.data.nom.type))+3,
     vertex.shape = V(g.data.nom.type)$shape, edge.width = E(g.data.nom.type)$width/4,
     vertex.color = V(g.data.nom.type)$color,
     vertex.label = V(g.data.nom.type)$degree.lables.code, vertex.label.cex = 1,
     vertex.label.color = V(g.data.nom.type)$label.color, edge.curved = 0.1, main="")
legend(-1.3,0.95,pch = c(21,22), col = c("dodgerblue4","grey"), pt.bg = c("light blue","lightsalmon"), text.col = c("dodgerblue4","black"),
       legend = c("PM","DOM"), bty = "n", pt.cex = 2)
mtext("a", side = 3, line = -2, adj = 0.025, cex = 1, font = 2)#subfigure numbering
screen(2)
#legend for the PM
legend("top", pch = V(g.data.nom.type)$degree.lables.code[V(g.data.nom.type)$degree.lables.code != "" & V(g.data.nom.type)$type == TRUE],
       legend = V(g.data.nom.type)$name[V(g.data.nom.type)$degree.lables.code != "" & V(g.data.nom.type)$type == TRUE], col = "black", cex = 0.7, y.intersp = 0.15, bty = "n")
screen(3)
#legend for the DOM
legend("top", pch = V(g.data.nom.type)$degree.lables.code[V(g.data.nom.type)$degree.lables.code != "" & V(g.data.nom.type)$type == FALSE],
       legend = V(g.data.nom.type)$name[V(g.data.nom.type)$degree.lables.code != "" & V(g.data.nom.type)$type == FALSE], col = "black",cex = 0.7, y.intersp = 0.15, bty = "n")
screen(4)
par(mar = rep(0.1, 4), cex = 0.7)
hist(degree(g.data.nom.type), main = "", border = NA, col = "grey", xlim = c(0, 50), ann =FALSE, xaxt="n", yaxt = "n", ylim = c(0, 200))
axis(1, at = seq(0, 50, by = 10))#x axis ticks
axis(2, at = seq(0, 200,by = 100))#y axis ticks
mtext("Count", side = 2, line = 2.2, adj = 0.5, cex = 0.6)#y axis annotation
mtext("Node degree", side = 1, line = 2.2, adj = 0.5, cex = 0.6)#x axis annotation
mtext("b", side = 3, line = -0.01, adj = 0.01, cex = 1, font = 2)#subfigure numbering
screen(5)
par(new=TRUE)
par(mar = rep(0.1, 4), cex = 0.7)
hist(E(g.data.nom.type)$weight,main="", border = NA, col = "grey70", ann =FALSE, xaxt="n", yaxt = "n")
axis(1, at = seq(1, 14, by = 3))#x axis ticks
axis(2, at = seq(0, 400, by = 200),labels = seq(0, 400, by = 200))#y axis ticks
mtext("Count", side = 2, line = 2.2, adj = 0.5, cex = 0.6)#y axis annotation
mtext("Link weight", side = 1, line = 2.2, adj = 0.5, cex = 0.6)#x axis annotation
mtext("c", side = 3, line = -0.01, adj = 0.01, cex = 1, font = 2)#subfigure numbering
close.screen(all.screens = TRUE)
dev.off()

# Network evolution from 1990-2015 ----------------------------------------
##In order to simulate the high- and low-diversity networks, a citation network is created where the DOM-PM combinations studied in each publication
#are assigned as attributes to the respective node (that correspond to the DOI of that publication). Uncomment the lines below if you'd like to 
#recreate the citation network (but it requires the references list of each publication in a seperate file, these files should be located in the
#folder: citation_network). Else, the already required citation network will be uploaded from the working directory as: "g.cit.only" variable.

  #Creating the database for the network, where all origin (citing paper) is on the first column and targets (cited papers)
  #are on the second column
  # d.data.cit <- data.frame()#initializing an empty data frame
  # doi.exper.relevant  <- unique(d.data.nom.type$DOI[d.data.nom.type$DOI != ""])
  # all.files <- list.files(path = "../citation_network", pattern = "^10\\..+?[^~]$")#read in all files, excluding the copies that end with "~"
  # all.files <- sapply(all.files, function(x) paste("../citation_network", x, sep = "/"))#appending the file name to the path so it can be accessed from the current working directory
  # for(i in 1:length(all.files)){
  #   file <- read.delim(all.files[i], header = FALSE, quote = "")
  #   #parent reference details:
  #   origin.DOI <- tolower(gsub("^\\s+|\\s+$", "", unlist(strsplit(as.character(file[1, ]), split = " - "))[2]))#obtaining the parent reference DOI and removing heading white spaces and making all lower case
  #   origin.title <- sub("^\\s+", "", unlist(strsplit(as.character(file[2, ]), split = " - "))[2])#obtaining the parent reference title and removing heading white spaces
  #   #make a doi list:
  #   list.doi <- sapply(1:nrow(file), function(x) {if(grepl("^DO.*?10.*", as.character(file[x, ]))){tolower(gsub("^\\s+|\\s+$", "", unlist(strsplit(as.character(file[x, ]), "DO\\s+-\\s+"))[2]))#obtaining the DOI and removing heading and trailing white spaces...
  #   }else{NA}})
  #   list.doi <- list.doi[!is.na(list.doi)]#remove NA entries
  #   #print(length(list.doi))
  #   d.data.cit <- rbind(d.data.cit, data.frame(rep(origin.DOI, length(list.doi)), list.doi))
  #   #print(all.files[i])
  # }
  # colnames(d.data.cit) <- c("origin","target")
  # #Building the citation network, directed, network: for only the DOIs present in the experiemntal network
  # g.cit <- graph.data.frame(d.data.cit, directed = TRUE)
  # #extract the DOIs that exist in the citation network:
  # doi.exper.relevant = as.character(doi.exper.relevant)[as.character(doi.exper.relevant)%in%V(g.cit)$name]
  # g.cit.only <- induced_subgraph(g.cit, vids = doi.exper.relevant)
  # #attaching to each publication node the list of PM and DOM types it studied
  # for(doi in doi.exper.relevant){
  #   curr.reference <- as.character(unique(d.data.nom.type$reference[d.data.nom.type$DOI == as.character(doi)]))#taking the refernec name for the given doi, and later extract the materials of that reference name, this is done because of the structure of the database where each material combination is has the reference
  #   #(mendely key identifier but for all experiemnts from the same reference there is only one DOI).
  #   curr.comb <- d.data.nom.type[d.data.nom.type$reference == curr.reference[1],c("ENP",nom.data.column)]#there could be more than oen combinations used, to have an NOM type extended instead of NOM type
  #   V(g.cit.only)$material.employ[V(g.cit.only)$name == as.character(doi)] <- list(curr.comb)
  # }
  # save(g.cit.only, file = "gcit.Rdata")
  # #creates a list of all possible DOM-PM combinations and save it to the working directory as Rdata file, uncomment to rerun if the database as changed,
  # #else, it will be loaded as "all.comb" variable
  # all.comb <- c()
  # for(i in unique(d.data.nom.type$ENP)){
  #   all.comb <- rbind(all.comb,data.frame(rep(i,length(unique(d.data.nom.type[,nom.data.column]))),unique(d.data.nom.type[,nom.data.column])))
  # }
  # all.comb <- unique(all.comb)#accumulate all the unique tested combinations
  # colnames(all.comb) <- c("ENP",nom.data.column)
  # all.comb$comb <- paste(all.comb[[1]],all.comb[[2]],sep = "-")
  # save(all.comb,file = "allComb.Rdata")

load("gcit.Rdata")#the data is stored in the variable: "g.cit.only"
load("allComb.Rdata")#the data is stored in the variable: "all.comb"
###
# The following code calculates the combinations diversity index over time for the empirical network and for the high- and low-diversity simulated networks
###
start.year <- 1990#the year from which the simulations/analysis start
#empirical network:
g.empir.year <- adj.to.graph(d.data.nom.type[d.data.nom.type$year<= start.year,c("ENP",nom.data.column)], names = unique(d.data.nom.type$ENP))#the empirical network at the year the analysis starts
empir.diversity.time <- diversity.over.time(d.data.nom.type[ , c("ENP",nom.data.column,"year")], start.year)#diversity over time for the detailed DOM column
empir.diversity.DOMbroad <- diversity.over.time(d.data.nom.type[ , c("ENP","NOM.type","year")], start.year)#diversity values when considering the broad definition of DOM
#for the DOM by location, make every non standard humic substance unique:

d.data.nom.type$NOM.by.location  <- as.character(d.data.nom.type$NOM.by.location)
modify.to.unique  <- c("sea water", "ground water", "surface water", "STP influent", "STP effluent", "river water",
                       "swamp total DOC", "soil humic acid", "EPS", "drinking water total DOC", "peat soil total DOC")#a list of DOM names that should be made unique under the assumption that each nonstandard DOM is unique
#go over all the DOM and if are int he ebove list from different DOI then paste the given DOI as a unique identifier, since all other are already specificed by location in this column
for(i in unique(as.character(d.data.nom.type$reference))){
  dom.to.modify  <- d.data.nom.type[d.data.nom.type$reference == i, "NOM.by.location"]#all the DOM that reported in a given reference
  for(j in unique(as.character(dom.to.modify))){#for each DOM specified for the given reference
    if(j %in% modify.to.unique){
      d.data.nom.type[as.character(d.data.nom.type[ , "NOM.by.location"]) == j & as.character(d.data.nom.type$reference) == i, "NOM.by.location"]  <- paste(j, i, sep = "__")
}}}
length(d.data.nom.type$NOM.by.location)#How many DOM types are there, when non standard humic substances are treated as unique by location?

empir.diversity.DOMbyloc <- diversity.over.time(d.data.nom.type[ , c("ENP","NOM.by.location","year")], start.year)#diversity values when considering the DOM definition by location (for humic substances)

#Linear regression of the empir.diversity values:
d.diversity  <- data.frame(diversity = empir.diversity.time, time = 1:length(empir.diversity.time))
lm.diversity  <- lm(diversity ~ time, data = d.diversity)
summary(lm.diversity)
tsdisplay(resid(lm.diversity))
tsdisplay(resid(ar.burg(ts(resid(lm.diversity)), order.max = 1)))#can an AR(1) process explain the dependency structrue in the residuals? yes, they form white noise
corStruc  <- corARMA(form = ~time, p = 1)#define the correlation structure of the residuals
gls.data <- gls(model = diversity ~ time, data = d.diversity, correlation = corStruc)#fit a generalized least square using an AR(1) correlation structure of the residuals.
signif(confint(gls.data), digits = 3)#the 95% confidence interval
signif(coef(gls.data)[1], digits = 3)#the intercept in 3 significant digits
signif(coef(gls.data)[2], digits = 3)#the slope in 3 significant digits
  
#simulated networks:
##one simulation: obtaining network structures
g.high.list <- change.graph.trend(g.cit.only, start.year, d.data.nom.type, 1, seed = 1)#high repetition from the 90th (low diversity)
g.high.rep  <- diversity.over.time(g.high.list,start.year)[length(empir.diversity.time)]#diversity of the high repetition (low diversity) network at the last time step
g.high <- adj.to.graph(g.high.list[ , 1:2], names = unique(g.high.list$ENP))#the high repetition network (low diversity) at the end of the simulation (the latest year of published experiment in the database)
#g.high <- g.high.list[[1]]#the high repetition network (low diversity) at the end of the simulation (the latest year of published experiment in the database)
g.low.list <- change.graph.trend(g.cit.only, start.year, d.data.nom.type, 0, seed = 1)#low repetition from the 90th (high diversity)
g.low.rep <- diversity.over.time(g.low.list,start.year)[length(empir.diversity.time)]#diversity of low repetition (high diversity) at the last time step
g.low <- adj.to.graph(g.low.list[ , 1:2], names = unique(g.low.list$ENP))#the low repetition network (high diversity) at the end of the simulation (the latest year of published experiment in the database)

#   #Uncomment in case you'd like to simulate again the high and low diversity networks, else the ready csv files will be read from the working directories for ploting the results
#   s.sample <- 1000#number of simulated network for high and low repetition (low and high diversities, respectively)
#   #printing the results into file as the loop runs, also the dataframe will be need to be transposed once it is loaded again into R:
#   sapply(774:1000, function(x) {print(x); g.low.list <- change.graph.trend(g.cit.only, start.year, d.data.nom.type, 0, seed = x);
#                              write.table(t(diversity.over.time(g.low.list, start.year)),
#                              file = "high_diversity.csv", append = TRUE, sep=",", col.names = FALSE, row.names = FALSE)})
#   sapply(97:1000,function(x) {print(x);g.high.list <- change.graph.trend(g.cit.only, start.year, d.data.nom.type, 1, seed = x);
#                               write.table(t(diversity.over.time(g.high.list, start.year)),
#                               file = "low_diversity.csv", append = TRUE, sep=",", col.names = FALSE, row.names = FALSE)})

high.diversity  <- t(read.csv("high_diversity.csv", header = FALSE))#read in the high diversity simulated values
low.diversity  <- t(read.csv("low_diversity.csv", header = FALSE))#read in the low diversity simulated values

#Diversity of ENP and DOM individualy (materials' diverstity):
pm.diversity <- sapply(1990:2015,function(x) length(unique(d.data.nom.type$ENP[d.data.nom.type$year <= x]))/nrow(d.data.nom.type[d.data.nom.type$year <= x,]))
dom.diversity <- sapply(1990:2015,function(x) length(unique(d.data.nom.type[d.data.nom.type$year <= x, nom.data.column]))/nrow(d.data.nom.type[d.data.nom.type$year <= x,]))
pm.diversity.type <- sapply(1990:2015,function(x) length(unique(d.data.nom.type$ENP.type[d.data.nom.type$year <= x]))/nrow(d.data.nom.type[d.data.nom.type$year <= x,]))
dom.diversity.type <- sapply(1990:2015,function(x) length(unique(d.data.nom.type$NOM.type[d.data.nom.type$year <= x]))/nrow(d.data.nom.type[d.data.nom.type$year <= x,]))
#the number of experiments in the last decrease period 2012-2015:
nrow(d.data.nom.type[d.data.nom.type$year >= 2012, ])
signif(nrow(d.data.nom.type[d.data.nom.type$year >= 2012, ])/nrow(d.data.nom.type), digits = 2)#the fraction of experiments in the last period

# Systemic trends in the experimental diversity - variables prep. ---------
#all contingency tables for the analysis are be between the years 1990-2015
d.data.nom.type.1990 <- d.data.nom.type[d.data.nom.type$year >= 1990, ]
required.years <- 1990:2015

# Systemic trends in experimental designs - DOM types temporal trends ---------------
#unaggregated DOM:
dom.type2.freq <- table(d.data.nom.type.1990$NOM.type.detailed, d.data.nom.type.1990$year)
dom.type2.freq <- dom.type2.freq[!rownames(dom.type2.freq) == "", ]#remove empty rows
dom.type2.freq <- table.check(dom.type2.freq, required.years)
#order the DOM types according to group-1 DOM, group-2 DOM and water samples:
dom.type.match <- unique(d.data.nom.type[ , c("NOM.type.detailed", "NOM.type")])#match each dom to its type
order.dom.type2.freq <- unlist(sapply(rownames(dom.type2.freq), function(x) {
  t <- as.character(dom.type.match$NOM.type[as.character(dom.type.match$NOM.type.detailed) == x][1])
  if(t%in%group1.dom){return(1)}
  if(t%in%water.samples){return(0)}
  if(t%in%group2.dom){return(2)}
  }))
order.dom.type2.freq <- sort(order.dom.type2.freq)
o <- match(names(order.dom.type2.freq), rownames(dom.type2.freq))

dom.type2.freq <- dom.type2.freq[o, ]
#color based on group1, group2 and water samples:
#group1.dom.colors <- colorRampPalette(c("blue", "darkorchid1"))(sum(order.dom.type2.freq == 1))
group1.dom.colors <- colorRampPalette(c("red2", "gray27"))(sum(order.dom.type2.freq == 1))
#for group 1 color differently things coming out of river and differently things coming out of soil and sigma aldrich humic acid:
group1.dom.river.hs <- which(rownames(dom.type2.freq)[order.dom.type2.freq == 1]%in% c("river humic acid", "river fulvic acid", "river total DOM", "river humic substance"))
group1.dom.river.hs.col <- colorRampPalette(c("lightskyblue", "lightskyblue2"))(length(group1.dom.river.hs))
for(i in 1:length(group1.dom.river.hs)){
  group1.dom.colors[group1.dom.river.hs[i]] <- group1.dom.river.hs.col[i]
}

#group2.dom.colors <- colorRampPalette(c("green", "firebrick1"))(sum(order.dom.type2.freq == 2))
group2.dom.colors <- colorRampPalette(c("seagreen3", "seagreen3"))(sum(order.dom.type2.freq == 2))
#water.samples.colors <- colorRampPalette(c("orange", "yellow"))(sum(order.dom.type2.freq == 0))
water.samples.colors <- colorRampPalette(c("orange1", "orange1"))(sum(order.dom.type2.freq == 0))#all water samples in one color
dom.type2.freq.col <- c(water.samples.colors, group1.dom.colors, group2.dom.colors)
#dom.type2.freq.col <- rainbow(nrow(dom.type2.freq), s = 0.7)
dom.max.use.2012.2015 <- sort(apply(dom.type2.freq, 1, sum), decreasing = T)[1:5]#names of DOM most used between 2012:2015
dom.max.use.2012.2015.col <- unlist(sapply(names(dom.max.use.2012.2015), function(x) dom.type2.freq.col[rownames(dom.type2.freq) == x]))
dom.max.use.2012.2015.names <- c(names(dom.max.use.2012.2015), "water samples","group-2 DOM")#adding gruop 2 DOM
dom.max.use.2012.2015.col <- c(dom.max.use.2012.2015.col, "orange1","seagreen3")


#aggregate the DOM data into larger groups:
saccharides <- c("polysaccharide", "gum arabic", "disaccharide")
d.data.nom.type.1990$NOM.type.3 <- sapply(as.character(d.data.nom.type.1990$NOM.type), function(x) {
  if(x%in%water.samples){return("water sample")}
  if(x%in%saccharides){return("saccharides")}
  if(x%in%group1.dom){return("group-1 DOM")}
  else{return(x)}
})

dom.type3.freq <- table(d.data.nom.type.1990$NOM.type.3, d.data.nom.type.1990$year)
#anacdotal DOM to remove from figure:
#dom.type3.freq <- dom.type3.freq[!row.names(dom.type3.freq)%in%c("chelating agent", 
 #                                                                "surfactant", "vitamine", "fatty acid", "glycolipid", "polyphenol", "polymer"), ]
dom.types3.freq.col <- rainbow(length(rownames(dom.type3.freq)), s = 0.7)
dom.type3.freq <- table.check(dom.type3.freq, required.years)

#
#only humic substances:
d.data.nom.type.1990$NOM.type2 <- sapply(as.character(d.data.nom.type.1990$NOM.type), 
                                         function(x) if(x %in% group1.dom){return("group-1 DOM")}else{return(x)})

#devision only of the humic substances:
d.data.humic.sub <- d.data.nom.type.1990[d.data.nom.type.1990$NOM.type2 == "group-1 DOM", ]
d.data.humic.sub$NOM.type.detailed <- as.character(d.data.humic.sub$NOM.type.detailed)
#hs.type.freq <- table(d.data.humic.sub$NOM.type.detailed, d.data.humic.sub$year)

#divide the materials according to PCA analysis:
fresh.water.fulvic.acid <- c("ground water fulvic acid", "lake fulvic acid", "river fulvic acid", "small stream fulvic acid",
                             "surface water fulvic acid", "ground water fulvic acid", "wetland fulvic acid")
fresh.water.humic.acid <- c("lake humic acid", "river humic acid", "ground water humic acid", "wetland humic acid")
fresh.water.nom <- c("swamp total DOM", "river total DOM", "lake total DOM",
                     "small stream total DOM", "drinking water total DOM", "surface water DOM",
                     "wetland total DOM")
fresh.water.humic.substance <- c("river humic substance", "wetland humic substance", "lake humic substance")
fresh.water.nonhumic.fraction <- c("lake hydrophilic acid", "wetland hydrophilic neutral", "wetland hydrophilic acid", "hydrophobic neutral",
                                   "river hydrophilic acid")
d.data.humic.sub$NOM.type.agg <- sapply(d.data.humic.sub$NOM.type.detailed, function(x){
  if(x%in%fresh.water.fulvic.acid){return("fresh water fulvic acid")}
  if(x%in%fresh.water.humic.acid){return("fresh water humic acid")}
  if(x%in%fresh.water.nom){return("fresh water NOM")}
  if(x%in%fresh.water.humic.substance){return("humic substance")}
  if(x%in%fresh.water.nonhumic.fraction){return("non humic fraction")}
  #if(x%in%wetland.nom){return("wetland NOM")}
  else{return(x)}
})
#anacdotal DOM to remove:
dom.types.remove <- c("unspecific humic acid", "unspecific fulvic acid", "mixed fulvic and humic aldrich acids"
                      ,"synthetic humic acid", "synthetic fulvic acid", "plant humic acid", "DOM coating")#anacdotal DOM to remove
hs.type.freq <- table(d.data.humic.sub$NOM.type.agg, d.data.humic.sub$year)
hs.type.freq <- hs.type.freq[!rownames(hs.type.freq)%in%dom.types.remove, ]

hs.type.freq <- table.check(hs.type.freq, required.years)#make sure all years 1990-2015 are represented

hs.type.freq.col <- c("blueviolet", "darkorchid1", "darkorchid4", "gray36", "gold", "darkslategray1", "darkslategray3", "lightskyblue",
                      "khaki2", "khaki3", "lightsalmon", "lightsalmon3", "mistyrose1", "mistyrose3", "goldenrod3", "olivedrab3",
                      "navajowhite3", "navajowhite4", "palegoldenrod", "palevioletred1", "palevioletred3", "pink3", "darkorange")
#hs.type.colors <- data.frame(type = row.names(hs.type.freq), col = hs.type.colors.list)

#now break down the identity of the fresh water humic substances:
fresh.water.hs.freq <- table(d.data.humic.sub$NOM.type.detailed[d.data.humic.sub$NOM.type.agg%in%c("fresh water fulvic acid", "fresh water humic acid",
                                                                                                   "fresh water NOM")], 
                             d.data.humic.sub$year[d.data.humic.sub$NOM.type.agg%in%c("fresh water fulvic acid", "fresh water humic acid",
                                                                                      "fresh water NOM")])
#anacdotal fresh water DOM to remove:
#fresh.water.hs.freq <- fresh.water.hs.freq[!rownames(fresh.water.hs.freq)%in% c("surface water DOM", "surface water fulvic acid"), ]
fresh.water.hs.freq <- table.check(fresh.water.hs.freq, required.years)#make sure all required columns are present
fresh.water.hs.freq.col <- rainbow(nrow(fresh.water.hs.freq), s = 0.7)

# Systemic trends in experimental design - PM types temporal trends ----------------
#
#temporal trend in the PM types used by core material:
pm.core.type.freq <- table(d.data.nom.type.1990$ENP.type, d.data.nom.type.1990$year)
pm.core.type.freq <- pm.core.type.freq[!row.names(pm.core.type.freq) == "", ]#remove rows with name ""
#type.order <- rownames(pm.core.type.freq)#before removing anacdotal instances since they will show up in the detailed description
pm.core.type.freq <- pm.core.type.freq[!row.names(pm.core.type.freq)%in%c("boron", "silicon"), ]#remove anacdotal PM types
pm.core.type.freq <- table.check(pm.core.type.freq, required.years)
pm.core.type.freq.col <- c("cornflowerblue", "darkgoldenrod2", "brown2", "darkolivegreen3", "hotpink3", 
                           "sienna2", "navajowhite2", "seagreen4", "orangered")#, "lightsteelblue3", "darkorange1")
#pm.types.colors <- data.frame(type = rownames(pm.type.freq), col = colors)
#How many of the materials in the last years are the same material but coated one? how many are unique core?
#
#temporal trend in the PM types used by coating and chemical compound:
pm.specif.freq <- table(d.data.nom.type.1990$ENP, d.data.nom.type.1990$year)
pm.specif.freq <- table.check(pm.specif.freq, required.years)
pm.specif.freq <- pm.specif.freq[!rownames(pm.specif.freq) == "", ]#remove first row since its empty
#reorder pm.specif.freq according to the types order in pm.core.type.freq:
#mapping materials to their types:
enp.type.order <- 1:length(unique(d.data.nom.type$ENP.type))
names(enp.type.order) <- as.character(sort(unique(d.data.nom.type$ENP.type)))
material.type.map <- unlist(sapply(rownames(pm.specif.freq), function(x){
  t <- as.character(d.data.nom.type.1990$ENP.type[d.data.nom.type.1990$ENP == x])[1]
  return(as.numeric(enp.type.order[names(enp.type.order) == t][1]))
  }))
length(material.type.map)
material.type.map <- sort(material.type.map)
rows.order <- match(names(material.type.map), rownames(pm.specif.freq))
#material.type.map <- material.type.map[rows.order]
pm.specif.freq <- pm.specif.freq[rows.order, ]
#pm.specif.freq.col <- rainbow(nrow(pm.specif.freq), s = 0.7)
# pm.type.col.ranges <- data.frame(start = c("orchid", "orangered", "orange4", "gray1", "lightblue",
#                                   "palegreen", "rosybrown", "gold", "deeppink", "darksalmon", "green4"),
#                end = c("orchid4", "orangered4", "orange", "gray5", 
#                       "royalblue4", "palegreen4", "rosybrown4", "goldenrod4", "deeppink4", "indianred2", "green4"))
#matching the colors to the types, such that all non metal and non metal oxide are the maroon2 color. See ''enp.type.order'' variable for the matches
pm.type.col.ranges <- data.frame(start = c(rep("maroon2", 3), "lightblue", rep("maroon2", 2), "gold", rep("maroon2", 4)),#monochromatic color for non metal and non metal oxides
                                 end = c(rep("maroon2", 3), "royalblue4", rep("maroon2", 2), "goldenrod4", rep("maroon2", 4)))

pm.specif.freq.col <- c()
for(i in 1:(length(enp.type.order))){
  pm.specif.freq.col <- c(pm.specif.freq.col, colorRampPalette(c(as.character(pm.type.col.ranges$start[i]), 
                                               as.character(pm.type.col.ranges$end[i])))(sum(as.numeric(material.type.map) == i)))
}
#color all the non metal and non metal oxide in one color

#d.ps.specif.freq <- data.frame(name = rownames(pm.specif.freq), col = pm.specif.freq.col)
#display only the top 10 most used materials:
names(pm.specif.freq) <- rownames(pm.specif.freq)
max.names.specif.pm <- names(sort(apply(pm.specif.freq[, as.character(2012:2015)], 1,  sum), decreasing = T)[1:5])
max.names.specif.pm.type <- sapply(max.names.specif.pm, function(x) material.type.map[names(material.type.map) == x])
max.names.specif.pm <- max.names.specif.pm[order(max.names.specif.pm.type, decreasing = T)]
#max.names.specif.pm <- c("ZnO", "TiO2", "Fe2O3", "pvp-Ag", "cit-Ag", "cit-Au", "pvp-Ag")#manual ordering, because I am just tired!
max.colors.spec.pm <- unlist(sapply(max.names.specif.pm, function(x) pm.specif.freq.col[rownames(pm.specif.freq) == x]))
max.names.specif.pm <- c(max.names.specif.pm, "non metal and non metal oxides")
max.colors.spec.pm <- c(max.colors.spec.pm, "maroon1")


# Figure 3 new- diversity evolution 1990-2015 --------------------------------
min.years <- c(TRUE,sapply(2:(length(empir.diversity.time)-1), function(x) (empir.diversity.time[x] < empir.diversity.time[x-1] 
                                                                            & empir.diversity.time[x] < empir.diversity.time[x +1])),FALSE)#years of min repetition (local min)
max.years <- c(FALSE,sapply(2:(length(empir.diversity.time)-1), function(x) (empir.diversity.time[x] > empir.diversity.time[x-1] 
                                                                             & empir.diversity.time[x] > empir.diversity.time[x +1])),TRUE)#years of max repetition (local max)
time.axis <- start.year:max(d.data.nom.type$year)#x axis values

#par(xpd = NA)

m  <- rbind(c(0.1, 0.9, 0.55, 0.95),#main plot combination diversity
            c(0.08, 0.84, 0.3, 0.5),#temporal trends in DOM frequency
            c(0.08, 0.84, 0.1, 0.29))#temporal trends in the PM frequncy
m  <- rbind(c(0.1, 0.9, 0.5, 0.99),#main plot combination diversity
            c(0.08, 0.84, 0.3, 0.45),#temporal trends in DOM frequency
            c(0.08, 0.84, 0.1, 0.29))#temporal trends in the PM frequncy

pdf("Figure3.pdf")
split.screen(m)
screen(1)
par(mar = c(0, 0, 0, 0), cex.axis = 0.7)
plot(time.axis, empir.diversity.time, xlab = "year", ylab = "H(W)", col = "tomato", pch = 8,
     ylim = c(0,1.1), xlim = c(start.year,2018), bty = "l", ann=FALSE, xaxt="n")
axis(1, at = seq(1990, 2015, by =5), labels = FALSE)
mtext(expression("Combinations diversity index D"[comb]), side = 2, line = 2.2, adj = 0.5, cex = 1)#y axis annotation

rect(time.axis[min.years], xright = time.axis[max.years], ybottom = 0, col = "#80808005", ytop = max(high.diversity), border = NA)#periods of repetition, in col the last digits are for transperancy

polygon(c(rev(time.axis), time.axis), c(rev(empir.diversity.DOMbyloc), empir.diversity.DOMbroad), col = '#FA585810',ylim = c(0.2,1),border = NA)#variation due to the different definition of DOM in the dataset, too detialed definition and coarse definitions results in a substantial vairation in the diversity index
par(new=T)
boxplot(t(high.diversity), axes=F, ylab="", xlab="", col = NA, border = "#80808070", boxwex = 0.7, at = time.axis, add = TRUE, cex = 0.7)
par(new=T)
boxplot(t(low.diversity), axes=F, ylab="", xlab="", col = NA, border = "#80808070", boxwex = 0.7, at = time.axis, add = TRUE, cex = 0.7)
#panel numbering:
text(x = start.year, y = 1.1, label = "a", font = 2)#subfigure numbering for the start year network
#subfigures numbering:
text(x = start.year -0.5, y = 1, label = "1", font = 2)#subfigure numbering for the start year network
text(x = 2015, y = g.low.rep + 0.07, label = "2", font = 2)#subfigure numbering for high diversity network
text(x = 2015, y = empir.diversity.time[length(time.axis)] + 0.1, label = "4", font = 2)#subfigure numbering for the empirical network
text(x = 2015, y = g.high.rep + 0.1, label = "3",font = 2)#subfigure numbering for the low diversity network
par(fig = c(0.13, 0.22, 0.80, 0.98), new=T)
plot(g.empir.year, vertex.frame.color = NA, vertex.size = 5, vertex.shape = V(g.empir.year)$shape, vertex.label ="", edge.width = sqrt(E(g.empir.year)$width), vertex.color = "black")
par(fig = c(0.75, 1, 0.85, 0.98), new=T)
plot(g.low, vertex.frame.color = NA, vertex.size = 5, vertex.shape = V(g.high)$shape, vertex.label ="", edge.width = sqrt(sqrt(E(g.high)$width)), vertex.color = "black")#high diversity network (low repetition)
par(fig = c(0.75, 1, 0.7, 0.83), new=T)
plot(g.data.nom.type, vertex.frame.color = NA, vertex.size = 5, vertex.shape = V(g.data.nom.type)$shape,vertex.label = "", edge.width = sqrt(E(g.data.nom.type)$width), vertex.color =  "black")
par(fig = c(0.75, 1, 0.55, 0.68),new=T)
plot(g.high, vertex.frame.color = NA, vertex.size = 5, vertex.shape = V(g.low)$shape,vertex.label = "", edge.width = sqrt(E(g.low)$width)+1, vertex.color = "black")#low diversity network
screen(2)
par(mar = c(0.1, 0.1, 0.1, 0.1), cex.axis = 0.7)
#br <- barplot(dom.type3.freq, las = 2, col = as.character(dom.types3.freq.col), border = NA, xaxt = "n", ann=FALSE, bty = "u", axes=F, xpd = TRUE)
br <- barplot(dom.type2.freq, las = 2, col = dom.type2.freq.col, border = NA, xaxt = "n", ann=FALSE, bty = "u", axes=F, xpd = TRUE)
axis(4, line = -0.9)
#mtext("Number of experiments", side = 4, line = 1, adj = 0.5, cex = 0.8)#y axis annotation
#legend(x = br[11], y = 130, legend = rownames(dom.type3.freq), fill = as.character(dom.types3.freq.col),
 #      border = NA, bty = "n", x.intersp = 0.1, y.intersp = 0.7, title = "DOM", cex = 0.6)
  legend(x = br[11], y = 130, legend = dom.max.use.2012.2015.names, fill = dom.max.use.2012.2015.col,
       border = NA, bty = "n", x.intersp = 0.1, y.intersp = 0.7, title = "DOM", cex = 0.6)
par(new = T)
plot(dom.diversity, type = "l", bty = "n", ylim = c(0.01, 0.5), ylab = "diversity", ann=FALSE, axes =FALSE)
axis(2, line = -0.9)
#mtext(expression("Material diversity index D"[mat]), side = 2, line = 1.7, adj = 0.5, cex = 0.8)#y axis annotation
text(x = 1.7, y = 0.5, label = "b", font = 2)#subfigure numbering for the start year network
screen(3)
par(mar = c(0.1, 0.1, 0.1, 0.1), cex.axis = 0.7)
#br <- barplot(pm.core.type.freq, ann=FALSE, ylab = "", xlab = "", col = pm.core.type.freq.col, axes=F, 
 #             xpd = TRUE, border =NA, bty = "l", xaxt = "n")
br <- barplot(pm.specif.freq, ann=FALSE, ylab = "", xlab = "", col = pm.specif.freq.col, axes=F, 
              xpd = TRUE, border =NA, bty = "l", xaxt = "n")
axis(4, line = -0.9)
mtext("Number of experiments", side = 4, line = 1, adj = -4.5, cex = 0.8)#y axis annotation
#using in legends the most studied ones
#legend(x = br[11], y = 130, legend = rownames(pm.core.type.freq), fill =pm.core.type.freq.col, border = NA, bty = "n", x.intersp = 0.1, y.intersp = 0.7,
 #      title = "PM core material", cex = 0.6)
legend(x = br[11], y = 120, legend = max.names.specif.pm, fill = max.colors.spec.pm, border = NA, bty = "n", x.intersp = 0.1, y.intersp = 0.7,
       title = "PM", cex = 0.6)
axis(1, at = br[seq(1, 26, by = 5)], labels = c("1990", "1995", "2000", "2005", "2010", "2015"))
par(new = T)
plot(pm.diversity, type = "l", bty = "n", ylim = c(0.01, 0.5), ylab = "diversity", ann=FALSE, axes =FALSE)
axis(2, line = -0.9)
mtext(expression("Material diversity index D"[mat]), side = 2, line = 1.7, adj = -1.2, cex = 0.8)#y axis annotation
text(x = 1.7, y = 0.5, label = "c", font = 2)#subfigure numbering for the start year network
mtext("Year", side = 1, line = 2.2, adj = 0.5, cex = 1)#y axis annotation
close.screen(all.screens = TRUE)
dev.off()

# Systematic trends in experimental designs - PM coating effect-------------------------------
###
# The following compares the degree and weight distributions of coated vs. uncoated PM nodes
###

pm.names  <- V(g.data.nom.type)$name[V(g.data.nom.type)$type]
pm.with.coating <- sapply(pm.names, function(x) {if(grepl(x, pattern = "-")){return("orange")}else{return("seagreen3")}})#these are all PM with a "-" symbol, later on make sure that it is not part of the name of an extended PM that has no coating, for now assume "-" seperates coating from core

V(g.data.nom.type)$coating.color <- "#58585830"#the default color for nodes 
V(g.data.nom.type)$coating.color[V(g.data.nom.type)$type] <- pm.with.coating#modify the colors of PM nodes according to their coating/lack of

#number of coated PM:
sum(V(g.data.nom.type)$coating.color == "orange")
#number of uncoated PM:
sum(V(g.data.nom.type)$coating.color == "seagreen3")

coating.analysis  <- data.frame(deg = degree(g.data.nom.type)[V(g.data.nom.type)$type], coat = V(g.data.nom.type)$coating.color[V(g.data.nom.type)$type])
coating.analysis$coat  <- factor(coating.analysis$coat,levels=c("seagreen3", "orange"))#the order of the factor level "coating factor" should be rearranged such that the uncoated (seegreen3) material will be first
#nrow(coating.analysis) == sum(V(g.data.nom.type)$type)#just to make sure only the PM nodes are analyzed
degree.coated.stat  <- wilcox.test(coating.analysis$deg~coating.analysis$coat, paired = FALSE, alternative = "greater", exact = F)# degree is significanlty larger for the uncoated materials, think of performing a non parametric test for two samples, a non parameteric equivalent for two sample t test
#max degree of coated and uncoated PM:
max(coating.analysis$deg[coating.analysis$coat == "seagreen3"])
max(coating.analysis$deg[coating.analysis$coat == "orange"])
####calculate the number of experiments (sum of all links weights) for coated and uncoated materials
#omitting the highly used coated PM:
link.coat.nocitAg  <- sapply(get.edgelist(g.data.nom.type)[ ,1] , function(x) {if(x%in%V(g.data.nom.type)$name[V(g.data.nom.type)$coating.color == "orange"] & (x != "cit-Ag" &
                                                                                                                                                                x != "cit-Au" & x != "pvp-Ag")){return(TRUE)}else{return(FALSE)}})
edge.list = get.edgelist(g.data.nom.type)
n.expr.coated  <- unlist(sapply(V(g.data.nom.type)$name[V(g.data.nom.type)$type == TRUE],function(x) {rel.edges  <- sum(E(g.data.nom.type)$weight[which(edge.list[,1] == x)])
if(x%in%V(g.data.nom.type)$name[V(g.data.nom.type)$coating.color == "orange"] & (x != "cit-Ag" & x != "cit-Au" & x != "pvp-Ag")){return(rel.edges)}}))
n.expr.uncoated  <- unlist(sapply(V(g.data.nom.type)$name[V(g.data.nom.type)$type == TRUE], function(x) {rel.edges  <- sum(E(g.data.nom.type)$weight[which(edge.list[,1] == x)])
if(V(g.data.nom.type)$coating.color[V(g.data.nom.type)$name == x] == "seagreen3"){return(rel.edges)}}))

weight.coated.stat.nocit  <- wilcox.test(n.expr.uncoated, n.expr.coated, paired = FALSE, alternative = "greater", exact = F)
d.coated.nocit  <- data.frame(link = E(g.data.nom.type)$weight, coated = link.coat.nocitAg)

# Figure 4 - PM coating effect --------------------------------------------
#to plot the figure execute first:
par(xpd = NA)
m  <- rbind(c(0.1, 0.9, 0.29, 0.9),#for the network
            c(0.2, 0.5, 0.1, 0.25),#for the degree barplot
            c(0.6, 0.9, 0.1, 0.25))#for the edge weight barplot
#and then the following lines:
split.screen(m)
pdf("Figure4.pdf")
split.screen(m)
screen(1)
par(mar = rep(0,4))
plot(g.data.nom.type, vertex.frame.color = NA, vertex.size = sqrt(degree(g.data.nom.type))+3, vertex.label = "", vertex.shape = V(g.data.nom.type)$shape, edge.width = E(g.data.nom.type)$width/4, vertex.color = V(g.data.nom.type)$coating.color)
legend("topleft", legend =c("coated PM", "PM with\nunspecified\ncoating", "DOM"), col = c("orange", "seagreen3", "#58585850"), pch = c(16, 16, 15), cex = 0.7, bty = "n", border = "white", pt.cex = 2, y.intersp = 1.5)#is the vertical space between lines in the legend
mtext("a", side = 3, line = -0.5, adj = 0.01, cex = 1,font = 2)#subfigure numbering
screen(2)
par(mar = rep(0, 4),cex = 0.7)
boxplot(coating.analysis$deg ~ coating.analysis$coat, col = c("seagreen3", "orange"), names = c("Coating not\nspecified", "Coated"), ylab = "Number of DOM types", log = "y", boxwex = 0.4)
legend("topright", paste("p value = ", as.character(signif(digits = 3, degree.coated.stat$p.value)), sep = " "), bg="white", cex = 0.9)
mtext("Number of DOM types", side = 2, line = 2.2, adj = 0.5, cex = 0.8)#y axis annotation
mtext("b", side = 3, line = -1.2, adj = 0.02, cex = 1, font = 2)#subfigure numbering
screen(3)
par(mar = rep(0, 4), cex = 0.7)

boxplot(n.expr.uncoated, n.expr.coated, col = c("seagreen3", "orange"), names = c("Coating not\nspecified", "Coated"), ylab = "NUmber of experiments", log = "y", boxwex = 0.4)

legend("topright",paste("p value = ", as.character(signif(digits = 3,weight.coated.stat.nocit$p.value)),sep = " "),bg="white",border = "n",cex = 0.9)
mtext("Number of experiments", side = 2, line = 2.2, adj = 0.5, cex = 0.8)#y axis annotation
mtext("c", side = 3, line = -1., adj = 0.02, cex = 1, font = 2)#subfigure numbering
dev.off()
close.screen(all.screens = TRUE)


# Supporing Information ---------------------------------------------------
################################################


# Consistency of reported DOM parameters ------------------------------------


#DOM purchased from the International Humic substances Society:
  # d.papers.char <- read.csv("../newdata/papers_summary.csv", header = T)
  # d.papers.char <- data.frame(apply(d.papers.char, 2, endofline.remover))#remove potentialy leading and trailing "\n" symbols:
  # save(d.papers.char, file = "summary_papers.Rdata")

load("summary_papers.Rdata")#the variable is loaded as: "d.papers.char"
#papers not in the papers characterization database:
d.data.nom.type$DOI[!d.data.nom.type$DOI == ""][!unique(tolower(as.character(d.data.nom.type$DOI[!d.data.nom.type$DOI == ""])))%in%tolower(as.character(d.papers.char$DOI))]
#papers in the characterization database but not in the dom-pm database:
tolower(as.character(d.papers.char$DOI))[!tolower(as.character(d.papers.char$DOI))%in%unique(tolower(as.character(d.data.nom.type$DOI[!d.data.nom.type$DOI == ""])))]

d.papers.char <- data.frame(apply(d.papers.char, 2, na.zero)) #replace NA with zeros)
#elements for the flow chart:
nrow(d.papers.char)#number of papers
sum(d.papers.char$humic.substance == 1)#how many papers employ at least 1 humic acid? but we are interested in the number of publications wmploying at
#least one group-1 DOM:
length(unique(as.character(d.data.nom.type$Full.paper.title[d.data.nom.type$NOM.type%in%group1.dom])))

#not using at least 1 group-1 DOM:
length(unique(as.character(d.data.nom.type$DOI[!d.data.nom.type$DOI == ""]))) - length(unique(as.character(d.data.nom.type$Full.paper.title[d.data.nom.type$NOM.type%in%group1.dom])))

sum(d.papers.char$IHSS ==1)#employing humic substance from the IHSS
sum(d.papers.char$IHSS ==0 & d.papers.char$humic.substance == 1)#employing humic substance not from the IHSS (humic substance and not group1 DOM
#in general)
length(unique(as.character(d.data.nom.type$Full.paper.title[d.data.nom.type$NOM.type%in%group1.dom]))) - sum(d.papers.char$IHSS ==1) #not employing humic substances from the IHSS
sum(d.papers.char$IHSS ==1 & d.papers.char$IHSS.standar.specified ==1)#report standard of the humic substance employed
sum(d.papers.char$IHSS ==1 & d.papers.char$IHSS.standar.specified ==0)#don't report standard of the IHSS sample Cat. no


#Summary of the experiment types:
agg <- sum(sapply(d.papers.char$Experiment.type, function(x) grepl(x = x, pattern = "^1$")))#only aggregation
diss.precip <- sum(sapply(d.papers.char$Experiment.type, function(x) grepl(x = x, pattern = "(2)|(5)")))# dissolution and or precipitation
adsrp <- sum(sapply(d.papers.char$Experiment.type, function(x) grepl(x = x, pattern = "^4$")))
some.diss <- sum(sapply(d.papers.char$Experiment.type, function(x) grepl(x = x, pattern = "2")))#dissoluton
some.agg <- sum(sapply(d.papers.char$Experiment.type, function(x) grepl(x = x, pattern = "1")))#some aggregation

#Availability of chemical properties of group-1 DOM
d.dom.avail <- read.csv("../newdata/dom_propr_avail.csv", header = T)
head(d.dom.avail)
colnames(d.dom.avail)
d.dom.avail <- d.dom.avail[, c("X", "IHSS.", "suva", "Mw", "elemental.comp", "X13C.NMR..full.spectra.")]


#remove empty rows
d.dom.avail <- d.dom.avail[!d.dom.avail$X == "", ]
d.dom.avail <- d.dom.avail[order(d.dom.avail$IHSS), ]
d.dom.avail$X <- tolower(endofline.remover(d.dom.avail$X)) #replace NA with zeros)
100*nrow(d.dom.avail[complete.cases(d.dom.avail),])/nrow(d.dom.avail)#how many complete cases are there? in percent
100*sum(!is.na(d.dom.avail$Mw))/nrow(d.dom.avail)#how many MW is reported?
100*sum(!is.na(d.dom.avail$suva))/nrow(d.dom.avail)#how many MW is reported?
100*sum(!is.na(d.dom.avail$elemental.comp))/nrow(d.dom.avail)#how many MW is reported?
100*sum(!is.na(d.dom.avail$X13C.NMR..full.spectra.))/nrow(d.dom.avail)#how many MW is reported?

#sanity check:  are there names mismatches between the databases? besides unspecific fulvic and humic acids are not considered
unique(tolower(d.data.nom.type$NOM.by.location)[d.data.nom.type$NOM.type%in%group1.dom][!tolower(d.data.nom.type$NOM.by.location[d.data.nom.type$NOM.type%in%group1.dom])%in%d.dom.avail$X])
d.dom.avail$X[!d.dom.avail$X%in%tolower(d.data.nom.type$NOM.by.location[d.data.nom.type$NOM.type%in%group1.dom])]#one entry is not found

d.dom.avail.nona <- data.frame(apply(d.dom.avail[, c("suva", "Mw", "elemental.comp", "X13C.NMR..full.spectra.")], 2, na.zero)) #replace NA with zeros)
#d.dom.avail.nona <- d.dom.avail.nona[, c("suva", "Mw", "elemental.comp", "X13C.NMR..full.spectra.")]
rownames(d.dom.avail.nona) <- rownames(d.dom.avail$X)
d.dom.avail.nona$suva <- as.numeric(as.character(d.dom.avail.nona$suva))
d.dom.avail.nona$elemental.comp <- as.numeric(as.character(d.dom.avail.nona$elemental.comp))
d.dom.avail.nona$Mw <- as.numeric(as.character(d.dom.avail.nona$Mw))
d.dom.avail.nona$X13C.NMR..full.spectra. <- as.numeric(as.character(d.dom.avail.nona$X13C.NMR..full.spectra.))
str(d.dom.avail)

#number of papers that use DOM without one of the parameters:
#number of papers with DOM in group 1:
length(unique(as.character(d.data.nom.type$Full.paper.title[d.data.nom.type$NOM.type%in%group1.dom])))
#no molecular weight
length(unique(as.character(d.data.nom.type$Full.paper.title[tolower(d.data.nom.type$NOM.by.location)%in%as.character(d.dom.avail$X[d.dom.avail.nona$Mw == 0])])))
#no suva
length(unique(as.character(d.data.nom.type$Full.paper.title[tolower(d.data.nom.type$NOM.by.location)%in%as.character(d.dom.avail$X[d.dom.avail.nona$suva == 0])])))
#no elemental composition
length(unique(as.character(d.data.nom.type$Full.paper.title[tolower(d.data.nom.type$NOM.by.location)%in%as.character(d.dom.avail$X[d.dom.avail.nona$X13C.NMR..full.spectra. == 0])])))
length(unique(as.character(d.data.nom.type$Full.paper.title[tolower(d.data.nom.type$NOM.by.location)%in%as.character(d.dom.avail$X[d.dom.avail.nona$elemental.comp == 0])])))
#number of papers that use group1 DOM for which MW, elemental composition or nmr spectra are unknown:
length(unique(as.character(d.data.nom.type$Full.paper.title[
  d.data.nom.type$NOM.by.location%in%as.character(d.dom.avail$X[d.dom.avail.nona$Mw == 0 & 
  d.dom.avail.nona$elemental.comp == 0 & d.dom.avail.nona$X13C.NMR..full.spectra. == 0])])))
#number of publications which use at least 1 grouo1-DOM for which all parameters are known:
length(unique(as.character(c(d.data.nom.type$Full.paper.title[tolower(d.data.nom.type$NOM.by.location)%in%as.character(d.dom.avail$X[d.dom.avail.nona$Mw == 1])],
  d.data.nom.type$Full.paper.title[tolower(d.data.nom.type$NOM.by.location)%in%as.character(d.dom.avail$X[d.dom.avail.nona$X13C.NMR..full.spectra. == 1])],
  d.data.nom.type$Full.paper.title[tolower(d.data.nom.type$NOM.by.location)%in%as.character(d.dom.avail$X[d.dom.avail.nona$elemental.comp == 1])]))))

length(unique(as.character(d.data.nom.type$Full.paper.title[tolower(d.data.nom.type$NOM.by.location)%in%as.character(d.dom.avail$X[d.dom.avail.nona$Mw == 1]) |
               tolower(d.data.nom.type$NOM.by.location)%in%as.character(d.dom.avail$X[d.dom.avail.nona$X13C.NMR..full.spectra. == 1]) |
                 tolower(d.data.nom.type$NOM.by.location)%in%as.character(d.dom.avail$X[d.dom.avail.nona$elemental.comp == 1])])))

#number of publications which use at least 1 group1-DOM for which at least 1 value is unknown: (they can use multiple DOM for which the a value is known)
length(unique(as.character(c(d.data.nom.type$Full.paper.title[tolower(d.data.nom.type$NOM.by.location)%in%as.character(d.dom.avail$X[d.dom.avail.nona$Mw == 0])],
       d.data.nom.type$Full.paper.title[tolower(d.data.nom.type$NOM.by.location)%in%as.character(d.dom.avail$X[d.dom.avail.nona$X13C.NMR..full.spectra. == 0])],
       d.data.nom.type$Full.paper.title[tolower(d.data.nom.type$NOM.by.location)%in%as.character(d.dom.avail$X[d.dom.avail.nona$elemental.comp == 0])]))))

length(unique(as.character(d.data.nom.type$Full.paper.title[d.data.nom.type$NOM.by.location%in%as.character(d.dom.avail$X[d.dom.avail.nona$Mw == 0]) &
d.data.nom.type$NOM.by.location%in%as.character(d.dom.avail$X[d.dom.avail.nona$suva == 0])])))
length(unique(as.character(d.data.nom.type$Full.paper.title[d.data.nom.type$NOM.by.location%in%as.character(d.dom.avail$X[d.dom.avail.nona$suva == 0])])))

#number of experiments that use DOM with full categorization available:
length(as.character(d.data.nom.type$Full.paper.title[d.data.nom.type$NOM.type%in%group1.dom]))
nrow(d.data.nom.type[tolower(d.data.nom.type$NOM.by.location)%in%as.character(d.dom.avail$X[d.dom.avail.nona$Mw == 1]), ])
nrow(d.data.nom.type[tolower(d.data.nom.type$NOM.by.location)%in%as.character(d.dom.avail$X[d.dom.avail.nona$Mw == 0]), ])
nrow(d.data.nom.type[tolower(d.data.nom.type$NOM.by.location)%in%as.character(d.dom.avail$X[d.dom.avail.nona$suva == 1]), ])
nrow(d.data.nom.type[tolower(d.data.nom.type$NOM.by.location)%in%as.character(d.dom.avail$X[d.dom.avail.nona$suva == 0]), ])
nrow(d.data.nom.type[d.data.nom.type$NOM.by.location%in%as.character(d.dom.avail$X[d.dom.avail.nona$Mw == 0 & d.dom.avail.nona$suva == 0]), ])
nrow(d.data.nom.type[d.data.nom.type$NOM.by.location%in%as.character(d.dom.avail$X[d.dom.avail.nona$elemental.comp == 1]), ])
nrow(d.data.nom.type[d.data.nom.type$NOM.by.location%in%as.character(d.dom.avail$X[d.dom.avail.nona$elemental.comp == 0]), ])
nrow(d.data.nom.type[tolower(d.data.nom.type$NOM.by.location)%in%as.character(d.dom.avail$X[d.dom.avail.nona$X13C.NMR..full.spectra. == 1]), ])
nrow(d.data.nom.type[tolower(d.data.nom.type$NOM.by.location)%in%as.character(d.dom.avail$X[d.dom.avail.nona$X13C.NMR..full.spectra. == 0]), ])

heatmap.data.prop <- matrix(nrow = 1, ncol = 3)
for(i in 1:nrow(d.dom.avail.nona)){
  heatmap.data.prop <- rbind(heatmap.data.prop, t(rbind(rep(as.character(d.dom.avail$X[i]), 4), names(d.dom.avail.nona[i,]) ,as.numeric(d.dom.avail.nona[i,]))))
}

heatmap.data.prop <- data.frame(heatmap.data.prop)
heatmap.data.prop[, 3] <- heatmap.data.prop[, 3]
heatmap.data.prop <- heatmap.data.prop[-1, ]
colnames(heatmap.data.prop) <- c("material", "test", "carried.on")
#heatmap.data.prop[, 1] <- as.character(1:nrow(heatmap.data.prop))
is.IHSS <- sapply(as.character(unique(heatmap.data.prop[,1])), function(x) 
  if(d.dom.avail$IHSS.[d.dom.avail$X == x] == 1){return("seagreen3")}else{return("black")})
pdf("FigureProp.pdf", width = 12, height = 7)
p <- (ggplot(heatmap.data.prop, aes(material, test)) + 
        geom_tile(aes(fill = carried.on), colour = "white") +
        scale_fill_manual(drop=FALSE, values=colorRampPalette(c("white","grey50"))(2), na.value="#EEEEEE", name="done?") + 
        theme(axis.ticks = element_blank(), 
              axis.text.x = element_text(size = 6, angle = -90, hjust = 0),
              legend.title = element_text(colour = "grey30"),
              plot.margin = unit(c(0, 0, 0, 0), "cm"),
              axis.title.x=element_blank(),
              axis.title.y=element_blank()
        ))
p
dev.off()

#is there correlation between papers that employ group-1 DOM that all it's dom are characterized in terms of Mw, and seperately in terms of SUVA?
#for this we extract all papers that employ more than 1 group-1 DOM and see if the all materials in this papers are fully caracterized with either
#of the properties: MW or SUVA or elemental composition or NMR spectra (for NMR spectra we expect that there will be many papers with full characterizatio
#if they employ humic substances from the IHSS only).
d.data.compat <- matrix(nrow = 1, ncol = 3)
for(i in unique(d.data.nom.type$DOI[!d.data.nom.type$DOI == ""])){
  #go over all publications in the database,
  used.dom <- d.data.nom.type[
    d.data.nom.type$Full.paper.title == d.data.nom.type$Full.paper.title[d.data.nom.type$DOI == i], c("NOM.by.location", "NOM.type")]#all the dom of this publications
  used.dom.group1 <- tolower(unique(as.character(used.dom$NOM.by.location[used.dom$NOM.type%in%group1.dom])))#the used dom names of group1
  used.dom.group1 <- used.dom.group1[!used.dom.group1%in%c("unfa", "unha")]
  dom.len <- length(used.dom.group1)#the number of froup-1 DOM used in this publication
  if(dom.len > 1){#if there is more than 1 DOM from group 1 studied by this publication, the continue
    info.array <- c()
    if(sum(d.dom.avail.nona$Mw[d.dom.avail$X%in%used.dom.group1]) < dom.len & sum(d.dom.avail.nona$Mw[d.dom.avail$X%in%used.dom.group1]) > 0){#is there mixed information about their mw?
      info.array <- c(info.array, 1)#mixed information
    }
    else if(sum(d.dom.avail.nona$Mw[d.dom.avail$X%in%used.dom.group1]) == dom.len){
        info.array <- c(info.array, 2)
    }else{
      info.array <- c(info.array, 0)
    }#conssistently missing information
    if(sum(d.dom.avail.nona$suva[d.dom.avail$X%in%used.dom.group1]) < dom.len & sum(d.dom.avail.nona$suva[d.dom.avail$X%in%used.dom.group1]) > 0){#is there full information about their suva?
      info.array <- c(info.array, 1)#mixed information
    }
    else if(sum(d.dom.avail.nona$suva[d.dom.avail$X%in%used.dom.group1]) == dom.len){#is there full information about their suva?
      info.array <- c(info.array, 2)#conssitent and full information
    }else{
        info.array <- c(info.array, 0)#conssistently missing information
    }
    if(sum(d.dom.avail.nona$elemental.comp[d.dom.avail$X%in%used.dom.group1]) < dom.len & sum(d.dom.avail.nona$elemental.comp[d.dom.avail$X%in%used.dom.group1]) > 0){#is there full information about their suva?
      info.array <- c(info.array, 1)
    }
    else if(sum(d.dom.avail.nona$elemental.comp[d.dom.avail$X%in%used.dom.group1]) == dom.len){#is there full information about their elemental comp?
      info.array <- c(info.array, 2)#conssitent and full information
    }else{
      info.array <- c(info.array, 0)#conssistently missing information
    }
    if(sum(d.dom.avail.nona$X13C.NMR..full.spectra.[d.dom.avail$X%in%used.dom.group1]) < dom.len & sum(d.dom.avail.nona$X13C.NMR..full.spectra.[d.dom.avail$X%in%used.dom.group1]) >0){#is there full information about their suva?
      info.array <- c(info.array, 1)
    }
    else if(sum(d.dom.avail.nona$X13C.NMR..full.spectra.[d.dom.avail$X%in%used.dom.group1]) == dom.len){#is there full information about their suva?
      info.array <- c(info.array, 2)#conssitent and full information
    }
    else{#is there full information about their suva?
      info.array <- c(info.array, 0)#conssistently missing information
    }
    d.data.compat <- rbind(d.data.compat, t(rbind(rep(as.character(i), 4), 
                      c("molecular weight", "SUVA", "elmental composition", "13C NMR"), info.array)))
    }
}

nrow(d.data.compat)/4
d.data.compat <- data.frame(d.data.compat)
#d.data.compat[, 3] <- d.data.compat[, 3]#conver column to numeric
d.data.compat <- d.data.compat[-1, ]#remove first NA row
colnames(d.data.compat) <- c("publication", "test", "full.information")
#heatmap.data.prop[, 1] <- as.character(1:nrow(heatmap.data.prop))
pdf("FigureCompatibility.pdf", width = 12, height = 7)
p <- (ggplot(d.data.compat, aes(publication, test)) + 
        geom_tile(aes(fill = full.information), colour = "white") +
        scale_fill_manual(drop=FALSE, values= c("grey87", "grey50", "grey20"), na.value="#EEEEEE", name="DOM characteization",
                          labels = c("completely missing", "mixed", "full information")) + 
        theme(axis.ticks = element_blank(), 
              axis.text.x = element_text(size = 6, angle = -90, hjust = 0),
              legend.title = element_text(colour = "grey30"),
              plot.margin = unit(c(0, 0, 0, 0), "cm"),
              axis.title.x=element_blank(),
              axis.title.y=element_blank()
        ))
p
dev.off()
#summary:
l <- length(d.data.compat$publication[d.data.compat$test == "SUVA"])#number of papers that employ at least two group-1 DOM types:
table(d.data.compat[d.data.compat$test == "SUVA", "full.information"])/l
table(d.data.compat[d.data.compat$test == "molecular weight", "full.information"])/l
table(d.data.compat[d.data.compat$test == "elmental composition", "full.information"])/l
table(d.data.compat[d.data.compat$test == "13C NMR", "full.information"])/l
sum(aggregate(as.numeric(as.character(d.data.compat$full.information)), by = list(d.data.compat$publication), FUN = sum) == 8)/l #full information in all DOM parameters is 2*times number of tests
sum(aggregate(as.numeric(as.character(d.data.compat$full.information)), by = list(d.data.compat$publication), FUN = sum) == 2)/l 

# Frequency and prevalence of materials in the experiments ----------------


# abbreviation tables ----------------------------
#the following creates the abbreviation tables for the PM and DOM specified in the database
  # d.si.all <- read.csv("../newdata/newdata.csv")#the table of the database for the supporting information
  # #(contains lines that were removed from d.data.nom.type)
  # save(d.si.all, file = "siData.Rdata")

#creates an abbveriation table for the particulate matter employed in the experiments
load("siData.Rdata")#the read in data will be stored in d.di.all object (a dataframe)
ENP.abbreviations <- d.si.all[d.si.all$ENP.abbreviations != "", "ENP.abbreviations"]
ENP.abbreviations <- data.frame(colsplit(ENP.abbreviations,pattern = "\\s\\|\\s", names = c("abbreviation", "meaning","ref")))
ENP.abbreviations <- ENP.abbreviations[ ,c("abbreviation", "meaning")]
ENP.abbreviations$meaning <- sapply(ENP.abbreviations$meaning,function(x) gsub("^\\s+|\\s+$","", x))#replace any leading or trailing spaces
ENP.abbreviations$abbreviation <- sapply(ENP.abbreviations$abbreviation,function(x) gsub("^\\s+|\\s+$", "", x))#replace any leading or trailing spaces
ENP.abbreviations <- ENP.abbreviations[order(ENP.abbreviations$abbreviation), ]
#ENP.abbreviations$meaning <- sapply(ENP.abbreviations$meaning,function(x) gsub(pattern = "\\|\\s",replacement = "",x))
ENP.abbrev.table <- xtable(ENP.abbreviations, caption = "Particulate matter abbreviations", label = "table:PM:abbreviation", align = "lll", digits = rep(0, 3))
print(ENP.abbrev.table, include.rownames = FALSE, caption.placement = "top", tabular.environment = "longtable", sanitize.text.function = function(x) {x}, floating= FALSE)

#creates an abbveriation table for the DOM employed in the experiments
DOM.abbreviations <- d.si.all[d.si.all$NOM.type.abbreviations != "", "NOM.type.abbreviations"]
DOM.abbreviations <- data.frame(colsplit(DOM.abbreviations,pattern = "\\s\\|\\s", names = c("abbreviation", "meaning")))#split the line into abbreviation and the explanation
DOM.abbreviations$meaning <- sapply(DOM.abbreviations$meaning,function(x) gsub("^\\s+|\\s+$", "", x))#replace any leading or trailing spaces
DOM.abbreviations$abbreviation <- sapply(DOM.abbreviations$abbreviation,function(x) gsub("^\\s+|\\s+$", "", x))#replace any leading or trailing spaces
DOM.abbreviations <- DOM.abbreviations[order(DOM.abbreviations$abbreviation), ]
DOM.abbrev.table <- xtable(DOM.abbreviations, caption = "Dissolved organic matter abbreviations", label = "table:DOM:abbreviation", align = "lll", digits = rep(0, 3))
print(DOM.abbrev.table, include.rownames = FALSE, caption.placement = "top", tabular.environment = "longtable", sanitize.text.function = function(x) {x}, floating = FALSE)

# DOM temporal evolution Figure S3 ---------------
par(xpd = NA)
m  <- rbind(c(0.1, 0.9, 0.8, 0.99),#for the DOM overall
            c(0.1, 0.9, 0.4, 0.8),#for the NOM
            c(0.1,0.9,0.07,0.39))#for the fresh water NOM
#split.screen(m)
pdf("FigureS3.pdf")
split.screen(m)
screen(1)
par(mar = rep(0, 4), cex = 0.7)
barplot(dom.type3.freq, las = 2, col = as.character(dom.types3.freq.col), border = NA, xaxt = "n")
legend("topleft", legend = rownames(dom.type3.freq), fill = as.character(dom.types3.freq.col),
       border = NA, bty = "n", x.intersp = 0.1, y.intersp = 0.7, title = "DOM")
mtext(text = "a", font = 2, side = 2, line = 2, las = 2, adj = 1, at = 120)#subfigure numbering for the start year network
screen(2)
par(mar = rep(0, 4), cex = 0.7)
barplot(hs.type.freq, col = as.character(hs.type.freq.col), las = 2, border = NA, xaxt = "n")
legend("topleft", legend = rownames(hs.type.freq), fill = as.character(hs.type.freq.col),
       border = NA, bty = "n", x.intersp = 0.1, y.intersp = 0.8, cex = 0.9, title = "humic substances/hydrophilic macromolecules")
mtext("number of experiments", side = 2, line = 2.2, adj = 0.5, cex = 0.7)#y axis annotation
mtext(text = "b", font = 2, side = 2, line = 2, las = 2, adj = 1, at  = 60)#subfigure numbering for the start year network
screen(3)
par(mar = rep(0, 4), cex = 0.7)
barplot(fresh.water.hs.freq, las = 2, col = as.character(fresh.water.hs.freq.col), border = NA, ylab = "number of experiments", xlab = "year")
legend("topleft", legend = rownames(fresh.water.hs.freq), fill = as.character(fresh.water.hs.freq.col), 
       bty = "n", x.intersp = 0.1, y.intersp = 0.7, border = NA, title = "fresh water NOM / humic substances", cex = 1)
mtext(text = "c", font = 2, side = 2, line = 2, las = 2, adj = 1, at  = 50)#subfigure numbering for the start year network
mtext(text = "year", side = 1, line = 2.1, adj = 0.5)#
close.screen(all.screens = TRUE)
dev.off()
# PM temporal evolution Figure S4 --------------------------
par(xpd = NA)
m  <- rbind(c(0.1, 0.9, 0.5, 0.9),#for the PM by core
            c(0.1, 0.9, 0.1, 0.5))#for the specific PM
#split.screen(m)
pdf("FigureS4.pdf")
split.screen(m)
screen(1)
par(mar = rep(0, 4), cex = 0.7)
barplot(pm.core.type.freq, col = pm.core.type.freq.col, xaxt = "n", border = NA)#here each material with unique coating is unique.
legend("topleft", legend = rownames(pm.core.type.freq), x.intersp = 0.1, y.intersp = 0.7, fill = pm.core.type.freq.col,
       bty = "n", title = "PM chemcial compound type", border = NA)
mtext(text = "a", font = 2, side = 2, line = 2.9, las = 2, adj = 1, at = 120)
mtext(text = "Number of experiments", side = 2, line = 2, adj = 0.5, at = 0)#
screen(2)
par(mar = rep(0, 4), cex = 0.7)
barplot(pm.specif.freq, col = pm.specif.freq.col, las = 2, border = NA)
legend("topleft", legend = max.names.specif.pm, x.intersp = 0.1, y.intersp = 0.7, fill = max.colors.spec.pm,
       bty = "n", title = "PM chemical compound", border = NA)
mtext(text = "b", font = 2, side = 2, line = 2.9, las = 2, adj = 1, at = 120)#
mtext(text = "Year", side = 1, line = 2.1, adj = 0.5)#
close.screen(all.screens = TRUE)
dev.off()
# Systematic trends in experimental designs - isolated DOM vs. water samples --------
#water.sample  <- c("lake water", "lotic mesocosm water", "sea water", "lagoon water", "ground water", "river water", "STP effluent", "seawater mesocosm",
 #                  "storm water", "fresh water mesocosm", "surface water", "wetland water", "drinking water", "industrial wastewater", "STP influent",
  #                 "pond water", "small stream water", "swamp water")
water.sample <- as.character(unique(d.data.nom.type$NOM.type.extended[d.data.nom.type$NOM.type%in%water.samples]))
#PM is grey and the DOM is colored based on water sample vs. not water sample...
V(g.data.nom.type)$DOM.water  <- ""
V(g.data.nom.type)$DOM.water[V(g.data.nom.type)$name %in% water.sample]  <- "#81DAF5" 
V(g.data.nom.type)$DOM.water[!V(g.data.nom.type)$name %in% water.sample & V(g.data.nom.type)$type == FALSE]  <- "#BE81F7"
V(g.data.nom.type)$DOM.water[V(g.data.nom.type)$type] <-  "#58585830"
#number of water samples:
sum(V(g.data.nom.type)$DOM.water == "#81DAF5")
#number of isolated DOM
sum(V(g.data.nom.type)$DOM.water == "#BE81F7")
#V(g.data.nom.type)$pm.label  <- ""
#V(g.data.nom.type)$pm.label[V(g.data.nom.type)$type]  <- V(g.data.nom.type)$name[V(g.data.nom.type)$type]
#V(g.data.nom.type)$label.cex  <- 0.4

link.water  <- sapply(edge.list[,2], function(x) {if(x%in%water.sample){return(TRUE)}else{return(FALSE)}})#find the links that have water sample for DOM

#statistics regarding the differences in link weights and degree of isolated DOM vs. water samples:
dom.isol.weight.stat <- wilcox.test(E(g.data.nom.type)$weight[link.water],E(g.data.nom.type)$weight[!link.water],paired = FALSE,alternative = "greater")#no significant difference between the number of times a given water sample was studied vs. isolated DOM.
#what about the number of combinations of PM?
degree.waterdom <- degree(g.data.nom.type)[V(g.data.nom.type)$name %in% water.sample]#degree of isolated DOM nodes only degree.waterdom
degree.isoldom <- degree(g.data.nom.type)[!V(g.data.nom.type)$name %in% water.sample & V(g.data.nom.type)$type == FALSE]#degree of DOM water sample nodes only
dom.isol.degree.stat <- wilcox.test(degree.waterdom,degree.isoldom,paired = FALSE,alternative = "greater")#there is a significant difference betweent he diversity of PM used for water samples vs. this used for isolated DOM

#what is the overlap between PM tested with isolated DOM and water sample?
#extract the a list of 1 and 0 for all  PM that was studied with water sample:
#if the edge list of the water sample contains a given PM it will be given 0.5 else it'll be given 0: same for isolated DOM (here the PM is sorted based on its type, given by the column ENP.type)
pm.list.sorted.type  <- unique(d.data.nom.type$ENP[order(d.data.nom.type$ENP.type)])
pm.span.water.sample  <- sapply(pm.list.sorted.type, function(x) {if(x%in%edge.list[link.water,1]){return(0.5)}else{return(0)}})
sum(pm.span.water.sample == 0.5)#number of PM tested with water samples
pm.span.isolated.dom  <- sapply(pm.list.sorted.type, function(x) {if(x%in%edge.list[!link.water,1]){return(0.7)}else{return(0)}})
sum(pm.span.isolated.dom == 0.7)#number of PM tested with isolated DOM, of course there are some that were tested with both therefore the number up to more than the ttoal number of PM nodes

#the dimensions for the regtangles of the isolated DOM (for figure 4)
pm.isolated.dom  <- which(pm.span.isolated.dom == 0.7)#the PM that were tested with isolated DOM
xright.isolated <-pm.isolated.dom[unlist(sapply(1:(length(pm.isolated.dom) - 1),
                                                function(x) {if(pm.isolated.dom[x + 1] == pm.isolated.dom[x] + 1){if(x == length(pm.isolated.dom) - 1){return(c(FALSE,TRUE))}else{
                                                  return(FALSE)}}else{if(x == length(pm.isolated.dom) - 1){return(c(TRUE,TRUE))}else{return(TRUE)}}}))]#checks that the last index if of the one before last entry, if so, always
#set this entry to a right side boundry. If it will also be a left side boundry then this is an isolated entry that was tested
xleft.isolated <-pm.isolated.dom[unlist(sapply(2:length(pm.isolated.dom),
                                               function(x) {if(pm.isolated.dom[x] == pm.isolated.dom[x - 1] + 1){if(x == 2){return(c(TRUE,FALSE))}else{
                                                 return(FALSE)}}else{if(x == 2){return(c(TRUE,TRUE))}else{return(TRUE)}}}))]#Always define the first index as the left boundry.

#the dimensions for the regtangles of the water sample DOM (for figure 4):
pm.water.dom  <- which(pm.span.water.sample == 0.5)#the PM that was studied with water samples
xright.water <- pm.water.dom[unlist(sapply(1:(length(pm.water.dom) - 1), function(x) {if(pm.water.dom[x + 1] == pm.water.dom[x] + 1){
  if(x == length(pm.water.dom) - 1){return(c(FALSE,TRUE))}else{
    return(FALSE)}}else{if(x == length(pm.water.dom) - 1){return(c(TRUE,TRUE))}else{return(TRUE)}}}))]#checks that the last index if of the one before last entry, if so, always
#set this entry to a right side boundry. If it will also be a left side boundry then this is an isolated entry that was tested
xleft.water <- pm.water.dom[unlist(sapply(2:length(pm.water.dom), function(x) {if(pm.water.dom[x] == pm.water.dom[x - 1] + 1){
  if(x == 2){return(c(TRUE,FALSE))}else{
    return(FALSE)}}else{if(x == 2){return(c(TRUE,TRUE))}else{return(TRUE)}}}))]#Always define the first index as the left boundry.

# DOM isolated and water samples - Figure S5 -------------------------------
#The following figure compares the empirical and the subset networks:
par(xpd = NA)
m  <- rbind(c(0.1, 0.9, 0.4, 0.9),#for the network
            c(0.4, 0.6, 0.25, 0.4),#for the degree barplot
            c(0.05, 0.95, 0.15, 0.25))#for the material overlap

#split.screen(m)
pdf("FigureS5.pdf", width = 11, height = 7)
split.screen(m)
screen(1)
par(mar = rep(0., 4))
plot(g.data.nom.type, vertex.frame.color = NA, vertex.size = sqrt(degree(g.data.nom.type))+3, vertex.label = "", vertex.shape = V(g.data.nom.type)$shape, edge.width = E(g.data.nom.type)$width/4,vertex.color = V(g.data.nom.type)$DOM.water)
legend("topleft", legend = c("water sample", "isolated DOM", "PM"), col = c("#81DAF5", "#BE81F7", "grey"), pch = c(15, 15, 16), cex = 0.7, bty = "n", border = "white", pt.cex = 2, y.intersp = 1.5)# is the vertical space between lines in the legend
mtext("a", side = 3, line = -0.5, adj = 0.01, cex = 1, font = 2)#subfigure numbering

screen(2)
par(mar = rep(0.5, 4), cex = 0.7)
boxplot(degree.isoldom,degree.waterdom, col =c("#BE81F7", "#81DAF5"), las = 1, log = "y", names = c("Isolated DOM", "Water sample"), boxwex = 0.2)
legend("topright", paste("p value = ", as.character(signif(digits = 3, dom.isol.degree.stat$p.value)), sep = " "), bg="white", cex = 0.8)
mtext("Number of PM types", side = 2, line = 2.2, adj = 0.5, cex = 0.7)#y axis annotation
mtext("b", side = 3, line = -1.1, adj = -1.1, cex = 1, font = 2)#subfigure numbering

screen(3)
par(mar = rep(0, 4), cex = 0.45)
plot(pm.span.isolated.dom, ylim = c(0.1, 1.5), ann = FALSE, xaxt = "n", yaxt = "n", bty = "n", col = "#BE81F7", pch = "")
rect(xleft = xleft.isolated-0.5, ybottom = rep(0.7-0.03, length(xleft.isolated)), xright = xright.isolated+0.5, ytop = rep(0.7+0.03, length(xright.isolated)), col = "#BE81F7", border = NA)

rect(xleft = xleft.water-0.5, ybottom = rep(0.5-0.03, length(xleft.water)), xright = xright.water+0.5, ytop = rep(0.5+0.03, length(xright.water)), col = "#81DAF5", border = NA)
label.pm  <- as.character(pm.list.sorted.type)
#label.pm[seq(1, length(pm.list.sorted.type), by = 2)]  <- ""
label.pm <- sapply(label.pm,function(x) gsub("-|_", replacement = "", x = x))
axis(1, at = 1:length(label.pm), labels = label.pm, line = -1.5, tick = FALSE, las = 2, col = "grey", cex.axis = 1.4, font = 1)#x axis ticks
mtext("c", side = 3, line = -2, adj = 0.01, cex = 1, font = 2)#subfigure numbering
dev.off()
close.screen(all.screens = TRUE)


# heatmap materials distribution Figure S2 --------------------------------------
#The following code generates the heatmap description of the network, it needs the code for the coating analysis and the water sample analysis

heatmap.matrix  <- sapply(unique(as.character(d.data.nom.type$ENP)), function(x) {sapply(unique(as.character(d.data.nom.type[ ,nom.data.column])),
                                 function(y) {link = paste(x, y, sep ="--"); if(are.connected(g.data.nom.type, x, y)){
                                 return(E(g.data.nom.type)$weight[get.edge.ids(g.data.nom.type, c(x, y))])}else{return(0)}})})#returns a matrix of |PM|x|DOM| dimension

#heatmap.matrix <- heatmap.matrix[order(apply(heatmap.matrix,2,sum)),]#sort according to link weights

heatmap.matrix.melt  <- melt(heatmap.matrix)
colnames(heatmap.matrix.melt) <- c("DOM", "PM", "n.experiments")
heatmap.matrix.melt$DOM  <- factor(heatmap.matrix.melt$DOM, levels = unique(heatmap.matrix.melt$DOM))#the levels must be perserved in order
#for the coloring scheme below to work
heatmap.matrix.melt$PM  <- factor(heatmap.matrix.melt$PM, levels = unique(heatmap.matrix.melt$PM))
#Color DOM based on water sampled or isolated DOM:
pm.color.heatmap  <- sapply(unique(heatmap.matrix.melt$PM), function(x) V(g.data.nom.type)$coating.color[x == V(g.data.nom.type)$name])
#Color PM based on coated or not coated PM:
dom.color.heatmap  <- sapply(unique(heatmap.matrix.melt$DOM), function(x) V(g.data.nom.type)$DOM.water[x == V(g.data.nom.type)$name])
#levels(heatmap.matrix.melt$DOM)[levels(heatmap.matrix.melt$DOM) == "surface water NOM"]  <- "surface water DOM"
pdf("FigureS2.pdf", width = 10, height = 10)
(p <- ggplot(heatmap.matrix.melt, aes(PM, DOM)) + 
  geom_tile(aes(fill = n.experiments), colour = "white") + 
  #scale_x_discrete(expand = c(0.2,0))+
  scale_fill_gradient(low = "white", high = "black", name  = "#experiments", breaks=c(0, 7, max(E(g.data.nom.type)$weight))) + 
  scale_x_discrete(expand = c(0, -1)) +
  theme(axis.ticks = element_blank(), 
        axis.text.x = element_text(size = 4.1, angle = -90, hjust = 0, colour = pm.color.heatmap),
        axis.text.y = element_text(size = 4.1, colour = dom.color.heatmap),
        legend.title = element_text(colour = "grey30"),
        plot.margin = unit(c(0, 0, 0, 0), "cm")
  ))
dev.off()

# Supporting information - NOM source and its chemical composition --------
# NOM distribution in PC space Figure S6 --------------------------------------
#Distribution of labled humic substances in the first two priciple components space:

ggplot(d.pca.summary) +
  geom_vline(xintercept = 0, size = 0.1) +
  geom_hline(yintercept = 0, size = 0.1) +
  geom_point(aes(PC1, PC2, color = type), size = 3)+#, color = 'grey') +
  geom_label_repel(force = 10,
    aes(PC1, PC2, fill = type, label = source),
    fontface = 'bold', color = 'white',
    label.size = 1.5,
    size = 2,
    box.padding = unit(0.25, "lines"),
    point.padding = unit(0.5, "lines")
  ) +
  scale_fill_discrete(guide = guide_legend(ncol =  1, title = "NOM source")) +
  scale_color_discrete(guide = F) +#omit the redundant point legend
  #geom_text(aes(d.pca.summary[, 1], d.pca.summary[, 2], label = ""), show.legend  = F) +
  theme_classic(base_size = 7) +
  theme(axis.line.x = element_line(color="black", size = 0.1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 8),
        axis.line.y = element_line(color="black", size = 0.1))
ggsave("FigureS6.pdf", width = 15, height = 7.5)

#PCA with highlighted regions - Figure S7 --------------------------------
#regions for annotations:
#podzol:
podzol.region <- unlist(sapply(1:nrow(d.pca.summary), function(x) if(grepl(pattern = "podzol", x = d.pca.summary$type[x])){
                return(unlist(d.pca.summary[x, c("PC1", "PC2")]))}else{return(NA)}))
podzol.region.x <- range(podzol.region[names(podzol.region) == "PC1"]) + c(-0.05, 0.05)
podzol.region.y <- range(podzol.region[names(podzol.region) == "PC2"]) + c(-0.05, 0.05)
#peat humic acid:
peat.ha.region <- unlist(sapply(1:nrow(d.pca.summary), function(x) if(grepl(pattern = "^peat humic", x = d.pca.summary$type[x])){
  return(unlist(d.pca.summary[x, c("PC1", "PC2")]))}else{return(NA)}))
peat.ha.region.x <- range(peat.ha.region[names(peat.ha.region) == "PC1"]) + c(-0.05, 0.05)
peat.ha.region.y <- range(peat.ha.region[names(peat.ha.region) == "PC2"]) + c(-0.05, 0.05)
#river fulvic acid:
river.fa.region <- unlist(sapply(1:nrow(d.pca.summary), function(x) if(grepl(pattern = "^river fulvic", x = d.pca.summary$type[x])){
  return(unlist(d.pca.summary[x, c("PC1", "PC2")]))}else{return(NA)}))
river.fa.region.x <- range(river.fa.region[names(river.fa.region) == "PC1"]) + c(-0.05, 0.05)
river.fa.region.y <- range(river.fa.region[names(river.fa.region) == "PC2"]) + c(-0.05, 0.05)
#river humic acid region:
river.ha.region <- unlist(sapply(1:nrow(d.pca.summary), function(x) if(grepl(pattern = "^river humic", x = d.pca.summary$type[x]) &
                                                                       !grepl(pattern = "^ogeechee", x = d.pca.summary$source[x])){
  return(unlist(d.pca.summary[x, c("PC1", "PC2")]))}else{return(NA)}))
river.ha.region.x <- range(river.ha.region[names(river.ha.region) == "PC1"]) + c(-0.05, 0.05)
river.ha.region.y <- range(river.ha.region[names(river.ha.region) == "PC2"]) + c(-0.05, 0.05)
#marine:
marine.region <- unlist(sapply(1:nrow(d.pca.summary), function(x) if(grepl(pattern = "^ocean", x = d.pca.summary$type[x])){
  return(unlist(d.pca.summary[x, c("PC1", "PC2")]))}else{return(NA)}))
marine.region.x <- range(marine.region[names(marine.region) == "PC1"]) + c(-0.05, 0.05)
marine.region.y <- range(marine.region[names(marine.region) == "PC2"]) + c(-0.05, 0.05)
#estuarine sediment:
estuarine.sediment.region <- unlist(sapply(1:nrow(d.pca.summary), 
                             function(x) if(grepl(pattern = "^estuarine sediment", x = d.pca.summary$type[x])
                                            & !grepl(pattern = "chesapeake", x = d.pca.summary$source[x])){
  return(unlist(d.pca.summary[x, c("PC1", "PC2")]))}else{return(NA)}))
estuarine.sediment.region.x <- range(estuarine.sediment.region[names(estuarine.sediment.region) == "PC1"]) + c(-0.05, 0.05)
estuarine.sediment.region.y <- range(estuarine.sediment.region[names(estuarine.sediment.region) == "PC2"]) + c(-0.05, 0.05)
#soil humic acid1:
soil.1.region <- unlist(sapply(1:nrow(d.pca.summary), function(x) if((grepl(pattern = "^soil humic", x = d.pca.summary$type[x])) &
                                                                   !(grepl(pattern = "(summit)|(melbourne)", x = d.pca.summary$source[x]))){
  return(unlist(d.pca.summary[x, c("PC1", "PC2")]))}else{return(NA)}))
soil.1.region.x <- range(soil.1.region[names(soil.1.region) == "PC1"]) + c(-0.05, 0.05)
soil.1.region.y <- range(soil.1.region[names(soil.1.region) == "PC2"]) + c(-0.05, 0.05)

#soil fulvic acid fraction:
soil.faf.region <- unlist(sapply(1:nrow(d.pca.summary), function(x) if(grepl(pattern = "^soil fulvic acid fraction", x = d.pca.summary$type[x])){
  return(unlist(d.pca.summary[x, c("PC1", "PC2")]))}else{return(NA)}))
soil.faf.region.x <- range(soil.faf.region[names(soil.faf.region) == "PC1"]) + c(-0.05, 0.05)
soil.faf.region.y <- range(soil.faf.region[names(soil.faf.region) == "PC2"]) + c(-0.05, 0.05)

ggplot(d.pca.summary) +
  geom_vline(xintercept = 0, size = 0.1) +
  geom_hline(yintercept = 0, size = 0.1) +
  geom_point(aes(PC1, PC2, color = type), size = 3)+#, color = 'grey') +
  scale_fill_discrete(guide = guide_legend(ncol =  1, title = "NOM source")) +
  scale_color_discrete(guide = F) +#omit the redundant point legend
  #geom_text(aes(d.pca.summary[, 1], d.pca.summary[, 2], label = ""), show.legend  = F) +
  theme_classic(base_size = 7) +
  theme(axis.line.x = element_line(color="black", size = 0.1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 8),
        axis.line.y = element_line(color="black", size = 0.1)) +
  annotate("rect", xmin = river.fa.region.x[1] , xmax = river.fa.region.x[2], ymin = river.fa.region.y[1], ymax = river.fa.region.y[2],
           alpha = .2, fill = "skyblue") + 
  annotate("text", x = mean(river.fa.region.x) , y = river.fa.region.y[2] + 0.05, label = "1. river fulvic acids", fontface = 2) +
  annotate("rect", xmin = soil.faf.region.x[1] , xmax = soil.faf.region.x[2], ymin = soil.faf.region.y[1], ymax = soil.faf.region.y[2],
           alpha = .2, fill = "hotpink") + 
  annotate("text", x = mean(soil.faf.region.x) +0.5 , y = soil.faf.region.y[2] + 0.05, label = "2. soil fulvic acids fractions", fontface = 2) +
  annotate("rect", xmin = river.ha.region.x[1] , xmax = river.ha.region.x[2], ymin = river.ha.region.y[1], ymax = river.ha.region.y[2],
           alpha = .2, fill = "darkorchid1") + 
  annotate("text", x = mean(river.ha.region.x) , y = river.ha.region.y[2] + 0.05, label = "3. river humic acids", fontface = 2) +
  annotate("rect", xmin = marine.region.x[1] , xmax = marine.region.x[2], ymin = marine.region.y[1], ymax = marine.region.y[2],
           alpha = .2, fill = "seagreen3") + 
  annotate("text", x = mean(marine.region.x) , y = marine.region.y[2] + 0.05, label = "4. marine DOM", fontface = 2) +
  annotate("rect", xmin = estuarine.sediment.region.x[1] , xmax = estuarine.sediment.region.x[2], ymin = estuarine.sediment.region.y[1], ymax = estuarine.sediment.region.y[2],
           alpha = .2, fill = "darkgoldenrod1") + 
  annotate("text", x = mean(estuarine.sediment.region.x), y = estuarine.sediment.region.y[2] + 0.05, label = "5. estuarine sediment humic acids", fontface = 2) +
  annotate("rect", xmin = podzol.region.x[1] , xmax = podzol.region.x[2], ymin = podzol.region.y[1], ymax = podzol.region.y[2],
           alpha = .1, fill = "cyan3") + 
  # geom_segment(aes(x = d.pca.summary$PC1[d.pca.summary$source == "chesapeake bay sediment humic acid"]+0.01, 
  #                  y = d.pca.summary$PC2[d.pca.summary$source == "chesapeake bay sediment humic acid"]+0.01, 
  #                  xend = d.pca.summary$PC1[d.pca.summary$source == "kings bay sediment humic acid2"], 
  #                  yend = d.pca.summary$PC2[d.pca.summary$source == "kings bay sediment humic acid2"]),
  #              arrow = arrow(length = unit(0.1, "cm"))) + 
  annotate("text", x = podzol.region.x[1] , y = mean(podzol.region.y), label = "6. podzol soil humic acids", fontface = 2, angle = 90) +
  annotate("rect", xmin = peat.ha.region.x[1] , xmax = peat.ha.region.x[2], ymin = peat.ha.region.y[1], ymax = peat.ha.region.y[2],
           alpha = .2, fill = "seagreen3") + 
  annotate("text", x = mean(peat.ha.region.x) , y = peat.ha.region.y[2] + 0.05, label = "7. peat humic acids", fontface = 2) +
  annotate("text", x = d.pca.summary$PC1[d.pca.summary$source == "suwannee river fulvic acid (IHSS 1S101F)"] + 0.05, #annotation of the SRFA humic acid
           y = d.pca.summary$PC2[d.pca.summary$source == "suwannee river fulvic acid (IHSS 1S101F)"]+ 0.05, label = "A1", fontface = 2) + 
  annotate("text", x = d.pca.summary$PC1[d.pca.summary$source == "suwannee river fulvic acid (IHSS 2S101F)"] + 0.05, #annotation of the SRFA humic acid
           y = d.pca.summary$PC2[d.pca.summary$source == "suwannee river fulvic acid (IHSS 2S101F)"]+ 0.05, label = "A2", fontface = 2) +
  annotate("text", x = d.pca.summary$PC1[d.pca.summary$source == "suwannee river humic acid (IHSS 1S101H)"] + 0.05, #annotation of the SRHA humic acid
           y = d.pca.summary$PC2[d.pca.summary$source == "suwannee river humic acid (IHSS 1S101H)"]+ 0.05, label = "B1", fontface = 2) + 
  annotate("text", x = d.pca.summary$PC1[d.pca.summary$source == "suwannee river humic acid (IHSS 2S101H)"] + 0.05, #annotation of the SRHA humic acid
           y = d.pca.summary$PC2[d.pca.summary$source == "suwannee river humic acid (IHSS 2S101H)"]+ 0.05, label = "B2", fontface = 2) + 
  annotate("text", x = d.pca.summary$PC1[d.pca.summary$source == "ogeechee river humic acid"] + 0.05, #annotation of the Ogeechee humic acid
           y = d.pca.summary$PC2[d.pca.summary$source == "ogeechee river humic acid"]+ 0.05, label = "C", fontface = 2) + 
  annotate("text", x = d.pca.summary$PC1[d.pca.summary$source == "sigma aldrich humic acid"] + 0.05, #annotation of the aldrich humic acid
           y = d.pca.summary$PC2[d.pca.summary$source == "sigma aldrich humic acid"]+ 0.05, label = "D", fontface = 2) + 
  annotate("text", x = d.pca.summary$PC1[d.pca.summary$source == "elliot soil humic acid (IHSS 1S102H)"] + 0.05,#annotate the elliot soil humic acid
           y = d.pca.summary$PC2[d.pca.summary$source == "elliot soil humic acid (IHSS 1S102H)"]+ 0.05, label = "E1", fontface = 2) + 
  annotate("text", x = d.pca.summary$PC1[d.pca.summary$source == "elliot soil humic acid (IHSS 4S102H)"] + 0.05, #annotate the elliot soil humic acid
           y = d.pca.summary$PC2[d.pca.summary$source == "elliot soil humic acid (IHSS 4S102H)"]+ 0.05, label = "E2", fontface = 2) + 
  annotate("text", x = d.pca.summary$PC1[d.pca.summary$source == "chesapeake bay sediment humic acid"] + 0.05, #annotate the chesapeake sediement HA outlier
           y = d.pca.summary$PC2[d.pca.summary$source == "chesapeake bay sediment humic acid"]+ 0.05, label = "F", fontface = 2) +
ggsave("FigureS7.pdf", width = 10, height = 7.5)

# Euclidean distances in PC space Figure S8 ------------------------------------
#The following script analyzes the differences in distances of the materials inthe 2D space, spanned by the first PCs.

#First the distance of materials from similar environment is compared to distances of materials from different environments:
inter.type.dist <- c()
intra.type.dist <- c()
dist.matrix <- data.frame(as.matrix(dist(d.pca.summary[ ,1:2], upper = T, diag = T)))
for(i in 1:(nrow(dist.matrix)-1)){
  current.type <- as.character(d.pca.summary$type[i])
  for(j in (i+1):nrow(dist.matrix)){#going over only the upper half of the matrix without the diagonal
    if(as.character(d.pca.summary$type[j]) == current.type){
      intra.type.dist <- c(intra.type.dist, dist.matrix[i, j])
    }else{inter.type.dist <- c(inter.type.dist, dist.matrix[i, j])}
  }
}

#same analysis without soil humic acid:
# inter.type.dist.no.soil <- c()
# intra.type.dist.no.soil <- c()
# d.pca.summary.no.soil <- d.pca.summary[!d.pca.summary$type %in% c("soil humic acid"), ]
# dist.matrix.no.soil <- data.frame(as.matrix(dist(d.pca.summary.no.soil[, 1:2], upper = T, diag = T)))
# 
# for(i in 1:(nrow(dist.matrix.no.soil)-1)){
#   current.type <- as.character(d.pca.summary.no.soil$type[i])
#   for(j in (i+1):nrow(dist.matrix.no.soil)){#going over only the upper half of the matrix without the diagonal
#     if(as.character(d.pca.summary.no.soil$type[j]) == current.type){
#       intra.type.dist.no.soil <- c(intra.type.dist.no.soil, dist.matrix.no.soil[i, j])
#     }else{inter.type.dist.no.soil <- c(inter.type.dist.no.soil, dist.matrix.no.soil[i, j])}
#   }
# }
#second the same analysis is done only by shuffeling the lables of the materials to observe if the different distributions of distances
#is significant or not (y randomization)
inter.type.dist.yrand <- c()
intra.type.dist.yrand <- c()
d.pca.summary.yran <- d.pca.summary
for(y in 1:100){
  d.pca.summary.yran$type <- d.pca.summary.yran$type[sample(1:nrow(d.pca.summary.yran), replace = F)]#randomize the labels of the humic substances
  for(i in 1:(nrow(dist.matrix)-1)){#distances should not be recalculates only the seperation of same label vs. different label
    current.type <- as.character(d.pca.summary.yran$type[i])
    for(j in (i+1):nrow(dist.matrix)){#going over only the upper half of the matrix without the diagonal
      if(as.character(d.pca.summary.yran$type[j]) == current.type){
        intra.type.dist.yrand <- c(intra.type.dist.yrand, dist.matrix[i, j])
      }else{inter.type.dist.yrand <- c(inter.type.dist.yrand, dist.matrix[i, j])}
    }
  }
}

par(xpd = NA)
# m  <- rbind(c(0.0, 0.45, 0.1, 0.9),#original differences in distances
#             c(0.55, 1, 0.1, 0.9))#distances in ''y-randomized'' dataset

#m  <- rbind(c(0.0, 0.35, 0.3, 0.8),#original differences in distances
 #           c(0.31, 0.65, 0.3, 0.8),#distances without soil
  #          c(0.61, 0.95, 0.3, 0.8))#distances in ''y-randomized'' dataset

pdf("FigureS8.pdf")
split.screen(m)
screen(1)
hist(inter.type.dist, xlim = range(c(intra.type.dist, inter.type.dist)), 
     freq = F, ylim = c(0, 1), col = "#F39C1250", border = "grey", xlab = "euclidean distance", main = "", ylab ="")
hist(intra.type.dist, add = T, xlim = range(c(intra.type.dist, inter.type.dist)), freq = F, col = "#45993C70", border = "grey")
text(x = 0.1, y = 1., label = "a", font = 2)
#screen(2)
#hist(inter.type.dist.no.soil, xlim = range(c(intra.type.dist.no.soil, inter.type.dist.no.soil)), 
 #    freq = F, ylim = c(0, 1), col = "#F39C1250", border = "grey", xlab = "euclidean distance", main = "", ylab = "")
#hist(intra.type.dist.no.soil, add = T, xlim = range(c(intra.type.dist.no.soil, inter.type.dist.no.soil)), freq = F, col = "#45993C70", border = "grey")
#text(x = 0.1, y = 1., label = "b", font = 2)
screen(2)
hist(inter.type.dist.yrand, xlim = range(c(intra.type.dist.yrand, inter.type.dist.yrand)), 
     freq = F, ylim = c(0, 1), col = "#F39C1250", border = "grey", xlab = "euclidean distance", main = "")
hist(intra.type.dist.yrand, add = T, xlim = range(c(intra.type.dist.yrand, inter.type.dist.yrand)), freq = F, col = "#45993C70", border = "grey")
text(x = 0.1, y = 1., label = "b", font = 2)
legend("topright", legend = c("different environments", "similar environments"), fill  = c("#F39C1250", "#45993C70"), cex = 0.5, border = "grey", bty = "n")
close.screen(all.screens = T)
dev.off()


# Supporting information - Comparison to a random network -----------------

#number of random network realizations to create:
n.random <- 1000
#Simple graph measures, in order to compare the structural features only we need to remove the weights of the edges such that the graph is binary,
#since the random graph is binary
g.data.nom.type.binary  <- adj.to.graph(d.data.nom.type[c("ENP", nom.data.column)], d.data.nom.type$ENP)
E(g.data.nom.type.binary)$weight  <- 1
E(g.data.nom.type.binary)$width  <- 1
#diameter(g.data.nom.type.binary)
n.nodes <- vcount(g.data.nom.type.binary)
n.links  <- ecount(g.data.nom.type.binary)
p.degree <- ecount(g.data.nom.type.binary)/(sum(V(g.data.nom.type.binary)$type) * sum(V(g.data.nom.type.binary)$type == FALSE))#the probability of linking in a G(n,p) model that corresponds to 
#the empirical network. is adjusted to a bipartite network, since the average number of linksin for the nom is |ENP|*p which is also ecount(nom)/|nom| which is 
#p = ecount(nom)/(|nom|*|enp|) => ecount(network)/(|NOM|*|ENP|)
g.mean.degree <- 2*ecount(g.data.nom.type.binary)/vcount(g.data.nom.type.binary)#mean degree:
g.assor  <- assortativity_degree(g.data.nom.type.binary,directed = FALSE)
#is the graph sparse?
g.rho <- ecount(g.data.nom.type.binary)/(sum(V(g.data.nom.type.binary)$type) * sum(V(g.data.nom.type.binary)$type == FALSE))
g.diameter <- diameter(g.data.nom.type.binary)
g.path <- average.path.length(g.data.nom.type.binary)
#
random.values <- sapply(1:n.random,function(x) {
  set.seed(x)#set seed number for reproducibility
  g.data.nom.type.rand <- bipartite.random.game(sum(V(g.data.nom.type)$type == TRUE), sum(V(g.data.nom.type)$type != TRUE), type  = "gnp", p = p.degree)
  n.nodes <- vcount(g.data.nom.type.rand)
  n.links <- ecount(g.data.nom.type.rand)
  dim  <-  diameter(g.data.nom.type.rand)
  av.path  <-  average.path.length(g.data.nom.type.rand)
  deg.assor  <-  assortativity_degree(g.data.nom.type.rand)
  dens  <- n.links/(sum(V(g.data.nom.type.rand)$type)*sum(V(g.data.nom.type.rand)$type == FALSE))
  mean.degree <- 2*n.links/n.nodes
  return(data.frame(dim,av.path, deg.assor, dens, mean.degree, n.nodes, n.links))})
random.values  <- data.frame(t(apply(random.values, 2, unlist)))
#Calculate the empirical quantiles of each parameter (95\% two sides confidence interval will be the range of values between the 2.5% and 97.5% empirical quantiles)
random.values.95  <- data.frame(sapply(1:ncol(random.values), function(x) return(list(quantile(random.values[ ,x], probs = 0.025),
                     quantile(random.values[ ,x], probs = 0.975)))))
colnames(random.values.95)  <- names(random.values)#give meaningful names to each column

#if the observed parameter is strictly less than 2.5% quantile of strictly more than the 97.5% quantile the value is outside of the 95% confidence interval -> significant
if(g.diameter < as.numeric(random.values.95$dim[1]) || g.diameter > as.numeric(random.values.95$dim[2])){diam.sig <- "*"}else{diam.sig<- as.character()}
if(g.path < as.numeric(random.values.95$av.path[1]) || g.path > as.numeric(random.values.95$av.path[2])){path.sig <- "*"}else{path.sig<- as.character()}
if(g.rho < as.numeric(random.values.95$dens[1]) || g.rho > as.numeric(random.values.95$dens[2])){dens.sig <- "*"}else{dens.sig<- as.character()}
if(g.assor < as.numeric(random.values.95$deg.assor[1]) || g.assor > as.numeric(random.values.95$deg.assor[2])){assor.sig <- "*"}else{assor.sig <- as.character()}
if(g.mean.degree < as.numeric(random.values.95$mean.degree[1]) || g.mean.degree > as.numeric(random.values.95$mean.degree[2])){mean.deg.sig <- "*"}else{mean.deg.sig <- as.character()}

#Ensemble the summery table for the properties:
prop.summary.table <- data.frame("Parameter" = c("mean degree","diameter","average shortest path","density","degree assortativity"),
                      "Empirical network" = c(paste(signif(g.mean.degree,digits = 2),mean.deg.sig),paste(signif(g.diameter,digits = 2),diam.sig),paste(signif(g.path,digits = 2),path.sig),paste(signif(g.rho,digits = 2),dens.sig),paste(signif(g.assor,digits = 2),assor.sig)),
                      "Random network" = c(paste("[",as.character(signif(as.numeric(random.values.95$mean.degree[1]),digits = 2)),",",as.character(signif(as.numeric(random.values.95$mean.degree[2]),digits = 2)),"]",sep = ""),paste("[",as.character(signif(as.numeric(random.values.95$dim[1]),digits = 2)),",",as.character(signif(as.numeric(random.values.95$dim[2]),digits = 2)),"]",sep = ""),
                       paste("[",as.character(signif(as.numeric(random.values.95$av.path[1]),digits = 2)),",",as.character(signif(as.numeric(random.values.95$av.path[2]),digits = 2)),"]",sep = ""),
                       paste("[",as.character(signif(as.numeric(random.values.95$dens[1]),digits = 2)),",",as.character(signif(as.numeric(random.values.95$dens[2]),digits = 2)),"]",sep = ""),
                       paste("[",as.character(signif(as.numeric(random.values.95$deg.assor[1]),digits = 2)),",",as.character(signif(as.numeric(random.values.95$deg.assor[2]),digits = 2)),"]",sep = "")))
# Supporting information - Table S3 ---------------------------------------
print(
  xtable(prop.summary.table,caption = "properties",label = "table:basicnetworkprop",align = "llcc",digits = rep(0,4))
  ,include.rownames = FALSE,caption.placement ="top",table.placement = "H",floating = TRUE,sanitize.text.function = function(x) {x}, scalebox = 0.7)


# Resilience of the experimental network ----------------------------------

# Supporting information - bootstrap---------------
#The following section provides the code for the bootstrap analysis.
d.data.nom.type$comb  <- sapply(1:nrow(d.data.nom.type), function(x) paste(d.data.nom.type$ENP[x], d.data.nom.type[x, nom.data.column], sep = "-"))
#Label each DOM-PM as "new" or "old" based on the definition of diversity:
d.data.nom.type$label.diversity <- rep("old", nrow(d.data.nom.type))
for(i in unique(d.data.nom.type$year)){
  #this loop assigns the label "new" or "old" to each DOM-ENP combination,
  #based on the first definition of diversity
  #this year combinations:
  current.comb  <- unique(d.data.nom.type$comb[d.data.nom.type$year == i])
  n.new.comb <- sum(!current.comb%in%d.data.nom.type$comb[d.data.nom.type$year < i])#the number of unique and new comb added in the current year
  #now assign the label: "new" to n.new.comb and label: "old" to the rest of the combination in this year
  relv.rows  <- which(d.data.nom.type$year == i)#the rows the correspond to experiments that were published in year i
  #assign "new" label to the first n.new.comb rows from relv.rows and all the rest assign the label: "old"
  d.data.nom.type$label.diversity[relv.rows]  <-  c(rep("new", n.new.comb), rep("old", length(relv.rows) - n.new.comb))
}

R  <- 9999#number of bootstrap samples to generate
set.seed(1)
publication.boot  <- boot(d.data.nom.type, bootstrap.publications, R = R, strata = factor(d.data.nom.type$year))
publication.boot.ci.normal  <- sapply(1:length(publication.boot$t0), #95% confidence interval for the diversity index according to the bootstrap samples
                               function(x) {boot.ci(boot.out = publication.boot, type = c("norm"), index = x, conf = 0.95)$normal[2:3]})

# Bootstrap analysis - Figure S9 --------------------------------------
#The following figure depicts the diversity index trend and the 95% bootstrap normal confidence interval.
pdf("FigureS9.pdf")
plot(publication.boot$t0, ylim = c(0,1), col = "black", ylab = expression("D"[comb]), xlab = "Year", xaxt = "n", bty = "n", pch = 8)
axis(1, at = seq(0, 25, by =5), labels = seq(1990, 2015, by =5))
#add the 95% confidence interval based on the bootstrap samples:
publication.boot.ci.normal  <- sapply(1:length(publication.boot$t0), function(x) {boot.ci(boot.out = publication.boot, type = c("norm"), index = x, conf = 0.95)$normal[2:3]})
polygon(c(rev(1:26), 1:26), c(rev(publication.boot.ci.normal[1, ]), publication.boot.ci.normal[2, ]), col = '#2E9AFE30',ylim = c(0.2,1),border = NA)
legend("topleft", legend = c(expression("95% CI (normal) D"[comb])), bty = "n",fill = c('#2E9AFE30'),border = NA)
dev.off()

# The effect of considering only experiments with NOM as their DOM constituent --------
#The following code compares the empirical network with the subset network comprising of 
#only the experiments that employ humic substances as their DOM constituent.
is.humic.subs <- sapply(V(g.data.nom.type)$name[V(g.data.nom.type)$type == FALSE], function(x) {
  found  <- grep(pattern = "(humic)|(fulvic)|(total DOM)|(humic substance)",x = x);if(length(found) != 0){return(TRUE)}else{return(FALSE)}})
humic.nom.names <- V(g.data.nom.type)$name[V(g.data.nom.type)$type == FALSE][is.humic.subs]
pm.names <- V(g.data.nom.type)$name[V(g.data.nom.type)$type == TRUE]
#inducing subgraph to obtian only the PM and humic substances DOM nodes:
g.humic.dom.only <- induced_subgraph(g.data.nom.type,vids = c(pm.names,humic.nom.names))
#Eliminating the nodes that are not connected to anything:
g.components  <- clusters(g.humic.dom.only)
singeltons <- which(g.components$csize == 1)#what components are comprised of singeltons only?
humic.dom.exp.nodes <- names(g.components$membership)[!g.components$membership%in%singeltons]#which pm is are not member of these singletons components? I don't take only the lcc because there can be disconnected nodes that still have DOM of humic substance type and I would like to include these components as well in the dom humic only network
g.humic.dom.only <- induced_subgraph(g.humic.dom.only,vids = humic.dom.exp.nodes)
V(g.humic.dom.only)$degree.lables.code <- label.node.degree(g.humic.dom.only)#label only the 5 PM with the highest degree and the 5 DOM with the highest degree

#comparison to the original network:
#the differences between the general network and the humic only experiments network:
netwroks.compare  <- data.frame("all DOM" = c(as.character(signif(vcount(g.data.nom.type), digits = 3)), as.character(signif(sum(V(g.data.nom.type)$type), digits = 3)),
                                              as.character(signif(sum(V(g.data.nom.type)$type == FALSE), digits = 3)),
                                              as.character(signif(ecount(g.data.nom.type), digits = 3)),
                                              as.character(signif(sum(E(g.data.nom.type)$weight), digits = 3)),
                                              as.character(signif(assortativity_degree(g.data.nom.type), digits = 3)),
                                              as.character(signif(ecount(g.data.nom.type)
                                                                  /(sum(V(g.data.nom.type)$type)*sum(V(g.data.nom.type)$type == FALSE)), digits = 3)),
                                              as.character(signif(ecount(g.data.nom.type)/sum(E(g.data.nom.type)$weight), digits = 3))),
                                "humic DOM only" = c(as.character(signif(vcount(g.humic.dom.only), digits = 3)),as.character(signif(sum(V(g.humic.dom.only)$type), digits = 3)),
                                                     as.character(signif(sum(V(g.humic.dom.only)$type == FALSE), digits = 3)),
                                                     as.character(signif(ecount(g.humic.dom.only), digits = 3)),
                                                     as.character(signif(sum(E(g.humic.dom.only)$weight), digits = 3)),
                                                     as.character(signif(assortativity_degree(g.humic.dom.only), digits = 3)),
                                                     as.character(signif(ecount(g.humic.dom.only)
                                                                         /(sum(V(g.humic.dom.only)$type)*sum(V(g.humic.dom.only)$type == FALSE)), digits = 3)),
                                                     as.character(signif(ecount(g.humic.dom.only)/sum(E(g.humic.dom.only)$weight), digits = 3))))
row.names(netwroks.compare) <- c("number of nodes", "number of PM nodes", "number of DOM nodes", "number of combinations", "number of experiments",
                                 "degree assortativity", "density", "combinations diversity")
# NOM only - Figure S10 --------------------------------------

par(xpd = NA)
m  <- rbind(c(0.0, 0.45, 0.1, 0.9),#empirical network
            c(0.55, 1, 0.1, 0.9))#subset network
pdf("FigureS10.pdf")
split.screen(m)
screen(1)
par(mar = rep(0, 4))
plot(g.data.nom.type, vertex.frame.color = V(g.data.nom.type)$frame.color, vertex.size = log(degree(g.data.nom.type))+4, vertex.shape = V(g.data.nom.type)$shape,
     vertex.color = V(g.data.nom.type)$color, vertex.label = V(g.data.nom.type)$degree.lables.code, vertex.label.cex = 1, 
     vertex.label.color = "black", edge.width = E(g.data.nom.type)$width/4, main="", edge.curved = 0.1)
legend(-1.1, 1.5 , pch = c(21, 22), col = c("dodgerblue4", "grey"), pt.bg = c("light blue", "lightsalmon"), text.col = c("dodgerblue4", "black"),
       legend = c("PM","DOM"), bty = "n", pt.cex = 2, cex = 0.7, y.intersp = 1.5)
#legend for the PM
legend("bottomleft", pch = V(g.data.nom.type)$degree.lables.code[V(g.data.nom.type)$degree.lables.code != "" & V(g.data.nom.type)$type == TRUE],
       legend = V(g.data.nom.type)$name[V(g.data.nom.type)$degree.lables.code != "" & V(g.data.nom.type)$type == TRUE], col = "black", cex = 0.7,
       y.intersp = 1, bty = "n")
#legend for the DOM
legend("bottomright", pch = V(g.data.nom.type)$degree.lables.code[V(g.data.nom.type)$degree.lables.code != "" & V(g.data.nom.type)$type == FALSE],
       legend = V(g.data.nom.type)$name[V(g.data.nom.type)$degree.lables.code != "" & V(g.data.nom.type)$type == FALSE], col = "black", cex = 0.7,
       y.intersp = 1, bty = "n")
screen(2)
par(mar = rep(0, 4))
plot(g.humic.dom.only, vertex.frame.color = V(g.humic.dom.only)$frame.color, vertex.size = log(degree(g.humic.dom.only))+4, vertex.shape = V(g.humic.dom.only)$shape,
     vertex.color = V(g.humic.dom.only)$color, vertex.label = V(g.humic.dom.only)$degree.lables.code, vertex.label.cex = 1, 
     vertex.label.color = "black", edge.width = E(g.humic.dom.only)$width/4, main="", edge.curved = 0.1)
#legend for the PM
legend("bottomleft", pch = V(g.humic.dom.only)$degree.lables.code[V(g.humic.dom.only)$degree.lables.code != "" & V(g.humic.dom.only)$type == TRUE],
       legend = V(g.humic.dom.only)$name[V(g.humic.dom.only)$degree.lables.code != "" & V(g.humic.dom.only)$type == TRUE], col = "black", cex = 0.7,
       y.intersp = 1, bty = "n")
#legend for the DOM
legend("bottomright", pch = V(g.humic.dom.only)$degree.lables.code[V(g.humic.dom.only)$degree.lables.code != "" & V(g.humic.dom.only)$type == FALSE],
       legend = V(g.humic.dom.only)$name[V(g.humic.dom.only)$degree.lables.code != "" & V(g.humic.dom.only)$type == FALSE], col = "black", cex = 0.7,
       y.intersp = 1, bty = "n")
close.screen(all.screens = TRUE)
dev.off()

# Supporting infromation - Table S4 ---------------------------------------
print(
  xtable(netwroks.compare ,caption = "comparison",label = "table:comp",align = "llc",digits = rep(0,3))
  ,include.rownames = TRUE,caption.placement ="top",table.placement = "H",floating = TRUE,sanitize.text.function = function(x) {x}, scalebox = 0.7)

# Supporting information - Temporal dimension of the network's structure --------

#There is no code used besides the one for Figure S5 (see next section)

# Supporting information - Figure S11 --------------------------------------
#color the edges based on year the last experiments employing the two connected materials was performed:
#first convert all lines on the database that they will look like the edge notation in the g.data.nom.type network.
database.edges  <- sapply(1:nrow(d.data.nom.type), function(x) paste(d.data.nom.type$ENP[x], d.data.nom.type$NOM.type.detailed[x], sep = "-"))
network.edges  <- sapply(1:nrow(edge.list), function(x) paste(edge.list[x,1], edge.list[x,2], sep = "-"))
#now match the link to the line
for(i in 1:length(network.edges)){
  rows  <- which(database.edges == network.edges[i])#returns the rows in the database that corresponds to the edge comprising the DOM and PM combination
  #print(max(d.data.nom.type$year[rows]))
  E(g.data.nom.type)$year.latest[i] <- max(d.data.nom.type$year[rows])#the latest year the given combination was studied
}
#now convert these to colors
color.list.year <- colorRampPalette(c("blue", "skyblue", "green", "orange", "yellow"))(length(unique(E(g.data.nom.type)$year.latest)))#create color map to the year published
names(color.list.year) <- sort(unique(E(g.data.nom.type)$year.latest), decreasing = TRUE)#sort the colors
for(i in 1:ecount(g.data.nom.type)){
  #go over the link and assign to them the color for that latest year of publications
  E(g.data.nom.type)$year.color[i]  <- color.list.year[names(color.list.year) == E(g.data.nom.type)$year.latest[i]]
}
pdf("FigureS11.pdf")
plot(g.data.nom.type, vertex.frame.color = "grey", vertex.size = sqrt(degree(g.data.nom.type))+2, vertex.shape = V(g.data.nom.type)$shape,
     edge.width = E(g.data.nom.type)$width/4,main="", vertex.color = ifelse(V(g.data.nom.type)$type == TRUE, "light blue","lightsalmon"), 
     vertex.label = "", edge.color = E(g.data.nom.type)$year.color)
pnts <- cbind(x =c(-0.9, -1., -1, -0.9), y =c(1.0, 1.0, 0.5, 0.5))#points to poisition the legend on the plot
legend.gradient(pnts,rev(color.list.year), limits = range(E(g.data.nom.type)$year.latest), title = "")
legend(-1.1,1.2, pch = c(21,22), col = "grey", pt.bg = c("light blue", "lightsalmon"), text.col = c("blue", "black"), legend = c("PM", "DOM"), bty = "n",
       pt.cex = 2, horiz = TRUE)
dev.off()



