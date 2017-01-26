
# Install required packages -----------------------------------------------
load.required.packages <- function(){
  inst.pack <- installed.packages()[ ,1]#installed packages
  required.packages <- c("igraph", "ggplot2", "ggrepel", "reshape2",
                         "bio3d", "MASS", "SDMTools", "TeachingDemos",
                         "xtable", "boot", "gplots", "nlme", "forecast")
  for(i in required.packages){
    if(!is.element(i, inst.pack)){
    install.packages(i)#if not installed then install the package
    }
    library(i,character.only = TRUE)#load the package
  }
}
#load.required.packages()

# Read in the dataframe and convert the data into a graph object ----------
adj.to.graph <- function(d.data,names){
  #this function converts an adjacency matrix to a graph with weights attributes
  #it takes it takes as an argument an adjacency list (data frame) and a names list to seperte between the types of nodes -> to create a bipartite network
  g.data <- graph.adjacency(get.adjacency(graph.data.frame(d.data,directed =FALSE), sparse=FALSE), mode ="undirected", weighted = TRUE)
  #Basic vertex properties
  V(g.data)$type <- V(g.data)$name %in% names #convert to a bipartite network by checking of the value is in the ENP (true) or nom column (false)
  V(g.data)$color  <- ifelse(V(g.data)$type == TRUE, "light blue", "lightsalmon")#color by the node type, PM blue and DOM red
  V(g.data)$label.color[V(g.data)$type]  <- "darkslateblue"#for PM node type
  V(g.data)$frame.color[V(g.data)$type]  <- "dodgerblue4"
  V(g.data)$label.color[!V(g.data)$type]  <- "black"#for DOM type
  V(g.data)$frame.color[!V(g.data)$type]  <- "grey"
  V(g.data)$shape  <- ifelse(V(g.data)$type == TRUE, "circle","square")
  V(g.data)$label.cex = 0.7  
  E(g.data)$width <- E(g.data)$weight
  return(g.data)
}

graph.prop <- function(g.graph){
  #this function takes a graph and extract the basic information
  data.frame(n.nodes = vcount(g.graph),n.links = ecount(g.graph),
  mean.degree = 2*ecount(g.graph)/vcount(g.graph), shortest.path = average.path.length(g.graph), 
  assort = assortativity.degree(g.graph), diameter = diameter(g.graph),density = ecount(g.graph)/(sum(V(g.graph)$type) * sum(V(g.graph)$type == FALSE)))
}

# Add attributes to a graph's vertices ------------------------------------
label.node.degree <- function(g, top = 5){
  #this function calcualtes the top 5 nodes with the highest degree (of each node type -> assuming bipatite network) and returns a the graph instance with the corresponding vertices attributes such as label and label color and frame
  #for this nodes and "" label for all other nodes
  #sort the degree of each node type
  label.by.degree <- rep("", vcount(g))
  #for PM nodes:
  pm.degree <- sort(degree(g)[V(g)$type], decreasing = T)
  pm.degree.top <- names(pm.degree[1:top])#pm nodes with the highest degree
  label.by.degree[V(g)$name %in% pm.degree.top]  <- 1:top#assign numbers as lables to the PM nodes with the highest degree
  #for DOM nodes:
  dom.degree <- sort(degree(g)[V(g)$type == FALSE], decreasing = T)
  dom.degree.top <- names(dom.degree[1:top])#pm nodes with the highest degree
  label.by.degree[V(g)$name %in% dom.degree.top]  <- letters[1:top]#assign numbers as lables to the PM nodes with the highest degree
  return(label.by.degree)
}

# Analysis of diversity trend ---------------------------------------------
modify.subst <- function(doi,g.cit.only.copy,tendency,all.comb,material.in.network){
  #the doi is the currently analyzed node, tendency is the fraction to be imitated from the cited referecnes and all.comb are all the possible combinations avaailable in the\
  #original dataset
  #the number of combinations required by the analyzed node:
  n.exper <- nrow(V(g.cit.only.copy)$material.employ[V(g.cit.only.copy)$name == doi][[1]])
  #print(doi)
  #print(n.exper)
  #the nodes that the given node cites:
  out.neig <- unlist(neighborhood(graph = g.cit.only.copy,order = 1,nodes = doi,mode = "out"))#the first node is always the one that is being analyzed
  #the combinations of materials these nodes employ
  if(length(out.neig) == 1){#if it cites no other references it just returns the smae materials used in the analyzed node without modifications
    #print("hi")
    return(V(g.cit.only.copy)$material.employ[V(g.cit.only.copy)$name == doi])
  }else{out.neig <- out.neig[-1]}
  neig.comb <- data.frame()
  for(neig in names(out.neig)){
    this.neig.comb <- unique(V(g.cit.only.copy)$material.employ[V(g.cit.only.copy)$name == neig][[1]])
    this.neig.comb$comb <- paste(this.neig.comb[[1]],this.neig.comb[[2]],sep = "-")
    #accumulate all the combinations employed by the neighbours in a list (after converting each combination to a string where ENP and NOM are seperated by "-")
    neig.comb <- rbind(neig.comb,this.neig.comb)
    #print(this.neig.comb)
  }
  #modify the substances they employ according to repetition/innovation from their original neighbours int he original graph
  neig.comb.freq <- sort(table(neig.comb$comb),decreasing = TRUE)
  take.from.cited <- round(tendency*n.exper)#the number of combinations to take from the cited references
  take.from.all.comb <- n.exper - take.from.cited#the number of combinations to take from all other combinations that are not in the cited references
  com.not.neig <- all.comb$comb[!all.comb$comb %in% neig.comb$comb]#all the combinations that are not in the neighbouring cited referecnes
  if(take.from.cited == 0){
    #if no combination should be taken from the cited references:
    repeat.comb <- ""
  }else{
    if(length(unique(neig.comb.freq)) == 1){#if all combinations are unique and not repeated in the cited references,sample randomly, so there won't be bias to some letter
      if(length(neig.comb.freq) < take.from.cited){#if there are more comb required than available, sample randomly from the material that was already tested in the network
        repeat.comb <- c(names(neig.comb.freq),sample(material.in.network,take.from.cited-length(neig.comb.freq),replace=TRUE))
      }else{
        repeat.comb <- names(neig.comb.freq[sample(1:nrow(neig.comb),take.from.cited,replace =TRUE)])
      }
    }else{
      #if the cited combinations have unique frequencies
      if(length(neig.comb.freq) < take.from.cited){#if the length of the required combinations is less than the one available from its neighbours add more randomly combinations but ones that were already tested in the network
        repeat.comb <- c(names(neig.comb.freq)[1:take.from.cited],sample(material.in.network,take.from.cited-length(neig.comb.freq),replace=TRUE))
      }else{
        repeat.comb <- names(neig.comb.freq)[1:take.from.cited]
      }
    }}
  #replace to this combination in the currently analyzed node:
  #print(rbind(unique(neig.comb[neig.comb$comb %in% repeat.comb,1:2]),unique(all.comb[all.comb$comb %in% sample(com.not.neig,take.from.all.comb,replace =FALSE),1:2])))
  #return(list(rbind(unique(neig.comb[neig.comb$comb %in% repeat.comb,1:2]),unique(all.comb[all.comb$comb %in% sample(com.not.neig,take.from.all.comb,replace =FALSE),1:2]))))
  return(list(rbind(unique(all.comb[all.comb$comb %in% repeat.comb,1:2]),unique(all.comb[all.comb$comb %in% sample(com.not.neig,take.from.all.comb,replace =TRUE),1:2]))))
}

change.graph.trend <- function(g.graph,year.threshold,d.data,tendency,seed){
  #takes as arguments the original citation graph (that has the materials tested attributes), the year threshold from which it starts to modify researchers behaviour and d.data frame that contains all the relevant data regarding DOI, material used etc... tendency is the balance between innovation and repetition required:
  #first create a copy of the graph
  seed <- 1000*seed#set seed for reproducibilty, is multiplied by 1000 so the next call for the function will get different seed
  g.cit.only.copy <- g.graph
  #define a threshold of a year from which the modification will start to take place
  dois.threshold <- tolower(d.data$DOI[d.data$year >= year.threshold])
  dois.threshold <- dois.threshold[dois.threshold != ""]
  dois.threshold <- as.character(dois.threshold[dois.threshold %in% V(g.cit.only.copy)$name])
  #subset the graph to include only DOI with publication year above thereshold
  sub.g.cit.only <- induced_subgraph(g.cit.only.copy,as.character(dois.threshold))
  material.in.network <- c()# a growing list of the material that were already employed in the network at agiven stage
  while(vcount(sub.g.cit.only) > 1){
    #find among those verteces with outdegree == 0
    old.nodes <- names(which(degree(sub.g.cit.only,mode = "out") == 0))
    for(old in old.nodes){#accumulate all the combinations of materials that were alreay tested in the network
      material <- unique(V(g.cit.only.copy)$material.employ[V(g.cit.only.copy)$name == old][[1]])
      material.in.network <- c(material.in.network,paste(material[[1]],material[[2]],sep = "-"))
    }
    for(old in old.nodes){
      set.seed(seed)
      #iterate over all nodes in the old nodes and modify their employed content accroding to the repetition/innovation rule 
      V(g.cit.only.copy)$material.employ[V(g.cit.only.copy)$name == old] <- modify.subst(old,g.cit.only.copy,tendency,all.comb,material.in.network)
      #print(V(g.cit.only.copy)$material.employ[V(g.cit.only.copy)$name == old])
      seed <- seed + 1
    }
    #delete those vertexes from the subgraph/subset again and find again the ones with 0 outdegree
    sub.g.cit.only <- delete_vertices(sub.g.cit.only,old.nodes) 
    #continue untill all verteces were modified
  }
  #use the modified attrbutes in the g.cit.only.copy as a basis to create new experimental dataframe:
  #materials <- data.frame(ENP=c(),DOM = c(),year = c(),check.rows = FALSE)
  materials <- data.frame()
  for(doi in V(g.cit.only.copy)$name){
    mat.list <- V(g.cit.only.copy)$material.employ[V(g.cit.only.copy)$name == doi][[1]]
    #row.names(mat.list) <- NA
    #year.expr <- rep(unique(d.data.nom.type$year[d.data.nom.type$DOI == doi]),n = nrow(mat.list))
    if(nrow(mat.list) > 0){
      year.expr <- rep(unique(d.data.nom.type$year[d.data.nom.type$DOI == doi]), n = nrow(mat.list))
      mat.list$year <- year.expr
      #materials <- rbind(materials,cbind(mat.list,year.expr))
      materials <- rbind(materials, mat.list)
  }}
  colnames(materials) <- c("ENP","DOM","year")
  #return(list(convert.attrib.to.graph(g.cit.only.copy,d.data),g.cit.only.copy))
  return(materials)
}

diversity.over.time <- function(d.data,start.year){
  #given a dataframe of three columns: PM, DOM and year, this function returns the array of diversity values from the year given by star.year argument (it
  #also considers the experiments prior this year of course)
  diversity.array <- c()
  d.data$comb <- sapply(1:nrow(d.data), function(x) paste(d.data[x,1],d.data[x,2],sep = "-"))#create a forth column with the combination of DOM-PM as single strings
  for(i in start.year:max(d.data$year)){
    d.sliced <- d.data$comb[d.data$year <= i]
    diversity.val <- length(unique(d.sliced))/length(d.sliced)
    diversity.array <- c(diversity.array,diversity.val)
  }
  return(diversity.array)
}


# Analysis of DOM source and chemical composition -------------------------

analysis.DOM.source  <- function(d.data, all.substances){
  #This function performs pca on the d.data, assuming the first column is the identifier of the materials. 
  #all.substances that specifies for each DOM materail its source. The function returns a dataframe comprising of the first two pc scores and the type and identifier of each material
  #now eliminate the columns that have only missing values:
  na.whole.col  <- sapply(1:ncol(d.data), function(x) {sum(is.na(d.data[,x]))==nrow(d.data)})
  #eliminate the rows that have only missing values:
  na.whole.row <- sapply(1:nrow(d.data), function(x) {sum(is.na(d.data[x,])) == (ncol(d.data)-1)})#ncol()-1 since the first column is the identifier and it always has
  #a value
  d.data <- d.data[ ,!na.whole.col]
  d.data <- d.data[!na.whole.row, ]
  #if any value has the sign "<" for example, we have that "<0.05" is represented as "&lt;0.05" therefore there is a need to 
  #before further analysis, convert all columns to characters: because if they are factors, regex subst. won't work
  for(i in 1:ncol(d.data)){d.data[ ,i]  <- as.character(d.data[ ,i])}
  for(row in 1:nrow(d.data)){
    for(col in 3:length(d.data[row, ])){#skip the first three columns that are the col
      x  <- as.character(d.data[row, col])
      if(length(grep(x = x, pattern = "^&lt"))){
        d.data[row,col]  <- gsub(x = as.character(d.data[row, col]), pattern = "&lt\\;", replacement = "")
      }}}
  for(i in 2:ncol(d.data)){d.data[ ,i]  <- as.numeric(d.data[ ,i])}#convert back to numeric type
  
  #ommiting highly correlated variables (one of a pair)
  rownames(d.data) <- d.data[ ,1]
  d.data <- d.data[-1]
  cor.data <- cor(d.data[-1],use = "complete.obs")
  cor.data[upper.tri(cor.data)] <- 0#converting all the upper triangular valeus to zero
  diag(cor.data) <- 0#comverting the diagonal to zer, in order not to consider cor between a variable and itself
  data.new <- d.data[,!apply(cor.data, 2, function(x) any(x > 0.99))]
  
  #now take only complete cases: only now remove the missing values, because some column might be deleted due to high correlation and therefore existance
  #of missing values there should not affect the extraction of missing values
  data.new <- data.new[complete.cases(data.new), ]
  pca.data  <- princomp(data.new,cor = T)
  print(summary(pca.data,loadings = T))
  #print(pca.data$scores[,1])
  xlim  <- range(c(pca.data$scores[ ,1]))*1.5
  ylim  <- range(c(pca.data$scores[ ,2]))*1.5
  d.pca <- data.frame(pca1 = pca.data$scores[ ,1], pca2 = pca.data$scores[ ,2])
  rownames(d.pca) <- sapply(rownames(data.new), function(x) as.character(all.substances[which(all.substances[2] == x), 1]))
  d.pca$type  <- sapply(rownames(data.new), function(x) as.character(unlist(all.substances[which(all.substances[2] == x), 3])[1]))
  d.pca$type  <- factor(d.pca$type)
  #rownames(data.new) <- rownames(d.pca)
  return(d.pca)
}

percentage.in.quadrant <- function(data, rows){
  #this function calculates the percentage of materials environmental source in a given quadrant and returns the one that is most prevelant
  #along with its percentage of occurance
  perc.occur <- sort(table(data$type[rows])/sum(table(data$type[rows])), decreasing = T)#the percentage of occurance of each type
  fresh.water <- signif(sum(perc.occur[sapply(names(perc.occur), 
                                              function(x) {if(grepl(x, pattern = "(lake)|(river)|(groundwater)|(wetland)") 
                                                              & !grepl(x, pattern = "(sediment)")){return(TRUE)}else{return(FALSE)}})]), digits = 2)
  soil <- signif(sum(perc.occur[sapply(names(perc.occur), 
                                       function(x) {if(grepl(x, pattern = "(soil)|(podzol)") 
                                                       & !grepl(x, pattern = "(sediment)")){return(TRUE)}else{return(FALSE)}})]),digits = 2)
  peat <- signif(sum(perc.occur[sapply(names(perc.occur), 
                                       function(x) {if(grepl(x, pattern = "(peat)") 
                                                       & !grepl(x, pattern = "(sediment)")){return(TRUE)}else{return(FALSE)}})]), digits = 2)
  sediment <- signif(sum(perc.occur[sapply(names(perc.occur), 
                                           function(x) {if(grepl(x, pattern = "(sediment)")) {return(TRUE)}else{return(FALSE)}})]), digits = 2)
  marine <- signif(sum(perc.occur[sapply(names(perc.occur), 
                                         function(x) {if(grepl(x, pattern = "(ocean)") 
                                                         & !grepl(x, pattern = "(sediment)")){return(TRUE)}else{return(FALSE)}})]), digits = 2)
  coal <- signif(sum(perc.occur[sapply(names(perc.occur), 
                                       function(x) {if(grepl(x, pattern = "(coal)")){return(TRUE)}else{return(FALSE)}})]), digits = 2)
  aldrich <- signif(sum(perc.occur[sapply(names(perc.occur), 
                                          function(x) {if(grepl(x, pattern = "(aldrich)")){return(TRUE)}else{return(FALSE)}})]), digits = 2)
  sum.list <- c(fresh.water, soil, peat, sediment, marine, coal, aldrich)
  names(sum.list) <- c("fresh water", "soil", "peat", "sediment", "marine", "coal", "aldrich")
  #max.sum.list <- sum.list[sum.list == max(sum.list)]
  sum.list <- sum.list*100#convert to percentages
  max.sum.list <- sum.list[sum.list > 20]
  #concatenate all into a string
  answer <- ""
  for(i in 1:length(max.sum.list)){answer <- paste(paste(answer, names(max.sum.list[i]), sep = ""), paste(max.sum.list[i], "%\n ", sep = ""), sep = " ")}
  return(answer)
}


# Check contingency table -------------------------------------------------
table.check <- function(table, cols){
  #this function checks that a given contingency table (i.e. 'table') contains all the columns listed in 'col'
  #if it doesn't it adds this column in the right position. It assums that the variable 'col' is a numeric list
  if(length(colnames(table)) == length(cols)){return(table)}
  #else add the required columns
  n <- nrow(table)
  for(i in cols[!cols%in%colnames(table)]){#go over all the instances that are not in the table
    table <- cbind(table, rep(0, n))#append the last column to to account for the missing column
    #rename the column
    colnames(table)[ncol(table)] <- i
    #print(i)
  }
  #reorder to columns in an increasing way
  table <- table[ , order(colnames(table))]
  return(table)
}

endofline.remover <- function(x){
  #removes end of line and trailing and leading spaces
  x <- sapply(x, function(y) tolower(gsub(x = gsub(pattern = "^\n|\n$", x = y, replacement = ""), pattern = "^\\s|\\s$",replacement = "")))
}

na.zero <- function (x) {
  #this function replaces NA with 0 in a given array
  x[is.na(x)] <- 0
  return(x)
}

# Bootstrap analysis ------------------------------------------------------
diversity.func.boot <- function(d.data, year.start){
  #This function calculates the diversity of a given experimental dataset. It takes as arguments the dataset (assuming it contains the column: "label.diversity,
  # which contains the labels: "new" or "old" for the DOM-PM combinations, where new/old is given by the defintion of diversity). Year.start is the year from which to
  #start the calculations of diversity and the year from which to start the calculations
  diversity.trend <- c() # place holder for the trend values to be computed below
  for(year in year.start:max(d.data$year)){
    d.data.slice  <- d.data[d.data$year <= year, "label.diversity"]
    n.experiments  <- length(d.data.slice)#the number of experiments done up to the given year
    n.com <- sum(d.data.slice == "new")#the number of unique combinations studied up to the given year
    diversity.trend  <- c(diversity.trend, n.com/n.experiments)
  }
  return(diversity.trend)
}
bootstrap.publications <- function(data, i, year = 1990){
  #This function takes the data set of all DOIs and the index of the lines that are being resampled, where resampling is stratefied by years (so resampling is done from each year 
  #and not indiscriminately from the whole dataset) and calculates the diversity index over the years. Arguments: data = dataset, i = resampled rows, year = the year from which diversity should be calculated, by default is set to 1990
  data  <- data[i,]
  diversity.trend  <- diversity.func.boot(d.data = data, year.start = year)
}
