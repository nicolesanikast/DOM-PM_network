#the DOIs in the database:
load("dataBase.Rdata")
doi.exist <- as.character(tolower(d.data.nom.type$DOI[!d.data.nom.type$DOI == ""]))

extract.doi.list <- function(file, doi.exist){
  #this function extracts the DOI list from a ris file of SciFinder
  file <- read.delim(file, header = FALSE)
  #parent reference details:
  list.doi <- sapply(1:nrow(file), function(x) {if(grepl("^DO.*?10.*", as.character(file[x, ]))){
    d <- tolower(gsub("^\\s+|\\s+$", "", unlist(strsplit(as.character(file[x, ]), "DO\\s+-\\s+"))[2]));
    if(!d%in%doi.exist){return(d)}else{return(NA)}#obtaining the DOI and removing heading and trailing white spaces...
  }else{NA}})
  list.doi <- list.doi[!is.na(list.doi)]#remove NA entries
  #print(length(list.doi))
}
file <- "~/Downloads/dom_colloids_precipitation.ris"
strsplit(file, split = ".ris")[[1]]
t <- extract.doi.list(file, doi.exist)
write.csv(t, file = paste(strsplit(file, split = ".ris")[[1]], ".csv", sep = ""))
