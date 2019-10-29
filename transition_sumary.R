source('EndMTBooleanNetworkFunctions.R')
if(!require(igraph)){
  print("Installing igraph")
  install.packages("igraph")
}
library(igraph)
if(!require(xtable)){
  print("Installing xtable")
  install.packages("xtable")
}
library(xtable)
net_path <- "/Users/natha/OneDrive/Desktop/red/EndMT.net"
endmtnet <- loadNetwork(net_path)

# Sort the attractors
atr <- getAttractors(endmtnet, method="sat.exhaustive")
typesDF <- sortAttractorsDF(atr, sortAttractor)
groups <- form_attractor_groups(typesDF)
sizes <- get_group_sizes(typesDF, form_attractor_groups)
attractor_df <- getAttractorsDF(typesDF, atr)

tfolder <- "C:/Users/natha/OneDrive/Desktop/red/T_Results"
files <- list.files(path=tfolder, pattern="*.RData", full.names=TRUE, recursive=FALSE)
short_files <- list.files(path=tfolder, pattern="*.RData", full.names=FALSE, recursive=FALSE)
load(files[1])
types <- names(sizes)

get_ij <- function(file_name){
  spfn <- unlist(strsplit(x=file_name, split=" "))
  from <- spfn[1]
  to <- unlist(strsplit(x=spfn[3], split=".RData"))
  i <- which(types == from)
  j <- which(types == to)
  ij <- list("i"=i,"j"=j)
  return(ij)
}

get_group_characteristics <- function(analysisij){
  if(nrow(analysisij) == 0){
    return(list("active" = NULL, "inactive" = NULL))
  }
  active <- c()
  inactive <- c()
  for(i in 1:ncol(analysisij)){
    s <- sum(analysisij[, i])
    if(s == 0){
      inactive <- c(inactive, colnames(analysisij)[i])
    } else if (s == nrow(analysisij)){
      active <- c(active, colnames(analysisij)[i])
      #print("Hello world")
    }
  }
  return(list("active" = active, "inactive" = inactive))
}

pert_nrows <- matrix(0, nrow = length(types), ncol = length(types))
colnames(pert_nrows) <- types
rownames(pert_nrows) <- types
p_s_df <- data.frame(pert_nrows)
for(i in 1:length(files)){
  load(files[i])
  ij <- get_ij(short_files[i])
  nr_ij <- nrow(analysisij)
  pert_nrows[ij$i,ij$j] <- nr_ij
  characteristics <- get_group_characteristics(analysisij)
  a_ij <- paste(characteristics$active, collapse = ', ')
  i_ij <- paste(characteristics$inactive, collapse = ', ')
  if(!is.null(characteristics$inactive) && !is.null(characteristics$active)){
    c_ij <- paste(nr_ij,
                  "+(", a_ij, ")",
                  "-(", i_ij, ")", 
                  collapse = '')
  }else if(!is.null(characteristics$active)){
    c_ij <- paste(nr_ij,
                  "+(", a_ij, ")", 
                  collapse = '')
  }else if(!is.null(characteristics$inactive)){
    c_ij <- paste(nr_ij,
                  "-(", i_ij, ")", 
                  collapse = '')
  }else{
    c_ij <- nr_ij
  }
  p_s_df[ij$i,ij$j] <- c_ij
}

latex_p_s_df <- xtable(p_s_df)

n_pert_nrows <- pert_nrows
for(r in 1:nrow(n_pert_nrows)){
  s_row <- sum(n_pert_nrows[r,])
  for(c in 1:ncol(n_pert_nrows)){
    n_pert_nrows[r,c] <- n_pert_nrows[r,c]/s_row
  }
}

g1 <- graph_from_adjacency_matrix(adjmatrix=n_pert_nrows*100, 
                            mode = "directed", 
                            weighted = TRUE,
                            diag = TRUE, 
                            add.colnames = NULL, 
                            add.rownames = NA)

jpeg(file = "transition_graph.jpeg", width = 1800, height = 1800)
plot(g1,
     vertex.shape	= "none",
     vertex.label.cex = 5,
     vertex.label.color = "dark blue",
     vertex.size = 10 * (sizes**(1/4)),
     vertex.color="dark blue",
     edge.width = E(g1)$weight,
     edge.color = heat.colors(81),
     edge.arrow.size = 4,
     edge.arrow.width = 1,
     layout=layout.circle)
dev.off()
