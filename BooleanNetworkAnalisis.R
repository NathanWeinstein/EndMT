# This tools are intended to automate the analisis of a boolean network using BoolNet
# First check that BoolNet has been properly installed
if(!require(BoolNet)){
  print("Installing BoolNet")
  install.packages("BoolNet")
}else{
  print("BoolNet is properly installed")
}
library(BoolNet)
if(!require(deSolve)){
  print("Installing deSolve")
  install.packages("deSolve")
}
library(deSolve)
if(!require(doParallel)){
  print("Installing doParallel")
  install.packages("doParallel")
}
library(doParallel)
# General helper functions

# A function to get only the rows that correspond to the attractor reached from a state
getAttractor <- function(net, state, plot_path_file = "none") {
  path <- getPathToAttractor(net, state, includeAttractorStates = "all")
  dim1 <- dim(path)
  r1 <- dim1[1]
  r2 <- dim1[2]
  path0 <- getPathToAttractor(net, state, includeAttractorStates = "none")
  if(plot_path_file != "none"){
    w <- (r1 / 2) + 5
    h <- r2 / 4
    png(file = plot_path_file, width = w, height = h, units = "in", res = 100)
    plotSequence("sequence" = path)
    dev.off()
  }
  dim0 <- dim(path0)
  r0 <- dim0[1]
  return(path[(r0 + 1):r1,])
}
# The main method is analize_boolean_model
# net_path is the pathway to the network to analize
# required_attractor_types is a list of the expected types of attractors
# required_transitions is a list of n-tuples [initial_state, expected_attractor_type]
# attractor_sorter is a function written b the user that classifies the attractors
analize_boolean_model <- function(net_path, results_path, required_transitions = "none", required_attractor_types = "none", attractor_sorter, form_attractor_groups){
  # 1 load the network
  net <- loadNetwork(net_path)
  # 2 get the attractors and sort them
  atri <- getAttractors(net, method="sat.exhaustive")
  typesDF <- sortAttractorsDF(atri, attractor_sorter)
  groups <- form_attractor_groups(typesDF)
  # 3 verify the existence of he required transitions
  transition_results <- verify_transitions(net, required_transitions, attractor_sorter)
  # 4 analyze the conditions that foster the transition from one cell type to other cell types
  transition_analisis <- analyze_transitions(net, groups, atri, attractor_sorter)
  # 5 analyze all loss and gain of function mutantants
  mutant_results <- analyze_mutant_behavior(net, attractor_sorter, form_attractor_groups, length(groups))
  behavior <- list("network" = net, "attractors" = atri, "atractor_classes" = typesDF, "transition_results" = transition_results, "groups" = groups, "mutant_behavior" = mutant_behavior)
}

# veryfy_transitions is a fnction that checks if list of tranitions occures as expected
verify_transitions <- function(net, required_transitions, attractor_sorter){
  transition_results <- c()
  tr_names <- c()
  for(transition in required_transitions){
    attractor <- getAttractor(net, transition$initial_state)
    kind <- attractor_sorter(attractor)
    tr_names <- c(tr_names, transition$name)
    if(compare_logical(kind$celltype, transition$expected_behavior)){
      transition_results <- c(transition_results, TRUE)
    }else{
      transition_results <- c(transition_results, FALSE)
    }
  }
  names(transition_results) <- tr_names
  return(transition_results)
}

compare_logical <- function(l1, l2){
  for(i in 1:length(l1)){
    if(l1[i] != l2[i]){
      return(FALSE)
    }
  }
  return(TRUE)
}

# A function to analyze the transitions and finding the conditions that might cause
# a change in cell identity or behavior
# Each row in each group must include a column with the attractor numbers called "Attractor"

analyze_transitions <- function(net, groups, attractors, attractor_sorter, parameters){
  no_cores <- detectCores() - 1
  print("Number of processes:") 
  print(no_cores)
  registerDoParallel(cores=no_cores)  
  cl <- makeCluster(no_cores, type="FORK")
  groupnames <- names(groups)
  lg <- length(groups) 
  foreach(k = 0:(lg^2 - 1)) %dopar% {
    i <- (k %/% lg) + 1
    j <- (k %% lg) + 1
    print(c("i"=i, "j"=j))
    typej <- groups[[j]][1, 3:ncol(groups[[j]])]
    analysisij <- probe_transition(groups[[i]], typej, net, attractors, attractor_sorter, parameters)
    nameij <- paste(groupnames[i], groupnames[j], sep = " to ")
    print(nameij)
    fnij <-  paste(nameij, ".RData", sep = "")
    save(analysisij, file = fnij)
  }
  stopCluster(cl)
  return(1)
}

# This function probes the transition from one group of attractors to another
probe_transition <- function(groupi, typej, net, attractors, attractor_sorter, parameters){
  combinations <- data.frame() 
  attractorsgi <- groupi[, "Attractor"]
  for(atr in attractorsgi){
    atrseq <- getAttractorSequence(attractors, atr)
    for(r in 1:nrow(atrseq)) {
      row <- atrseq[r,]
      for(p in 0:(2**length(parameters) - 1)){
        pattern <- generate_binary_pattern(p, parameters)
        nstart <- apply_pattern_to_state(row, pattern)
        natr <- getAttractor(net, unlist(nstart))
        nsort <- attractor_sorter(natr)
        if(compare_logical(nsort$celltype, typej)){
          combinations <- rbind(combinations, pattern)
        }
      }
    }
  }
  if(nrow(combinations) > 0){
    colnames(combinations) <- parameters
  }
  # Return the patterns that lead to an attractor in group j.
  return(unique(combinations))
}

# This function probes the transition from one group of attractors to another
probe_transition_0 <- function(groupi, groupj, net, attractors, attractor_sorter, parameters){
  combinations <- data.frame()
  # First find the attractors that are present in group i and not in group j.
  attractorsgi <- groupi[, "Attractor"]
  attractorsgj <- groupj[, "Attractor"]
  dif_ij <- setdiff(attractorsgi, attractorsgj)
  # For each attractor in group i that is not in group j, 
  for(atr in dif_ij){
    atrseq <- getAttractorSequence(attractors, atr)
    for(r in 1:nrow(atrseq)) {
      row <- atrseq[r,]
      for(p in 0:(2**length(parameters) - 1)){
        pattern <- generate_binary_pattern(p, parameters)
        nstart <- apply_pattern_to_state(row, pattern)
        natr <- getAttractor(net, unlist(nstart))
        nsort <- attractor_sorter(natr)
        if(row_in_data_frame(nsort$celltype, groupj[3:ncol(groupj)])){
          combinations <- rbind(combinations, pattern)
        }
      }
    }
  }
  if(nrow(combinations) > 0){
    colnames(combinations) <- parameters
  }
  # Return the patterns that lead to an attractor in group j.
  return(unique(combinations))
}

row_in_data_frame <- function(row, data_frame){
  df <- t(data.frame(row))
  return(nrow(merge(df,data_frame))>0) 
}

apply_pattern_to_state <- function(state, pattern){
  nstate <- state
  for(i in 1:length(pattern)){
    nstate[names(pattern)[i]] <- pattern[i]
  }
  return(nstate)
}

# Returns a numeric type with one number for each variable and asigns the variables as names
generate_binary_pattern <- function(n, variables){
  if(0 <= n & n < 2**length(variables)){
    len <- length(variables)
    bin <- (as.integer(intToBits(n)))
    bins <- bin[1:len]
    bins <- rev(bins)
    names(bins) <- variables
    return(bins)
  }else{
    return("NONE")
  }
}

# A function to clasify all the attractors and return a data frame
sortAttractorsDF <- function(atri, attractor_sorter){
  attrs <- atri$attractors
  atrSeq1 <- getAttractorSequence(atri, 1)
  type1 <- sortAttractor(atrSeq1)
  colnames <- c("Attractor", "Order", names(type1$celltype))
  typesDF = data.frame(matrix(NA, ncol = length(colnames), nrow = length(attrs)))
  names(typesDF) <- colnames
  for(i in 1:length(attrs)){
    atrSeq <- getAttractorSequence(atri, i)
    type <- sortAttractor(atrSeq)
    rowi <- c(i, nrow(atrSeq), type$celltype)
    names(rowi) <- colnames
    typesDF[i,] <- rowi
  }
  return(typesDF)
}

# A function that forms groups of attractors based on a list of logical expressions that describes the 
# characteristics of each group and returns a list of data.frames, one for each defined group
form_attractor_groups <- function(typesDF, group_characteristics){
  groups <- list()
  for(i in 1:length(group_characteristics)){
    groups[[i]] <- subset(typesDF, group_characteristics[i])
  }
  names(groups) <- names(group_characteristics)
  return(groups)
}

# A function that receives an attractor type data frame and an AtractorInfo and returns a data.frame
# containing all the states of all the attractors in the index list
getAttractorsDF <- function(type_data_frame, attractor_info){
  attractor_indexes <- type_data_frame[,"Attractor"]
  nstates <- sum(type_data_frame[,"Order"])
  colnames <- c("Attractor", "State", unlist(attractor_info$stateInfo$genes))
  cellDF = data.frame(matrix(NA, ncol = length(colnames), nrow = nstates))
  names(cellDF) <- colnames
  r <- 1
  for(i in 1:length(attractor_indexes)){
    atrSeq <- getAttractorSequence(attractor_info, attractor_indexes[i])
    for(j in 1:nrow(atrSeq)) {
      pattern <- atrSeq[j,]
      colr <- c(i, j, unlist(pattern))
      names(colr) <- colnames
      cellDF[r,] <- colr
      r <- r + 1
    }
  }
  return(cellDF)
}

# This function receives an atractor group of the same format as returned by sortAttractorsDF
get_attractor_group_characteristics <- function(type_data_frame, attractor_info){
  atractors_df <- getAttractorsDF(type_data_frame, attractor_info)
  active_in_all_states <- c()
  inactive_in_all_states <- c()
  for(i in 3:ncol(atractors_df)){
    s <- sum(atractors_df[, i])
    if(s == 0){
      inactive_in_all_states <- c(inactive_in_all_states, colnames(atractors_df)[i])
    } else if (s == nrow(atractors_df)){
      active_in_all_states <- c(active_in_all_states, colnames(atractors_df)[i])
    }
  }
  return(list("active_in_all_states" = active_in_all_states, "inactive_in_all_states" = inactive_in_all_states))
}

analyze_mutant_behavior <- function(net, attractor_sorter, form_attractor_groups, number_of_groups){
  genes <- net$genes
  mutantsDF = data.frame(matrix(NA, ncol = number_of_groups + 1, nrow = 2*length(genes)))
  i <- 1
  for (gene in genes){
    netl <- fixGenes(net, gene, 0)
    atrl <- getAttractors(netl, method="sat.exhaustive")
    typesDFl <- sortAttractorsDF(atrl, attractor_sorter)
    sizesl <- get_group_sizes(typesDFl, form_attractor_groups)
    rowl <- c("mutation" = paste(gene, "-"), sizesl)
    mutantsDF[i,] <- rowl
    i <- i + 1
    netg <- fixGenes(net, gene, 1)
    atrg <- getAttractors(netg, method="sat.exhaustive")
    typesDFg <- sortAttractorsDF(atrg, attractor_sorter)
    sizesg <- get_group_sizes(typesDFg, form_attractor_groups)
    rowg <- c("mutation" = paste(gene, "+"), sizesg)
    mutantsDF[i,] <- rowg
    i <- i + 1
  }
  names(mutantsDF) <- c( "mutants", names(sizesg))
  return(mutantsDF)
}

get_group_sizes <- function(typesDF, form_attractor_groups){
  groups <- form_attractor_groups(typesDF)
  sizes <- rep(0, length(groups))
  names(sizes) <- names(groups)
  for(i in 1:length(groups)){
    sizes[i] <- nrow(groups[[i]])
  }
  return(sizes)
}

# A function that builds and returns the SQUAD rate of change function based on a "BooleanNetwork"
suqad_rate_of_change <- function(bn, parameters){
  rate_of_change <- function(t, state, parameteers){
    
  }
  return(rate_of_change)
}

#this function generates a function based on a boolean rule in boolnet format
generate_w <- function(boolean_expression){
  
}

# This i the squad function
squad <- function(x, w, h, g){
  num <- -exp(0.5 * h) + exp(-h * (w - 0.5))
  den <- (1 - exp(0.5 * h))*(1 + exp(-h * (w - 0.5)))
  deg <- g*x
  dxdt <- num/den - deg
  return(dxdt)
}

# A function to check the attractors using odes
attractors_as_odes <- function(atrinfo, times, parameters, dif_f, threshold = 1e-10){
  conserved <- c()
  atrs <- atrinfo$attractors
  final_states = data.frame(matrix(NA, 
                                   ncol = length(atrinfo$stateInfo$genes), 
                                   nrow = length(atrs)))
  for(i in 1:length(atrs)){
    atri <- getAttractorSequence(atrinfo, i)
    atri_1_n <- unlist(atri[1,])
    odeSimi <- ode(y = atri_1_n, times = times, func = dif_f, parms = parameters)
    statei <- odeSimi[nrow(odeSimi),]
    final_states[i,] <- statei[2:30]
    difi <- dif_f(1203, statei, parameters)
    check <- (difi[[1]] < threshold) & (difi[[1]] > -threshold)
    print(i)
    print(difi)
    print(check)
    if(mean(check) == 1){
      conserved <- c(conserved, TRUE)
    }else{
      conserved <- c(conserved, FALSE)
    }
  }
  names(final_states) <- atrinfo$stateInfo$genes
  return(list("conserved" = conserved, "final_states" = final_states))
}

# A function to follow n random trajectories in a squad model and storing all steady states
follow_n_trajectories <- function(n, genes, times, parameters, dif_f, threshold = 1e-10){
  final_states = data.frame(matrix(NA, 
                                   ncol = length(genes), 
                                   nrow = n))
  j <- 1
  for(i in 1:n){
    # Generate the random initial value
    random_initial_state <- runif(length(genes))
    names(random_initial_state) <- genes
    odeSimi <- ode(y = random_initial_state, times = times, func = dif_f, parms = parameters)
    statei <- odeSimi[nrow(odeSimi),]
    difi <- dif_f(1203, statei, parameters)
    check <- (difi[[1]] < threshold) & (difi[[1]] > -threshold)
    if(mean(check) == 1){
      final_states[j,] <- statei[2:30]
      j <- j + 1
    }
  }
  names(final_states) <- genes
  return(final_states)
}
# A function to follow n random trajectories respecting certain parameters in a squad model and storing all steady states
follow_n_trajectories_p <- function(n, genes, constants, times, parameters, dif_f, threshold = 1e-10){
  final_states = data.frame(matrix(NA, 
                                   ncol = length(genes), 
                                   nrow = n))
  constant_names <- names(constants)
  j <- 1
  for(i in 1:n){
    # Generate the random initial value
    random_initial_state <- runif(length(genes))
    names(random_initial_state) <- genes
    for(k in 1:length(constants)){
      random_initial_state[constant_names[k]] <- constants[k]
    }
    odeSimi <- ode(y = random_initial_state, times = times, func = dif_f, parms = parameters)
    statei <- odeSimi[nrow(odeSimi),]
    difi <- dif_f(1203, statei, parameters)
    check <- (difi[[1]] < threshold) & (difi[[1]] > -threshold)
    if(mean(check) == 1){
      final_states[j,] <- statei[2:30]
      j <- j + 1
    }
  }
  names(final_states) <- genes
  return(final_states)
}

# The functions that will allow the analisis of the conditions that cause a transition 
# in cell type using SQUAD
ode_transition  <- function(group, cgj, variables, times, parameters, dif_f, sorting_f, threshold = 1e-10){
  functional_patterns <- data.frame()
  for(i in 1:nrow(group)){
    for(p in 0:2**length(variables) - 1){
      pattern <- generate_binary_pattern(p, variables)
      nstart <- unlist(apply_pattern_to_state(group[i,], pattern))
      odeSimi <- ode(y = nstart, times = times, func = dif_f, parms = parameters)
      statei <- odeSimi[nrow(odeSimi),]
      type <- sorting_f(statei[2:30], threshold)
      if(compare_logical(cgj, type)){
        combinations <- rbind(functional_patterns, pattern)
      }
    }
  }
  if(nrow(functional_patterns) > 0){
    colnames(functional_patterns) <- variables
  }
  return(functional_patterns)
}