source('EndMTBooleanNetworkFunctions.R')

net_path <- "/Users/natha/OneDrive/Desktop/red/EndMT.net"
endmtnet <- loadNetwork(net_path)

# Sort the attractors
atr <- getAttractors(endmtnet, method="sat.exhaustive")
typesDF <- sortAttractorsDF(atr, sortAttractor)
groups <- form_attractor_groups(typesDF)
sizes <- get_group_sizes(typesDF, form_attractor_groups)
attractor_df <- getAttractorsDF(typesDF, atr)

ts_hits <- sizes
for(b in 1:length(ts_hits)){
  ts_hits[b] <- 0
}


n_states <- 10000000
intstates <- sample.int(2^length(endmtnet$genes) - 1, n_states)
for(i in intstates){
  specs <- generate_binary_pattern(i, endmtnet$genes)
  state <- generateState(endmtnet, specs)
  atri <- getAttractor(endmtnet, state)
  sortatri <- sortAttractor(atri)$celltype
  if(sortatri["Phalanx"]){
    ts_hits["Phalanxes"] <- ts_hits["Phalanxes"] + 1
  }else if(sortatri["Stalk"]){
    if(sortatri["MC"]){
      ts_hits["MCStalks"] <- ts_hits["MCStalks"] + 1
    }else{
      ts_hits["nMCStalks"] <- ts_hits["nMCStalks"] + 1
    }
  }else if(sortatri["Tip"]){
    if(sortatri["MC"]){
      ts_hits["MCTips"] <- ts_hits["MCTips"] + 1
    }else{
      ts_hits["nMCTips"] <- ts_hits["nMCTips"] + 1
    }
  }else if(sortatri["EC"]){
    if(sortatri["MC"]){
      ts_hits["MCEConly"] <- ts_hits["MCEConly"] + 1
    }else{
      ts_hits["EConly"] <- ts_hits["EConly"] + 1
    }
  }else if(sortatri["MC"]){
    ts_hits["MCsnECs"] <- ts_hits["MCsnECs"] + 1
  }else{
    ts_hits["nECsnMCs"] <- ts_hits["nECsnMCs"] + 1
  }
}
