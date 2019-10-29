source('EndMTBooleanNetworkFunctions.R')


# Load the network from a filewith BoolNet format
net_path <- "/Users/natha/OneDrive/Desktop/red/EndMT.net"
plot_folder_path <- "/Users/natha/OneDrive/Desktop/red/plot_folder"
if(file.exists(plot_folder_path)){
}else{
  dir.create(plot_folder_path)
}


endmtnet <- loadNetwork(net_path)

# Sort the attractors
atr <- getAttractors(endmtnet, method="sat.exhaustive")
typesDF <- sortAttractorsDF(atr, sortAttractor)
groups <- form_attractor_groups(typesDF)
sizes <- get_group_sizes(typesDF, form_attractor_groups)
groups_to_venn(typesDF)

# Find the characteristics of certain groups
chnECsnMCs <- get_attractor_group_characteristics(groups$nECsnMCs, atr)
ECs <- subset(typesDF, EC == TRUE)
chECs <- get_attractor_group_characteristics(ECs, atr)
chEConly <- get_attractor_group_characteristics(groups$EConly, atr)
MCs <- subset(typesDF, MC == TRUE)
chMCs <- get_attractor_group_characteristics(MCs, atr)
chMCsnECs <- get_attractor_group_characteristics(groups$MCsnECs, atr)
ECsMCs <- subset(ECs, MC==TRUE)
chECsMCs <- get_attractor_group_characteristics(ECsMCs, atr)
chMCEConly <- get_attractor_group_characteristics(groups$MCEConly, atr)
chPhalanxes <- get_attractor_group_characteristics(groups$Phalanxes, atr)
Stalks <- subset(typesDF, Stalk == TRUE)
chStalks <- get_attractor_group_characteristics(Stalks, atr)
chMCStalks <- get_attractor_group_characteristics(groups$MCStalks, atr)
chnMCStalks <- get_attractor_group_characteristics(groups$nMCStalks, atr)
Tips <- subset(typesDF, Tip == TRUE)
chTips <- get_attractor_group_characteristics(Tips, atr)
chMCTips <- get_attractor_group_characteristics(groups$MCTips, atr)
chnMCTips <- get_attractor_group_characteristics(groups$nMCTips, atr)

# Save the network, the attractors and the classifications
save(endmtnet, atr, typesDF, groups, sizes, file = "attractors.RData")

# Test the required transitions
required_transitions <- list()
#First the Phalnx endothelial cell behavior
specs <- c("FLI1" = 1) #The other genes are inactive by default
initEC <- generateState(endmtnet, specs)
expectedEC <- c("EC" = TRUE, "MC" = FALSE, "Phalanx" = TRUE, "Stalk" = FALSE, "Tip" = FALSE)
phalanx_transition <- list("initial_state" = initEC, "expected_behavior" = expectedEC, "name" = "FLI1 Phalanx")
required_transitions <- append(required_transitions, list(phalanx_transition))
# Prepare the other transitions
ECpath <- getPathToAttractor(endmtnet, initEC, includeAttractorStates = "all")
EC <- ECpath[dim(ECpath)[1],]
# Now Phalanx to tip
initTip <- unlist(EC)
initTip["VEGFA"] <- 1
expectedTip <- c("EC" = TRUE, "MC" = FALSE, "Phalanx" = FALSE, "Stalk" = FALSE, "Tip" = TRUE)
ectip_transition <- list("initial_state" = initTip, "expected_behavior" = expectedTip, "name" = "VEGFA Tip")
required_transitions <- append(required_transitions, list(ectip_transition))
# Now Tip to Phalanx
Tippath <- getPathToAttractor(endmtnet, initTip, includeAttractorStates = "all")
Tip <- Tippath[dim(Tippath)[1],]
iTip <- unlist(Tip)
iTip["VEGFA"] <- 0
iTip["DLL4"] <- 1
expectediTip <- c("EC" = TRUE, "MC" = FALSE, "Phalanx" = FALSE, "Stalk" = FALSE, "Tip" = FALSE)
tip_to_phalanx <- list("initial_state" = iTip, "expected_behavior" = expectediTip, "name" = "Tip to EC")
required_transitions <- append(required_transitions, list(tip_to_phalanx))
# Next Phalanx to Stalk
initStalk <- unlist(EC)
initStalk["WNT7a"] <- 1
expectedStalk <- c("EC" = TRUE, "MC" = FALSE, "Phalanx" = FALSE, "Stalk" = TRUE, "Tip" = FALSE)
ecstalk_transition <- list("initial_state" = initStalk, "expected_behavior" = expectedStalk, "name" = "WNT7a Stalk")
required_transitions <- append(required_transitions, list(ecstalk_transition))
tr_results <- verify_transitions(endmtnet, required_transitions, sortAttractor)
save(tr_results, file = "tr_results.RData")

# Now we shall test the functions that help us explore the perturbations that lead from 
# one cell type to other cell types
pattern <- c("AP1" = 1, "VEGFA" = 1, "TGF" = 1)
state_with_pattern <- apply_pattern_to_state(EC, pattern)
Phalanxes <- groups[["Phalanxes"]]
isPhalanxinPhalanxes <- row_in_data_frame(expectedEC, Phalanxes[,3:ncol(Phalanxes)]) # Should be TRUE
MCTips <- groups[["MCTips"]]
isPhalanxinTips <- row_in_data_frame(expectedEC, MCTips[,3:ncol(MCTips)]) # Should be FALSE
parameters <- c("DLL4","FGF2","FLI1","GATA2", "HIF1a", "PDGF_AB", "TGFB", "VEGFA", "WNT5b", "WNT7a")
#combPhalanxToMCTip <- probe_transition(Phalanxes, MCTips[1, 3:ncol(MCTips)], endmtnet, atr, sortAttractor, parameters)
nMCTips <- groups[["nMCTips"]]
#combPhalanxTonMCTip <- probe_transition(Phalanxes, nMCTips[1, 3:ncol(nMCTips)], endmtnet, atr, sortAttractor, parameters)
EConly <- groups[["EConly"]]
#combPhalanxToECsnSnT <- probe_transition(Phalanxes, EConly[1, 3:ncol(EConly)], endmtnet, atr, sortAttractor, parameters)
print("Now all the transitions")
#alltransitions <- analyze_transitions(endmtnet, groups, atr, sortAttractor, parameters)
#save(alltransitions, file = "alltransitions.RData")

#Analyze mutant behavior
mutants <- analyze_mutant_behavior(endmtnet, sortAttractor, form_attractor_groups, length(groups))
wt_mutants <- subset(mutants, nECsnMCs > 0 & EConly > 0 & Phalanxes > 0 & nMCStalks > 0 & MCStalks > 0 & nMCTips > 0 & MCTips > 0 & MCEConly > 0 & MCsnECs > 0) # 24 
nECsnMCs_mutants <- subset(mutants, nECsnMCs == 0) # 
EConly_mutants <- subset(mutants, EConly == 0) #
Phalanxes_mutants <- subset(mutants, Phalanxes == 0) # 15
nMCStalks_mutants <- subset(mutants, nMCStalks == 0) # 
MCStalks_mutants <- subset(mutants, MCStalks == 0) #
nMCTips_mutants <- subset(mutants, nMCTips == 0) # 
MCTips_mutants <- subset(mutants, MCTips == 0) #
MCEConly_mutants <- subset(mutants, MCEConly == 0) # 
MCsnECs_mutants <- subset(mutants, MCsnECs == 0) #

# Classification of the mutations that prevent phalanx behavior
p_only <-subset(mutants, nECsnMCs > 0 & EConly > 0 & Phalanxes == 0 & nMCStalks > 0 & MCStalks > 0 & nMCTips > 0 & MCTips > 0 & MCEConly > 0 & MCsnECs > 0)

print("All done before SQUAD")



