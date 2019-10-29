# First check that BoolNet has been properly installed
if(!require(BoolNet)){
  print("Installing BoolNet")
  install.packages("BoolNet")
}
library(BoolNet)
if(!require(VennDiagram)){
  print("Installing VennDiagram")
  install.packages("VennDiagram")
}
library(VennDiagram)
if(!require(deSolve)){
  print("Installing deSolve")
  install.packages("deSolve")
}
library(deSolve)

source('BooleanNetworkAnalisis.R')

# A function to classify a state
sortingHat <- function(pattern) {
  GATA2 <- (pattern["GATA2"] == 1)
  FLI1 <- (pattern["FLI1"] == 1)
  SNAI1 <- (pattern["SNAI1"] == 1)
  SNAI2 <- (pattern["SNAI2"] == 1)
  ZEB1 <- (pattern["ZEB1"] == 1)
  ZEB2 <- (pattern["ZEB2"] == 1)
  TWIST1 <- (pattern["TWIST1"] == 1)
  NRP1 <- (pattern["NRP1"] == 1)
  NOTCH <- (pattern["NOTCH"] == 1)
  CTNNB <- (pattern["CTNNB"] == 1)
  ETS1 <- (pattern["ETS1"] == 1)
  if(GATA2 && FLI1) {
    EC <- TRUE
  }else{
    EC <- FALSE
  } 
  if((SNAI1 || SNAI2) && ZEB1 && ZEB2 && TWIST1) {
    MC <- TRUE
  }else{
    MC <- FALSE
  }
  if(!NRP1 && !CTNNB && !SNAI1 && !SNAI2 && EC) {
    Phalanx <- TRUE
  }else{
    Phalanx <- FALSE
  }
  if(!NRP1 && CTNNB && SNAI2 && EC) {
    Stalk <- TRUE
  }else{
    Stalk <- FALSE
  }
  if(NRP1 && (ETS1 || NOTCH) && EC) {
    Tip <- TRUE
  }else{
    Tip <- FALSE
  }
  celltype <- c(EC, MC, Phalanx, Stalk, Tip)
  names(celltype) <- c("EC", "MC", "Phalanx", "Stalk", "Tip")
  return(celltype)
}


# A function to classify a state in the continuos model
sortingHatc <- function(pattern, threshold) {
  GATA2 <- (pattern["GATA2"] > threshold)
  FLI1 <- (pattern["FLI1"] > threshold)
  SNAI1 <- (pattern["SNAI1"] > threshold)
  SNAI2 <- (pattern["SNAI2"] > threshold)
  ZEB1 <- (pattern["ZEB1"] > threshold)
  ZEB2 <- (pattern["ZEB2"] > threshold)
  TWIST1 <- (pattern["TWIST1"] > threshold)
  NRP1 <- (pattern["NRP1"] > threshold)
  NOTCH <- (pattern["NOTCH"] > threshold)
  CTNNB <- (pattern["CTNNB"] > threshold)
  ETS1 <- (pattern["ETS1"] > threshold)
  if(GATA2 && FLI1) {
    EC <- TRUE
  }else{
    EC <- FALSE
  } 
  if((SNAI1 || SNAI2) && ZEB1 && ZEB2 && TWIST1) {
    MC <- TRUE
  }else{
    MC <- FALSE
  }
  if(!NRP1 && !CTNNB && !SNAI1 && !SNAI2 && EC) {
    Phalanx <- TRUE
  }else{
    Phalanx <- FALSE
  }
  if(!NRP1 && CTNNB && SNAI2 && EC) {
    Stalk <- TRUE
  }else{
    Stalk <- FALSE
  }
  if(NRP1 && ETS1 && EC) {
    Tip <- TRUE
  }else{
    Tip <- FALSE
  }
  celltype <- c(EC, MC, Phalanx, Stalk, Tip)
  names(celltype) <- c("EC", "MC", "Phalanx", "Stalk", "Tip")
  return(celltype)
}

# A function to classify an attractor
sortAttractor <- function(atrSequence) {
  celltype <- c(TRUE, TRUE, TRUE, TRUE, TRUE)
  pcelltype <- c(FALSE, FALSE, FALSE, FALSE, FALSE)
  names(celltype) <- c("EC", "MC", "Phalanx", "Stalk", "Tip")
  names(pcelltype) <- c("pEC", "pMC", "pPhalanx", "pStalk", "pTip")
  for(i in 1:nrow(atrSequence)) {
    pattern <- atrSequence[i,]
    ptype <- sortingHat(pattern)
    if(ptype["EC"]) {
      pcelltype["pEC"] <- TRUE
    }else{
      celltype["EC"] <- FALSE
    }
    if(ptype["MC"]) {
      pcelltype["pMC"] <- TRUE
    }else{
      celltype["MC"] <- FALSE
    }
    if(ptype["Phalanx"]) {
      pcelltype["pPhalanx"] <- TRUE
    }else{
      celltype["Phalanx"] <- FALSE
    }
    if(ptype["Stalk"]) {
      pcelltype["pStalk"] <- TRUE
    }else{
      celltype["Stalk"] <- FALSE
    }
    if(ptype["Tip"]) {
      pcelltype["pTip"] <- TRUE
    }else{
      celltype["Tip"] <- FALSE
    }
  }
  type <- list(celltype, pcelltype)
  names(type) <- c("celltype", "pcelltype")
  return(type)
}

# This function groups the attractos based on certain characteristics
form_attractor_groups <- function(typesDF){
  ECs <- subset(typesDF, EC == TRUE) # Perhaps only ECs are relevant?
  MCs <- subset(typesDF, MC == TRUE)
  # From here the groups make sense
  MCsnECs <- subset(MCs, EC == FALSE) # i
  MCEConly <- subset(MCs, EC == TRUE & Tip == FALSE & Stalk == FALSE & Phalanx == FALSE) #h
  nECsnMCs <- subset(typesDF, MC == FALSE & EC == FALSE) # a
  EConly <- subset(typesDF, EC == TRUE & MC == FALSE & Tip == FALSE & Stalk == FALSE & Phalanx == FALSE) # b
  Phalanxes <- subset(typesDF, Phalanx == TRUE)
  Stalks <- subset(typesDF, Stalk == TRUE)
  MCStalks <- subset(Stalks, MC == TRUE)
  nMCStalks <- subset(Stalks, MC == FALSE) # d
  Tips <- subset(typesDF, Tip == TRUE)
  MCTips <- subset(Tips, MC == TRUE) # g
  nMCTips <- subset(Tips, MC == FALSE) # f
  min_groups <- list("Phalanxes" = Phalanxes,"MCStalks" = MCStalks,"nMCStalks" = nMCStalks,"MCTips" = MCTips,"nMCTips" = nMCTips,"EConly" = EConly)
  all_groups <- list("nECsnMCs" = nECsnMCs, "EConly" = EConly, "Phalanxes" = Phalanxes, "nMCStalks" = nMCStalks, "MCStalks" = MCStalks, "nMCTips" = nMCTips, "MCTips" = MCTips, "MCEConly" = MCEConly,"MCsnECs" = MCsnECs)
  return(all_groups)
}

# A function to create a venn diagram based on the sizes of the attractor groups 
groups_to_venn <- function(typesDF){
  venn.plot <- draw.quintuple.venn(
    area1 = nrow(subset(typesDF, EC == TRUE)), # ECs
    area2 = nrow(subset(typesDF, MC == TRUE)), # MCs
    area3 = nrow(subset(typesDF, Phalanx == TRUE)), # Phalanxes
    area4 = nrow(subset(typesDF, Stalk == TRUE)), # Stalks
    area5 = nrow(subset(typesDF, Tip == TRUE)), # Tips
    n12 = nrow(subset(typesDF, EC == TRUE & MC == TRUE)), # MCsECs,
    n13 = nrow(subset(typesDF, EC == TRUE & Phalanx == TRUE)), # EC Phalanxes
    n14 = nrow(subset(typesDF, EC == TRUE & Stalk == TRUE)), # EC Stalks
    n15 = nrow(subset(typesDF, EC == TRUE & Tip == TRUE)), # EC Tips
    n23 = nrow(subset(typesDF, MC == TRUE & Phalanx == TRUE)),
    n24 = nrow(subset(typesDF, MC == TRUE & Stalk == TRUE)),
    n25 = nrow(subset(typesDF, MC == TRUE & Tip == TRUE)),
    n34 = nrow(subset(typesDF, Phalanx == TRUE & Stalk == TRUE)),
    n35 = nrow(subset(typesDF, Phalanx == TRUE & Tip == TRUE)),
    n45 = nrow(subset(typesDF, Stalk == TRUE & Tip == TRUE)),
    n123 = nrow(subset(typesDF, EC == TRUE & MC == TRUE & Phalanx == TRUE)),
    n124 = nrow(subset(typesDF, EC == TRUE & MC == TRUE & Stalk == TRUE)),
    n125 = nrow(subset(typesDF, EC == TRUE & MC == TRUE & Tip == TRUE)),
    n134 = nrow(subset(typesDF, EC == TRUE & Phalanx == TRUE & Stalk == TRUE)),
    n135 = nrow(subset(typesDF, EC == TRUE & Phalanx == TRUE & Tip == TRUE)),
    n145 = nrow(subset(typesDF, EC == TRUE & Stalk == TRUE & Tip == TRUE)),
    n234 = nrow(subset(typesDF, MC == TRUE & Phalanx == TRUE & Stalk == TRUE)),
    n235 = nrow(subset(typesDF, MC == TRUE & Phalanx == TRUE & Tip == TRUE)),
    n245 = nrow(subset(typesDF, MC == TRUE & Stalk == TRUE & Tip == TRUE)),
    n345 = nrow(subset(typesDF, Phalanx == TRUE & Stalk == TRUE & Tip == TRUE)),
    n1234 = nrow(subset(typesDF, EC == TRUE & MC == TRUE & Phalanx == TRUE & Stalk == TRUE)),
    n1235 = nrow(subset(typesDF, EC == TRUE & MC == TRUE & Phalanx == TRUE & Tip == TRUE)),
    n1245 = nrow(subset(typesDF, EC == TRUE & MC == TRUE & Stalk == TRUE & Tip == TRUE)),
    n1345 = nrow(subset(typesDF, EC == TRUE & Phalanx == TRUE & Stalk == TRUE & Tip == TRUE)),
    n2345 = nrow(subset(typesDF, MC == TRUE & Phalanx == TRUE & Stalk == TRUE & Tip == TRUE)),
    n12345 = nrow(subset(typesDF, EC == TRUE & MC == TRUE & Phalanx == TRUE & Stalk == TRUE & Tip == TRUE)),
    category = c("Endothelial", "Mesenchymal", "Phalanx", "Stalk", "Tip"),
    fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
    cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
    cat.dist = c(0.2, 0.2, 0.2, 0.2, 0.2),
    cat.pos = c(0, 330, 216, 144, 72),
    cat.cex = 1,
    margin = 0.05,
    cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8,1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5),
    ind = TRUE);
  # Writing to file
  pdf("Quintuple_Venn_diagram.pdf");
  grid.draw(venn.plot);
  dev.off();
}

# A continuous version of each update rule
W_AP1 <- function(WNT5b, SMAD2){
  return(max(WNT5b, SMAD2))
}

W_CTNNB <- function(WNT5b, WNT7a){
  return(max(WNT5b, WNT7a))
}

W_DLL4 <- function(DLL4){
  return(DLL4)
}

W_ETS1 <- function(VEGFR2, FGF2, ZEB2, AP1){
  return(max(VEGFR2, FGF2, ZEB2, AP1))
}

W_FGF2 <- function(FGF2){
  return(FGF2)
}

W_FLI1 <- function(FLI1,GATA2,HIF1a){
  return(max(FLI1,GATA2,HIF1a))
}

W_GATA2 <- function(FLI1, GATA2, VEGFR2){
  return(max(FLI1, GATA2, VEGFR2))
}

W_HIF1a <- function(HIF1a){
  return(HIF1a)
}

W_LEF1 <- function(LEF1, NRARP, SMAD2, NFkB, CTNNB){
  return(max(LEF1, NRARP, SMAD2, NFkB, CTNNB))
} 

W_NFkB <- function(PDGF_AB){
  return(PDGF_AB)
}

W_NOTCH <- function(DLL4, NRARP){
  return(min(DLL4, (1 - NRARP)))
} 

W_NRP1 <- function(VEGFA, ETS1, GATA2, NOTCH){# VEGFA & (ETS1 | GATA2) & !NOTCH
  return(min(VEGFA, max(ETS1, GATA2), 1 - NOTCH))
}

W_NRARP <- function(NOTCH){
  return(NOTCH)
}

W_PDGF_AB <- function(PDGF_AB){ #PDGF_AB
  return(PDGF_AB)
}

W_SMAD1 <- function(TGFBR, SMAD6){# TGFB & !SMAD6
  return(min(TGFBR, 1 - SMAD6))
}

W_SMAD2 <- function(TGFBR, SMAD6, NRP1){ # TGFB & !SMAD6 & NRP1
  return(min(TGFBR, 1 - SMAD6, NRP1))
}

W_SMAD6 <- function(FLI1, GATA2, NOTCH, SMAD1){# FLI1 & GATA2 & (NOTCH | SMAD1)
  return(min(FLI1, GATA2, max(NOTCH, SMAD1)))
}

W_SNAI1 <- function(NOTCH, CTNNB, HIF1a, STAT3, AP1, NFkB, SMAD2, TWIST1, SNAI1, SNAI2){# (NOTCH | CTNNB) & (HIF1a | STAT3 | CTNNB | AP1 | NFkB | SMAD2) & !TWIST1 & !SNAI1
  return(min(max(HIF1a, STAT3, CTNNB, AP1, NFkB, SMAD2), max(NOTCH, CTNNB), 1 - TWIST1, 1 - SNAI1, 1- SNAI2)) 
}

W_SNAI2 <- function(SNAI1, GATA2, SMAD2, SNAI2, TWIST1, LEF1, NOTCH){# !(SNAI1 & GATA2) & (SMAD2 | SNAI2 | TWIST1 | LEF1 | NOTCH)
  return(min(1 - min(SNAI1, GATA2), max(SMAD2, SNAI2, TWIST1, LEF1, NOTCH)))
}

W_STAT3 <- function(VEGFR2){
  return(VEGFR2)
}

W_TGFB <- function(TGFB){
  return(TGFB)
}

W_TGFBR<- function(TGFB, FGF2){
  min(TGFB, 1 - FGF2)
}

W_TWIST1 <- function(CTNNB, HIF1a){
  return(max(CTNNB, HIF1a))
} 

W_VEGFA <- function(WNT5b, TGFBR, STAT3, HIF1a, VEGFA){ #(!WNT5b & !TGFBR & (STAT3 | HIF1a)) | VEGFA
  return(max(min(1 - WNT5b , 1 - TGFBR , max(STAT3, HIF1a)), VEGFA))
}

W_VEGFR2 <- function(VEGFA, SNAI1, NOTCH, FLI1, GATA2, ETS1){# VEGFA & !SNAI1 & !NOTCH & (FLI1 | GATA2 | ETS1)
  return(min(VEGFA, 1 - SNAI1, 1 - NOTCH, max(FLI1, GATA2, ETS1)))
}

W_WNT5b <- function(WNT5b){
  return(WNT5b)
}

W_WNT7a <- function(WNT7a){
  return(WNT7a)
}

W_ZEB1 <- function(SNAI2, CTNNB, LEF1){
  return(max(SNAI2, CTNNB, LEF1))
}

W_ZEB2 <- function(HIF1a, ETS1, SMAD2, AP1, SNAI1){
  return(max(HIF1a, ETS1, SMAD2, AP1, SNAI1))
}

squadEndMT <- function(t, state, parameters){
  with(as.list(c(state, parameters)),{
    dAP1 <- squad(AP1, W_AP1(WNT5b, SMAD2), h, g)
    dCTNNB<- squad(CTNNB, W_CTNNB(WNT5b, WNT7a), h, g)
    dDLL4 <- squad(DLL4, W_DLL4(DLL4), h, g)
    dETS1 <- squad(ETS1, W_ETS1(VEGFR2, FGF2, ZEB2, AP1), h, g)
    dFGF2 <- squad(FGF2, W_FGF2(FGF2), h, g)
    dFLI1 <- squad(FLI1, W_FLI1(FLI1,GATA2,HIF1a), h, g)
    dGATA2 <- squad(GATA2, W_GATA2(FLI1, GATA2, VEGFR2), h, g)
    dHIF1a <- squad(HIF1a,W_HIF1a(HIF1a), h, g)
    dLEF1 <- squad(LEF1, W_LEF1(LEF1, NRARP, SMAD2, NFkB, CTNNB), h, g)
    dNFkB <- squad(NFkB, W_NFkB(PDGF_AB), h, g)
    dNOTCH <- squad(NOTCH, W_NOTCH(DLL4, NRARP), h, g)
    dNRARP <- squad(NRARP, W_NRARP(NOTCH), h, g)
    dNRP1 <- squad(NRP1, W_NRP1(VEGFA, ETS1, GATA2, NOTCH), h, g)
    dPDGF_AB <- squad(PDGF_AB, W_PDGF_AB(PDGF_AB), h, g)
    dSMAD1 <- squad(SMAD1, W_SMAD1(TGFBR, SMAD6), h, g)
    dSMAD2 <- squad(SMAD2, W_SMAD2(TGFBR, SMAD6, NRP1), h, g)
    dSMAD6 <- squad(SMAD6, W_SMAD6(FLI1, GATA2, NOTCH, SMAD1), h, g)
    dSNAI1 <- squad(SNAI1, W_SNAI1(NOTCH, CTNNB, HIF1a, STAT3, AP1, NFkB, SMAD2, TWIST1, SNAI1, SNAI2), h, g)
    dSNAI2 <- squad(SNAI2, W_SNAI2(SNAI1, GATA2, SMAD2, SNAI2, TWIST1, LEF1, NOTCH), h, g)
    dSTAT3 <- squad(STAT3, W_STAT3(VEGFR2), h, g)
    dTGFB <- squad(TGFB, W_TGFB(TGFB), h, g)
    dTGFBR <- squad(TGFBR, W_TGFBR(TGFB, FGF2), h, g)
    dTWIST1 <- squad(TWIST1, W_TWIST1(CTNNB, HIF1a), h, g)
    dVEGFA <- squad(VEGFA, W_VEGFA(WNT5b, TGFBR, STAT3, HIF1a, VEGFA), h, g)
    dVEGFR2 <- squad(VEGFR2, W_VEGFR2(VEGFA, SNAI1, NOTCH, FLI1, GATA2, ETS1), h, g)
    dWNT5b <- squad(WNT5b, W_WNT5b(WNT5b), h, g)
    dWNT7a <- squad(WNT7a, W_WNT7a(WNT7a), h, g)
    dZEB1 <- squad(ZEB1, W_ZEB1(SNAI2, CTNNB, LEF1), h, g)
    dZEB2 <- squad(ZEB2, W_ZEB2(HIF1a, ETS1, SMAD2, AP1, SNAI1), h, g)
    rate <- list(c(dAP1,
                   dCTNNB,
                   dDLL4,
                   dETS1,
                   dFGF2,
                   dFLI1,
                   dGATA2,
                   dHIF1a,
                   dLEF1,
                   dNFkB,
                   dNOTCH,
                   dNRP1,
                   dNRARP,
                   dPDGF_AB,
                   dSMAD1,
                   dSMAD2,
                   dSMAD6,
                   dSNAI1,
                   dSNAI2,
                   dSTAT3,
                   dTGFB,
                   dTGFBR,
                   dTWIST1,
                   dVEGFA,
                   dVEGFR2,
                   dWNT5b,
                   dWNT7a,
                   dZEB1,
                   dZEB2))
    return(rate)
  })
}

mutants_to_numbers <- function(mutants){
  m_numbers <- list()
  m_names <- mutants[, "mutants"]
  for(i in 1:nrow(mutants)){
    m_number <- 0
    for(j in 2:ncol(mutants)){
      if(mutants[i,j] > 0){
        m_number <- m_number + 2**(j-1)
      }
    }
    m_numbers[[i]] <- m_number
  }
  names(m_numbers) <- m_names
  return(m_numbers)
}




