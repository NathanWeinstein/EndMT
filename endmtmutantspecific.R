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
attractor_df <- getAttractorsDF(typesDF, atr)

# AP1 mutations, gain of function not known
netAP1l <- fixGenes(endmtnet, "AP1", 0)
atrAP1l <- getAttractors(netAP1l, method="sat.exhaustive")
AP1l_typesDF <- sortAttractorsDF(atrAP1l, sortAttractor)
AP1l_attractor_df <- getAttractorsDF(AP1l_typesDF, atrAP1l)
SNAI1wt <- mean(attractor_df[,"SNAI1"])
SNAI1AP1l <- mean(AP1l_attractor_df[,"SNAI1"])
if(SNAI1wt > SNAI1AP1l){
  print("AP1l reproduced correctly")
}

# CTNNB mutations, CTNNBl multicelular effect
netCTNNBg <- fixGenes(endmtnet, "CTNNB", 1)
atrCTNNBg <- getAttractors(netCTNNBg, method="sat.exhaustive")
CTNNBg_typesDF <- sortAttractorsDF(atrCTNNBg, sortAttractor)
CTNNBg_attractor_df <- getAttractorsDF(CTNNBg_typesDF, atrCTNNBg)
SNAI2wt <- mean(attractor_df[,"SNAI2"])
SNAI2CTNNBg <- mean(CTNNBg_attractor_df[,"SNAI2"])
TWIST1wt <- mean(attractor_df[,"TWIST1"])
TWIST1CTNNBg <- mean(CTNNBg_attractor_df[,"TWIST1"])
ZEB1wt <- mean(attractor_df[,"ZEB1"])
ZEB1CTNNBg <- mean(CTNNBg_attractor_df[,"ZEB1"])
if((SNAI2wt < SNAI2CTNNBg) && (TWIST1wt < TWIST1CTNNBg) && (ZEB1wt < ZEB1CTNNBg)){
  print("CTNNBg reproduced correctly")
}

# DLL4g Reduces EC response to VEGF,
netDLL4g <- fixGenes(endmtnet, "DLL4", 1)
atrDLL4g <- getAttractors(netDLL4g, method="sat.exhaustive")
DLL4g_typesDF <- sortAttractorsDF(atrDLL4g, sortAttractor)
DLL4g_sizes <- get_group_sizes(DLL4g_typesDF, form_attractor_groups)
# No phalanx and no tip cell behaviors
DLL4g_attractor_df <- getAttractorsDF(DLL4g_typesDF, atrDLL4g)

# FGF2l Reduces VEGFR2 expression?
netFGF2l <- fixGenes(endmtnet, "FGF2", 0)
atrFGF2l <- getAttractors(netFGF2l, method="sat.exhaustive")
FGF2l_typesDF <- sortAttractorsDF(atrFGF2l, sortAttractor)
FGF2l_sizes <- get_group_sizes(FGF2l_typesDF, form_attractor_groups)
FGF2l_attractor_df <- getAttractorsDF(FGF2l_typesDF, atrFGF2l)
VEGFR2wt <- mean(attractor_df[,"VEGFR2"])
VEGFR2FGF2l <- mean(FGF2l_attractor_df[,"VEGFR2"])
if(VEGFR2wt > VEGFR2FGF2l){
  print("FGF2l Reduces VEGFR2 expression")
}

# FGF2g Supresses EndMT in some cases?
netFGF2g <- fixGenes(endmtnet, "FGF2", 1)
atrFGF2g <- getAttractors(netFGF2g, method="sat.exhaustive")
FGF2g_typesDF <- sortAttractorsDF(atrFGF2g, sortAttractor)
FGF2g_sizes <- get_group_sizes(FGF2g_typesDF, form_attractor_groups)
FGF2g_attractor_df <- getAttractorsDF(FGF2g_typesDF, atrFGF2g)
EndMT_FGF2g <- sum(FGF2g_sizes[c(5,7,8,9)])/sum(FGF2g_sizes)
EndMT_wt <- sum(sizes[c(5,7,8,9)])/sum(sizes)
if(EndMT_wt > EndMT_FGF2g){
  print("FGF2g Supresses EndMT in some cases")
}

# FLI1l iduces EndMT?
netFLI1l <- fixGenes(endmtnet, "FLI1", 0)
atrFLI1l <- getAttractors(netFLI1l, method="sat.exhaustive")
FLI1l_typesDF <- sortAttractorsDF(atrFLI1l, sortAttractor)
FLI1l_sizes <- get_group_sizes(FLI1l_typesDF, form_attractor_groups)
FLI1l_attractor_df <- getAttractorsDF(FLI1l_typesDF, atrFLI1l)
if(FLI1l_sizes["MCsnECs"]/sum(FLI1l_sizes) > sizes["MCsnECs"]/sum(sizes)){
  print("FLI1l induces EndMT")
}

# FLI1g Increases the expression of TIE2, GATA2, and VEGFR2?
netFLI1g <- fixGenes(endmtnet, "FLI1", 1)
atrFLI1g <- getAttractors(netFLI1g, method="sat.exhaustive")
FLI1g_typesDF <- sortAttractorsDF(atrFLI1g, sortAttractor)
FLI1g_sizes <- get_group_sizes(FLI1g_typesDF, form_attractor_groups)
FLI1g_attractor_df <- getAttractorsDF(FLI1g_typesDF, atrFLI1g)
GATA2wt <- mean(attractor_df[,"GATA2"])
GATA2_FLI1g <- mean(FLI1g_attractor_df[,"GATA2"])
VEGFR2_FLI1g <- mean(FLI1g_attractor_df[,"VEGFR2"])
if((GATA2wt < GATA2_FLI1g) && (VEGFR2wt < VEGFR2_FLI1g)) {
  print("FLI1g Increases the expression of GATA2, and VEGFR2")
}

# GATA2l triggers EndMT?
netGATA2l <- fixGenes(endmtnet, "GATA2", 0)
atrGATA2l <- getAttractors(netGATA2l, method="sat.exhaustive")
GATA2l_typesDF <- sortAttractorsDF(atrGATA2l, sortAttractor)
GATA2l_sizes <- get_group_sizes(GATA2l_typesDF, form_attractor_groups)
GATA2l_attractor_df <- getAttractorsDF(GATA2l_typesDF, atrGATA2l)
EndMT_GATA2l <- sum(GATA2l_sizes[c(5,7,8,9)])/sum(GATA2l_sizes)
if(EndMT_wt < EndMT_GATA2l){
  print("GATA2l induces EndMT")
}
wt_NRP1 <- mean(attractor_df[,"NRP1"])
GATA2l_NRP1 <- mean(GATA2l_attractor_df[,"NRP1"])
if(wt_NRP1 > GATA2l_NRP1){
  print("GATA2l inhibits NRP1 expression")
}
VEGFR2_GATA2l <- mean(GATA2l_attractor_df[,"VEGFR2"])
if(VEGFR2wt < VEGFR2_GATA2l){
  print("GATA2l inhibits VEGFR2 expression")
}

# HIF1l dencreases the expression of SNAI1 and TWIST1?
netHIF1l <- fixGenes(endmtnet, "HIF1a", 0)
atrHIF1l <- getAttractors(netHIF1l, method="sat.exhaustive")
HIF1l_typesDF <- sortAttractorsDF(atrHIF1l, sortAttractor)
HIF1l_sizes <- get_group_sizes(HIF1l_typesDF, form_attractor_groups)
HIF1l_attractor_df <- getAttractorsDF(HIF1l_typesDF, atrHIF1l)
SNAI1wt <- mean(attractor_df[,"SNAI1"])
SNAI1_HIF1l <- mean(HIF1l_attractor_df[,"SNAI1"])
TWIST1_HIF1l <- mean(HIF1l_attractor_df[,"TWIST1"])
EndMT_HIF1l <- sum(HIF1l_sizes[c(5,7,8,9)])/sum(HIF1l_sizes)
if(SNAI1wt > SNAI1_HIF1l){
  print("HIF1l dencreases the expression of SNAI1")
}
if(TWIST1wt > TWIST1_HIF1l){
  print("HIF1l dencreases the expression of TWIST1")
}
if(EndMT_wt > EndMT_HIF1l){
  print("HIF1l prevents EndMT")
}

# HIF1g inncreases the expression of SNAI1 and induces EndMT?
netHIF1g <- fixGenes(endmtnet, "HIF1a", 1)
atrHIF1g <- getAttractors(netHIF1g, method="sat.exhaustive")
HIF1g_typesDF <- sortAttractorsDF(atrHIF1g, sortAttractor)
HIF1g_sizes <- get_group_sizes(HIF1g_typesDF, form_attractor_groups)
HIF1g_attractor_df <- getAttractorsDF(HIF1g_typesDF, atrHIF1g)
SNAI1_HIF1g <- mean(HIF1g_attractor_df[,"SNAI1"])
EndMT_HIF1g <- sum(HIF1g_sizes[c(5,7,8,9)])/sum(HIF1g_sizes)
if(SNAI1wt < SNAI1_HIF1g){
  print("HIF1g increases the expression of SNAI1")
}
if(EndMT_wt < EndMT_HIF1g){
  print("HIF1g induces EndMT")
}

# LEF1g induces EndMT?
netLEF1g <- fixGenes(endmtnet, "LEF1", 1)
atrLEF1g <- getAttractors(netLEF1g, method="sat.exhaustive")
LEF1g_typesDF <- sortAttractorsDF(atrLEF1g, sortAttractor)
LEF1g_sizes <- get_group_sizes(LEF1g_typesDF, form_attractor_groups)
LEF1g_attractor_df <- getAttractorsDF(LEF1g_typesDF, atrLEF1g)
EndMT_LEF1g <- sum(LEF1g_sizes[c(5,7,8,9)])/sum(LEF1g_sizes)
if(EndMT_wt < EndMT_LEF1g){
  print("LEF1g induces EndMT")
}

# NFkBl inhibits SNAI1?
netNFkBl <- fixGenes(endmtnet, "NFkB", 0)
atrNFkBl <- getAttractors(netNFkBl, method="sat.exhaustive")
NFkBl_typesDF <- sortAttractorsDF(atrNFkBl, sortAttractor)
NFkBl_sizes <- get_group_sizes(NFkBl_typesDF, form_attractor_groups)
NFkBl_attractor_df <- getAttractorsDF(NFkBl_typesDF, atrNFkBl)
SNAI1_NFkBl <- mean(NFkBl_attractor_df[,"SNAI1"])
if(SNAI1wt > SNAI1_NFkBl){
  print("NFkBl inhibits SNAI1")
}

# NFkBg?
netNFkBg <- fixGenes(endmtnet, "NFkB", 1)
atrNFkBg <- getAttractors(netNFkBg, method="sat.exhaustive")
NFkBg_typesDF <- sortAttractorsDF(atrNFkBg, sortAttractor)
NFkBg_sizes <- get_group_sizes(NFkBg_typesDF, form_attractor_groups)
NFkBg_attractor_df <- getAttractorsDF(NFkBg_typesDF, atrNFkBg)
SNAI1_NFkBg <- mean(NFkBg_attractor_df[,"SNAI1"])
ZEB2_NFkBg <- mean(NFkBg_attractor_df[,"ZEB2"])
ZEB2wt <- mean(attractor_df[,"ZEB2"])
EndMT_NFkBg <- sum(NFkBg_sizes[c(5,7,8,9)])/sum(NFkBg_sizes)
if(SNAI1wt < SNAI1_NFkBg){
  print("NFkBg induces SNAI1 expression")
}
if(ZEB2wt < ZEB2_NFkBg){
  print("NFkBg induces ZEB2 expression")
}
if(EndMT_wt < EndMT_NFkBg){
  print("NFkBg induces EndMT")
}

# NOTCHl?
netNOTCHl <- fixGenes(endmtnet, "NOTCH", 0)
atrNOTCHl <- getAttractors(netNOTCHl, method="sat.exhaustive")
NOTCHl_typesDF <- sortAttractorsDF(atrNOTCHl, sortAttractor)
NOTCHl_sizes <- get_group_sizes(NOTCHl_typesDF, form_attractor_groups)
NOTCHl_attractor_df <- getAttractorsDF(NOTCHl_typesDF, atrNOTCHl)
wt_Tip <- sum(sizes[c(6,7)])/sum(sizes)
NOTCHl_Tip <- sum(NOTCHl_sizes[c(6,7)])/sum(NOTCHl_sizes)
if(wt_Tip < NOTCHl_Tip){
  print("NOTCHl promotes tip cell behavior")
}
NRP1_wt <- mean(attractor_df[,"NRP1"])
NRP1_NOTCHl <- mean(NOTCHl_attractor_df[,"NRP1"])
if(NRP1_wt < NRP1_NOTCHl){
  print("NOTCHl increases NRP1 expression")
}

# NOTCHg?
netNOTCHg <- fixGenes(endmtnet, "NOTCH", 1)
atrNOTCHg <- getAttractors(netNOTCHg, method="sat.exhaustive")
NOTCHg_typesDF <- sortAttractorsDF(atrNOTCHg, sortAttractor)
NOTCHg_sizes <- get_group_sizes(NOTCHg_typesDF, form_attractor_groups)
NOTCHg_attractor_df <- getAttractorsDF(NOTCHg_typesDF, atrNOTCHg)
SNAI1_NOTCHg <- mean(NOTCHg_attractor_df[,"SNAI1"])
if(SNAI1wt < SNAI1_NOTCHg){
  print("NOTCHg increases SNAI1 expression")
}
SNAI2_NOTCHg <- mean(NOTCHg_attractor_df[,"SNAI2"])
if(SNAI2wt < SNAI2_NOTCHg){
  print("NOTCHg increases SNAI2 expression")
}
wt_Stalk <- sum(sizes[c(4,5)])/sum(sizes)
NOTCHg_Stalk <- sum(NOTCHg_sizes[c(4,5)])/sum(NOTCHg_sizes)
if(wt_Stalk < NOTCHg_Stalk){
  print("NOTCHg promotes stalk cell behavior")
}

# NRP1l?
netNRP1l <- fixGenes(endmtnet, "NRP1", 0)
atrNRP1l <- getAttractors(netNRP1l, method="sat.exhaustive")
NRP1l_typesDF <- sortAttractorsDF(atrNRP1l, sortAttractor)
NRP1l_sizes <- get_group_sizes(NRP1l_typesDF, form_attractor_groups)
NRP1l_attractor_df <- getAttractorsDF(NRP1l_typesDF, atrNRP1l)
SMAD1wt <- mean(attractor_df[,"SMAD1"])
SMAD2wt <- mean(attractor_df[,"SMAD2"])
SMAD1NRP1l <- mean(NRP1l_attractor_df[,"SMAD1"])
SMAD2NRP1l <- mean(NRP1l_attractor_df[,"SMAD2"])
if(SMAD1wt < SMAD1NRP1l){
  print("NRP1l upregulates SMAD1")
}
if(SMAD2wt < SMAD2NRP1l){
  print("NRP1l downregulates SMAD2")
}

# NRP1g?
netNRP1g <- fixGenes(endmtnet, "NRP1", 1)
atrNRP1g <- getAttractors(netNRP1g, method="sat.exhaustive")
NRP1g_typesDF <- sortAttractorsDF(atrNRP1g, sortAttractor)
NRP1g_sizes <- get_group_sizes(NRP1g_typesDF, form_attractor_groups)
NRP1g_attractor_df <- getAttractorsDF(NRP1g_typesDF, atrNRP1g)
SMAD1NRP1g <- mean(NRP1g_attractor_df[,"SMAD1"])
SMAD2NRP1g <- mean(NRP1g_attractor_df[,"SMAD2"])
if(SMAD1wt > SMAD1NRP1g){
  print("NRP1g downregulates SMAD1")
}
if(SMAD2wt < SMAD2NRP1g){
  print("NRP1g upregulates SMAD2")
}
EndMT_NRP1g <- sum(NRP1g_sizes[c(5,7,8,9)])/sum(NRP1g_sizes)
if(EndMT_wt < EndMT_NRP1g){
  print("NRP1g induces EndMT")
}

# NRARPg inhibits NOTCH?
netNRARPg <- fixGenes(endmtnet, "NRARP", 1)
atrNRARPg <- getAttractors(netNRARPg, method="sat.exhaustive")
NRARPg_typesDF <- sortAttractorsDF(atrNRARPg, sortAttractor)
NRARPg_sizes <- get_group_sizes(NRARPg_typesDF, form_attractor_groups)
NRARPg_attractor_df <- getAttractorsDF(NRARPg_typesDF, atrNRARPg)
NOTCHwt <- mean(attractor_df[,"NOTCH"])
NOTCHNRARPg <- mean(NRARPg_attractor_df[,"NOTCH"])
if(NOTCHwt > NOTCHNRARPg){
  print("NRARPg inhibits NOTCH")
}

# PDGF_ABl inhibits EndMT?
netPDGFl <- fixGenes(endmtnet, "PDGF_AB", 0)
atrPDGFl <- getAttractors(netPDGFl, method="sat.exhaustive")
PDGFl_typesDF <- sortAttractorsDF(atrPDGFl, sortAttractor)
PDGFl_sizes <- get_group_sizes(PDGFl_typesDF, form_attractor_groups)
PDGFl_attractor_df <- getAttractorsDF(PDGFl_typesDF, atrPDGFl)
EndMT_PDGFl <- sum(PDGFl_sizes[c(5,7,8,9)])/sum(PDGFl_sizes)
if(EndMT_wt > EndMT_PDGFl){
  print("PDGF_ABl inhibits EndMT")
}

# PDGF_ABg induces EndMT?
netPDGFg <- fixGenes(endmtnet, "PDGF_AB", 1)
atrPDGFg <- getAttractors(netPDGFg, method="sat.exhaustive")
PDGFg_typesDF <- sortAttractorsDF(atrPDGFg, sortAttractor)
PDGFg_sizes <- get_group_sizes(PDGFg_typesDF, form_attractor_groups)
PDGFg_attractor_df <- getAttractorsDF(PDGFg_typesDF, atrPDGFg)
EndMT_PDGFg <- sum(PDGFg_sizes[c(5,7,8,9)])/sum(PDGFg_sizes)
if(EndMT_wt < EndMT_PDGFg){
  print("PDGF_ABg induces EndMT")
}
NFkBPDGFg <- mean(PDGFg_attractor_df[,"NFkB"])
NFkBwt <- mean(attractor_df[,"NFkB"])
if(NFkBPDGFg > NFkBwt){
  print("PDGFg upregulated NFkB")
}
SNAI1PDGFg <- mean(PDGFg_attractor_df[,"SNAI1"])
if(SNAI1PDGFg > SNAI1wt){
  print("PDGFg induces SNAI1 expression")
}

# SMAD1l inhibits SNAI1 expression?
netSMAD1l <- fixGenes(endmtnet, "SMAD1", 0)
atrSMAD1l <- getAttractors(netSMAD1l, method="sat.exhaustive")
SMAD1l_typesDF <- sortAttractorsDF(atrSMAD1l, sortAttractor)
SMAD1l_sizes <- get_group_sizes(SMAD1l_typesDF, form_attractor_groups)
SMAD1l_attractor_df <- getAttractorsDF(SMAD1l_typesDF, atrSMAD1l)
SNAI1SMAD1l <- mean(SMAD1l_attractor_df[,"SNAI1"])
if(SNAI1SMAD1l < SNAI1wt){
  print("SMAD1l inhibits SNAI1 expression")
}

# SMAD2l inhibits SNAI1 expression?
netSMAD2l <- fixGenes(endmtnet, "SMAD2", 0)
atrSMAD2l <- getAttractors(netSMAD2l, method="sat.exhaustive")
SMAD2l_typesDF <- sortAttractorsDF(atrSMAD2l, sortAttractor)
SMAD2l_sizes <- get_group_sizes(SMAD2l_typesDF, form_attractor_groups)
SMAD2l_attractor_df <- getAttractorsDF(SMAD2l_typesDF, atrSMAD2l)
SNAI1SMAD2l <- mean(SMAD2l_attractor_df[,"SNAI1"])
if(SNAI1SMAD2l < SNAI1wt){
  print("SMAD2l inhibits SNAI1 expression")
}
TGFSNAI1SMAD2l <- nrow(subset(SMAD2l_attractor_df, TGFB == 1 & SNAI1 == 1))/nrow(SMAD2l_attractor_df)
TGFSNAI1wt <- nrow(subset(attractor_df, TGFB == 1 & SNAI1 == 1))/nrow(attractor_df)
if(TGFSNAI1SMAD2l < TGFSNAI1wt){
  print("SMAD2l inhibits TGFB-mediated SNAI1 expression")
}

# SMAD1_2l inhibits SNAI1 expression?
netSMAD1_2l <- fixGenes(endmtnet, c("SMAD1","SMAD2"), c(0,0))
atrSMAD1_2l <- getAttractors(netSMAD1_2l, method="sat.exhaustive")
SMAD1_2l_typesDF <- sortAttractorsDF(atrSMAD1_2l, sortAttractor)
SMAD1_2l_sizes <- get_group_sizes(SMAD1_2l_typesDF, form_attractor_groups)
SMAD1_2l_attractor_df <- getAttractorsDF(SMAD1_2l_typesDF, atrSMAD1_2l)
SNAI1SMAD1_2l <- mean(SMAD1_2l_attractor_df[,"SNAI1"])
if(SNAI1SMAD1_2l < SNAI1wt){
  print("SMAD1_2l inhibits SNAI1 expression")
}
TGFSNAI1SMAD1_2l <- nrow(subset(SMAD1_2l_attractor_df, TGFB == 1 & SNAI1 == 1))/nrow(SMAD1_2l_attractor_df)
TGFSNAI1wt <- nrow(subset(attractor_df, TGFB == 1 & SNAI1 == 1))/nrow(attractor_df)
if(TGFSNAI1SMAD1_2l < TGFSNAI1wt){
  print("SMAD1_2l inhibits TGFB-mediated SNAI1 expression")
}

# SMAD2g inhibits SNAI1 expression?
netSMAD2g <- fixGenes(endmtnet, "SMAD2", 1)
atrSMAD2g <- getAttractors(netSMAD2g, method="sat.exhaustive")
SMAD2g_typesDF <- sortAttractorsDF(atrSMAD2g, sortAttractor)
SMAD2g_sizes <- get_group_sizes(SMAD2g_typesDF, form_attractor_groups)
SMAD2g_attractor_df <- getAttractorsDF(SMAD2g_typesDF, atrSMAD2g)
SNAI1SMAD2g <- mean(SMAD2g_attractor_df[,"SNAI2"])
if(SNAI1SMAD2g < SNAI1wt){
  print("SMAD2g inhibits SNAI1 expression")
}

# SMAD6l induces SMAD1 activation?
netSMAD6l <- fixGenes(endmtnet, "SMAD6", 0)
atrSMAD6l <- getAttractors(netSMAD6l, method="sat.exhaustive")
SMAD6l_typesDF <- sortAttractorsDF(atrSMAD6l, sortAttractor)
SMAD6l_sizes <- get_group_sizes(SMAD6l_typesDF, form_attractor_groups)
SMAD6l_attractor_df <- getAttractorsDF(SMAD6l_typesDF, atrSMAD6l)
SMAD1SMAD6l <- mean(SMAD6l_attractor_df[,"SMAD1"])
SMAD1wt <- mean(attractor_df[,"SMAD1"])
if(SMAD1SMAD6l > SMAD1wt){
  print("SMAD6l induces SMAD1 activation")
}

# SMAD6g induces SMAD1 activation?
netSMAD6g <- fixGenes(endmtnet, "SMAD6", 1)
atrSMAD6g <- getAttractors(netSMAD6g, method="sat.exhaustive")
SMAD6g_typesDF <- sortAttractorsDF(atrSMAD6g, sortAttractor)
SMAD6g_sizes <- get_group_sizes(SMAD6g_typesDF, form_attractor_groups)
SMAD6g_attractor_df <- getAttractorsDF(SMAD6g_typesDF, atrSMAD6g)
SMAD1SMAD6g <- mean(SMAD6g_attractor_df[,"SMAD1"])
SMAD1wt <- mean(attractor_df[,"SMAD1"])
if(SMAD1SMAD6g < SMAD1wt){
  print("SMAD6g inhibits SMAD1 activation")
}

# SNAI1g induces EndMT?
netSNAI1g <- fixGenes(endmtnet, "SNAI1", 1)
atrSNAI1g <- getAttractors(netSNAI1g, method="sat.exhaustive")
SNAI1g_typesDF <- sortAttractorsDF(atrSNAI1g, sortAttractor)
SNAI1g_sizes <- get_group_sizes(SNAI1g_typesDF, form_attractor_groups)
SNAI1g_attractor_df <- getAttractorsDF(SNAI1g_typesDF, atrSNAI1g)
EndMT_SNAI1g <- sum(SNAI1g_sizes[c(5,7,8,9)])/sum(SNAI1g_sizes)
if(EndMT_wt < EndMT_SNAI1g){
  print("SNAI1g induces EndMT")
}

# SNAI2g induces EndMT?
netSNAI2g <- fixGenes(endmtnet, "SNAI2", 1)
atrSNAI2g <- getAttractors(netSNAI2g, method="sat.exhaustive")
SNAI2g_typesDF <- sortAttractorsDF(atrSNAI2g, sortAttractor)
SNAI2g_sizes <- get_group_sizes(SNAI2g_typesDF, form_attractor_groups)
SNAI2g_attractor_df <- getAttractorsDF(SNAI2g_typesDF, atrSNAI2g)
EndMT_SNAI2g <- sum(SNAI2g_sizes[c(5,7,8,9)])/sum(SNAI2g_sizes)
if(EndMT_wt < EndMT_SNAI2g){
  print("SNAI2g induces EndMT")
}

# STAT3l indhibits SNAI1?
netSTAT3l <- fixGenes(endmtnet, "STAT3", 0)
atrSTAT3l <- getAttractors(netSTAT3l, method="sat.exhaustive")
STAT3l_typesDF <- sortAttractorsDF(atrSTAT3l, sortAttractor)
STAT3l_sizes <- get_group_sizes(STAT3l_typesDF, form_attractor_groups)
STAT3l_attractor_df <- getAttractorsDF(STAT3l_typesDF, atrSTAT3l)
SNAI1STAT3l <- mean(STAT3l_attractor_df[,"SNAI1"])
if(SNAI1STAT3l < SNAI1wt){
  print("STAT3l inhibits SNAI1")
}

# STAT3g upregulates SNAI1?
netSTAT3g <- fixGenes(endmtnet, "STAT3", 1)
atrSTAT3g <- getAttractors(netSTAT3g, method="sat.exhaustive")
STAT3g_typesDF <- sortAttractorsDF(atrSTAT3g, sortAttractor)
STAT3g_sizes <- get_group_sizes(STAT3g_typesDF, form_attractor_groups)
STAT3g_attractor_df <- getAttractorsDF(STAT3g_typesDF, atrSTAT3g)
SNAI1STAT3g <- mean(STAT3g_attractor_df[,"SNAI1"])
if(SNAI1STAT3g > SNAI1wt){
  print("STAT3g upregulates SNAI1")
}

# TGFBl inhibits EndMT?
netTGFBl <- fixGenes(endmtnet, "TGFB", 0)
atrTGFBl <- getAttractors(netTGFBl, method="sat.exhaustive")
TGFBl_typesDF <- sortAttractorsDF(atrTGFBl, sortAttractor)
TGFBl_sizes <- get_group_sizes(TGFBl_typesDF, form_attractor_groups)
TGFBl_attractor_df <- getAttractorsDF(TGFBl_typesDF, atrTGFBl)
EndMT_TGFBl <- sum(TGFBl_sizes[c(5,7,8,9)])/sum(TGFBl_sizes)
if(EndMT_wt > EndMT_TGFBl){
  print("TGFBl inhibits EndMT")
}

# TGFBg induces SNAI1 expression?
netTGFBg <- fixGenes(endmtnet, "TGFB", 1)
atrTGFBg <- getAttractors(netTGFBg, method="sat.exhaustive")
TGFBg_typesDF <- sortAttractorsDF(atrTGFBg, sortAttractor)
TGFBg_sizes <- get_group_sizes(TGFBg_typesDF, form_attractor_groups)
TGFBg_attractor_df <- getAttractorsDF(TGFBg_typesDF, atrTGFBg)
SNAI1TGFBg <- mean(TGFBg_attractor_df[,"SNAI1"])
if(SNAI1TGFBg > SNAI1wt){
  print("TGFBg upregulates SNAI1")
}
EndMT_TGFBg <- sum(TGFBg_sizes[c(5,7,8,9)])/sum(TGFBg_sizes)
if(EndMT_wt < EndMT_TGFBg){
  print("TGFBg induces EndMT")
}

# TGFBRl inhibits EndMT?
netTGFBRl <- fixGenes(endmtnet, "TGFBR", 0)
atrTGFBRl <- getAttractors(netTGFBRl, method="sat.exhaustive")
TGFBRl_typesDF <- sortAttractorsDF(atrTGFBRl, sortAttractor)
TGFBRl_sizes <- get_group_sizes(TGFBRl_typesDF, form_attractor_groups)
TGFBRl_attractor_df <- getAttractorsDF(TGFBRl_typesDF, atrTGFBRl)
EndMT_TGFBRl <- sum(TGFBRl_sizes[c(5,7,8,9)])/sum(TGFBRl_sizes)
if(EndMT_wt > EndMT_TGFBRl){
  print("TGFBRl inhibits EndMT")
}

# TGFBRg induces EndMT?
netTGFBRg <- fixGenes(endmtnet, "TGFBR", 1)
atrTGFBRg <- getAttractors(netTGFBRg, method="sat.exhaustive")
TGFBRg_typesDF <- sortAttractorsDF(atrTGFBRg, sortAttractor)
TGFBRg_sizes <- get_group_sizes(TGFBRg_typesDF, form_attractor_groups)
TGFBRg_attractor_df <- getAttractorsDF(TGFBRg_typesDF, atrTGFBRg)
EndMT_TGFBRg <- sum(TGFBRg_sizes[c(5,7,8,9)])/sum(TGFBRg_sizes)
if(EndMT_wt < EndMT_TGFBRg){
  print("TGFBRg inducess EndMT")
}
SMAD2_TGFBRg <- mean(TGFBRg_attractor_df[,"SMAD2"])
if(SMAD2_TGFBRg > SMAD2wt){
  print("TGFBRg upregulates SMAD2")
}

# TWIST1l inhibits EndMT?
netTWIST1l <- fixGenes(endmtnet, "TWIST1", 0)
atrTWIST1l <- getAttractors(netTWIST1l, method="sat.exhaustive")
TWIST1l_typesDF <- sortAttractorsDF(atrTWIST1l, sortAttractor)
TWIST1l_sizes <- get_group_sizes(TWIST1l_typesDF, form_attractor_groups)
TWIST1l_attractor_df <- getAttractorsDF(TWIST1l_typesDF, atrTWIST1l)
EndMT_TWIST1l <- sum(TWIST1l_sizes[c(5,7,8,9)])/sum(TWIST1l_sizes)
if(EndMT_wt > EndMT_TWIST1l){
  print("TWIST1l inhibits EndMT")
}

# TWIST1g induceds EndMT?
netTWIST1g <- fixGenes(endmtnet, "TWIST1", 1)
atrTWIST1g <- getAttractors(netTWIST1g, method="sat.exhaustive")
TWIST1g_typesDF <- sortAttractorsDF(atrTWIST1g, sortAttractor)
TWIST1g_sizes <- get_group_sizes(TWIST1g_typesDF, form_attractor_groups)
TWIST1g_attractor_df <- getAttractorsDF(TWIST1g_typesDF, atrTWIST1g)
EndMT_TWIST1g <- sum(TWIST1g_sizes[c(5,7,8,9)])/sum(TWIST1g_sizes)
if(EndMT_wt < EndMT_TWIST1g){
  print("TWIST1g induces EndMT")
}
SMAD2_TWIST1g <- mean(TWIST1g_attractor_df[,"SMAD2"])
if(SMAD2_TWIST1g > SMAD2wt){
  print("TWIST1g upregulates SMAD2")
}
TGFBR_TWIST1g <- mean(TWIST1g_attractor_df[,"TGFBR"])
TGFBRwt <- mean(attractor_df[,"TGFBR"])
if(TGFBR_TWIST1g > TGFBRwt){
  print("TWIST1g upregulates TGFBR")
}

# VEGFAl?
netVEGFAl <- fixGenes(endmtnet, "VEGFA", 0)
atrVEGFAl <- getAttractors(netVEGFAl, method="sat.exhaustive")
VEGFAl_typesDF <- sortAttractorsDF(atrVEGFAl, sortAttractor)
VEGFAl_sizes <- get_group_sizes(VEGFAl_typesDF, form_attractor_groups)
VEGFAl_attractor_df <- getAttractorsDF(VEGFAl_typesDF, atrVEGFAl)
VEGFR2VEGFAl <- mean(VEGFAl_attractor_df[,"VEGFR2"])
VEGFR2wt <- mean(attractor_df[,"VEGFR2"])
if(VEGFR2VEGFAl < VEGFR2wt){
  print("VEGFRAl downregulates VEGFR2")
}

# VEGFAg
netVEGFAg <- fixGenes(endmtnet, "VEGFA", 1)
atrVEGFAg <- getAttractors(netVEGFAg, method="sat.exhaustive")
VEGFAg_typesDF <- sortAttractorsDF(atrVEGFAg, sortAttractor)
VEGFAg_sizes <- get_group_sizes(VEGFAg_typesDF, form_attractor_groups)
VEGFAg_attractor_df <- getAttractorsDF(VEGFAg_typesDF, atrVEGFAg)
VEGFR2VEGFAg <- mean(VEGFAl_attractor_df[,"VEGFR2"])
EndMT_VEGFAg <- sum(VEGFAg_sizes[c(5,7,8,9)])/sum(VEGFAg_sizes)
if(EndMT_wt > EndMT_VEGFAg){
  print("VEGFAg indhibits EndMT")
}

# WNT5bl
netWNT5bl <- fixGenes(endmtnet, "WNT5b", 0)
atrWNT5bl <- getAttractors(netWNT5bl, method="sat.exhaustive")
WNT5bl_typesDF <- sortAttractorsDF(atrWNT5bl, sortAttractor)
WNT5bl_sizes <- get_group_sizes(WNT5bl_typesDF, form_attractor_groups)
WNT5bl_attractor_df <- getAttractorsDF(WNT5bl_typesDF, atrWNT5bl)
EndMT_WNT5bl <- sum(WNT5bl_sizes[c(5,7,8,9)])/sum(WNT5bl_sizes)
if(EndMT_wt > EndMT_WNT5bl){
  print("WNT5bl inhibits EndMT")
}

# WNT5bg
netWNT5bg <- fixGenes(endmtnet, "WNT5b", 1)
atrWNT5bg <- getAttractors(netWNT5bg, method="sat.exhaustive")
WNT5bg_typesDF <- sortAttractorsDF(atrWNT5bg, sortAttractor)
WNT5bg_sizes <- get_group_sizes(WNT5bg_typesDF, form_attractor_groups)
WNT5bg_attractor_df <- getAttractorsDF(WNT5bg_typesDF, atrWNT5bg)
EndMT_WNT5bg <- sum(WNT5bg_sizes[c(5,7,8,9)])/sum(WNT5bg_sizes)
if(EndMT_wt < EndMT_WNT5bg){
  print("WNT5bg induces EndMT")
}

# WNT7ag
netWNT7ag <- fixGenes(endmtnet, "WNT7a", 1)
atrWNT7ag <- getAttractors(netWNT7ag, method="sat.exhaustive")
WNT7ag_typesDF <- sortAttractorsDF(atrWNT7ag, sortAttractor)
WNT7ag_sizes <- get_group_sizes(WNT7ag_typesDF, form_attractor_groups)
WNT7ag_attractor_df <- getAttractorsDF(WNT7ag_typesDF, atrWNT7ag)
TWIST1_WNT5bg <- mean(WNT7ag_attractor_df[,"TWIST1"])
if(TWIST1wt < TWIST1_WNT5bg){
  print("WNT7ag induces TWIST1 expression")
}

# ZEB1l
netZEB1l <- fixGenes(endmtnet, "ZEB1", 0)
atrZEB1l <- getAttractors(netZEB1l, method="sat.exhaustive")
ZEB1l_typesDF <- sortAttractorsDF(atrZEB1l, sortAttractor)
ZEB1l_sizes <- get_group_sizes(ZEB1l_typesDF, form_attractor_groups)
ZEB1l_attractor_df <- getAttractorsDF(ZEB1l_typesDF, atrZEB1l)
EndMT_ZEB1l <- sum(ZEB1l_sizes[c(5,7,8,9)])/sum(ZEB1l_sizes)
if(EndMT_wt > EndMT_ZEB1l){
  print("ZEB1l inhibits EndMT")
}

# ZEB1g
netZEB1g <- fixGenes(endmtnet, "ZEB1", 1)
atrZEB1g <- getAttractors(netZEB1g, method="sat.exhaustive")
ZEB1g_typesDF <- sortAttractorsDF(atrZEB1g, sortAttractor)
ZEB1g_sizes <- get_group_sizes(ZEB1g_typesDF, form_attractor_groups)
ZEB1g_attractor_df <- getAttractorsDF(ZEB1g_typesDF, atrZEB1g)
EndMT_ZEB1g <- sum(ZEB1g_sizes[c(5,7,8,9)])/sum(ZEB1g_sizes)
if(EndMT_wt < EndMT_ZEB1g){
  print("ZEB1g induces EndMT")
}

# ZEB2l
netZEB2l <- fixGenes(endmtnet, "ZEB2", 0)
atrZEB2l <- getAttractors(netZEB2l, method="sat.exhaustive")
ZEB2l_typesDF <- sortAttractorsDF(atrZEB2l, sortAttractor)
ZEB2l_sizes <- get_group_sizes(ZEB2l_typesDF, form_attractor_groups)
ZEB2l_attractor_df <- getAttractorsDF(ZEB2l_typesDF, atrZEB2l)
EndMT_ZEB2l <- sum(ZEB2l_sizes[c(5,7,8,9)])/sum(ZEB2l_sizes)
if(EndMT_wt > EndMT_ZEB2l){
  print("ZEB2l inhibits EndMT")
}

# ZEB2g
netZEB2g <- fixGenes(endmtnet, "ZEB2", 1)
atrZEB2g <- getAttractors(netZEB2g, method="sat.exhaustive")
ZEB2g_typesDF <- sortAttractorsDF(atrZEB2g, sortAttractor)
ZEB2g_sizes <- get_group_sizes(ZEB2g_typesDF, form_attractor_groups)
ZEB2g_attractor_df <- getAttractorsDF(ZEB2g_typesDF, atrZEB1g)
EndMT_ZEB2g <- sum(ZEB2g_sizes[c(5,7,8,9)])/sum(ZEB2g_sizes)
if(EndMT_wt < EndMT_ZEB2g){
  print("ZEB2g induces EndMT")
}

# ETS1gFLI1g, is equivalent to uniform laminar shear stress
netETS1gFLI1g <- fixGenes(endmtnet, c("ETS1", "FLI1"), c(1,1))
atrETS1gFLI1g <- getAttractors(netETS1gFLI1g, method="sat.exhaustive")
ETS1gFLI1g_typesDF <- sortAttractorsDF(atrETS1gFLI1g, sortAttractor)
ETS1gFLI1g_sizes <- get_group_sizes(ETS1gFLI1g_typesDF, form_attractor_groups)
ETS1gFLI1g_attractor_df <- getAttractorsDF(ETS1gFLI1g_typesDF, atrETS1gFLI1g)
EndMT_ETS1gFLI1g <- sum(ETS1gFLI1g_sizes[c(5,7,8,9)])/sum(ETS1gFLI1g_sizes)
if(ETS1gFLI1g_sizes["MCsnECs"] == 0){
  print("Uniform laminar shear stress prevents full EndMT")
  if(EndMT_wt < EndMT_ETS1gFLI1g){
    print("Uniform laminar shear stress induces partial EndMT")
  }
}


# ETS1lFLI1l, is equivalent to low, unstable shear stress
netETS1lFLI1l <- fixGenes(endmtnet, c("ETS1", "FLI1"), c(0,0))
atrETS1lFLI1l <- getAttractors(netETS1lFLI1l, method="sat.exhaustive")
ETS1lFLI1l_typesDF <- sortAttractorsDF(atrETS1lFLI1l, sortAttractor)
ETS1lFLI1l_sizes <- get_group_sizes(ETS1lFLI1l_typesDF, form_attractor_groups)
ETS1lFLI1l_attractor_df <- getAttractorsDF(ETS1lFLI1l_typesDF, atrETS1lFLI1l)
if(ETS1lFLI1l_sizes["MCsnECs"] > sizes["MCsnECs"]){
  print("Low shear stress induces full EndMT")
}
EndMT_ETS1lFLI1l <- sum(ETS1lFLI1l_sizes[c(5,7,8,9)])/sum(ETS1lFLI1l_sizes)
if(EndMT_wt < EndMT_ETS1lFLI1l){
  print("Low laminar shear stress induces EndMT")
}