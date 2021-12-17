here::i_am("scripts/create_model.R")

library(sismonr)
library(tidyverse)
library(here)

## ----------------------------- ##
## Creating the in silico system ##
## ----------------------------- ##

## Gene ID - name correspondence
genes_name2id <- tibble("ID" = as.character(1:9),
                        "name" = c("MYB", ## 1
                                   "bHLH1", ## 2
                                   "WDR", ## 3
                                   "MBW1", ## 4
                                   "bHLH2", ## 5
                                   "MBW2", ## 6
                                   "MYBrep", ## 7
                                   "R3-MYB", ## 8
                                   "DFR"), ## 9
                        "colour" = c("#70AD47", ## 1
                                     "#FFC000", ## 2
                                     "#C38649", ## 3
                                     "#E37D00", ## 4
                                     "#9DC3E6", ## 5
                                     "#0015B0", ## 6
                                     "#92D050", ## 7 
                                     "#FFCDD8", ## 8
                                     "#FF0000")) ## 9

id2names <- genes_name2id$name
names(id2names) <- genes_name2id$ID

colours <- genes_name2id$colour
names(colours) <- genes_name2id$ID

## ----------------------------- ##
## Creating the in silico system ##
## ----------------------------- ##

## We create a system with 9 genes, and no regulatory
## interactions (they will be added manually)
colsystem <- createInSilicoSystem(empty = T,
                                  G = 9,
                                  PC.p = 1,
                                  ## all genes are regulators of transcription:
                                  PC.TC.p = 1, 
                                  PC.TL.p = 0,
                                  PC.RD.p = 0,
                                  PC.PD.p = 0,
                                  PC.PTM.p = 0,
                                  PC.MR.p = 0,
                                  ploidy = 1)

## Changing the kinetic parameters of the genes
colsystem$genes$TCrate <- c(5, 0.5, 0.5, 0.01, 0.1, 0.01, 0.01, 0.1, 0.5)
colsystem$genes$TLrate <- c(0.1, 0.002, 0.01, 0.01, 0.001, 0.01, 0.01, 0.01, 0.001)
colsystem$genes$RDrate <- c(0.1, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01)
colsystem$genes$PDrate <- c(0.01, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001)

## Adding  regulatory interactions in the system
interactions <- list(list("edge" = c(1, 4),
                          "regsign" = "1",
                          "kinetics" = list("TCbindingrate" = 0.1,
                                            "TCunbindingrate" = 2,
                                            "TCfoldchange" = 10)),
                     list("edge" = c(2, 4),
                          "regsign" = "1",
                          "kinetics" = list("TCbindingrate" = 0.1,
                                            "TCunbindingrate" = 2,
                                            "TCfoldchange" = 10)),
                     list("edge" = c(3, 4),
                          "regsign" = "1",
                          "kinetics" = list("TCbindingrate" = 0.1,
                                            "TCunbindingrate" = 2,
                                            "TCfoldchange" = 10)),
                     list("edge" = c(4, 5),
                          "regsign" = "1",
                          "kinetics" = list("TCbindingrate" = 0.1,
                                            "TCunbindingrate" = 2,
                                            "TCfoldchange" = 25)),
                     list("edge" = c(1, 6),
                          "regsign" = "1",
                          "kinetics" = list("TCbindingrate" = 0.1,
                                            "TCunbindingrate" = 2,
                                            "TCfoldchange" = 10)),
                     list("edge" = c(3, 6),
                          "regsign" = "1",
                          "kinetics" = list("TCbindingrate" = 0.1,
                                            "TCunbindingrate" = 2,
                                            "TCfoldchange" = 10)),
                     list("edge" = c(5, 6),
                          "regsign" = "1",
                          "kinetics" = list("TCbindingrate" = 0.1,
                                            "TCunbindingrate" = 2,
                                            "TCfoldchange" = 10)),
                     list("edge" = c(6, 7),
                          "regsign" = "1",
                          "kinetics" = list("TCbindingrate" = 0.1,
                                            "TCunbindingrate" = 2,
                                            "TCfoldchange" = 50)),
                     list("edge" = c(6, 8),
                          "regsign" = "1",
                          "kinetics" = list("TCbindingrate" = 0.1,
                                            "TCunbindingrate" = 2,
                                            "TCfoldchange" = 50)),
                     list("edge" = c(6, 9),
                          "regsign" = "1",
                          "kinetics" = list("TCbindingrate" = 0.1,
                                            "TCunbindingrate" = 2,
                                            "TCfoldchange" = 15)),
                     list("edge" = c(7, 5),
                          "regsign" = "-1",
                          "kinetics" = list("TCbindingrate" = 0.05,
                                            "TCunbindingrate" = 2)),
                     list("edge" = c(8, 6),
                          "regsign" = "-1",
                          "kinetics" = list("TCbindingrate" = 0.01,
                                            "TCunbindingrate" = 2)),
                     list("edge" = c(8, 4),
                          "regsign" = "-1",
                          "kinetics" = list("TCbindingrate" = 0.01,
                                            "TCunbindingrate" = 2)))

for(inter in interactions){
  colsystem <- addEdge(colsystem,
                       inter$edge[1],
                       inter$edge[2],
                       regsign = inter$regsign,
                       kinetics = inter$kinetics)
}


plotGRN(colsystem, nodes_label = id2names)



## ---------------------------------- ##
## Creating the in silico individuals ##
## ---------------------------------- ##

colpop <- createInSilicoPopulation(3000,
                                   colsystem,
                                   initialNoise = F,
                                   ngenevariants = 5)


## Changing the initial conditions:
## only bHLH1 and WDR are constitutively expressed.
for(i in names(colpop$individualsList)){
  colpop$individualsList[[i]]$InitAbundance$GCN1$R =
    colpop$individualsList[[i]]$InitAbundance$GCN1$R * 
    c(0, 1, 1, 0, 0, 0, 0, 0, 0)
  colpop$individualsList[[i]]$InitAbundance$GCN1$P =
    colpop$individualsList[[i]]$InitAbundance$GCN1$P * 
    c(0, 1, 1, 0, 0, 0, 0, 0, 0)
}


colpop_1ind <- colpop
colpop_1ind$individualsList <- colpop_1ind$individualsList[1]

## --------------------------------- ##
## Saving the system and individuals ##
## --------------------------------- ##

save(colsystem, colpop, id2names, colours, file = here("output/sismonr_anthocyanin_model.RData"))

## ------------------------------------------------------------------ ##
## Simulating the expression profiles of the genes for one individual ##
## ------------------------------------------------------------------ ##

set.seed(123)
sim <- simulateInSilicoSystem(colsystem,
                              colpop_1ind, 
                              simtime = 3000,
                              ntrials = 1)
sim$runningtime / 60
sum(sim$runningtime)

plotSimulation(sim$Simulation, 
               labels = id2names,
               colours = colours)

ggsave(here("output/plot_one_simulation.png"), width = 12, height = 8)
