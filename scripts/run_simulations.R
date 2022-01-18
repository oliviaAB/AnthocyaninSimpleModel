here::i_am("scripts/run_simulations.R")

library(sismonr)
library(here)

#load(here("output/sismonr_anthocyanin_model.RData"))
load(here("output/sismonr_anthocyanin_model_no_cycles.RData"))

set.seed(123)
sim <- simulateParallelInSilicoSystem(colsystem,
                                      colpop, 
                                      simtime = 2000,
                                      ntrials = 1,
                                      no_cores = 4)

#save(sim, file = here("output/summer_project_rupert_2021/simulations.RData"))             
save(sim, file = here("output/summer_project_rupert_2021/simulations_no_cycles.RData"))             