here::i_am("scripts/get_simulated_data.R")

library(sismonr)
library(tidyverse)
library(here)

load(here("output/simulations.RData"))

## For the pcalg package, input data must be in the form: rows = observations, columns = variables

simulated_data <- as_tibble(sim$Simulation) %>% 
  filter(time == max(time)) %>%  ## select the last time-point of the simulation
  select(starts_with("R")) %>% ## select the RNA abundances only
  rename_with(function(x){paste0("G", str_extract(x, "(?<=R)\\d+"))}) ## rename columns as G[gene ID]

write_csv(simulated_data, file = here("output/simulated_data.csv"))


## Making a plot of the network
load(here("output/sismonr_anthocyanin_model.RData"))

nodes_labels <- paste0("G", 1:ncol(simulated_data))
names(nodes_labels) <- paste0(1:ncol(simulated_data))

plotGRN(colsystem, nodes_label = nodes_labels) ## saved a screen capture as output/model_grn.png

## Getting the adjacency matrix of the network
true_adjacency_matrix <- matrix(data = 0,
                                nrow = ncol(simulated_data),
                                ncol = ncol(simulated_data))

rownames(true_adjacency_matrix) <- colnames(true_adjacency_matrix) <- nodes_labels

edges_list <- getEdges(colsystem) %>% 
  select(from, to) %>% 
  mutate_all(as.numeric)

for(i in 1:nrow(edges_list)) true_adjacency_matrix[edges_list$from[i], edges_list$to[i]] <- 1

write.csv(true_adjacency_matrix,
          file = here("output/true_adjacency_matrix.csv"), 
          row.names = TRUE)

## Read this adjacency matrix with:
## as.matrix(read.csv(here("output/true_adjacency_matrix.csv"), row.names = 1, header = TRUE))
