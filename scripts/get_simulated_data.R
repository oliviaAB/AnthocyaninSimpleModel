here::i_am("scripts/get_simulated_data.R")

library(sismonr)
library(tidyverse)
library(here)

load(here("output/simulations.RData"))

## For the PC algorithm, input data must be in the form: rows = observations, columns = variables

simulated_data <- as_tibble(sim$Simulation) %>% 
  filter(time == max(time)) %>%  ## select the last time-point of the simulation
  select(starts_with("R")) %>% ## select the RNA abundances only
  rename_with(function(x){paste0("G", str_extract(x, "(?<=R)\\d+"))}) ## rename columns as G[gene ID]

write_csv(simulated_data, file = here("output/simulated_data.csv"))

