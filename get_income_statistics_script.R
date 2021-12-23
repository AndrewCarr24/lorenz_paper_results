library(tidyverse)
library(readr)

# Reading PUMS data from 2011-2015 ACS and dropping incomes outside PUMAs
pums_data <- read_csv("pums_data.csv") %>% 
  rename(area = PUMID) %>% select(area, income, wt) %>% filter(area != 0)

# Reading get_income_statistics function
source("get_income_statistics.R")

# Reading MCIB function
source("mcib_func.R")

# Reading lorenz interpolation function 
source("lorenz_int.R")

# Get income statistics estimated by Lorenz interpolation, CDF interpolation, and MCIB
income_statistics <- get_income_statistics(input_data = pums_data, input_func = lorenz_int, num_cores = 11)














