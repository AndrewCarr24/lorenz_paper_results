library(tidyverse)
library(readr)

# Reading PUMS data from 2011-2015 ACS and dropping incomes outside PUMAs
setwd("~/lorenz_revision3/code")
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

























# results_test %>% mutate(row_idx = 1:n()) %>% select(row_idx, everything()) %>% 
#   mutate(theil_abs_error = abs(lor_theil - true_theil), theil_error = lor_theil - true_theil) %>% 
#   arrange(desc(theil_abs_error)) %>% select(row_idx, true_theil, lor_theil, theil_error, lor_mean_errors, mcib_mean_errors) %>% 
#   .[1,c("lor_mean_errors", "mcib_mean_errors")] %>% map(., ~.x) %>% .[[1]] %>% .[[1]] %>% .[1:15] %>% sum
# 
# 
# results_test$mcib_mean_errors[[1]] %>% .[1:15] %>% sum
# results_test$lor_mean_errors[[1]] %>% .[1:15] %>% sum
# 
# 
# Reduce(`+`, map2(results_test$true_bin_means, results_test$lor_bin_means, ~.y-.x))/200
# Reduce(`+`, map2(results_test$true_bin_means, results_test$mcib_bin_means, ~.y-.x))/200
# 
# 
# 


###################################################################################################
#### Modify lorenz_v5 (changing lorenz interpolation method) and test again #####
##########################################################################################



# # Changing interpolation method 
# lorenz_to_coefs_cont2 <- function(lorenz_df, bounds, G){
#   
#   coefs_lst <- list()
#   
#   for(idx in 1:(nrow(lorenz_df)-2)){
#     
#     x_1 <- lorenz_df$x[idx]
#     x_2 <- lorenz_df$x[idx+1]
#     
#     y_1 <- lorenz_df$y[idx]
#     y_2 <- lorenz_df$y[idx+1]
#     
#     A <- as.matrix(cbind(c(x_1^3, x_2^3, 3*G*x_1^2, 3*G*x_2^2), c(x_1^2, x_2^2, 2*G*x_1, 2*G*x_2), c(x_1, x_2, G, G), c(1, 1, 0, 0)))
#     y <- as.matrix(c(y_1, y_2, bounds[idx], bounds[idx+1]))
#     
#     coefs_lst[[idx]] <- solve(A, y)
#     
#   }
#   
#   return(coefs_lst)
#   
# }
# 
# fit_poly_to_top_cont <- function(lorenz_df, lorenz_coefs, slope_factor){
#   
#   
#   # Fitting quadratic function with start slope to the top category
#   x_1 = lorenz_df$x[(length(lorenz_df$x)-1)]
#   x_2 = lorenz_df$x[(length(lorenz_df$x))]
#   start_slope = 3*lorenz_coefs[[nrow(lorenz_df)-2]][1]*x_1^2 + 2*lorenz_coefs[[nrow(lorenz_df)-2]][2]*x_1 + lorenz_coefs[[nrow(lorenz_df)-2]][3]
#   
#   A = matrix(c(x_1^2, x_2^2, 2*x_1,
#                x_1, x_2, 1,
#                1, 1, 0), 3, 3)
#   
#   y_1 <- lorenz_df$y[(length(lorenz_df$y)-1)]
#   y_2 <- lorenz_df$y[(length(lorenz_df$y))]
#   y = matrix(c(y_1, y_2, start_slope), 3, 1)
#   quad_coefs <- solve(A, y)
#   
#   mid_x <- (lorenz_df$x[(length(lorenz_df$x)-1)] + lorenz_df$x[(length(lorenz_df$x))])/2
#   slope_mid_x <- 2*quad_coefs[1]*mid_x + quad_coefs[2]
#   
#   # Cubic to two points and two slopes
#   x_1 = lorenz_df$x[(length(lorenz_df$x)-1)]
#   x_2 = lorenz_df$x[(length(lorenz_df$x))]
#   
#   A = matrix(c(x_1^3, x_2^3, 3*x_1^2, 3*mid_x^2,
#                x_1^2, x_2^2, 2*x_1, 2*mid_x,
#                x_1, x_2, 1, 1,
#                1, 1, 0, 0), 4, 4)
#   
#   y = matrix(c(lorenz_df$y[(length(lorenz_df$y)-1)],
#                lorenz_df$y[length(lorenz_df$y)],
#                start_slope,
#                slope_factor*slope_mid_x), 4, 1)
#   
#   return(solve(A, y))
# }
# 
# # Changing method of getting incs
# get_unweighted_incs_mod <- function(lorenz_df, lorenz_coefs, num_cats, num_splits, G){
#   
#   map(1:num_cats, function(idx_main){
#     
#     xs <- seq(lorenz_df$x[idx_main], lorenz_df$x[idx_main+1], (lorenz_df$x[idx_main+1]-lorenz_df$x[idx_main])/num_splits)
#     ys <- lorenz_coefs[[idx_main]][1]*xs^3 + lorenz_coefs[[idx_main]][2]*xs^2 + lorenz_coefs[[idx_main]][3]*xs + lorenz_coefs[[idx_main]][4]
#     slopes <- map(1:num_splits, function(idx){ (ys[idx+1]-ys[idx])/(xs[idx+1]-xs[idx]) }) %>% unlist
#     
#     return(slopes*G)
#     
#   })
#   
# }
# 
# lorenz_v5_mod <- function(freqs, bounds, G, inc_share = NULL, stat = "gini", slope_parm = .9){
#   
#   mcib_means <- get_mcib_means(freqs, bounds)
#   mcib_means[length(mcib_means)+1] <- if(freqs[16]>0){ (G*sum(freqs) - sum(freqs[1:15]*mcib_means))/freqs[16] }else{0}
#   
#   mcib_means <- mean_correct(freqs, bounds, mcib_means, G)
#   lorenz_df <- tibble(x = c(0, cumsum(freqs)/sum(freqs)), y = c(0, cumsum(freqs*mcib_means)/sum(freqs*mcib_means))) %>% distinct()
#   
#   # Adding income shares (if provided)
#   if(!is.null(inc_share[[1]][1])){
#     markers <- inc_share[[1]]
#     shares <- inc_share[[2]]
#     lorenz_df <- bind_rows(lorenz_df %>% tibble::add_column(type='est'), 
#                            tibble(x = markers, y = shares) %>% tibble::add_column(type='true')) %>% 
#       arrange(x) %>% 
#       filter(!((x %in% markers) & (type == "est"))) %>% 
#       filter(!ranges_func(y, inc_shares_for_lorenz) | type == "true") 
#   }
#   
#   lorenz_coefs <- lorenz_to_coefs_cont2(lorenz_df, bounds, G)
#   
#   lorenz_coefs[[(length(lorenz_coefs)+1)]] <- fit_poly_to_top_cont(lorenz_df, lorenz_coefs, slope_parm)
#   
#   # ## Interpolate to get ineq stats
#   num_splits <- 100
#   
#   inc_lst <- get_unweighted_incs_mod(lorenz_df, lorenz_coefs, length(lorenz_coefs), num_splits, G)
#   # wts <- freqs[freqs!=0]/sum(freqs[freqs!=0])
#   wts <- diff(lorenz_df$x)
#   inc_lst <- map(inc_lst, function(incs){ map(incs, function(inc){if(inc<0){0}else{inc}}) %>% unlist})
#   
#   if(stat == "gini"){
#     
#     return(dineq::gini.wtd(inc_lst %>% unlist, rep(wts, each = num_splits)))
#     
#   }else if(stat == "theil"){
#     
#     return(dineq::theil.wtd(inc_lst %>% unlist, rep(wts, each = num_splits)))
#     
#   }else if(stat == "sd"){
#     
#     wts_lng <- rep(wts, each = num_splits) 
#     
#     Hmisc::wtd.var(inc_lst %>% unlist, weights=wts_lng, method=c('unbiased', 'ML')) %>% sqrt %>% return
#     
#   }else if(stat == "atkinson"){
#     
#     wts_lng <- rep(wts, each = num_splits) 
#     
#     atkinson_wgt(inc_lst %>% unlist, wts_lng, epsilon = .2, wscale = 1000, G = G)
#     
#   }else if(stat == "inc_share"){
#     
#     get_inc_shares(inc_lst, wts, num_splits)
#     
#   }else if(stat == "quantile"){
#     
#     incs <- rep(1, length(inc_lst %>% unlist)) * rep(wts,each=100)
#     marker <- cumsum(incs)/sum(incs)
#     map(c(.2, .4, .6, .8, .95), ~(inc_lst %>% unlist)[which(marker > .x)[1]]) %>% unlist
#     
#   }
#   
# }
# 
# 
# test_lorenz(pums_data_ab, lorenz_v5_mod)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# as_tibble(list(c(1, 2, 3), c("a", "b", "c")), .name_repair = "unique")
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 