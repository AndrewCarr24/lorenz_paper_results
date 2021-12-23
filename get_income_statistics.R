
# Takes PUMS data, Lorenz int function, number of computer cores / Returns statistics (rel perc bias/RMSE) of income dist
get_income_statistics <- function(input_data, input_func, num_cores){
  
  # Helper - Creates income category column based on exact incomes
  group_incomes <- function(df){
    df %>% mutate(income_cat = ifelse(income < 10000, 1,
                                      ifelse(income >= 10000 & income < 15000, 2,
                                             ifelse(income >= 15000 & income < 20000, 3,
                                                    ifelse(income >= 20000 & income < 25000, 4,
                                                           ifelse(income >= 25000 & income < 30000, 5,
                                                                  ifelse(income >= 30000 & income < 35000, 6,
                                                                         ifelse(income >= 35000 & income < 40000, 7,
                                                                                ifelse(income >= 40000 & income < 45000, 8,
                                                                                       ifelse(income >= 45000 & income < 50000, 9,
                                                                                              ifelse(income >= 50000 & income < 60000, 10,
                                                                                                     ifelse(income >= 60000 & income < 75000, 11,
                                                                                                            ifelse(income >= 75000 & income < 100000, 12,
                                                                                                                   ifelse(income >= 100000 & income < 125000, 13,
                                                                                                                          ifelse(income >= 125000 & income < 150000, 14,
                                                                                                                                 ifelse(income >= 150000 & income < 200000, 15, 
                                                                                                                                        16)))))))))))))))) 
  }
  
  
  # Takes incomes and bounds / Returns bin means
  get_bin_means <- function(input_incs, input_bounds){
    
    closed_bin_means <- map(1:(length(input_bounds)-1), function(idx){
      mean(input_incs[input_incs >= input_bounds[idx] & input_incs < input_bounds[idx+1]])
    }
    ) %>% unlist
    
    bin_means <- c(closed_bin_means, mean(input_incs[input_incs >= input_bounds[length(input_bounds)]]))
    
    return(bin_means)
    
  }
  
  # Takes df columns / Returns relative percent bias
  rel_perc_bias <- function(col1, col2){
    
    e_j <- map2(col1, col2, function(c1, c2){
      
      100*((col1 - col2)/col2)
      
    }) %>% unlist 
    
    return( mean(e_j, na.rm = T) )
    
  }
  
  # Takes df columns / Returns relative percent RMSE
  rel_perc_rmse <- function(col1, col2){
    
    e_j <- map2(col1, col2, function(c1, c2){
      
      100*((col1 - col2)/col2)
      
    }) %>% unlist 
    
    return( sqrt(mean( (e_j)^2, na.rm = T )) )
    
  }
  
  future::plan(future::multisession, workers = num_cores)
  
  # Getting income statistics (Storing as results_tbl)
  results_tbl <- input_data %>% group_split(area) %>%  
    
    furrr::future_map(., function(sample_puma){
      
      # Get freqs, bounds, slope parm, distribution mean from sample_puma
      sample_puma <- sample_puma[rep(seq_len(nrow(sample_puma)), sample_puma$wt),] 
      ex_freqs <- sample_puma %>% group_incomes() %>% group_by(income_cat) %>% summarise(n = n(), .groups = 'drop_last') %>% .[["n"]]
      ex_bounds <- c(0, 10000, 15000, 20000, 25000, 30000, 35000, 40000, 45000, 50000, 60000, 75000, 100000, 125000, 150000, 200000)
      ex_slope_parm <- .9
      
      # When mean is known
      ex_G <- sum(sample_puma$income)/length(sample_puma$income)
      # When mean is unknown
      # bin_mids <- c((ex_bounds[1:15]+ex_bounds[2:16])/2, 300000)
      # ex_G <- sum(bin_mids*ex_freqs)/nrow(sample_puma)
      
      # True Income Statistics (based on exact incomes rather than income categories) #
      true_gini <- dineq::gini.wtd(sample_puma$income)
      true_theil <- dineq::theil.wtd(sample_puma$income)
      true_sd <- sd(sample_puma$income)
      true_at <- DescTools::Atkinson(sample_puma$income, parameter = .2)
      true_quantiles <- unname(quantile(sample_puma$income, c(.2, .4, .6, .8, .95)))
      true_inc_shares <- map(1:(length(true_quantiles)+1), function(idx){
        incs <- sample_puma$income
        if(idx == 1)
          sum(incs[incs < true_quantiles[idx]])/sum(incs) 
        else if(idx == (length(true_quantiles)+1))
          sum(incs[incs >= true_quantiles[idx-1]])/sum(incs)
        else
          sum(incs[incs >= true_quantiles[idx-1] & incs < true_quantiles[idx]])/sum(incs) 
      }) %>% unlist
      true_bin_means <- get_bin_means(sample_puma$income, ex_bounds)
      true_closed_bin_agg <- sum(sample_puma$income[sample_puma$income < 200000])
      
      # Lorenz interpolation estimates #
      # Applying arguments to lorenz interpolation function to get inequality stats, bin means, quantiles, income shares
      gini_args <-    list(freqs = ex_freqs, bounds = ex_bounds, G = ex_G, slope_parm = ex_slope_parm, stat = "gini")
      theil_args <-  list(freqs = ex_freqs, bounds = ex_bounds, G = ex_G, slope_parm = ex_slope_parm, stat = "theil")
      at_args <-  list(freqs = ex_freqs, bounds = ex_bounds, G = ex_G, slope_parm = ex_slope_parm, stat = "atkinson")
      sd_args <-        list(freqs = ex_freqs, bounds = ex_bounds, G = ex_G, slope_parm = ex_slope_parm, stat = "sd")
      bin_mean_args <- list(freqs = ex_freqs, bounds = ex_bounds, G = ex_G, slope_parm = ex_slope_parm, stat = "bin_means")
      quantile_args <- list(freqs = ex_freqs, bounds = ex_bounds, G = ex_G, slope_parm = ex_slope_parm, stat = "quantile")
      inc_share_args <- list(freqs = ex_freqs, bounds = ex_bounds, G = ex_G, slope_parm = ex_slope_parm, stat = "inc_share")
      lor_gini <- do.call(input_func, gini_args)
      lor_theil <- do.call(input_func, theil_args)
      lor_at <- do.call(input_func, at_args)
      lor_sd <- do.call(input_func, sd_args)
      lor_quantiles <-  do.call(input_func, quantile_args)
      lor_inc_shares <-  do.call(input_func, inc_share_args)
      lor_bin_means <- do.call(input_func, bin_mean_args)
      lor_closed_bin_agg <- sum(lor_bin_means[1:15]*ex_freqs[1:15])
      
      # CDF interpolation estimates # 
      # When mean is known
      bin_fit <- binsmooth::splinebins(c(ex_bounds[2:length(ex_bounds)], NA), ex_freqs, ex_G)
      # When mean is unknown 
      # bin_fit <- binsmooth::splinebins(c(ex_bounds[2:length(ex_bounds)], NA), ex_freqs)
      cdf_quantiles <- binsmooth::sb_percentiles(bin_fit, c(20, 40, 60, 80, 95))
      cdf_gini <- binsmooth::gini(bin_fit)
      cdf_theil <- binsmooth::theil(bin_fit)
      cdf_incs <- bin_fit$splineInvCDF(runif(2000))
      cdf_sd <- sd(cdf_incs)
      cdf_at <- DescTools::Atkinson(cdf_incs, parameter = .2)
      cdf_inc_tot <- sum(cdf_incs)
      cdf_is_20 <- sum(cdf_incs[cdf_incs < cdf_quantiles[1]])/cdf_inc_tot
      cdf_is_40 <- sum(cdf_incs[cdf_incs < cdf_quantiles[2] & cdf_incs >= cdf_quantiles[1]])/cdf_inc_tot
      cdf_is_60 <- sum(cdf_incs[cdf_incs < cdf_quantiles[3] & cdf_incs >= cdf_quantiles[2]])/cdf_inc_tot
      cdf_is_80 <- sum(cdf_incs[cdf_incs < cdf_quantiles[4] & cdf_incs >= cdf_quantiles[3]])/cdf_inc_tot
      cdf_is_100 <- sum(cdf_incs[cdf_incs >= cdf_quantiles[4]])/cdf_inc_tot
      cdf_inc_shares <- c(cdf_is_20, cdf_is_40, cdf_is_60, cdf_is_80, cdf_is_100)
      cdf_bin_means <- get_bin_means(cdf_incs, ex_bounds)
      cdf_closed_bin_agg <- sum(cdf_bin_means[1:15]*ex_freqs[1:15])
      
      # MCIB Atkinson (Atkinson not implemented in MCIB STATA module)
      mcib_at <- mcib_func(ex_freqs, ex_bounds, ex_G, eta = .2)
      mcib_bin_means <- mcib_func(ex_freqs, ex_bounds, ex_G)
      mcib_closed_bin_agg <- sum(mcib_bin_means[1:15]*ex_freqs[1:15])
      
      # Results # 
      income_results <- tibble(
        # Inequality stats
        true_gini = true_gini, true_theil = true_theil, true_at = true_at, true_sd = true_sd, 
        lor_gini = lor_gini, lor_theil = lor_theil, lor_at = lor_at, lor_sd = lor_sd, 
        cdf_gini = cdf_gini, cdf_theil = cdf_theil, cdf_at = cdf_at, cdf_sd = cdf_sd, 
        # MCIB Atkinson
        mcib_at = mcib_at,
        # Quantiles
        true_q20 = true_quantiles[1], true_q40 = true_quantiles[2], true_q60 = true_quantiles[3],
        true_q80 = true_quantiles[4], true_q95 = true_quantiles[5], 
        lor_q20 = lor_quantiles[1], lor_q40 = lor_quantiles[2], lor_q60 = lor_quantiles[3], 
        lor_q80 = lor_quantiles[4], lor_q95 = lor_quantiles[5], 
        cdf_q20 = cdf_quantiles[1], cdf_q40 = cdf_quantiles[2], cdf_q60 = cdf_quantiles[3],
        cdf_q80 = cdf_quantiles[4], cdf_q95 = cdf_quantiles[5],
        # Income shares
        true_is20 = true_inc_shares[1], true_is40 = true_inc_shares[2], true_is60 = true_inc_shares[3], 
        true_is80 = true_inc_shares[4], true_is100 = true_inc_shares[5] + true_inc_shares[6], 
        lor_is20 = lor_inc_shares[1], lor_is40 = lor_inc_shares[2], lor_is60 = lor_inc_shares[3], 
        lor_is80 = lor_inc_shares[4], lor_is100 = lor_inc_shares[5], 
        cdf_is20 = cdf_is_20, cdf_is40 = cdf_is_40, cdf_is60 = cdf_is_60,  
        cdf_is80 = cdf_is_80, cdf_is100 = cdf_is_100,
        # Bin Means
        true_b1 = true_bin_means[1], true_b2 = true_bin_means[2], true_b3 = true_bin_means[3], true_b4 = true_bin_means[4],
        true_b5 = true_bin_means[5], true_b6 = true_bin_means[6], true_b7 = true_bin_means[7], true_b8 = true_bin_means[8],
        true_b9 = true_bin_means[9], true_b10 = true_bin_means[10], true_b11 = true_bin_means[11], true_b12 = true_bin_means[12],
        true_b13 = true_bin_means[13], true_b14 = true_bin_means[14], true_b15 = true_bin_means[15], true_b16 = true_bin_means[16],
        true_closed_bin_agg = true_closed_bin_agg,
        lor_b1 =  lor_bin_means[1], lor_b2 = lor_bin_means[2], lor_b3 = lor_bin_means[3], lor_b4 = lor_bin_means[4],
        lor_b5 =  lor_bin_means[5], lor_b6 = lor_bin_means[6], lor_b7 = lor_bin_means[7], lor_b8 = lor_bin_means[8],
        lor_b9 =  lor_bin_means[9], lor_b10 = lor_bin_means[10], lor_b11 = lor_bin_means[11], lor_b12 = lor_bin_means[12],
        lor_b13 = lor_bin_means[13], lor_b14 = lor_bin_means[14], lor_b15 = lor_bin_means[15], lor_b16 = lor_bin_means[16],
        lor_closed_bin_agg = lor_closed_bin_agg,
        cdf_b1 = cdf_bin_means[1], cdf_b2 = cdf_bin_means[2], cdf_b3 = cdf_bin_means[3], cdf_b4 = cdf_bin_means[4],
        cdf_b5 = cdf_bin_means[5], cdf_b6 = cdf_bin_means[6], cdf_b7 = cdf_bin_means[7], cdf_b8 = cdf_bin_means[8],
        cdf_b9 = cdf_bin_means[9], cdf_b10 = cdf_bin_means[10], cdf_b11 = cdf_bin_means[11], cdf_b12 = cdf_bin_means[12],
        cdf_b13 = cdf_bin_means[13], cdf_b14 = cdf_bin_means[14], cdf_b15 = cdf_bin_means[15], cdf_b16 = cdf_bin_means[16],
        cdf_closed_bin_agg = cdf_closed_bin_agg,
        mcib_b1 = mcib_bin_means[1], mcib_b2 = mcib_bin_means[2], mcib_b3 = mcib_bin_means[3], mcib_b4 = mcib_bin_means[4],
        mcib_b5 = mcib_bin_means[5], mcib_b6 = mcib_bin_means[6], mcib_b7 = mcib_bin_means[7], mcib_b8 = mcib_bin_means[8],
        mcib_b9 = mcib_bin_means[9], mcib_b10 = mcib_bin_means[10], mcib_b11 = mcib_bin_means[11], mcib_b12 = mcib_bin_means[12],
        mcib_b13 = mcib_bin_means[13], mcib_b14 = mcib_bin_means[14], mcib_b15 = mcib_bin_means[15], mcib_b16 = mcib_bin_means[16],
        mcib_closed_bin_agg = mcib_closed_bin_agg
      )
      
      return(income_results)
      
    }) %>% do.call('rbind', .)
  
  ## MCIB Estimates - From output of MCIB STATA module 
  # When mean is known
  mcib_results <- readstata13::read.dta13("~/Downloads/puma_results_1_15.dta")
  # When mean is unknown
  # mcib_results <- readstata13::read.dta13("~/Downloads/puma_results_no_mean.dta")
  mcib_results_fin <- as_tibble(mcib_results) %>% select(gini, theil, sd, p20, p40, p60, p80, p95, shrq1, shrq2, shrq3, shrq4, shrq5) %>% 
    mutate(across(shrq1:shrq5, ~.x/100)) %>% 
    rename(mcib_gini = gini, mcib_theil = theil, mcib_sd = sd, mcib_q20 = p20, mcib_q40 = p40, mcib_q60 = p60, mcib_q80 = p80, 
           mcib_q95 = p95, mcib_is_20 = shrq1, mcib_is_40 = shrq2, mcib_is_60 = shrq3, mcib_is_80 = shrq4, mcib_is_100 = shrq5)
  
  # Adding MCIB results to results_tbl
  results_tbl <- bind_cols(results_tbl, mcib_results_fin)
  
  # Getting ineq_results - Relative Percent Bias/RMSE of all income statistics
  ineq_results <- results_tbl %>% 
    summarise(
      # LI Results
      lor_gini_bias = rel_perc_bias(lor_gini, true_gini), lor_gini_rmse = rel_perc_rmse(lor_gini, true_gini),
      lor_theil_bias = rel_perc_bias(lor_theil, true_theil), lor_theil_rmse = rel_perc_rmse(lor_theil, true_theil),
      lor_at_bias = rel_perc_bias(lor_at, true_at), lor_at_rmse = rel_perc_rmse(lor_at, true_at),
      lor_sd_bias = rel_perc_bias(lor_sd, true_sd), lor_sd_rmse = rel_perc_rmse(lor_sd, true_sd), 
      lor_q20_bias =   rel_perc_bias(lor_q20, true_q20), lor_q20_rmse = rel_perc_rmse(lor_q20, true_q20),
      lor_q40_bias =   rel_perc_bias(lor_q40, true_q40), lor_q40_rmse = rel_perc_rmse(lor_q40, true_q40),
      lor_q60_bias =   rel_perc_bias(lor_q60, true_q60), lor_q60_rmse = rel_perc_rmse(lor_q60, true_q60),
      lor_q80_bias =   rel_perc_bias(lor_q80, true_q80), lor_q80_rmse = rel_perc_rmse(lor_q80, true_q80),
      lor_q95_bias =   rel_perc_bias(lor_q95, true_q95), lor_q95_rmse = rel_perc_rmse(lor_q95, true_q95),
      lor_is20_bias =  rel_perc_bias(lor_is20, true_is20), lor_is20_rmse = rel_perc_rmse(lor_is20, true_is20),
      lor_is40_bias =  rel_perc_bias(lor_is40, true_is40), lor_is40_rmse = rel_perc_rmse(lor_is40, true_is40),
      lor_is60_bias =  rel_perc_bias(lor_is60, true_is60), lor_is60_rmse = rel_perc_rmse(lor_is60, true_is60),
      lor_is80_bias =  rel_perc_bias(lor_is80, true_is80), lor_is80_rmse = rel_perc_rmse(lor_is80, true_is80),
      lor_is100_bias = rel_perc_bias(lor_is100, true_is100), lor_is100_rmse = rel_perc_rmse(lor_is100, true_is100),
      lor_b1_bias = rel_perc_bias(lor_b1, true_b1), lor_b1_rmse = rel_perc_rmse(lor_b1, true_b1),
      lor_b2_bias = rel_perc_bias(lor_b2, true_b2), lor_b2_rmse = rel_perc_rmse(lor_b2, true_b2),
      lor_b3_bias = rel_perc_bias(lor_b3, true_b3), lor_b3_rmse = rel_perc_rmse(lor_b3, true_b3),
      lor_b4_bias = rel_perc_bias(lor_b4, true_b4), lor_b4_rmse = rel_perc_rmse(lor_b4, true_b4),
      lor_b5_bias = rel_perc_bias(lor_b5, true_b5), lor_b5_rmse = rel_perc_rmse(lor_b5, true_b5),
      lor_b6_bias = rel_perc_bias(lor_b6, true_b6), lor_b6_rmse = rel_perc_rmse(lor_b6, true_b6),
      lor_b7_bias = rel_perc_bias(lor_b7, true_b7), lor_b7_rmse = rel_perc_rmse(lor_b7, true_b7),
      lor_b8_bias = rel_perc_bias(lor_b8, true_b8), lor_b8_rmse = rel_perc_rmse(lor_b8, true_b8),
      lor_b9_bias = rel_perc_bias(lor_b9, true_b9), lor_b9_rmse = rel_perc_rmse(lor_b9, true_b9),
      lor_b10_bias = rel_perc_bias(lor_b10, true_b10), lor_b10_rmse = rel_perc_rmse(lor_b10, true_b10),
      lor_b11_bias = rel_perc_bias(lor_b11, true_b11), lor_b11_rmse = rel_perc_rmse(lor_b11, true_b11),
      lor_b12_bias = rel_perc_bias(lor_b12, true_b12), lor_b12_rmse = rel_perc_rmse(lor_b12, true_b12),
      lor_b13_bias = rel_perc_bias(lor_b13, true_b13), lor_b13_rmse = rel_perc_rmse(lor_b13, true_b13),
      lor_b14_bias = rel_perc_bias(lor_b14, true_b14), lor_b14_rmse = rel_perc_rmse(lor_b14, true_b14),
      lor_b15_bias = rel_perc_bias(lor_b15, true_b15), lor_b15_rmse = rel_perc_rmse(lor_b15, true_b15),
      lor_b16_bias = rel_perc_bias(lor_b16, true_b16), lor_b16_rmse = rel_perc_rmse(lor_b16, true_b16),
      lor_closed_bin_agg_bias = rel_perc_bias(lor_closed_bin_agg, true_closed_bin_agg),
      lor_closed_bin_agg_rmse = rel_perc_rmse(lor_closed_bin_agg, true_closed_bin_agg),
      
      # CDF-I Results
      cdf_gini_bias = rel_perc_bias(cdf_gini, true_gini), cdf_gini_rmse = rel_perc_rmse(cdf_gini, true_gini),
      cdf_theil_bias = rel_perc_bias(cdf_theil, true_theil), cdf_theil_rmse = rel_perc_rmse(cdf_theil, true_theil),
      cdf_at_bias = rel_perc_bias(cdf_at, true_at), cdf_at_rmse = rel_perc_rmse(cdf_at, true_at),
      cdf_sd_bias = rel_perc_bias(cdf_sd, true_sd), cdf_sd_rmse = rel_perc_rmse(cdf_sd, true_sd), 
      cdf_q20_bias = rel_perc_bias(cdf_q20, true_q20), cdf_q20_rmse = rel_perc_rmse(cdf_q20, true_q20),
      cdf_q40_bias = rel_perc_bias(cdf_q40, true_q40), cdf_q40_rmse = rel_perc_rmse(cdf_q40, true_q40),
      cdf_q60_bias = rel_perc_bias(cdf_q60, true_q60), cdf_q60_rmse = rel_perc_rmse(cdf_q60, true_q60),
      cdf_q80_bias = rel_perc_bias(cdf_q80, true_q80), cdf_q80_rmse = rel_perc_rmse(cdf_q80, true_q80),
      cdf_q95_bias = rel_perc_bias(cdf_q95, true_q95), cdf_q95_rmse = rel_perc_rmse(cdf_q95, true_q95),
      cdf_is20_bias = rel_perc_bias(cdf_is20, true_is20), cdf_is20_rmse = rel_perc_rmse(cdf_is20, true_is20),
      cdf_is40_bias = rel_perc_bias(cdf_is40, true_is40), cdf_is40_rmse = rel_perc_rmse(cdf_is40, true_is40),
      cdf_is60_bias = rel_perc_bias(cdf_is60, true_is60), cdf_is60_rmse = rel_perc_rmse(cdf_is60, true_is60),
      cdf_is80_bias = rel_perc_bias(cdf_is80, true_is80), cdf_is80_rmse = rel_perc_rmse(cdf_is80, true_is80),
      cdf_is100_bias = rel_perc_bias(cdf_is100, true_is100), cdf_is100_rmse = rel_perc_rmse(cdf_is100, true_is100),
      cdf_b1_bias = rel_perc_bias(cdf_b1, true_b1), cdf_b1_rmse = rel_perc_rmse(cdf_b1, true_b1),
      cdf_b2_bias = rel_perc_bias(cdf_b2, true_b2), cdf_b2_rmse = rel_perc_rmse(cdf_b2, true_b2),
      cdf_b3_bias = rel_perc_bias(cdf_b3, true_b3), cdf_b3_rmse = rel_perc_rmse(cdf_b3, true_b3),
      cdf_b4_bias = rel_perc_bias(cdf_b4, true_b4), cdf_b4_rmse = rel_perc_rmse(cdf_b4, true_b4),
      cdf_b5_bias = rel_perc_bias(cdf_b5, true_b5), cdf_b5_rmse = rel_perc_rmse(cdf_b5, true_b5),
      cdf_b6_bias = rel_perc_bias(cdf_b6, true_b6), cdf_b6_rmse = rel_perc_rmse(cdf_b6, true_b6),
      cdf_b7_bias = rel_perc_bias(cdf_b7, true_b7), cdf_b7_rmse = rel_perc_rmse(cdf_b7, true_b7),
      cdf_b8_bias = rel_perc_bias(cdf_b8, true_b8), cdf_b8_rmse = rel_perc_rmse(cdf_b8, true_b8),
      cdf_b9_bias = rel_perc_bias(cdf_b9, true_b9), cdf_b9_rmse = rel_perc_rmse(cdf_b9, true_b9),
      cdf_b10_bias = rel_perc_bias(cdf_b10, true_b10), cdf_b10_rmse = rel_perc_rmse(cdf_b10, true_b10),
      cdf_b11_bias = rel_perc_bias(cdf_b11, true_b11), cdf_b11_rmse = rel_perc_rmse(cdf_b11, true_b11),
      cdf_b12_bias = rel_perc_bias(cdf_b12, true_b12), cdf_b12_rmse = rel_perc_rmse(cdf_b12, true_b12),
      cdf_b13_bias = rel_perc_bias(cdf_b13, true_b13), cdf_b13_rmse = rel_perc_rmse(cdf_b13, true_b13),
      cdf_b14_bias = rel_perc_bias(cdf_b14, true_b14), cdf_b14_rmse = rel_perc_rmse(cdf_b14, true_b14),
      cdf_b15_bias = rel_perc_bias(cdf_b15, true_b15), cdf_b15_rmse = rel_perc_rmse(cdf_b15, true_b15),
      cdf_b16_bias = rel_perc_bias(cdf_b16, true_b16), cdf_b16_rmse = rel_perc_rmse(cdf_b16, true_b16),
      cdf_closed_bin_agg_bias = rel_perc_bias(cdf_closed_bin_agg, true_closed_bin_agg),
      cdf_closed_bin_agg_rmse = rel_perc_rmse(cdf_closed_bin_agg, true_closed_bin_agg),
      
      # MCIB Results
      mcib_gini_bias = rel_perc_bias(mcib_gini, true_gini), mcib_gini_rmse = rel_perc_rmse(mcib_gini, true_gini),
      mcib_theil_bias = rel_perc_bias(mcib_theil, true_theil), mcib_theil_rmse = rel_perc_rmse(mcib_theil, true_theil),
      mcib_at_bias = rel_perc_bias(mcib_at, true_at), mcib_at_rmse = rel_perc_rmse(mcib_at, true_at),
      mcib_sd_bias = rel_perc_bias(mcib_sd, true_sd), mcib_sd_rmse = rel_perc_rmse(mcib_sd, true_sd), 
      mcib_q20_bias = rel_perc_bias(mcib_q20, true_q20), mcib_q20_rmse = rel_perc_rmse(mcib_q20, true_q20),
      mcib_q40_bias = rel_perc_bias(mcib_q40, true_q40), mcib_q40_rmse = rel_perc_rmse(mcib_q40, true_q40),
      mcib_q60_bias = rel_perc_bias(mcib_q60, true_q60), mcib_q60_rmse = rel_perc_rmse(mcib_q60, true_q60),
      mcib_q80_bias = rel_perc_bias(mcib_q80, true_q80), mcib_q80_rmse = rel_perc_rmse(mcib_q80, true_q80),
      mcib_q95_bias = rel_perc_bias(mcib_q95, true_q95), mcib_q95_rmse = rel_perc_rmse(mcib_q95, true_q95),
      mcib_is20_bias = rel_perc_bias(mcib_is_20, true_is20), mcib_is20_rmse = rel_perc_rmse(mcib_is_20, true_is20),
      mcib_is40_bias = rel_perc_bias(mcib_is_40, true_is40), mcib_is40_rmse = rel_perc_rmse(mcib_is_40, true_is40),
      mcib_is60_bias = rel_perc_bias(mcib_is_60, true_is60), mcib_is60_rmse = rel_perc_rmse(mcib_is_60, true_is60),
      mcib_is80_bias = rel_perc_bias(mcib_is_80, true_is80), mcib_is80_rmse = rel_perc_rmse(mcib_is_80, true_is80),
      mcib_is100_bias = rel_perc_bias(mcib_is_100, true_is100), mcib_is100_rmse = rel_perc_rmse(mcib_is_100, true_is100),
      mcib_b1_bias = rel_perc_bias(mcib_b1, true_b1), mcib_b1_rmse = rel_perc_rmse(mcib_b1, true_b1),
      mcib_b2_bias = rel_perc_bias(mcib_b2, true_b2), mcib_b2_rmse = rel_perc_rmse(mcib_b2, true_b2),
      mcib_b3_bias = rel_perc_bias(mcib_b3, true_b3), mcib_b3_rmse = rel_perc_rmse(mcib_b3, true_b3),
      mcib_b4_bias = rel_perc_bias(mcib_b4, true_b4), mcib_b4_rmse = rel_perc_rmse(mcib_b4, true_b4),
      mcib_b5_bias = rel_perc_bias(mcib_b5, true_b5), mcib_b5_rmse = rel_perc_rmse(mcib_b5, true_b5),
      mcib_b6_bias = rel_perc_bias(mcib_b6, true_b6), mcib_b6_rmse = rel_perc_rmse(mcib_b6, true_b6),
      mcib_b7_bias = rel_perc_bias(mcib_b7, true_b7), mcib_b7_rmse = rel_perc_rmse(mcib_b7, true_b7),
      mcib_b8_bias = rel_perc_bias(mcib_b8, true_b8), mcib_b8_rmse = rel_perc_rmse(mcib_b8, true_b8),
      mcib_b9_bias = rel_perc_bias(mcib_b9, true_b9), mcib_b9_rmse = rel_perc_rmse(mcib_b9, true_b9),
      mcib_b10_bias = rel_perc_bias(mcib_b10, true_b10), mcib_b10_rmse = rel_perc_rmse(mcib_b10, true_b10),
      mcib_b11_bias = rel_perc_bias(mcib_b11, true_b11), mcib_b11_rmse = rel_perc_rmse(mcib_b11, true_b11),
      mcib_b12_bias = rel_perc_bias(mcib_b12, true_b12), mcib_b12_rmse = rel_perc_rmse(mcib_b12, true_b12),
      mcib_b13_bias = rel_perc_bias(mcib_b13, true_b13), mcib_b13_rmse = rel_perc_rmse(mcib_b13, true_b13),
      mcib_b14_bias = rel_perc_bias(mcib_b14, true_b14), mcib_b14_rmse = rel_perc_rmse(mcib_b14, true_b14),
      mcib_b15_bias = rel_perc_bias(mcib_b15, true_b15), mcib_b15_rmse = rel_perc_rmse(mcib_b15, true_b15),
      mcib_b16_bias = rel_perc_bias(mcib_b16, true_b16), mcib_b16_rmse = rel_perc_rmse(mcib_b16, true_b16),
      mcib_closed_bin_agg_bias = rel_perc_bias(mcib_closed_bin_agg, true_closed_bin_agg),
      mcib_closed_bin_agg_rmse = rel_perc_rmse(mcib_closed_bin_agg, true_closed_bin_agg)
    )  
  
  return(ineq_results)
  
}

