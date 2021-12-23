
# Implementation of MCIB to get statistics not available in STATA module - 

# Takes income distribution arguments (freqs,bounds,mean) and optional Atkinson eta / Returns either exact incomes or 
# Atkinson estimate depending on whether eta is specified
mcib_func <- function(freqs, bounds, mean, eta = NA){
  
  # Takes incomes and bounds / Returns bin means
  get_bin_means <- function(input_incs, input_bounds){
    
    closed_bin_means <- map(1:(length(input_bounds)-1), function(idx){
      
      mean(input_incs[input_incs >= input_bounds[idx] & input_incs < input_bounds[idx+1]])
      
    }
    ) %>% unlist
    
    bin_means <- c(closed_bin_means, mean(input_incs[input_incs >= input_bounds[length(input_bounds)]]))
    
    return(bin_means)
    
  }
  
  # Takes numeric freq vec (just closed bracket freqs) and lower bound vec and returns mcib coords as tibble
  freq_and_bounds_to_mcib_coords <- function(freqs, bounds){
    return(as.data.frame(cbind(x = (bounds[1:(length(bounds)-1)]+bounds[2:length(bounds)])/2, y = freqs/diff(bounds))))
  }
  
  
  # Takes mcib coords tibble (midpoint and rel_freq) and returns slopes and ints in list
  mcib_coords_to_slopes_ints <- function(points){
    
    points_lst <- split(points %>% as.matrix() %>% unname, seq(nrow(points))) %>% unname
    
    slopes <- lapply(seq_along(points_lst), function(idx){
      
      if(idx != length(points_lst)){
        (points_lst[[idx+1]][2]-points_lst[[idx]][2])/(points_lst[[idx+1]][1]-points_lst[[idx]][1])
      }else{
        (points_lst[[idx-1]][2]-points_lst[[idx]][2])/(points_lst[[idx-1]][1]-points_lst[[idx]][1])
      }
    })
    
    slopes <- sapply(seq_along(slopes), function(idx){
      if(idx != 1 & idx != length(slopes)){
        return((slopes[[idx-1]]+slopes[[idx]])/2)
      }else{
        return(slopes[[idx]])
      }
    }) %>% unlist
    
    slopes <- lapply(seq_along(slopes), function(idx){
      
      if(idx == 1){return(0)}else if(idx == length(slopes)){
        return(slopes[idx])
      }else{
        if(points_lst[[idx-1]][2] > points_lst[[idx]][2] & points_lst[[idx+1]][2] > points_lst[[idx]][2]){
          return(0)
        }else if(points_lst[[idx-1]][2] < points_lst[[idx]][2] & points_lst[[idx+1]][2] < points_lst[[idx]][2]){
          return(0)
        }else{
          return(slopes[idx])
        }
      }
      
    }) %>% unlist
    
    ints <- Map(function(slope, coords){
      -slope*coords[1] + coords[2]
    }, slopes, points_lst) %>% unlist
    
    return(list(slopes, ints))
    
  }
  
  # Takes lines that cross the x-axis and rotates them at the rel_freq point so that they remain positive over the bracket (new x-intercept is bracket boundary)
  x_adj <- function(slopes_ints, bounds, mcib_coords){
    
    bad_idx <- which(Map(function(m, b, idx){
      
      LB <- bounds[idx]
      UB <- bounds[idx+1]
      x_int <- -b/m
      x_int >= LB & x_int <= UB
      
    }, slopes_ints[[1]], slopes_ints[[2]], seq_along(slopes_ints[[1]])) %>% unlist)
    
    map_output <- Map(function(m, b, idx){
      
      if(idx %in% bad_idx){
        
        if(m < 0){
          
          m_new <- (as.numeric(mcib_coords[idx, 2]) - 0)/(as.numeric(mcib_coords[idx, 1]) - bounds[idx+1])
          b_new <- (-1)*m_new*bounds[idx+1]
          return(c(m_new, b_new))
          
        }else if(m > 0){
          
          m_new <- (as.numeric(mcib_coords[idx, 2]) - 0)/(as.numeric(mcib_coords[idx, 1]) - bounds[idx])
          b_new <- (-1)*m_new*bounds[idx+1]
          return(c(m_new, b_new))
          
        }
      }else{
        
        return(c(m, b))
        
      }
    }, slopes_ints[[1]], slopes_ints[[2]], seq_along(slopes_ints[[1]]))
    
    lapply(do.call('rbind', map_output) %>% as.data.frame(), function(x){x}) %>% unname
  }
  
  # Takes slopes and ints and lower bounds and returns MCIB means of closed intervals
  MCIB_closed_bracket_means <- function(slopes_ints, lower_bounds10, freqs){
    
    mean_integral <- function(m, b, x){
      return((m/3)*x^3 + (b/2)*x^2)}
    
    sapply(1:(length(lower_bounds10)-1), function(bracket_idx){
      
      LB <- (lower_bounds10)[bracket_idx]
      UB <- (lower_bounds10)[bracket_idx+1]
      
      m <- slopes_ints[[1]][bracket_idx]
      b <- slopes_ints[[2]][bracket_idx]
      
      mean_b <- mean_integral(m, b, UB) - mean_integral(m, b, LB)
      
      return(mean_b/freqs[bracket_idx])
    })
  }
  
  # Takes rescaled slopes_ints, rescaled bounds (pdf scale), and F_x and returns polynomial coeficients of CDF
  # of each bracket
  slopes_ints_rescaled_to_poly_coefs <- function(slopes_ints_rescaled, bounds, F_x){
    
    Map(function(a, b, x_0, y_0){
      
      c <- y_0 - (a*x_0^2)/2 - b*x_0
      
      return(c(a,b,c))
    }, slopes_ints_rescaled[[1]], slopes_ints_rescaled[[2]], bounds, F_x)
  }
  
  # Takes top bracket mean and returns rescaled parameters of the pareto distribution
  top_bracket_mean_to_rescaled_pareto_parms <- function(top_bracket_mean, bounds, N){
    
    beta <- (bounds/sqrt(N))[length(bounds)]
    alpha <- top_bracket_mean/(top_bracket_mean-bounds[length(bounds)])
   
    return(c(alpha, beta))
    
  }
  
  # Takes a number between 0 and 1 and returns draw from the MCIB-estimated pdf
  inverse_cdf_full <- function(inputs, bounds, poly_coefs, pareto_parms, F_x){
    
    sapply(inputs, function(y){
      
      idx <- which(sapply(seq_along(bounds[1:(length(bounds)-1)]), function(z){
        
        y > bounds[z] & y < bounds[z+1]
        
      }))
      
      if(idx < (length(bounds)-1)){
        
        a_b_c <- poly_coefs[[idx]]
        
        a <- a_b_c[1]
        b <- a_b_c[2]
        c <- a_b_c[3]
        
        if(a == 0){
          
          return((y - c)/b)
          
        }else{
          
          return(inverse_quadratic(a,b,c,y)[1])
        }
        
      }else{
        
        alpha <- pareto_parms[1]
        beta <- pareto_parms[2]
        
        return(inverse_pareto_cdf(alpha, beta, 1-F_x[length(F_x)], y))
        
      }
      
    })
    
  }
  
  # Helper - inverse of the quadratic formula
  inverse_quadratic <- function(a,b,c, x_0){
    a <- a/2
    c <- c - x_0
    return(c((-b+sqrt(b^2-4*a*c))/(2*a), (-b-sqrt(b^2-4*a*c))/(2*a)))
  }
  
  # Inverse cdf of pareto used for transform sampling top bracket
  inverse_pareto_cdf <- function(alpha, beta, top_prop, input){
    
    # Scaling input
    input <- (input-(1-top_prop))/(1-(1-top_prop))
    
    # Capping input at .995 (following MCIB)
    input <- min(input, .995)
    
    return(beta/(1-input)^(1/alpha))
  }
  
  # Auxiliary Functions for lorenz_interp #
  freqs_to_mcib_means <- function(freqs, bounds, mean){
    
    N <- sum(freqs)
    agg <- N*mean
    F_x <- (cumsum(freqs)/sum(freqs))[1:(length(freqs)-1)]
    
    mcib_coords <- freq_and_bounds_to_mcib_coords(freqs[1:(length(freqs)-1)], bounds)
    slopes_ints <- mcib_coords_to_slopes_ints(mcib_coords)
    
    # Fixing lines that cross x-axis
    slopes_ints <-  x_adj(slopes_ints, bounds, mcib_coords)
    
    # Storing closed bracket means
    closed_bracket_means <- MCIB_closed_bracket_means(slopes_ints, bounds, freqs)
    
    # Storing top bracket mean estimate
    top_bracket_mean <- (agg - sum(freqs[1:(length(freqs)-1)]*closed_bracket_means, na.rm = T))/
      freqs[length(freqs)]
    
    return(c(closed_bracket_means, top_bracket_mean))
    
  }

  #################################
  
  N <- sum(freqs)
  agg <- mean*N
  rescaled_bounds <- bounds[2:length(bounds)]/sqrt(N)
  F_x <- (cumsum(freqs)/sum(freqs))[1:(length(freqs)-1)]
  
  mcib_coords <- freq_and_bounds_to_mcib_coords(freqs[1:(length(freqs)-1)], bounds)
  slopes_ints <- mcib_coords_to_slopes_ints(mcib_coords)
  slopes_ints_rescaled <- mcib_coords_to_slopes_ints(mcib_coords/sqrt(N))
  slopes_ints_rescaled <- x_adj(slopes_ints_rescaled, bounds/sqrt(N), mcib_coords/sqrt(N))
  
  # Polynomial coefficients of CDF brackets
  poly_coefs <- slopes_ints_rescaled_to_poly_coefs(slopes_ints_rescaled, rescaled_bounds, F_x)
  
  # Storing closed bracket means
  closed_bracket_means <- MCIB_closed_bracket_means(slopes_ints, bounds, freqs)
  
  # Storing top bracket mean estimate
  top_bracket_mean <- (agg - sum(freqs[1:(length(freqs)-1)]*closed_bracket_means))/
    freqs[length(freqs)]
  
  # Storing rescaled parms of pareto distribution for top bracket
  pareto_parms <- top_bracket_mean_to_rescaled_pareto_parms(top_bracket_mean, bounds, N)
  
  # Using inverse-CDF method to take 3000 samples from MCIB-based pdf
  samples <- inverse_cdf_full(runif(3000), c(0, F_x, 1), poly_coefs, pareto_parms, F_x)*sqrt(N)
  
  if(is.na(eta)){
    return(get_bin_means(samples, bounds))
  }else{
    return(DescTools::Atkinson(samples, parameter = eta))
  }
  
}






