
# Takes coordinates of an LC / Returns coordinates that define quadratic-cubic function to approximate LC
lorenz_to_coefs <- function(lorenz_df, slope_parm, bounds, G){
  
  ## Helper - Adjusts cubic applied to top bin based to ensure cubic does not violate conditions of a Lorenz curve
  cubic_adjust <- function(bin_idx, x_1, x_2, y_1, bounds, G){
    
    A <- matrix(c(x_1^2, 2*x_1, 2*x_2,
                  x_1,   1,     1,
                  1,     0,     0), 3, 3)
    y <- matrix(c(y_1, bounds[bin_idx]/G, bounds[bin_idx+1]/G), 3, 1)
    
    coefs <- solve(A, y)
    a <- coefs[1]
    b <- coefs[2]
    c <- coefs[3]
    bin_mid <- (x_1+x_2)/2
    mid_slope <- 2*a*bin_mid + b
    
    # Cubic  
    A = matrix(c(x_1^3, 3*x_1^2, 3*x_2^2, 3*bin_mid^2,
                 x_1^2, 2*x_1,   2*x_2,   2*bin_mid,
                 x_1,   1,       1,       1,
                 1,     0,       0,       0), 4, 4)
    
    # Incrementally reducing slope at midpoint of top bin
    y <- matrix(c(y_1, bounds[bin_idx]/G, bounds[bin_idx+1]/G, mid_slope*.9999), 4, 1)
    
    coef_test <- solve(A, y)
    
    while(6*coef_test[1]*x_1 + 2*coef_test[2] > 0 & 6*coef_test[1]*x_2 + 2*coef_test[2] > 0){
      
      coef_orig <- coef_test
      mid_slope <- mid_slope*.9999
      y <- matrix(c(y_1, bounds[bin_idx]/G, bounds[bin_idx+1]/G, mid_slope), 4, 1)
      coef_test <- solve(A, y)
      
    }
    
    return(coef_orig)
    
  }
  
  
  coefs_lst <- list()
  for(i in 1:15){
    
    bin_idx <- i
    
    x_1 <- lorenz_df$x[bin_idx]
    x_2 <- lorenz_df$x[bin_idx+1]
    x_3 <- lorenz_df$x[bin_idx+2]
    
    A = matrix(c(x_1^3, 3*x_1^2, 3*x_2^2, 3*x_3^2,
                 x_1^2, 2*x_1,   2*x_2,   2*x_3,
                 x_1,   1,       1,       1,
                 1,     0,       0,       0), 4, 4)
    
    if(bin_idx == 1){first_y <- 0}else{first_y <- a*x_1^3 + b*x_1^2 + c*x_1 + d}
    if(bin_idx == 15){
      last_y <- 1  
      A[4,] <- c(1, 1, 1, 1)
    }else{
      last_y <- bounds[bin_idx+2]/G
    }
    
    y = matrix(c(first_y, bounds[bin_idx]/G, bounds[bin_idx+1]/G, last_y), 4, 1)
    
    coefs <- solve(A, y)
    
    coefs_lst[[i]] <- coefs
    
    a <- coefs[1] 
    b <- coefs[2]
    c <- coefs[3]
    d <- coefs[4]
    
    # Checking valid lorenz
    while(6*a*x_1 + 2*b < 0 | 6*a*x_2 + 2*b < 0){
      
      # print(paste("up ping", i))
      x_3 <- x_3*1.0001
      
      A = matrix(c(x_1^3, 3*x_1^2, 3*x_2^2, 3*x_3^2,
                   x_1^2, 2*x_1,   2*x_2,   2*x_3,
                   x_1,   1,       1,       1,
                   1,     0,       0,       0), 4, 4)
      
      if(bin_idx == 1){first_y <- 0}else{first_y <- coefs_lst[[bin_idx-1]][1]*x_1^3 + coefs_lst[[bin_idx-1]][2]*x_1^2 + 
        coefs_lst[[bin_idx-1]][3]*x_1 + coefs_lst[[bin_idx-1]][4]}
      
      if(bin_idx == 15){
        last_y <- last_y*.9999
        A[4,] <- c(1, 1, 1, 1)
      }else{
        last_y <- bounds[bin_idx+2]/G
      }
      
      y = matrix(c(first_y, bounds[bin_idx]/G, bounds[bin_idx+1]/G, last_y), 4, 1)
      
      coefs <- solve(A, y)
      
      coefs_lst[[i]] <- coefs
      
      a <- coefs[1] 
      b <- coefs[2]
      c <- coefs[3]
      d <- coefs[4]
      
      if(bin_idx == 15){
        
        coefs_lst[[i]] <- coefs <- cubic_adjust(bin_idx, x_1, x_2, first_y, bounds, G)
        a <- coefs[1] 
        b <- coefs[2]
        c <- coefs[3]
        d <- coefs[4]
        
      }
      
    }
    
  }
  
  
  ##########
  
  # Fitting quadratic function with start slope to the top category
  x_1 = lorenz_df$x[16]
  x_2 = 1
  start_slope = 200000/G 
  
  A = matrix(c(x_1^2, x_2^2, 2*x_1,
               x_1, x_2, 1,
               1, 1, 0), 3, 3)
  
  if(length(coefs_lst[[15]]) == 4){
    y_1 <- coefs_lst[[15]][1]*x_1^3 + coefs_lst[[15]][2]*x_1^2 + coefs_lst[[15]][3]*x_1 + coefs_lst[[15]][4]
  }else if(length(coefs_lst[[15]]) == 3){
    y_1 <- coefs_lst[[15]][1]*x_1^2 + coefs_lst[[15]][2]*x_1 + coefs_lst[[15]][3]
  }
  y_2 <- 1
  y = matrix(c(y_1, y_2, start_slope), 3, 1)
  quad_coefs <- solve(A, y)
  
  mid_x <- (lorenz_df$x[(length(lorenz_df$x)-1)] + lorenz_df$x[(length(lorenz_df$x))])/2
  slope_mid_x <- 2*quad_coefs[1]*mid_x + quad_coefs[2]
  
  # Cubic to two points and two slopes
  A = matrix(c(x_1^3, x_2^3, 3*x_1^2, 3*mid_x^2,
               x_1^2, x_2^2, 2*x_1, 2*mid_x,
               x_1, x_2, 1, 1,
               1, 1, 0, 0), 4, 4)
  
  y = matrix(c(y_1,
               y_2,
               start_slope,
               slope_parm*slope_mid_x), 4, 1)
  
  final_coefs <- solve(A, y)
  a <- final_coefs[1]; b <- final_coefs[2]
  while((6*a*x_1 + 2*b < 0 | 6*a*x_2 + 2*b < 0) & slope_parm < 1){
    
    slope_parm <- slope_parm*1.001
    y = matrix(c(y_1,
                 y_2,
                 start_slope,
                 slope_parm*slope_mid_x), 4, 1)
    final_coefs <- solve(A, y)
    a <- final_coefs[1]; b <- final_coefs[2]
  }
  
  print(slope_parm)
  
  coefs_lst[[length(coefs_lst) + 1]] <- final_coefs
  
  return(coefs_lst)
  
}






