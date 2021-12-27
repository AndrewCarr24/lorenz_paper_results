library(tidyverse)
library(readr)
library(ggrepel)
library(gridExtra)


# Reading get_income_statistics function
source("get_income_statistics.R")

# Reading MCIB function
source("mcib_func.R")

# Reading lorenz interpolation function 
source("lorenz_int.R")

bounds <- c(0, 10000, 15000, 20000, 25000, 30000, 35000, 40000, 45000, 50000, 60000, 75000, 100000, 125000, 150000, 200000)

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

# Takes x coordinates of an LC / Returns coordinates that define quadratic-cubic function to approximate LC
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
  
  coefs_lst[[length(coefs_lst) + 1]] <- final_coefs
  
  return(coefs_lst)
  
}


# Reading PUMS data from 2011-2015 ACS and dropping incomes outside PUMAs
pums_data <- read_csv("pums_data.csv") %>% 
  rename(area = PUMID) %>% select(area, income, wt) %>% filter(area != 0)




####################################
############ Plot 1 ################
####################################

# Helper for plot 1
get_segment <- function(point, slope, len){
  
  hyp <- sqrt(1^2 + slope^2)
  amt <- len/hyp
  
  upper_x <- point[1] + amt*1
  upper_y <- point[2] + amt*slope
  
  lower_x <- point[1] - amt*1
  lower_y <- point[2] - amt*slope
  
  if(lower_x < 0){lower_x <- 0; lower_y <- 0}
  
  return( tibble(x1 = lower_x, y1 = lower_y, x2 = upper_x, y2 = upper_y) )
  
}

# Try MSA 23
sample_puma <- pums_data %>% group_split(area) %>% .[[23]]
sample_puma <- sample_puma[rep(seq_len(nrow(sample_puma)), sample_puma$wt),] 
freqs <- sample_puma %>% group_incomes() %>% group_by(income_cat) %>% summarise(n = n(), .groups = 'drop_last') %>% .[["n"]]
G <- sum(sample_puma$income)/length(sample_puma$income)
lorenz_df <- tibble(x = c(0, cumsum(freqs)/sum(freqs))) %>% distinct()

# Getting slopes and coordinate data for plot 
slopes <- bounds[1:4]/G
x_1 <- lorenz_df$x[1] 
y_1 <- 0
x_2 <- lorenz_df$x[2]
x_3 <- lorenz_df$x[3]
x_4 <- lorenz_df$x[4]

# Interpolating bottom bin
coefs <- lorenz_to_coefs(lorenz_df, slope_parm = .9, bounds, G)[[1]]
a <- coefs[1]; b <- coefs[2]; c <- coefs[3]; d <- coefs[4]
y_2 <- a*x_2^3 + b*x_2^2 + c*x_2 + d
y_3 <- a*x_3^3 + b*x_3^2 + c*x_3 + d

# Fourth slope
coefs <- lorenz_to_coefs(lorenz_df, slope_parm = .9, bounds, G)[[2]]
a <- coefs[1]; b <- coefs[2]; c <- coefs[3]; d <- coefs[4]
y_4 <- a*x_4^3 + b*x_4^2 + c*x_4 + d

seg1 <- get_segment(c(x_1, y_1), slopes[1], len = .04)
seg2 <- get_segment(c(x_2, y_2), slopes[2], len = .03)
seg3 <- get_segment(c(x_3, y_3), slopes[3], len = .03)
seg4 <- get_segment(c(x_4, y_4), slopes[4], len = .03)

p1a <- ggplot(lorenz_df) + geom_vline(aes(xintercept = x), lty = 'dashed', size = .25) + 
  geom_hline(mapping = aes(yintercept = 0)) +
  geom_vline(mapping = aes(xintercept = 0)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_point(data = tibble(x = 0, y = 0), mapping = aes(x, y), size = .75) + 
  scale_x_continuous(limits = c(-.01,1), expand = c(0, .01)) +
  scale_y_continuous(limits = c(0,1), expand = c(0, .01)) + 
  ylab("") + xlab("") + 
  geom_segment(data = seg1, mapping = aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment"), lty = 'solid', size = .5) + 
  geom_segment(data = seg2, mapping = aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment"), lty = 'solid', size = .5) + 
  geom_segment(data = seg3, mapping = aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment"), lty = 'solid', size = .5) +
  scale_colour_discrete(guide=FALSE)

##

coefs_lst <- lorenz_to_coefs(lorenz_df, .9, bounds, G)
xs <- seq(lorenz_df$x[1], lorenz_df$x[2]-.00001, (lorenz_df$x[2]-lorenz_df$x[1])/1000)
ys <- map(xs, function(x){
  idx <- which(x < lorenz_df$x)[1]  -1
  a <- coefs_lst[[idx]][1]; b <- coefs_lst[[idx]][2]; c <- coefs_lst[[idx]][3]; d <- coefs_lst[[idx]][4]
  a*x^3 + b*x^2 + c*x + d
}) %>% unlist
first_seg_line <- tibble(x = xs, y = ys)

p1b <- ggplot(lorenz_df) + geom_vline(aes(xintercept = x), lty = 'dashed', size = .25) + 
  geom_hline(mapping = aes(yintercept = 0)) +
  geom_vline(mapping = aes(xintercept = 0)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_point(data = tibble(x = 0, y = 0), mapping = aes(x, y), size = .75) + 
  scale_x_continuous(limits = c(-.01,1), expand = c(0, .01)) +
  scale_y_continuous(limits = c(0,1), expand = c(0, .01)) + 
  ylab("") + xlab("") + 
  geom_segment(data = seg1, mapping = aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment"), lty = 'solid', size = .5) + 
  geom_segment(data = seg2, mapping = aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment"), lty = 'solid', size = .5) + 
  geom_segment(data = seg3, mapping = aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment"), lty = 'solid', size = .5) +
  geom_segment(data = seg4, mapping = aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment"), lty = 'solid', size = .5) +
  scale_colour_discrete(guide=FALSE) +
  geom_line(data = first_seg_line, mapping = aes(x = xs, y = ys)) + 
  geom_point(data = tibble(x = x_2, y = coefs_lst[[1]][1]*x_2^3 + coefs_lst[[1]][2]*x_2^2 + coefs_lst[[1]][3]*x_2 + coefs_lst[[1]][4]), 
             mapping = aes(x, y), size = .75)


##

# Getting line segment for partially interpolated Lorenz curve
coefs_lst <- lorenz_to_coefs(lorenz_df, .9, bounds, G)
xs2 <- seq(lorenz_df$x[1], 1, 1/10000)
ys2 <- map(xs2, function(x){
  
  idx <- which(x < lorenz_df$x)[1]  -1
  a <- coefs_lst[[idx]][1]; b <- coefs_lst[[idx]][2]; c <- coefs_lst[[idx]][3]; d <- coefs_lst[[idx]][4]
  
  a*x^3 + b*x^2 + c*x + d
  
}) %>% unlist
next_seg_line <- tibble(x = xs2[1:10000], y = ys2)

# Getting line segments for slopes (bins 3 - 12)
seg_lst <- list()
for(i in 1:15){
  x_val <- lorenz_df$x[i]
  a <- coefs_lst[[i]][1]; b <- coefs_lst[[i]][2]; c <- coefs_lst[[i]][3]; d <- coefs_lst[[i]][4]
  y_val <- a*x_val^3 + b*x_val^2 + c*x_val + d 
  seg_lst[[i]] <- get_segment(c(x_val, y_val), bounds[i]/G, len = .02)
}

# Top point
x_val <- lorenz_df$x[15]
a <- coefs_lst[[15]][1]; b <- coefs_lst[[15]][2]; c <- coefs_lst[[15]][3]; d <- coefs_lst[[15]][4]
y_val <- a*x_val^3 + b*x_val^2 + c*x_val + d 

p1c <- ggplot(lorenz_df) + geom_vline(aes(xintercept = x), lty = 'dashed', size = .25) + 
  geom_hline(mapping = aes(yintercept = 0)) +
  geom_vline(mapping = aes(xintercept = 0)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_point(data = tibble(x = 0, y = 0), mapping = aes(x, y), size = .75) + 
  scale_x_continuous(limits = c(-.01,1), expand = c(0, .01)) +
  scale_y_continuous(limits = c(0,1), expand = c(0, .01)) + 
  ylab("") + xlab("") + 
  geom_line(data = next_seg_line, mapping = aes(x = x, y = y)) + 
  geom_segment(data = seg1, mapping = aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment"), lty = 'solid', size = .5) + 
  geom_segment(data = seg2, mapping = aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment"), lty = 'solid', size = .5) + 
  geom_segment(data = seg_lst[[1]], mapping = aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment"), lty = 'solid', size = .5) + 
  geom_segment(data = seg_lst[[2]], mapping = aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment"), lty = 'solid', size = .5) + 
  geom_segment(data = seg_lst[[3]], mapping = aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment"), lty = 'solid', size = .5) + 
  geom_segment(data = seg_lst[[4]], mapping = aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment"), lty = 'solid', size = .5) + 
  geom_segment(data = seg_lst[[5]], mapping = aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment"), lty = 'solid', size = .5) + 
  geom_segment(data = seg_lst[[6]], mapping = aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment"), lty = 'solid', size = .5) + 
  geom_segment(data = seg_lst[[7]], mapping = aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment"), lty = 'solid', size = .5) + 
  geom_segment(data = seg_lst[[8]], mapping = aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment"), lty = 'solid', size = .5) + 
  geom_segment(data = seg_lst[[9]], mapping = aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment"), lty = 'solid', size = .5) + 
  geom_segment(data = seg_lst[[10]], mapping = aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment"), lty = 'solid', size = .5) + 
  geom_segment(data = seg_lst[[11]], mapping = aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment"), lty = 'solid', size = .5) + 
  geom_segment(data = seg_lst[[12]], mapping = aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment"), lty = 'solid', size = .5) + 
  geom_segment(data = seg_lst[[13]], mapping = aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment"), lty = 'solid', size = .5) + 
  geom_segment(data = seg_lst[[14]], mapping = aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment"), lty = 'solid', size = .5) + 
  geom_segment(data = seg_lst[[15]], mapping = aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment"), lty = 'solid', size = .5) + 
  geom_point(tibble(xval=1, yval=1), mapping = aes(xval, yval), size = .75) + 
  scale_colour_discrete(guide=FALSE) +
  geom_line(data = first_seg_line, mapping = aes(x = xs, y = ys)) + 
  geom_label_repel(data = tibble(xtext = .99, ytext = .99, lab = "Special Cubic For Top Bin"), mapping = aes(xtext, ytext, label = lab), 
                   nudge_x = -.4, nudge_y = -.15, size = 2.5)




p1 <- grid.arrange(p1a + labs(title = "", y= "Cumulative Income Share") + theme(plot.title = element_text(size = 9, hjust = .5), axis.text.x = element_text(size = 7),axis.text.y = element_text(size = 7), axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8, angle = 90)),
             p1b + labs(title = "", x = "F(x) - Cumulative Population Share")+ theme(plot.title = element_text(size = 9, hjust = .5),axis.text.x = element_text(size = 7),axis.text.y = element_text(size = 6), axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8)), 
             p1c + labs(title = "")+ theme(plot.title = element_text(size = 9, hjust = .5),axis.text.x = element_text(size = 7),axis.text.y = element_text(size = 7), axis.title.x = element_text(size = 8),axis.title.y = element_text(size = 8)),
             nrow = 1, ncol = 3)
p1



####################################
############ Plot 2 ################
####################################


# Pulling single PUMA from data -- change puma_num to get estimated PDF of a different PUMA
puma_num <- 130 
sample_puma <- pums_data %>% group_split(area) %>% .[[puma_num]]
freqs <- sample_puma %>% group_incomes() %>% group_by(income_cat) %>% summarise(n = n(), .groups = 'drop_last') %>% .[["n"]]
G <- sum(sample_puma$income)/length(sample_puma$income)
lorenz_df <- tibble(x = c(0, cumsum(freqs)/sum(freqs))) %>% distinct()
lorenz_coefs <- lorenz_to_coefs(lorenz_df, slope_parm = .9, bounds, G)
bounds_mod <- c(bounds, 400000)

# Getting tibble with x,y coordinates of PDF
pdf_tbl <- map2(seq_along(bounds_mod[1:(length(bounds_mod)-1)]), lorenz_coefs, function(idx, coefs){
  
  lb <- bounds_mod[idx]
  ub <- bounds_mod[idx+1]
  a <- coefs[1]; b <- coefs[2]; c <- coefs[3]
  
  xs <- seq(lb, ub, (ub-lb)/500)
  
  if(idx != 15){
    ys <- map(xs, function(x){ return( 1/(G*sqrt(4*b^2 - 12*a*(c - x/G))) ) }) %>% unlist
  }else{
    
    ys <- map(xs, function(x){ return( min(1/(G*sqrt(4*b^2 - 12*a*(c - x/G))), .000004) ) }) %>% unlist
    
  }
  
  return(tibble(x = xs, y = ys))
  
}) %>% do.call('rbind', .)


p2 <- ggplot(pdf_tbl, aes(x = x, y = y)) + geom_line() + ylim(0, .000015) + 
  xlab("x") + ylab("f(x)") + 
  scale_x_continuous(breaks = c(0, 100000, 200000, 300000, 400000), labels = c("0", "$100,000", "$200,000", "$300,000", "$400,000"))
p2






####################################
############ Plot 3 ################
####################################


# Get income statistics estimated by Lorenz interpolation, CDF interpolation, and MCIB
income_statistics <- get_income_statistics(input_data = pums_data, input_func = lorenz_int, sum_stats = FALSE, num_cores = 11)

# Creating tibble of Atkinson estimates of each method
data_for_plot <- income_statistics %>% 
  select(mcib_at, 
         cdf_at, 
         lor_at, 
         true_at) %>% 
  mutate(
    MCIB = mcib_at-true_at, 
    CDF =   cdf_at-true_at,
    Lorenz = lor_at-true_at) %>%
  select(MCIB, CDF, Lorenz, true_at) %>% 
  gather(metric, value,-true_at) %>% 
  mutate(metric = fct_relevel(metric, c("MCIB", "CDF", "Lorenz")))
hline_data <- data_for_plot %>% group_by(metric) %>% summarise(bias=mean(value), sd=sd(value))

# Atkinson plot
p3 <- data_for_plot %>% 
  ggplot(aes(true_at,value)) +
  facet_wrap(~metric) + 
  geom_point(size = .2) + 
  geom_hline(hline_data, mapping = aes(yintercept = bias)) + 
  geom_hline(hline_data, mapping = aes(yintercept = bias-sd), lty='dashed') + 
  geom_hline(hline_data, mapping = aes(yintercept = bias+sd), lty='dashed') + 
  expand_limits(y = c(-.02, .02)) + 
  labs(x="Atkinson Index", y="Residual") + 
  theme(plot.title = element_text(hjust = 0.5))
p3

