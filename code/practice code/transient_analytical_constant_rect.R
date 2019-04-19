source("C:/Users/Vic/Box/Quennessen_Thesis/practice code/Leslie.R")
source("C:/Users/Vic/Box/Quennessen_Thesis/practice code/do_Leslie2_constrect.R")

transient_analytical_constant_rect = function() {
  
  ### Run analytical (linear) transient model, with constant recruitment
  
  ##################### Life history parameters of species ###################

  ###### Black rockfish (Sebastes melanops)
  L_inf <- 44.2
  k <- 0.33
  t0 <- 0.075
  a2 <- 1.68e-5
  b2 <- 3
  a <- 1
  b <- 3
  a_mat <- 12 
  t_c <- 4
  max_age <- 35
  M <- 0.14
  f <- 0.3
  ############################################################################
  
  ### Create size distribution
  
  m_at_age <- 0:max_age
  
  # Define length 
  l_a <- exp(-1*M*m_at_age)
  len_a <- L_inf*(1-exp(-1*k*(m_at_age - t0)))
  len_a <- ifelse(len_a < 0, 0, len_a) # replace any negative len_a values with 0
  std_a <- len_a*0.1                   # standard deviation
  
  len_vec <- seq(0, L_inf*0.3, 1)      # length bins, real value 0.3 = 1.3
  # matrix of length bins
  len_mat <- matrix(rep(len_vec, length(m_at_age)), ncol = length(m_at_age))
  # matrix of mean lengths
  len_mat2 <- matrix(rep(len_a, each = length(len_vec)), nrow = length(len_vec))
  # matrix of difference from mean length
  len_mat3 <- len_mat - len_mat2
  # matrix of standard deviations
  std_mat <- matrix(rep(std_a, each = length(len_vec)), nrow = length(len_vec))
  
  # probability of being in each bin
  LM <- pnorm(len_mat3 + 0.5, 0, std_mat) - pnorm(len_mat3 - 0.5, 0, std_mat)
  LM[which(is.nan(LM))] <- 0
  
  ###########################################################################
  
  D_target <- 1
  tc_vec <- matrix(c(4, 8), ncol = 2)
  len <- length(tc_vec)
  Rfact <- 0.0
  t <- 200
  
  Z_start <- Nt <- Nt_i <- array(rep(0, (max_age + 1)*(t+1)*len), c(max_age + 1, t + 1, len))
  D_out <- D_out_i <- array(rep(0, (t + 1)*len), c(t + 1, len))
  L1_tc <- L2_tc <- L1_tc_i <- L2_tc_i <- array(rep(0, 2), dim(tc_vec))
  
  for (i in 1:length(tc_vec)) {
    Z_start[, , i] <- do_Leslie2_constrect(max_age, a_mat, tc_vec[i], L_inf, k, t0, 
                                   a2, b2, a, M, t, NaN, f)[[2]]
    
    obj <- do_Leslie2_constrect(max_age, a_mat, tc_vec[i], L_inf, k, t0, 
                                     a2, b2, a, M, t, Z_start[, t + 1, i], 0)
    Nt[, , i]    <- obj[[2]]
    L1_tc[i] <- obj[[3]]
    L2_tc[i] <- obj[[4]]
    D_out[, i] <- obj[[5]]
    
    obj_i <- do_Leslie2_constrect(max_age, a_mat, tc_vec[i], L_inf, k, t0, 
                               a2, b2, a, M, t, Z_start[, t + 1, i], 0, Rfact)
    Nt_i[, , i]    <- obj_i[[2]]
    L1_tc_i[i] <- obj_i[[3]]
    L2_tc_i[i] <- obj_i[[4]]
    D_out_i[, i] <- obj_i[[5]]
  }
  
  doFigs1 = T
  
  if (doFigs1 == T) {
    par(mfrow = c(1,2))
    x = seq(2, 10, by = 2)
    for (k in 1:len) {
      for (i in length(x)) {
        plot(0:max_age, Nt[, x[i], k])
      }
    }
  }
  
}
  