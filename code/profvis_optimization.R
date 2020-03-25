library(profvis)
profvis({
  
  Species = 'BR_OR_2015'
  R0 = 1e5
  A = 5
  MPA = 3
  Time1 = 50
  Time2 = 20
  Error = 0.05
  Transects = 24
  Plotting = FALSE
  Final_DRs = c(0.2, 0.4, 0.6, 0.8, 1.0)
  Control_rules = c(1:6)
  Plot_individual_runs <- F
  LDP = 0.1
  BM <- FALSE
  Output.Yield = TRUE
  Output.Biomass = TRUE
  Output.SSB = TRUE
  Output.N = TRUE
  Output.Density.Ratio = TRUE
  Output.Effort = TRUE
  Output.FM = FALSE
  Output.Catch = FALSE
  Output.Abundance = FALSE
  Surveys = TRUE
  Fishing = TRUE
  Stochasticity = TRUE
  Recruitment_mode = 'pool'
  Allocation = 'IFD'
  Adult_movement = TRUE
  Fishery_management = TRUE
  Years_sampled = 1
  Areas_sampled = 'all'
  Ind_sampled = 'all'
  
  ###### Error handling ########################################################
  
  # classes of variables
  if (!is.character(Species)) {
    stop('Study species must be a character string.')}
  if (R0 %% 1 != 0) {stop('R0 must be an integer value.')}
  if (A %% 1 != 0) {stop('A must be an integer value.')}
  if (MPA %% 1 != 0) {stop('MPA must be an integer value.')}
  if (Time1 %% 1 != 0) {stop('Time1 must be an integer value.')}
  if (Time2 %% 1 != 0) {stop('Time2 must be an integer value.')}
  if (!is.character(Recruitment_mode)) {
    stop('Recruitment mode must be a character value.')}
  if (!is.numeric(Error)) {stop('Error must be a numeric value.')}
  if (!is.logical(Stochasticity)) {
    stop('Stochasticity must be a logical value.')}
  if (!is.logical(Surveys)) {stop('Surveys must be a logical value.')}
  if (!is.logical(Fishery_management)) {
    stop('Fishery_management must be a logical value.')}
  if (!is.logical(Fishing)) {stop('Fishing must be a logical value.')}
  if (Transects %% 1 != 0) {stop('Transects must be an integer value.')}
  if (!is.logical(Adult_movement)) {
    stop('Adult_movement must be a logical value.')}
  if (!is.logical(Plotting)) {stop('Plotting must be a logical value.')}
  if (!is.numeric(Final_DRs)) {stop('Final_DRs must be a numeric vector.')}
  if (Years_sampled %% 1 != 0 && !is.null(Years_sampled)) {
    stop('Years_sampled must be an integer value or NULL.')}
  if (!is.character(Areas_sampled) && !is.null(Areas_sampled)) {
    stop('Areas_sampled must be a character value or NULL.')}
  if (!is.character(Ind_sampled) && !is.null(Ind_sampled)) {
    stop('Ind_sampled must be a character value or NULL.')}
  if (!is.character(Allocation)) {stop('Allocation must be a character value.')}
  if (!is.logical(BM)) {stop('BM must be a logical value.')}
  if (!is.numeric(LDP)) {stop('LDP must be a numeric value.')}
  if (sum(Control_rules %% 1 != 0) != 0) {
    stop('Control_rules must be a vector of integers.')}
  if (!is.logical(Output.FM)) {stop('Output.FM must be a logical value.')}
  if (!is.logical(Output.N)) {stop('Output.N must be a logical value.')}
  if (!is.logical(Output.Abundance)) {
    stop('Output.Abundance must be a logical value.')}
  if (!is.logical(Output.Biomass)) {
    stop('Output.Biomass must be a logical value.')}
  if (!is.logical(Output.SSB)) {stop('Output.SSB must be a logical value.')}
  if (!is.logical(Output.Catch)) {stop('Output.Catch must be a logical value.')}
  if (!is.logical(Output.Yield)) {stop('Output.Yield must be a logical value.')}
  if (!is.logical(Output.Effort)) {
    stop('Output.Effort must be a logical value.')}
  if (!is.logical(Output.Density.Ratio)) {
    stop('Output.Density.Ratio must be a logical value.')}
  
  # acceptable values
  if (R0 <= 0) {stop('R0 must be greater than 0.')}
  if (A <= 0) {stop('A must be greater than 0.')}
  if (MPA <= 0) {stop('MPA must be greater than 0.')}
  if (Time1 <= 0) {stop('Time1 must be greater than 0.')}
  if (Time2 <= 0) {stop('Time2 must be greater than 0.')}
  if (Recruitment_mode != 'pool' && Recruitment_mode != 'closed' &&
      Recruitment_mode != 'regional_DD' && Recruitment_mode != 'local_DD') {
    stop('Recruitment_mode must be either "pool", "closed", "regional_DD", or
         "local_DD".')}
  if (Error < 0) {stop('Error must be greater than or equal to 0.')}
  if (Transects <= 0) {stop('Transects must be greater than 0.')}
  if (sum(Final_DRs <= 0) > 0) {
    stop('All values in Final_DRs must be greater than 0.')}
  if (Years_sampled <= 0 && !is.null(Years_sampled)) {
    stop('Years_sampled must be greater than 0 or NULL.')}
  if (is.numeric(Years_sampled) && Years_sampled <= 0) {
    stop('Years_sampled must be greater than 0 or NULL.')}
  if (is.character(Areas_sampled) && Areas_sampled != 'far' &&
      Areas_sampled != 'all' ) {
    stop('Areas_sampled must be either "far" or "all" or NULL.')}
  if (is.character(Ind_sampled) && Ind_sampled != 'mature' &&
      Ind_sampled != 'all') {
    stop('Ind_sampled must be either "mature" or "all" or NULL.')}
  if (LDP < 0 || LDP > 1) {stop('LDP must be between 0 and 1.')}
  if (sum(Control_rules <= 0) > 0) {
    stop('All values in Control_rules must be greater than 0.')}
  if (is.null(Years_sampled) && BM == FALSE) {
    stop('BM must be TRUE for Years_sampled to be NULL.')}
  if (is.null(Areas_sampled) && BM == FALSE) {
    stop('BM must be TRUE for Areas_sampled to be NULL.')}
  if (is.null(Ind_sampled) && BM == FALSE) {
    stop('BM must be TRUE for Ind_sampled to be NULL.')}
  
  # relationtional values
  if (MPA > A) {stop('MPA must be less than or equal to A.')}
  
  ##############################################################################
  
  ##### Load life history characteristics for species ##########################
  
  par <- parameters(Species)
  
  Max_age                <- par[[1]]        # maximum age
  M                      <- par[[2]]        # natural mortality
  Rec_age                <- par[[3]]        # age at recruitment
  WA <- par[[4]];  WB    <- par[[5]]        # weight at length parameters (f)
  A1 <- par[[6]];  L1    <- par[[7]]        # growth parameters (f)
  A2 <- par[[8]];  L2    <- par[[9]]
  K                      <- par[[10]]
  L50                    <- par[[11]]       # length at 50% maturity
  K_mat                  <- par[[12]]       # slope of maturity curve
  H                      <- par[[13]]       # steepness
  Phi                    <- par[[14]]       # unfished recruits per spawner
  Sigma_R                <- par[[15]]       # recruitment standard deviation
  Rho_R                  <- par[[16]]       # recruitment autocorrelation
  AMP                    <- par[[17]]       # adult movement proportion
  D                      <- par[[18]]       # depletion
  Fb                     <- par[[19]]       # fishing mortality to cause D
  P                      <- par[[20]]       # proportion of positive transects
  X                      <- par[[21]]       # mean of positive transects
  SP                     <- par[[22]]       # std of positive transects
  Fleets                 <- par[[23]]       # fishery fleet names
  Alpha                  <- par[[24]]       # slope for upcurve
  Beta                   <- par[[25]]       # slope for downcurve
  F_fin                  <- par[[26]]       # F_fin for fishery, 0 if asymptotic
  A50_up                 <- par[[27]]       # A50 for upcurve
  A50_down               <- par[[28]]       # A50 for downcurve
  Cf                     <- par[[29]]       # fraction of fishery caught / fleet
  
  ##### Population Dynamics - Non-Time Varying #################################
  
  # Initialize arrays for time-varying dynamics
  IA <- initialize_arrays(A, MPA, Final_DRs, Time1, Time2, R0, Rec_age, Max_age,
                          A1, L1, A2, L2, K, WA, WB, K_mat, Fb, L50, Sigma_R,
                          Rho_R, Fleets, Alpha, A50_up, A50_down, F_fin, Beta,
                          Cf, P, X, SP, M, Control_rules, Phi, Stochasticity, D,
                          Transects, H, Surveys, Fishing, Error,
                          Recruitment_mode, LDP, Ind_sampled)
  
  Inside           <- IA[[1]]     # Area(s) in the marine reserve
  Outside          <- IA[[2]]     # Areas not in the marine reserve
  FDR              <- IA[[3]]
  TimeT            <- IA[[4]]     # total amount of timesteps (years)
  L                <- IA[[5]]     # Length at age, dim = 1*age
  W                <- IA[[6]]     # Weight at age, dim = 1*age
  S                <- IA[[7]]     # Selectivity at age
  Mat              <- IA[[8]]     # Fraction mature at age, dim = 1*age
  A50_mat          <- IA[[9]]     # Age at which fraction mature > 0.5
  CR               <- IA[[10]]    # Number of control rules
  Nat_mortality    <- IA[[11]]    # Range of potential natural mortality values
  NM               <- IA[[12]]    # Number of potential natural mortality values
  N                <- IA[[13]]    # Population size, dim = age*area*time
  SSB              <- IA[[14]]    # Spawning stock biomass, dim = area*time
  Biomass          <- IA[[15]]    # Biomass, dim = area*time
  Eps              <- IA[[16]]    # Epsilon vector, dim = area*time*CR
  B0               <- IA[[17]]    # Unfished spawning stock biomass
  Count            <- IA[[18]]    # Species count when sampling, dim = area*time
  Sigma_S          <- IA[[19]]    # Sampling normal standard deviation
  NuS              <- IA[[20]]    # Sampling normal variable, dim = area*time*CR
  Delta            <- IA[[21]]    # Constant of proportionality
  Gamma            <- IA[[22]]    # Gamma
  FM               <- IA[[23]]    # Fishing mortality rate, dim = age*area*time
  E                <- IA[[24]]    # nominal fishing effort in each area
  Catch            <- IA[[25]]    # Catch at age
  Yield            <- IA[[26]]    # Yield per area
  Density_ratio    <- IA[[27]]    # Density ratios
  ENM              <- IA[[28]]    # nm value that represents "true" M
  Abundance        <- IA[[29]]    # Abundance of all and/or mature individuals
  
  ##### Population Dynamics - Time Varying #####################################
  for (fdr in 1:FDR) {
    
    # for (t in Time1:TimeT) {
    for (t in (Rec_age + 1):TimeT) {
      
      for (cr in 1:CR) {
        
        for (nm in 1:NM) {
          
          # effort allocation
          E <- effort_allocation(t, cr, nm, fdr, Allocation, E, Yield, Time1,
                                 Inside, Outside)
          
          # If there is adult movement, add movement
          if (Adult_movement == TRUE) {N <- movement(t, cr, nm, fdr, N, A, AMP)}
          
          # Recruitment / larval movement (if applicable)
          R <- recruitment(t, cr, nm, fdr, SSB, A, R0, H, B0, Eps, Sigma_R,
                           Rec_age, Recruitment_mode, LDP)
          
          
          
          # biology
          PD <- pop_dynamics(t, cr, nm, fdr, Rec_age, Max_age, SSB,
                             N, W, Mat, A, Fb, E, S, NM, FM, A50_mat,
                             Biomass, Abundance, Fishing, Nat_mortality, R,
                             Ind_sampled)
          
          FM[, , t, cr, nm, fdr]               <- PD[[1]]
          N[, , t, cr, nm, fdr]                <- PD[[2]]
          Biomass[, t, cr, nm, fdr]            <- PD[[3]]
          SSB[, t, cr, nm, fdr]                <- PD[[4]]
          Abundance[, t, cr, nm, fdr, ]        <- PD[[5]]
          
          # sampling
          if (Surveys == TRUE) {
            Count[, t, , , cr, nm, fdr] <- sampling(t, cr, nm, fdr, Delta,
                                                    Gamma, Abundance, Transects,
                                                    X, Count, NuS, A,
                                                    Ind_sampled)
          }
          
          # fishing
          if (Fishing == TRUE) {
            Catch[, , t, cr, nm, fdr] <- catch(t, cr, nm, fdr, FM,
                                               Nat_mortality, N, A, Fb, E,
                                               Catch)
            Yield[, t, cr, nm, fdr] <- colSums(Catch[, , t, cr, nm, fdr]*W)
          }
          
        }
        
        # management
        if (Fishery_management == TRUE && t > Time1 && t < TimeT) {
          E[, t, cr, , fdr] <- control_rule(t, cr, nm, fdr, A, E, Count, Time1,
                                            TimeT, Transects, Nat_mortality,
                                            Final_DRs, Inside, Outside,
                                            Areas_sampled, Ind_sampled,
                                            Years_sampled, BM)
        }
        
        # calculate true density ratio
        Density_ratio[t, cr, fdr] <- true_DR(t, cr, fdr, Abundance, Inside,
                                             Outside, Density_ratio, ENM)
      }
    }
  }
  
  
})
