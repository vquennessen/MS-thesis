#' Victoria Quennessen
#' Age-Structured Population Dynamics Biological Sub-Model
#' 4/18/2019

popDyn = function(max_age,                  # maximum age
                  M,                        # natural mortality
                  rec_age,                  # age at recruitment
                  af, bf,                   # weight at length parameters (f)
                  am, bm,                   # weight at length parameters (m)
                  a1f, L1f, a2f, L2f, Kf,   # growth parameters (f)
                  a1m, L1m, a2m, L2m, Km,   # growth parameters (m)
                  L50,                      # length at 50% maturity
                  k_mat,                    # slope of maturity curve
                  ldp,                      # larval drift proportion
                  R0,                       # unfished recruitment
                  h,                        # steepness
                  phi,                      # unfished recruits per spawner
                  sigma_R,                  # recruitment standard deviation
                  rho_R,                    # recruitment autocorrelation
                  p,                        # adult movement proportion
                  D,                        # depletion
                  Fb,                       # fishing mortality to cause D
                  r,                        # Proportion of positive transects 
                                            #       in PISCO monitoring data
                  x,                        # mean of positive transects
                  sp) {                     # std of positive transects
  
  parameters("black rockfish")
}                     