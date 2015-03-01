# functions to evaluate model in Dowdy et al. (2006)
tbdx_neg <- function(
## population size
   pop=c(48727, 50984, 15, 44, 72, 158),      # starting size of 6 compartments
   allow.hiv=FALSE,                           # flag to allow HIV compartments
## model span
   year_start=0, year_stop=200, timestep=0.02,
   age_lower=15, age_upper=49, agespan=35,    # age_upper-age_lower+1
## TB transmission
   lam_0=1.8195e-4,                           # "calibrated to provide an overall annual risk of TB infection of 2.0%/year at equilibrium"
   i_smneg=0.22,                              # relative infectiousness of smear neg cases
   lat=0.2614309,                             # 1-(1-timestep)^age_lower; LTBI prevalence in lower age entrants; value at equilib when ARTI=0.02
   fp_neg=0.45, fp_hiv=0.45, fp_aids=0.35,    # fraction incident TB that is smear pos
## TB progression
   pp_neg=0.07, pp_hiv=0.07, pp_aids=0.56,    # proportion TB infections primary progressors
   lr_neg=0.28, lr_hiv=0.28, lr_aids=0.75,    # proportion TB reinfections in primary progressors
   er_neg=0.001, er_hiv=0.001, er_aids=0.079, # endogenous reactivation rate
   cure=0, # rate spontaneous cure active TB back to LTBI; sensitivity analysis only
## TB diagnosis and treatment
   acc_neg=0.734, acc_hiv=0.734, acc_aids=0.734, # fraction TB patients presenting to clinics
   tx=0.82,                                   # success of Tx once TB diagnosed
   dx=4,                                      # rate of TB detection after symptoms develop
   sn_smpos=0.80,                             # sensitivity diagnosis smear-pos; 0.80 standard, 0.98 molecular test or culture
   sn_smneg=0.25,                             # sensitivity diagnosis smear-neg; 0.20 standard, 0.50 molecular, 0.85 culture
   lfu=0,                                     # fraction of detected patients lost to f/u before Tx, zero at baseline
## HIV transmission and progression
   hivinc=0, #0.035,                          # HIV incidence; provides an equilibrium HIV prevalence of 16.8%
   hivprog=0, #1.739,                         # HIV progression to AIDS; inverse of mean time from HIV infection to AIDS
## mortality and aging
   mu_neg=0.009 + 1/agespan,                  # mortality: HIV-, no TB
   mu_smpos_neg=0.35 + 1/agespan,             # mortality: HIV-, smear+ TB
   mu_smneg_neg=0.1 + 1/agespan,              # mortality: HIV-, smear- TB
   mu_hiv=0.056 + 1/agespan,                  # mortality: HIV+, no TB
   mu_tb_hiv=0.777 + 1/agespan,               # mortality: HIV+, active TB  0.777 ~ 1-exp(-1.5)
   mu_aids=0.352 + 1/agespan,                 # mortality: AIDS, no TB
   mu_tb_aids=0.865 + 1/agespan               # mortality: AIDS, active TB  0.865 ~ 1-exp(-2.0)
) {

## set up objects
   # times over which model is evaluated
   times <- seq(from=year_start, to=year_stop, by=timestep) # more efficient to receive times rather than calculate each time
   ntimes <- length(times)
   # objects in which store lambda(t) and compartment sizes over time
   TSP <- TSN <- lam <- M <- rep(0, ntimes)
   N_neg <- matrix(0, nr=ntimes, nc=6, 
      dimnames=list(NULL, c("S", "L", "DSP", "USP", "DSN", "USN")))

   # initial conditions at first time
   N_neg[1, ] <- pop[1:6]
   #A. TSP
   TSP[1] <- sum(N_neg[1, 3:4])
   #B. TSN
   TSN[1] <- sum(N_neg[1, 5:6])
   #D. M
   M[1] <- mu_neg * sum(N_neg[1, 1:2]) + mu_smpos_neg * sum(N_neg[1, 3:4]) + 
      + mu_smneg_neg * sum(N_neg[1, 5:6])
   #C. lambda
   lam[1] <- lam_0 * (TSP[1] + i_smneg * TSN[1])

   for(tt in seq_len(ntimes - 1)) {
   # 1. S_neg
   N_neg[tt + 1, 1] <- N_neg[tt, 1] + timestep * 
      ( (1 - lat) * M[tt] +
      - (mu_neg + hivinc + lam[tt]) * N_neg[tt, 1] )
   # 2. L_neg
   N_neg[tt + 1, 2] <- N_neg[tt, 2] + timestep * 
      ( lat * M[tt] +
      + (1 - pp_neg) * lam[tt] * N_neg[tt, 1] +
      - (mu_neg + hivinc + er_neg + lr_neg * pp_neg * lam[tt]) * N_neg[tt, 2] +
      + (dx * tx * sn_smpos) * N_neg[tt, 3] +
      + (dx * tx * sn_smneg) * N_neg[tt, 5] +
      + cure * sum(N_neg[tt, 3:6]) )
   # 3. DSP_neg
   N_neg[tt + 1, 3] <- N_neg[tt, 3] + timestep * 
      ( fp_neg * acc_neg * pp_neg * lam[tt] * N_neg[tt, 1] +
      + fp_neg * acc_neg * (er_neg + lr_neg * pp_neg * lam[tt]) * N_neg[tt, 2] +
      - (mu_smpos_neg + hivinc + dx * tx * sn_smpos + cure) * N_neg[tt, 3] )
   # 4. USP_neg
   N_neg[tt + 1, 4] <- N_neg[tt, 4] + timestep * 
      ( fp_neg * (1 - acc_neg) * pp_neg * lam[tt] * N_neg[tt, 1] +
      + fp_neg * (1 - acc_neg) * (er_neg + lr_neg * pp_neg * lam[tt]) * N_neg[tt, 2] +
      - (mu_smpos_neg + hivinc + cure) * N_neg[tt, 4] )
   # 5. DSN_neg
   N_neg[tt + 1, 5] <- N_neg[tt, 5] + timestep * 
      ( (1 - fp_neg) * acc_neg * pp_neg * lam[tt] * N_neg[tt, 1] +
      + (1 - fp_neg) * acc_neg * (er_neg + lr_neg * pp_neg * lam[tt]) * N_neg[tt, 2] +
      - (mu_smneg_neg + hivinc + dx * tx * sn_smneg + cure) * N_neg[tt, 5] )
   # 6. USN_neg
   N_neg[tt + 1, 6] <- N_neg[tt, 6] + timestep * 
      ( (1 - fp_neg) * (1 - acc_neg) * pp_neg * lam[tt] * N_neg[tt, 1] +
      + (1 - fp_neg) * (1 - acc_neg) * (er_neg + lr_neg * pp_neg * lam[tt]) * N_neg[tt, 2] +
      - (mu_smneg_neg + hivinc + cure) * N_neg[tt, 6] )

   #A. TSP
   TSP[tt + 1] <- sum(N_neg[tt + 1, 3:4])
   #B. TSN
   TSN[tt + 1] <- sum(N_neg[tt + 1, 5:6])
   #D. M
   M[tt + 1] <- mu_neg * sum(N_neg[tt + 1, 1:2]) +
      + mu_smpos_neg * sum(N_neg[tt + 1, 3:4]) +
      + mu_smneg_neg * sum(N_neg[tt + 1, 5:6])
   #C. lam
   lam[tt + 1] <- lam_0 * (TSP[tt + 1] + i_smneg * TSN[tt + 1])
   } # end for(tt in seq_len(ntimes - 1))
   compart <- N_neg #cbind(N_neg, N_hiv, N_aids)

   #mortality
   mortality <- (mu_smpos_neg-1/35)*TSP+(mu_smneg_neg-1/35)*TSN

   #prevalence
   prevalence <- TSN + TSP

   #incidences from susceptible, latent, reinfections capable of primary progression
   inc1 <- (pp_neg*lam)*N_neg[tt+1, 1]
   inc2 <- er_neg*N_neg[tt+1, 2]
   inc3 <- lr_neg*pp_neg*lam*N_neg[tt+1, 2]
   incidence <- inc1+ inc2 + inc3

   dimnames(compart)[[2]] <- paste(c("S", "L", "DSP", "USP", "DSN", "USN"), "neg", sep="_")
   invisible(list(pop=compart, lam=lam, TSP=TSP, TSN=TSN, M=M, incidence=incidence, prevalence=prevalence, mortality=mortality))


}
   x_neg <- tbdx_neg()

   do.call("cbind", x_neg)[c(1:6, 9996:10001),]


   plot(x_neg$lam, type="l")
   plot(x_neg$incidence, type="l")
   plot(x_neg$mortality, type="l")
   plot(x_neg$prevalence, type="l")



