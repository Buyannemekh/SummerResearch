# functions to evaluate model in Dowdy et al. (2006)
# tbdx takes input parameters essentially as specified in manuscript
#    exception: acc is now acc_{neg,hiv,aids}
#    initializes objects and iterates over them
# tbdx_neg conveniently limits to HIV-neg compartments

tbdx <- function(
## population size
   pop=c(48195, 34421, 13, 33, 57, 103, 6762, 5799, 2, 3, 6, 6, 
      2446, 2014, 15, 22, 63, 40),            # starting size of 18 compartments
   allow.hiv=TRUE,                            # flag to allow HIV compartments
## model span
   year_start=0, year_stop=200, timestep=0.02,
   age_lower=15, age_upper=49, agespan=35,    # age_upper-age_lower+1
## TB transmission
   lam_0= 1.483185e-4,  #0.0001519506,    #                          # "calibrated to provide an overall annual risk of TB infection of 2.0%/year at equilibrium"
   i_smneg=0.22,                              # relative infectiousness of smear neg cases
   lat=1-(0.98)^15,                             # 1-(1-timestep)^age_lower; LTBI prevalence in lower age entrants; value at equilib when ARTI=0.02
   fp_neg=0.45, fp_hiv=0.45, fp_aids=0.35,    # fraction incident TB that is smear pos
## TB progression
   pp_neg=0.07, pp_hiv=0.07, pp_aids=0.56,    # proportion TB infections primary progressors
   lr_neg=0.28, lr_hiv=0.28, lr_aids=0.75,    # proportion TB reinfections in primary progressors
   er_neg=0.001, er_hiv=0.079/5.1, er_aids=0.079, # endogenous reactivation rate
   cure=0, # rate spontaneous cure active TB back to LTBI; sensitivity analysis only
## TB diagnosis and treatment
   acc_neg=0.734, acc_hiv=0.734, acc_aids=0.734, # fraction TB patients presenting to clinics
   tx=0.82,                                   # success of Tx once TB diagnosed
   dx_pos=4, dx_neg=4,                                      # rate of TB detection after symptoms develop
   sn_smpos=0.80,                             # sensitivity diagnosis smear-pos; 0.80 standard, 0.98 molecular test or culture
   sn_smneg=0.20,                             # sensitivity diagnosis smear-neg; 0.20 standard, 0.50 molecular, 0.85 culture
   lfu=0,                                     # fraction of detected patients lost to f/u before Tx, zero at baseline
## HIV transmission and progression
   hivinc=0.035,                              # HIV incidence; provides an equilibrium HIV prevalence of 16.8%
   hivprog=1/6.9,                            # HIV progression to AIDS; inverse of mean time from HIV infection to AIDS
## mortality and aging
   mu_neg=0.009 + 1/agespan,                  # mortality: HIV-, no TB
   mu_smpos_neg=0.35 + 1/agespan,             # mortality: HIV-, smear+ TB
   mu_smneg_neg=0.1 + 1/agespan,              # mortality: HIV-, smear- TB
   mu_hiv=0.056 + 1/agespan,                  # mortality: HIV+, no TB
   mu_tb_hiv=1.5 + 1/agespan,               # mortality: HIV+, active TB  0.777 ~ 1-exp(-1.5)
   mu_aids=0.352 + 1/agespan,                 # mortality: AIDS, no TB
   mu_tb_aids=2 + 1/agespan               # mortality: AIDS, active TB  0.865 ~ 1-exp(-2.0)
) {

## set up objects
   # times over which model is evaluated
   times <- seq(from=year_start, to=year_stop, by=timestep) # more efficient to receive times rather than calculate each time
   ntimes <- length(times)
   # objects in which store lambda(t) and compartment sizes over time
   TSP <- TSN <- lam <- M <- rep(0, ntimes)
   N_neg <- N_hiv <- N_aids <- matrix(0, nr=ntimes, nc=6, 
      dimnames=list(NULL, c("S", "L", "DSP", "USP", "DSN", "USN")))

   # initial conditions at first time
   N_neg[1, ] <- pop[1:6]
   N_hiv[1, ] <- pop[7:12]
   N_aids[1, ] <- pop[13:18]
   #A. TSP
   TSP[1] <- sum(N_neg[1, 3:4], N_hiv[1, 3:4], N_aids[1, 3:4])
   #B. TSN
   TSN[1] <- sum(N_neg[1, 5:6], N_hiv[1, 5:6], N_aids[1, 5:6])
   #D. M
   M[1] <- mu_neg * sum(N_neg[1, 1:2]) + mu_smpos_neg * sum(N_neg[1, 3:4]) + 
      + mu_smneg_neg * sum(N_neg[1, 5:6]) +
      + mu_hiv * sum(N_hiv[1, 1:2]) + mu_tb_hiv * sum(N_hiv[1, 3:6]) +
      + mu_aids * sum(N_aids[1, 1:2]) + mu_tb_aids * sum(N_aids[1, 3:6])
   #C. lambda
   lam[1] <- lam_0 * (TSP[1] + i_smneg * TSN[1])

## iterate over time, populating compartments by row
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
      + (dx_pos * tx * sn_smpos) * N_neg[tt, 3] +
      + (dx_neg * tx * sn_smneg * (1-lfu)) * N_neg[tt, 5] +
      + cure * sum(N_neg[tt, 3:6]) )
   # 3. DSP_neg
   N_neg[tt + 1, 3] <- N_neg[tt, 3] + timestep * 
      ( fp_neg * acc_neg * pp_neg * lam[tt] * N_neg[tt, 1] +
      + fp_neg * acc_neg * (er_neg + lr_neg * pp_neg * lam[tt]) * N_neg[tt, 2] +
      - (mu_smpos_neg + hivinc + dx_pos * tx * sn_smpos + cure) * N_neg[tt, 3] )
   # 4. USP_neg
   N_neg[tt + 1, 4] <- N_neg[tt, 4] + timestep * 
      ( fp_neg * (1 - acc_neg) * pp_neg * lam[tt] * N_neg[tt, 1] +
      + fp_neg * (1 - acc_neg) * (er_neg + lr_neg * pp_neg * lam[tt]) * N_neg[tt, 2] +
      - (mu_smpos_neg + hivinc + cure) * N_neg[tt, 4] )
   # 5. DSN_neg
   N_neg[tt + 1, 5] <- N_neg[tt, 5] + timestep * 
      ( (1 - fp_neg) * acc_neg * pp_neg * lam[tt] * N_neg[tt, 1] +
      + (1 - fp_neg) * acc_neg * (er_neg + lr_neg * pp_neg * lam[tt]) * N_neg[tt, 2] +
      - (mu_smneg_neg + hivinc + dx_neg * tx * sn_smneg * (1-lfu) + cure) * N_neg[tt, 5] )
   # 6. USN_neg
   N_neg[tt + 1, 6] <- N_neg[tt, 6] + timestep * 
      ( (1 - fp_neg) * (1 - acc_neg) * pp_neg * lam[tt] * N_neg[tt, 1] +
      + (1 - fp_neg) * (1 - acc_neg) * (er_neg + lr_neg * pp_neg * lam[tt]) * N_neg[tt, 2] +
      - (mu_smneg_neg + hivinc + cure) * N_neg[tt, 6] )

   if(allow.hiv) {
   # 7. S_hiv
   N_hiv[tt + 1, 1] <- N_hiv[tt, 1] + timestep * 
      ( hivinc * N_neg[tt, 1] +
      - (mu_hiv + hivprog + lam[tt]) * N_hiv[tt, 1] )
   # 8. L_hiv
   N_hiv[tt + 1, 2] <- N_hiv[tt, 2] + timestep * 
      ( hivinc * N_neg[tt, 2] + (1 - pp_hiv) * lam[tt] * N_hiv[tt, 1] +
      - (mu_hiv + hivprog + er_hiv + lr_hiv * pp_hiv * lam[tt]) * N_hiv[tt, 2] +
      + dx_pos * tx * sn_smpos * N_hiv[tt, 3] + 
      + dx_neg * tx * sn_smneg * (1-lfu) * N_hiv[tt, 5] )
   # 9. DSP_hiv
   N_hiv[tt + 1, 3] <- N_hiv[tt, 3] + timestep * 
      ( hivinc * N_neg[tt, 3] +
      + fp_hiv * acc_hiv * pp_hiv * lam[tt] * N_hiv[tt, 1] +
      + fp_hiv * acc_hiv * (er_hiv + lr_hiv * pp_hiv * lam[tt]) * N_hiv[tt, 2] +
      - (mu_tb_hiv + hivprog + dx_pos * tx * sn_smpos) * N_hiv[tt, 3] )
   #10. USP_hiv
   N_hiv[tt + 1, 4] <- N_hiv[tt, 4] + timestep * 
      ( hivinc * N_neg[tt, 4] +
      + fp_hiv * (1 - acc_hiv) * pp_hiv * lam[tt] * N_hiv[tt, 1] +
      + fp_hiv * (1 - acc_hiv) * (er_hiv + lr_hiv * pp_hiv * lam[tt]) * N_hiv[tt, 2] +
      - (mu_tb_hiv + hivprog) * N_hiv[tt, 4] )
   #11. DSN_hiv
   N_hiv[tt + 1, 5] <- N_hiv[tt, 5] + timestep * 
      ( hivinc * N_neg[tt, 5] +
      + (1 - fp_hiv) * acc_hiv * pp_hiv * lam[tt] * N_hiv[tt, 1] +
      + (1 - fp_hiv) * acc_hiv * (er_hiv + lr_hiv * pp_hiv * lam[tt]) * N_hiv[tt, 2] +
      - (mu_tb_hiv + hivprog + dx_neg * tx * sn_smneg * (1-lfu)) * N_hiv[tt, 5] )
   #12. USN_hiv
   N_hiv[tt + 1, 6] <- N_hiv[tt, 6] + timestep * 
      ( hivinc * N_neg[tt, 6] +
      + (1 - fp_hiv) * (1 - acc_hiv) * pp_hiv * lam[tt] * N_hiv[tt, 1] +
      + (1 - fp_hiv) * (1 - acc_hiv) * (er_hiv + lr_hiv * pp_hiv * lam[tt]) * N_hiv[tt, 2] +
      - (mu_tb_hiv + hivprog) * N_hiv[tt, 6] )

   #13. S_aids
   N_aids[tt + 1, 1] <- N_aids[tt, 1] + timestep * 
      ( hivprog * N_hiv[tt, 1] +
      - (mu_aids + lam[tt]) * N_aids[tt, 1] )
   #14. L_aids
   N_aids[tt + 1, 2] <- N_aids[tt, 2] + timestep * 
      ( hivprog * N_hiv[tt, 2] +
      + (1 - pp_aids) * lam[tt] * N_aids[tt, 1] +
      - (mu_aids + er_aids + lr_aids * pp_aids * lam[tt]) * N_aids[tt, 2] +
      + dx_pos * tx * sn_smpos * N_aids[tt, 3] +
      + dx_neg * tx * sn_smneg * (1-lfu) * N_aids[tt, 5] )
   #15. DSP_aids
   N_aids[tt + 1, 3] <- N_aids[tt, 3] + timestep * 
      ( hivprog * N_hiv[tt, 3] +
      + fp_aids * acc_aids * pp_aids * lam[tt] * N_aids[tt, 1] +
      + fp_aids * acc_aids * (er_aids + lr_aids * pp_aids * lam[tt]) * N_aids[tt, 2] +
      - (mu_tb_aids + dx_pos * tx * sn_smpos) * N_aids[tt, 3] )
   #16. USP_aids
   N_aids[tt + 1, 4] <- N_aids[tt, 4] + timestep * 
      ( hivprog * N_hiv[tt, 4] +
      + fp_aids * (1 - acc_aids) * pp_aids * lam[tt] * N_aids[tt, 1] +
      + fp_aids * (1 - acc_aids) * (er_aids + lr_aids * pp_aids * lam[tt]) * N_aids[tt, 2] +
      - (mu_tb_aids) * N_aids[tt, 4] )
   #17. DSN_aids
   N_aids[tt + 1, 5] <- N_aids[tt, 5] + timestep * 
      ( hivprog * N_hiv[tt, 5] +
      + (1 - fp_aids) * acc_aids * pp_aids * lam[tt] * N_aids[tt, 1] +
      + (1 - fp_aids) * acc_aids * (er_aids + lr_aids * pp_aids * lam[tt]) * N_aids[tt, 2] +
      - (mu_tb_aids + dx_neg * tx * sn_smneg * (1-lfu)) * N_aids[tt, 5] )
   #18. USN_aids
   N_aids[tt + 1, 6] <- N_aids[tt, 6] + timestep * 
      ( hivprog * N_hiv[tt, 6] +
      + (1 - fp_aids) * (1 - acc_aids) * pp_aids * lam[tt] * N_aids[tt, 1] +
      + (1 - fp_aids) * (1 - acc_aids) * (er_aids + lr_aids * pp_aids * lam[tt]) * N_aids[tt, 2] +
      - (mu_tb_aids) * N_aids[tt, 6] )
   #A. TSP
   TSP[tt + 1] <- sum(N_neg[tt + 1, 3:4], N_hiv[tt + 1, 3:4], N_aids[tt + 1, 3:4]) # {DSP,USP}_{neg,hiv,aids}
   #B. TSN
   TSN[tt + 1] <- sum(N_neg[tt + 1, 5:6], N_hiv[tt + 1, 5:6], N_aids[tt + 1, 5:6]) # {DSN,USN}_{neg,hiv,aids}
   #D. M
   M[tt + 1] <- mu_neg * sum(N_neg[tt + 1, 1:2]) + # mu_neg, {S,L}_neg
      + mu_smpos_neg * sum(N_neg[tt + 1, 3:4]) +   # mu_smpos_neg, {DSP,USP}_neg
      + mu_smneg_neg * sum(N_neg[tt + 1, 5:6]) +   # mu_smneg_neg, {DSN,USN}_neg
      + mu_hiv * sum(N_hiv[tt + 1, 1:2]) +         # mu_hiv, {S,L}_hiv
      + mu_tb_hiv * sum(N_hiv[tt + 1, 3:6]) +      # mu_tb_hiv, {DSP,USP,DSN,USN}_aids
      + mu_aids * sum(N_aids[tt + 1, 1:2]) +       # mu_aids, {S,L}_aids
      + mu_tb_aids * sum(N_aids[tt + 1, 3:6])      # mu_tb_aids, {DSP,USP,DSN,USN}_aids
   } else { # allow.hiv == FALSE
   #A. 
   TSP[tt + 1] <- sum(N_neg[tt + 1, 3:4])
   #B. 
   TSN[tt + 1] <- sum(N_neg[tt + 1, 5:6])
   #D. 
   M[tt + 1] <- mu_neg * sum(N_neg[tt + 1, 1:2]) +
      + mu_smpos_neg * sum(N_neg[tt + 1, 3:4]) +
      + mu_smneg_neg * sum(N_neg[tt + 1, 5:6])
   }
   #C. lam
   lam[tt + 1] <- lam_0 * (TSP[tt + 1] + i_smneg * TSN[tt + 1])
   } # end for(tt in seq_len(ntimes - 1))

   # incidence pp * lam * S + (er + lr * pp * lam) * L -- rearranged
   incid_tb <- lam * ( pp_neg * (N_neg[, 1] + lr_neg * N_neg[, 2]) + 
      pp_hiv * (N_hiv[, 1] + lr_hiv * N_hiv[, 2]) + 
      pp_aids * (N_aids[, 1] + lr_aids * N_aids[, 2]) ) +
      er_neg * N_neg[, 2] + er_hiv * N_hiv[, 2] + er_aids * N_aids[, 2]
   prev_tb <- TSP + TSN
   mort_tb <- (mu_smpos_neg - 1/agespan) * rowSums(N_neg[, 3:4]) +
      (mu_smneg_neg - 1/agespan) * rowSums(N_neg[, 5:6]) +
      (mu_tb_hiv - 1/agespan) * rowSums(N_hiv[, 3:6]) +
      (mu_tb_aids - 1/agespan) * rowSums(N_aids[, 3:6])
   
   compart <- cbind(N_neg, N_hiv, N_aids)
   dimnames(compart)[[2]] <- c(outer(c("S", "L", "DSP", "USP", "DSN", "USN"), 
      c("neg", "hiv", "aids"), paste, sep="_"))
   result <- list(times=times, pop=compart, lam=lam, TSP=TSP, TSN=TSN, M=M, 
      incid_tb=incid_tb, prev_tb=prev_tb, mort_tb=mort_tb)
   attr(result, "call") <- match.call()
   result
}

# wrapper for model that limits to 6 HIV-negative compartments; check against previous results
tbdx_neg <- function(pop=c(48727, 50984, 15, 44, 72, 158), lam_0=1.8195e-4, ...) {
   result <- tbdx(pop=c(pop, rep(0, 12)), lam_0=lam_0, hivinc=0, hivprog=0, allow.hiv=FALSE, ...)
   attr(result, "call") <- match.call()
   result
}

x_neg <- tbdx_neg()
do.call("cbind", x_neg)[c(1:6, 9996:10001), ]

x <- tbdx()
do.call("cbind", x)[c(1:6, 9996:10001), ]

## show how values of lam_0 were computed
uniroot(function(x) tbdx(lam_0=x)$lam[10001]-0.02, c(1e-5, 2e-4), trace=2, tol=1e-10)
# list(root = 0.00013471, f.root = -4.002136e-09, iter = 6L, init.it = NA_integer_, estim.prec = 5e-11)

uniroot(function(x) tbdx_neg(lam_0=x)$lam[10001]-0.02, c(1e-5, 2e-4), trace=2, tol=1e-10)
# list(root = 0.00018195, f.root = -5.526124e-11, iter = 6L, init.it = NA_integer_, estim.prec = 5e-11)

## assess convergence to equilibrium, different lam_0 and starting pop
y <- tbdx(lam_0 = 1e-4, pop=c(51932, 27996, 16, 16, 91, 125, 5799, 3990, 3, 2, 33, 5, 9226, 728, 2, 13, 21, 2))
y_popchange <- apply(y$pop, 2, diff)
y_abschange <- rowSums(abs(y_popchange))
y_maxchange <- apply(abs(y_popchange), 1, max)
plot(y$times[-1], y_abschange, type="l")
plot(y$times[-1], y_abschange, type="l", ylim=c(0,2))
lines(y$times[-1], y_maxchange, col="red")
plot(y$times[-1], y_abschange, type="l", ylim=c(1e-5,10), log="y")
lines(y$times[-1], y_maxchange, col="red")
min(seq_len(10000)[y_abschange<1])    #  750 ~15 yr
min(seq_len(10000)[y_abschange<0.15]) # 1646 ~33 yr
min(seq_len(10000)[y_abschange<0.01]) # 4255 ~85 yr
