pop <- 100000
timestep <- 0.02

N <- pop 



timespan <- seq(from=0, to=200, by=step)


#### Parameters
## TB transmission

lam_0 <- 1.9e-04          #annual risk of TB infection, per smear-positive active case
i_smneg  <- 0.22      #relative infectiousness of smear-negative cases
lat <- 1-(0.98)^15           #prevalence of latent TB in 15-year-old recruits
fp_neg <- 0 # 0.45     #Fraction of incident TB that is smear-positive
fp_hiv <- 0.45 	#.45
fp_aids <- 0     #0.35


###TB progression

pp_neg <- 0 #0.07     #proportion of TB infections causing immediate progression to active TB
pp_hiv <- 0.07        #0.07
pp_aids <- 0       #0.56

lr_neg <- 0  #0.28       #proportion of TB reinfections capable of primary progression
lr_hiv <- 0.28 #.28
lr_aids <- 0 #0.75

er_neg <- #0.1/100         #endogenous reactivation rate, per year
er_hiv <- 0.001
er_aids <- 0 #0.079

cure <- 0             #rate of spontaneous cure in HIV-negative patients from active to latent TB



####HIV transmission and progression #######################################

hivinc <- 16.8/100 # 0       #annual risk of HIV infection (HIV incidence)
hivprog  <- 0    #mean rate of progression from HIV to AIDS

####TB diagnosis and treatment

acc <- 0.734        #fraction of TB patients presenting to clinics
tx  <- 0.82         #success of treatment, once TB is diagnosed
dx   <- 4.0        #rate of TB detection after symptoms develop
sn_smpos   <- 0.98    #0.8   #sensitivity of diagnosis, smear-positive
sn_smneg  <- 0.5       #0.25   #sensitivity of diagnosis, smear-negative
lfu <- 0         #proportion of detected patients lost to follow-up before treatment

####Mortality and aging



mu_neg  <- 0.9/100+1/35      #rate of death or population exit, no active TB
mu_hiv  <- 0      
mu_aids <- 0

mu_smpos_neg <- 35/100+1/35   #rate of death or population exit, active TB, stratified by smear status among HIV-negatives
mu_smneg_neg <- 10/100+1/35
mu_tb_hiv <- 0
mu_tb_aids <- 0




##################



DSP_neg.init <- pop/24  #15.79977   smear-positive patients with access to TB diagnostic sevices
USP_neg.init <- pop/24 #45.41314    #smear-positive patients without access to TB diagnostic services
DSP_hiv.init <- 0
USP_hiv.init <- 0
DSP_aids.init <- 0
USP_aids.init <- 0

DSN_neg.init <-  pop/24  #61.12574     #smear-negative patients with access to TB diagnostic sevices
USN_neg.init <- pop/24  #163.4309      #smear-negative patients without access to TB diagnostic services
DSN_hiv.init <- 0
USN_hiv.init <- 0
DSN_aids.init <- 0
USN_aids.init <- 0


S_neg.init <- pop/2 #47881.8      #susceptible
S_hiv.init <- 0
S_aids.init <- 0
L_neg.init <- pop/3  #51832.43     #latently-infected
L_hiv.init <- 0
L_aids.init <- 0


run.time <- 200
end.step <- run.time/timestep +1


DSP_neg <- USP_neg <- DSP_hiv <- USP_hiv <- DSP_aids <- USP_aids <- numeric(end.step)
DSP_neg[1] <- DSP_neg.init 
USP_neg[1] <- USP_neg.init 
DSP_hiv[1] <- DSP_hiv.init 
USP_hiv[1] <- USP_hiv.init 
DSP_aids[1] <- DSP_aids.init 
USP_aids[1] <- USP_aids.init



DSN_neg <- USN_neg <- DSN_hiv <- USN_hiv <- DSN_aids <- USN_aids <- numeric(end.step)
DSN_neg[1] <- DSN_neg.init 
USN_neg[1] <- USN_neg.init 
DSN_hiv[1] <- DSN_hiv.init 
USN_hiv[1] <- USN_hiv.init 
DSN_aids[1] <- DSN_aids.init 
USN_aids[1] <- USN_aids.init



S_neg <- L_neg <- S_hiv <- L_hiv <- S_aids <- L_aids <- numeric(end.step)
S_neg[1] <- S_neg.init 
L_neg[1] <- L_neg.init 
S_hiv[1] <- S_hiv.init 
L_hiv[1] <- L_hiv.init 
S_aids[1] <- S_aids.init 
L_aids[1] <- L_aids.init



#########################

##############################

#Rate of change in compartment populations



TSP <- TSN <- lam <- M <- numeric(end.step)
#TSP[1] <- 0
#TSN[1] <- 0
#lam[1] <- lam_0*(TSP[1]+i_smneg*TSN[1])
#M[1] <- mu_neg*(S_neg[1]+L_neg[1])+mu_hiv*(S_hiv[1]+L_hiv[1])+mu_aids*(S_aids[1]+L_aids[1])+
#		mu_smpos_neg*(DSP_neg[1]+USP_neg[1])+mu_smneg_neg*(DSN_neg[1]+USN_neg[1])+mu_tb_hiv*(DSP_hiv[1]+
#		USP_hiv[1]+DSN_hiv[1]+USN_hiv[1])+mu_tb_aids*(DSP_aids[1]+USP_aids[1]+DSN_aids[1]+USN_aids[1])

for(j in 1:(end.step-1)){



#A
#total number of smear-positive cases at time t
	TSP[j] <- DSP_neg[j]+USP_neg[j]+DSP_hiv[j]+DSP_aids[j]+USP_hiv[j]+USP_aids[j]

#B
#total number of smear-negative cases at time t
	TSN[j] <- DSN_neg[j]+USN_neg[j]+DSN_hiv[j]+DSN_aids[j]+USN_hiv[j]+USN_aids[j]


#C
#the force of infection (annual rate of TB infection) at time t
	lam[j] <- lam_0*(TSP[j]+i_smneg*TSN[j])

#D
#total mortality plus aging at time t
	M[j]= mu_neg*(S_neg[j]+L_neg[j])+mu_hiv*(S_hiv[j]+L_hiv[j])+mu_aids*(S_aids[j]+L_aids[j])+
		mu_smpos_neg*(DSP_neg[j]+USP_neg[j])+mu_smneg_neg*(DSN_neg[j]+USN_neg[j])+mu_tb_hiv*(DSP_hiv[j]+
		USP_hiv[j]+DSN_hiv[j]+USN_hiv[j])+mu_tb_aids*(DSP_aids[j]+USP_aids[j]+DSN_aids[j]+USN_aids[j])



#1
#individuals who are HIV-negative and never previously infected with TB

	S_neg[j+1] <- S_neg[j] + timestep*((1-lat)*M[j] - (lam[j] + hivinc + mu_neg)*S_neg[j])

#2
	L_neg[j+1] <- L_neg[j] + timestep*(lat*M[j]+lam[j]*S_neg[j]*(1-pp_neg)-(er_neg +hivinc +mu_neg + lr_neg*pp_neg*lam[j])*L_neg[j]+(dx*tx*sn_smpos +cure)*DSP_neg[j]+cure*USP_neg[j]+(dx*tx*sn_smneg +cure)*DSN_neg[j] + cure*USN_neg[j])
 

#3
	DSP_neg[j+1] <- DSP_neg[j] + timestep*((fp_neg*acc*pp_neg*lam[j])*S_neg[j]+(fp_neg*acc*(er_neg +lr_neg*pp_neg*lam[j]))*L_neg[j]-(hivinc +mu_smpos_neg + dx*tx*sn_smpos + cure)*DSP_neg[j])


#4
	USP_neg[j+1] <- USP_neg[j] + timestep*((fp_neg*(1-acc)*pp_neg*lam[j])*S_neg[j]+(fp_neg*(1-acc)*(er_neg +lr_neg*pp_neg*lam[j]))*L_neg[j]-(hivinc +mu_smpos_neg + cure)*USP_neg[j])

#5
	DSN_neg[j+1] <- DSN_neg[j] + timestep*((1-fp_neg)*acc*pp_neg*lam[j]*S_neg[j]+(1-fp_neg)*acc*(er_neg +lr_neg*pp_neg*lam[j])*L_neg[j]-(hivinc + mu_smneg_neg + dx*tx*sn_smneg+cure)*DSN_neg[j])


#6
	USN_neg[j+1] <- USN_neg[j] + timestep*((1-fp_neg)*(1-acc)*pp_neg*lam[j]*S_neg[j]+(1-fp_neg)*(1-acc)*(er_neg +lr_neg*pp_neg*lam[j])*L_neg[j]-(hivinc+mu_smneg_neg+cure)*USN_neg[j])



}



#A
#total number of smear-positive cases at time t
	TSP[end.step] <- DSP_neg[end.step]+USP_neg[end.step]+DSP_hiv[end.step]+DSP_aids[end.step]+USP_hiv[end.step]+USP_aids[end.step]

#B
#total number of smear-negative cases at time t
	TSN[end.step] <- DSN_neg[end.step]+USN_neg[end.step]+DSN_hiv[end.step]+DSN_aids[end.step]+USN_hiv[end.step]+USN_aids[end.step]


#C
#the force of infection (annual rate of TB infection) at time t
	lam[end.step] <- lam_0*(TSP[end.step]+i_smneg*TSN[end.step])

#D
#total mortality plus aging at time t
	M[j]= mu_neg*(S_neg[end.step]+L_neg[end.step])+mu_hiv*(S_hiv[end.step]+L_hiv[end.step])+mu_aids*(S_aids[end.step]+L_aids[end.step])+
		mu_smpos_neg*(DSP_neg[end.step]+USP_neg[end.step])+mu_smneg_neg*(DSN_neg[end.step]+USN_neg[end.step])+mu_tb_hiv*(DSP_hiv[end.step]+
		USP_hiv[end.step]+DSN_hiv[end.step]+USN_hiv[end.step])+mu_tb_aids*(DSP_aids[end.step]+USP_aids[end.step]+DSN_aids[end.step]+USN_aids[end.step])



#values of each compartments at the end
lam[end.step]
S_neg[end.step]
L_neg[end.step]
USN_neg[end.step]
USP_neg[end.step]
DSP_neg[end.step]
DSN_neg[end.step]


#plotting
plot(lam, type="l", col="blue")
plot(S_neg, type="l", col = "blue", ylim=c(0,100000))
plot(L_neg, type="l", col = "blue", ylim=c(0,100000))
plot(S_neg, type="l", col = "blue", ylim=c(0,100000))
plot(DSN_neg, type="l", col = "blue")
plot(DSP_neg, type="l", col = "blue")
plot(USP_neg, type="l", col = "blue")
plot(USN_neg, type="l", col = "blue")


#prevalence
prevalence <- TSN + TSP

#incidences from susceptible, latent, reinfections capable of primary progression
inc1 <- (pp_neg*lam)*S_neg
inc2 <- er_neg*L_neg
inc3 <- lr_neg*pp_neg*lam*L_neg
incidence <- inc1+ inc2 + inc3


#values of prevalence and incidence at the end
prevalence[end.step]
inc1[end.step]
inc2[end.step]
inc3[end.step]
incidence[end.step]

#plotting
plot(prevalence, type="l", col = "blue" )
plot(incidence, type="l", col = "blue")

plot(inc1, type="l", col = "blue")
lines(inc2, type="l", lty=3, lwd=2, col = "red")
lines(inc3, type="l", lty=3, lwd=2, col = "green")

#plotting zooming in - looking at the equilibrium values
plot(inc1, type="l", col = "blue", ylim=c(0,100))
lines(inc2, type="l", lty=3, lwd=2, ylim=c(0,100), col = "red")
lines(inc3, type="l", lty=3, lwd=2, ylim=c(0,100), col = "green")






