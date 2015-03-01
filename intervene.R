## follows on from file tbdx-model-function-v2.R   
## source("tbdx-model-function-v2.R")  

## hack code to produce Table 2 
yr.st <- 200	## probably 60 would do, haven't checked 
#standard
x <- tbdx(year_stop=yr.st)
int_ipm <- do.call("cbind", x)[yr.st/0.02+1,c("incid_tb","prev_tb","mort_tb")]
# rapid molecular testing
x <- tbdx(year_stop=yr.st,sn_smpos = 0.98, sn_smneg = 0.5)
int_ipm <- rbind(int_ipm, do.call("cbind", x)[yr.st/0.02+1,c("incid_tb","prev_tb","mort_tb")])
# TB culture
x <- tbdx(year_stop=yr.st,dx_neg = 3, sn_smpos = 0.98, sn_smneg = 0.85, lfu = 0.2)
int_ipm <- rbind(year_stop=100,int_ipm, do.call("cbind", x)[yr.st/0.02+1,c("incid_tb","prev_tb","mort_tb")])
# Community wide ACF
x <- tbdx(year_stop=yr.st,acc_neg=0.823, acc_hiv=0.823, acc_aids=0.823)
int_ipm <- rbind(int_ipm, do.call("cbind", x)[yr.st/0.02+1,c("incid_tb","prev_tb","mort_tb")])
# HIV-targeted ACF
x <- tbdx(year_stop=yr.st,acc_hiv=0.911, acc_aids=0.911)
int_ipm <- rbind(int_ipm, do.call("cbind", x)[yr.st/0.02+1,c("incid_tb","prev_tb","mort_tb")])
# HAART
x <- tbdx(year_stop=yr.st,hivprog=0.5/6.9)
int_ipm <- rbind(int_ipm, do.call("cbind", x)[yr.st/0.02+1,c("incid_tb","prev_tb","mort_tb")])


int_ipm <- int_ipm[-1,]	## 1st row is year_stop
rownames(int_ipm) <- c("std", "RMT", "TBc", "Comm ACF", "HIV ACF", "HAART") 
reduce <- round(100*t(1-t(int_ipm[-1,])/int_ipm[1,]),1)

round(int_ipm)
round(100*t(1-t(int_ipm[-1,])/int_ipm[1,]),1)

