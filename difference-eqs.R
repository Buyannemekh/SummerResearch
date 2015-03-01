## graphs produced here were used in file Hill-talk.pdf 
## for PE Fellow Orientation 2013.12.05 
## Modeling Infectious Diseases
## Simulation, Uncertainly and Estimaiton of Parameters

## Reference: 
## Vynnycky & White text ## p.82, 1918 flu in Cumberland, MD
## model 4.2 online http://anintroductiontoinfectiousdiseasemodelling.com/  

 
## Data source: 
## Frost WH, Sydenstricker E.
## Influenza in Maryland. Preliminary statistics of certain localities.
## Public Health Reports 1919; 34(11): 491 - 504
## http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1996898/pdf/pubhealthreporig02712-0001.pdf
##weekly cases

## graphs produced here were used in file Hill-talk.pdf 
## for PE Fellow Orientation 2013.12.05 
## 

##set.seed(18400206)

windows(record=TRUE)

## influenza data for Cumberland, MD 1918.09.07 - 1918.11.30 
## weekly cases, week dates, total population  
cases <- c(6+3+0+0+0+2,4+2+3+4+2+4+0,4+8+6+4+3+11+5,9+31+19+24+23+24+47
,26+32+138+57+67+91+98,103+94+99+65+117+83+88,47+66+135+34+29+44+24
,42+12+18+26+12+19+11,12+14+7+4+3+25+3,6+5+2+3+10+6+4,7+6+4+2+1+3+2
,1+6+4+3+4+0+2,1+0+1+2+4+1+2)

weeks <- c("09/07","09/14","09/21","09/28","10/05","10/12","10/19"
,"10/26","11/02","11/09","11/16","11/23","11/30")

pop <- 5234

## graph these data (slide 5)
par(mar=c(5, 4, 1, 2) + 0.1) 
plot(1:length(cases),cases,type="b",xlab="week ending in 1918"
,xaxt="n",yaxt="n",bty="n",ylim=c(0,700),pch=16,lwd=2,ylab="incident cases/week"
     #,mar=c(5, 4, 1, 2) + 0.1
     )
axis(side=1,at=1:length(cases),labels=FALSE)
## slant the dates
text(x=1:length(cases),y=-70
,srt=45,adj=0.9,labels=weeks
,xpd=TRUE,cex=0.9)
axis(side=2,at=seq(0,700,by=100)
,labels=seq(0,700,by=100),las=1)


## paramters for SEIR differecne equations (slides 6-9) 
weeks <- c("08/31",weeks)  
## SEIR difference equations
# time units are days
step <- 0.5

N <- pop 
latent.pd <- 2
infcts.pd <- 2
R0 <- 2.1

## initial conditions for compartment sizes 
R.init <- round(0.3*N)
I.init <- 2
E.init <- 0
S.init <- N-E.init-I.init-R.init

## reporting fraction 
prop.rpt <- 0.8

## calculate compartment transfer probabilities 
## for the given time step 
# progression probability for units in 1/2 days
f <- step*(1/latent.pd)

# recovery probability for units in 1/2 days
r <- step*(1/infcts.pd)

# ECR (effectve contact rate) from R0 = beta/r
beta <- R0*r

## start things one week before at 08/31 
## this is 13 weeks from 08/31 until 11/30 

## specify how long to run the model 
timespan <- seq(from=0,to=13*7,by=step)

## set up S, E, I vectors and intialize
S <- E <- I <- numeric(length(timespan))
S[1] <- S.init
E[1] <- E.init
I[1] <- I.init
#R[1] <- R.init

## loop for the difference equations 
## updating iteratively 
for (j in 1:(length(timespan)-1)){
S[j+1] <- S[j] - beta*S[j]*I[j]/N
E[j+1] <- E[j] + beta*S[j]*I[j]/N - f*E[j]
I[j+1] <- I[j] + f*E[j] - r*I[j]
}
R <- N-S-E-I

## calculate the incidence (new cases) 
incid <- f*E  
report <- prop.rpt*incid #incid in time increments of 0.5 day
weekly.report <- numeric(14) ## entries 2 to 14 will be compared with the data 
for (i in 1:13) weekly.report[i+1] <- sum(report[(14*(i-1)+2):(14*i+1)])  

## graph data and model output (slide 10)     
par(mar=c(5, 4, 1, 2) + 0.1) 
plot(1:length(cases),cases,type="b",xlab="week ending in 1918"
     ,xaxt="n",yaxt="n",bty="n",xlim=c(0,length(cases)),ylim=c(0,700),pch=16,lwd=2
     ,ylab="incident cases/week"
     )
axis(side=1,at=0:length(cases),labels=FALSE)
## slant the dates
text(x=0:length(cases),y=-70
     ,srt=45,adj=0.9,labels=weeks
     ,xpd=TRUE,cex=0.9)
axis(side=2,at=seq(0,700,by=100)
     ,labels=seq(0,700,by=100),las=1)

lines(0:length(cases),weekly.report,type="b",lwd=2,pch=16,col="firebrick2")

legend("topright",c("observed","model"),cex=0.8,lty=c(1,1),lwd=c(2,2)
       ,pch=c(1,1),col=c("black","firebrick2"),bty="n")



#### different inputs
R0 <- 2.3 
# ECR from R0 = beta/r
beta <- R0*r

for (j in 1:(length(timespan)-1)){
  S[j+1] <- S[j] - beta*S[j]*I[j]/N
  E[j+1] <- E[j] + beta*S[j]*I[j]/N - f*E[j]
  I[j+1] <- I[j] + f*E[j] - r*I[j]
}
R <- N-S-E-I

incid <- f*E  
report <- prop.rpt*incid #incid in time increments of 0.5 day
weekly.report2 <- numeric(14) ## entries 2 to 14 will be compared with the data 
for (i in 1:13) weekly.report2[i+1] <- sum(report[(14*(i-1)+2):(14*i+1)]) 

### another 
R0 <- 2.1 
# ECR from R0 = beta/r
beta <- R0*r

latent.pd <- 3 
# progression probability for units in 1/2 days
f <- step*(1/latent.pd)

for (j in 1:(length(timespan)-1)){
  S[j+1] <- S[j] - beta*S[j]*I[j]/N
  E[j+1] <- E[j] + beta*S[j]*I[j]/N - f*E[j]
  I[j+1] <- I[j] + f*E[j] - r*I[j]
}
R <- N-S-E-I

incid <- f*E  
report <- prop.rpt*incid #incid in time increments of 0.5 day
weekly.report3 <- numeric(14) ## entries 2 to 14 will be compared with the data 
for (i in 1:13) weekly.report3[i+1] <- sum(report[(14*(i-1)+2):(14*i+1)]) 


## graph these together (slide 11) 
par(mar=c(5, 4, 1, 2) + 0.1) 
plot(1:length(cases),cases,type="b",xlab="week ending in 1918"
     ,xaxt="n",yaxt="n",bty="n",xlim=c(0,length(cases)),ylim=c(0,700),pch=16,lwd=2
     ,ylab="incident cases/week"
)
axis(side=1,at=0:length(cases),labels=FALSE)
## slant the dates
text(x=0:length(cases),y=-70
     ,srt=45,adj=0.9,labels=weeks
     ,xpd=TRUE,cex=0.9)
axis(side=2,at=seq(0,700,by=100)
     ,labels=seq(0,700,by=100),las=1)

lines(0:length(cases),weekly.report,type="b",lwd=2,pch=16,col="firebrick2")
lines(0:length(cases),weekly.report2,type="b",lwd=2,pch=16,col="olivedrab")
lines(0:length(cases),weekly.report3,type="b",lwd=2,pch=16,col="royalblue1")

legend("topright",c("observed","previous model",expression(R[0]==2.3),"latent period 3 days"),cex=0.8,lty=rep(1,4)
       ,lwd=rep(2,4),pch=rep(1,4),col=c("black","firebrick2","olivedrab","royalblue1"),bty="n")



################### fitting to observed data 
## back to original latent period 
latent.pd <- 2 
# progression probability for units in 1/2 days
f <- step*(1/latent.pd) 


## fitting R0 to weekly reported cases

## set up a function for the least squares difference 
## between model output and observed data 
## equivalent to normal errors (slide 13)
output.lsq <- function(x){
  # ECR from R0 = beta/r
  beta <- x*r 
  S <- E <- I <- numeric(length(timespan)) 
  S[1] <- S.init
  E[1] <- E.init
  I[1] <- I.init
  #R[1] <- R.init
  
  for (j in 1:(length(timespan)-1)){
    S[j+1] <- S[j] - beta*S[j]*I[j]/N
    E[j+1] <- E[j] + beta*S[j]*I[j]/N - f*E[j]
    I[j+1] <- I[j] + f*E[j] - r*I[j] 
  }
  R <- N-S-E-I
  
  incid <- f*E 
  report <- prop.rpt*incid
  weekly.report <- numeric(14)
  for (i in 1:13) weekly.report[i+1] <- sum(report[(14*(i-1)+2):(14*i+1)]) 
  LSQ <- sum((weekly.report[-1]-cases)^2) ## least squares 
  return(LSQ) 
  #NLL <- -sum(dpois(cases,lambda=weekly.report[-1],log = TRUE))
  #return(NLL) 
}## end of output.lsq 

 
## run the LSQ optimizer ( = normal errors MLE) 
lsq <- optimize(f=output.lsq,interval=c(1.5,3.5))
R0.lsq <- lsq$minimum
## run the differene eq model inputting this value 
beta <- R0.lsq*r 
S <- E <- I <- numeric(length(timespan)) 
S[1] <- S.init
E[1] <- E.init
I[1] <- I.init
#R[1] <- R.init

for (j in 1:(length(timespan)-1)){
  S[j+1] <- S[j] - beta*S[j]*I[j]/N
  E[j+1] <- E[j] + beta*S[j]*I[j]/N - f*E[j]
  I[j+1] <- I[j] + f*E[j] - r*I[j] 
}
R <- N-S-E-I

incid <- f*E
report <- prop.rpt*incid
weekly.report.lsq <- numeric(14)
for (i in 1:13) weekly.report.lsq[i+1] <- sum(report[(14*(i-1)+2):(14*i+1)])




## set up a function for negative Poisson log likelihood 
## based on observed data (slide 13) 
output.poi <- function(x){
  # ECR from R0 = beta/r
  beta <- x*r 
  S <- E <- I <- numeric(length(timespan)) 
  S[1] <- S.init
  E[1] <- E.init
  I[1] <- I.init
  #R[1] <- R.init
  
  for (j in 1:(length(timespan)-1)){
    S[j+1] <- S[j] - beta*S[j]*I[j]/N
    E[j+1] <- E[j] + beta*S[j]*I[j]/N - f*E[j]
    I[j+1] <- I[j] + f*E[j] - r*I[j] 
  }
  R <- N-S-E-I
  
  incid <- f*E 
  report <- prop.rpt*incid
  weekly.report <- numeric(14)
  for (i in 1:13) weekly.report[i+1] <- sum(report[(14*(i-1)+2):(14*i+1)]) 
  NLL <- -sum(dpois(cases,lambda=weekly.report[-1],log = TRUE))
  return(NLL) 
}#end of output.poi 


## solve for the Poisson MLE
poi <- optimize(f=output.poi,interval=c(1.5,3.5))
R0.poi <- poi$minimum
## run the differene eq model inputting this value  
beta <- R0.poi*r 
S <- E <- I <- numeric(length(timespan)) 
S[1] <- S.init
E[1] <- E.init
I[1] <- I.init
#R[1] <- R.init

for (j in 1:(length(timespan)-1)){
  S[j+1] <- S[j] - beta*S[j]*I[j]/N
  E[j+1] <- E[j] + beta*S[j]*I[j]/N - f*E[j]
  I[j+1] <- I[j] + f*E[j] - r*I[j] 
}
R <- N-S-E-I

incid <- f*E
report <- prop.rpt*incid
weekly.report.poi <- numeric(14)
for (i in 1:13) weekly.report.poi[i+1] <- sum(report[(14*(i-1)+2):(14*i+1)]) 


 

## graph data and the 2 dufferent fits (slide 14) 
par(mar =  c(5, 4, 4, 2) + 0.1)
plot(1:length(cases),cases,type="b",xlab="week ending in 1918"
     ,xaxt="n",yaxt="n",bty="n",xlim=c(0,length(cases)),ylim=c(0,700),pch=16,lwd=2
     ,ylab="incident cases/week"
)
axis(side=1,at=0:length(cases),labels=FALSE)
## slant the dates
text(x=0:length(cases),y=-70
     ,srt=45,adj=0.9,labels=weeks
     ,xpd=TRUE,cex=0.9)
axis(side=2,at=seq(0,700,by=100)
     ,labels=seq(0,700,by=100),las=1)

lines(0:length(cases),weekly.report.lsq,type="b",lwd=2,pch=16,col="olivedrab")
lines(0:length(cases),weekly.report.poi,type="b",lwd=2,pch=16,col="royalblue1")

legend("topright",c("observed","LSQ fit","Poisson fit"),cex=0.8,lty=c(1,1,1)
       ,lwd=c(2,2,2),pch=c(1,1,1),col=c("black","olivedrab","royalblue1"),bty="n")
title(main=bquote(paste("LSQ: ",R[0]==.(round(R0.lsq,2)),"; Poisson: ",R[0]==.(round(R0.poi,2))))) 


####################################
## using built in mle function in R
## starting estimate R0 = 2
model.lsq <- mle(output.lsq, start=list(x=2)) 
summary(model.lsq) 

model.poi <- mle(output.poi, start=list(x=2)) 
summary(model.poi)