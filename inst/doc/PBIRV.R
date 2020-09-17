## ----message=FALSE------------------------------------------------------------
library(PBIR)
library(survival)

## -----------------------------------------------------------------------------
set.seed(100)
n=100
error=rnorm(n)
tr=exp(rnorm(n)+error+0.5)
tp=exp(rnorm(n)+error)
tr[tp<tr]=Inf
tc=runif(n, 3, 8.5)
t2response=pmin(tr, tc)
delta_response=1*(tr<tc)
t2progression=pmin(tp, tc)
delta_progression=1*(tp<tc)


## -----------------------------------------------------------------------------
round(head(cbind(t2response, delta_response=1*(tr<tc), t2progression, delta_progression)), 4)


## -----------------------------------------------------------------------------
fit=PBIR1(t2PROGRESSION=t2progression, STATUS_PROGRESSION=delta_progression, 
          t2RESPONSE=t2response,       STATUS_RESPONSE=delta_response)


## -----------------------------------------------------------------------------

round(fit[95:100, ], 4)


## -----------------------------------------------------------------------------
tt=c(0, fit$time)
PBIR=c(0, fit$PBIR)
ci.up=c(0, fit$ci.up)
ci.low=c(0, fit$ci.low)

B=length(tt)
tt=rep(tt, rep(2, B))[-1]
PBIR=rep(PBIR, rep(2, B))[-(2*B)]
ci.up=rep(ci.up, rep(2, B))[-(2*B)]
ci.low=rep(ci.low, rep(2, B))[-(2*B)]

          
plot(range(tt), range(na.omit(c(ci.low, ci.up))), 
     xlab="time",  ylab="PBIR",  
     main="PBIR", type="n")
lines(tt, PBIR, col=1, lwd=2)
lines(tt, ci.low, col=2)
lines(tt, ci.up, col=2)




## -----------------------------------------------------------------------------
fit=PBIR1(t2PROGRESSION=t2progression, STATUS_PROGRESSION=delta_progression, 
          t2RESPONSE=t2response,       STATUS_RESPONSE=delta_response,
          time=c(2,4,6))

round(fit, 4)

## -----------------------------------------------------------------------------
fit=PBIR1(t2PROGRESSION=t2progression, STATUS_PROGRESSION=delta_progression, 
          t2RESPONSE=t2response,       STATUS_RESPONSE=delta_response,
          time=c(7,10))

round(fit, 4)

## -----------------------------------------------------------------------------
set.seed(100)
n=200
trt=rbinom(n, 1, 0.5)
error=rnorm(n)
tr=exp(rnorm(n)+error-trt*0.5+0.5)
tp=exp(rnorm(n)+error+trt*0.25)
tr[tp<tr]=Inf
tc=runif(n, 3, 8.5)
t2response=pmin(tr, tc)
delta_response=1*(tr<tc)
t2progression=pmin(tp, tc)
delta_progression=1*(tp<tc)



## -----------------------------------------------------------------------------
round(head(cbind(t2response, delta_response=1*(tr<tc), t2progression, delta_progression, trt)), 4)


## -----------------------------------------------------------------------------
fit2=PBIR2(t2PROGRESSION=t2progression, STATUS_PROGRESSION=delta_progression, 
           t2RESPONSE=t2response,       STATUS_RESPONSE=delta_response, 
           TRT=trt)


## -----------------------------------------------------------------------------

round(fit2[215:220, ], 4)


## -----------------------------------------------------------------------------
tt=fit2$time
diff=fit2$diff
low=fit2$ci.low
up=fit2$ci.up
B=length(tt)+1

tt=c(0, tt)
diff=c(0, diff)
low=c(0, low)
up=c(0, up)
tt=rep(tt, rep(2, B))[-1]
diff=rep(diff, rep(2, B))[-(2*B)]
low=rep(low, rep(2, B))[-(2*B)]
up=rep(up, rep(2, B))[-(2*B)]

plot(range(c(tt, 0)), range(c(low, up)), 
     xlab="time", ylab="difference in PBIR", 
     lwd=2, type="n")
lines(tt, diff, lwd=2, col=3)
lines(tt, low,  col=2)
lines(tt, up, col=2)
lines(range(fit$time), rep(0, 2), col=4, lty=4)



## -----------------------------------------------------------------------------
fit2=PBIR2(t2PROGRESSION=t2progression, STATUS_PROGRESSION=delta_progression, 
           t2RESPONSE=t2response,       STATUS_RESPONSE=delta_response, 
           TRT=trt,  time=c(2, 4, 6))

round(fit2, 4)


## -----------------------------------------------------------------------------
fit2=PBIR2(t2PROGRESSION=t2progression, STATUS_PROGRESSION=delta_progression, 
           t2RESPONSE=t2response,       STATUS_RESPONSE=delta_response, 
           TRT=trt,  time=c(2, 4, 6, 10))

round(fit2, 4)

## -----------------------------------------------------------------------------
set.seed(100)
n=200
trt=rbinom(n, 1, 0.5)
error=rnorm(n)
tr=exp(rnorm(n)+error-trt*0.5+0.5)
tp=exp(rnorm(n)+error+trt*0.25)
tr[tp<tr]=Inf
tc=runif(n, 3, 8.5)
t2response=pmin(tr, tc)
delta_response=1*(tr<tc)
t2progression=pmin(tp, tc)
delta_progression=1*(tp<tc)



## -----------------------------------------------------------------------------

fit1=mduration(t2PROGRESSION=t2progression[trt==1],
               STATUS_PROGRESSION=delta_progression[trt==1], 
               t2RESPONSE=t2response[trt==1],     
               STATUS_RESPONSE=delta_response[trt==1])

fit0=mduration(t2PROGRESSION=t2progression[trt==0],
               STATUS_PROGRESSION=delta_progression[trt==0], 
               t2RESPONSE=t2response[trt==0],     
               STATUS_RESPONSE=delta_response[trt==0])

## -----------------------------------------------------------------------------

fit1

fit0


## -----------------------------------------------------------------------------

fit1=mduration(t2PROGRESSION=t2progression[trt==1],
               STATUS_PROGRESSION=delta_progression[trt==1], 
               t2RESPONSE=t2response[trt==1],     
               STATUS_RESPONSE=delta_response[trt==1], 
               time.max=6.75)

fit0=mduration(t2PROGRESSION=t2progression[trt==0],
               STATUS_PROGRESSION=delta_progression[trt==0], 
               t2RESPONSE=t2response[trt==0],     
               STATUS_RESPONSE=delta_response[trt==0], 
               time.max=6.75)

diff=(fit1$meandor.est-fit0$meandor.est)
se=sqrt(fit1$meandor.se^2+fit0$meandor.se^2)
z=diff/se
ci=diff+c(-1,1)*qnorm(0.975)*se
pvalue=1-pchisq(z^2, 1)
result=cbind(diff, ci[1], ci[2], pvalue)
colnames(result)=c("difference in DOR", "CI.low", "CI.high", "p.value")
print(round(result, 4))


## -----------------------------------------------------------------------------
fit=CRR(t2PROGRESSION=t2progression, STATUS_PROGRESSION=delta_progression, 
         t2RESPONSE=t2response,       STATUS_RESPONSE=delta_response, 
         TRT=trt,  time=c(1,2,3,4,5))

fit

## -----------------------------------------------------------------------------

# Estimate the CRR over all the jump points of step functions
fit=CRR(t2PROGRESSION=t2progression,
         STATUS_PROGRESSION=delta_progression,
         t2RESPONSE=t2response,
         STATUS_RESPONSE=delta_response,
         TRT=trt)

# Plot the estimated PBIR by group
tt1=c(0, fit$result1$time)
crr1=c(0, fit$result1$CRR)
B1=length(tt1)
tt1=rep(tt1, rep(2, B1))[-1]
crr1=rep(crr1, rep(2, B1))[-(2*B1)]
tt0=c(0, fit$result0$time)
crr0=c(0, fit$result0$CRR)
B0=length(tt0)
tt0=rep(tt0, rep(2, B0))[-1]
crr0=rep(crr0, rep(2, B0))[-(2*B0)]
plot(range(c(fit$result1$time, fit$result0$time)),
     range(c(fit$result1$CRR, fit$result0$CRR)),
     xlab="time",  ylab="CRR",
     main="black: group 0; red: group 1", type="n")
lines(tt0, crr0, col=1, lwd=2)
lines(tt1, crr1, col=2, lwd=2)

