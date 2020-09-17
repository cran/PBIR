#' Estimate the PBIR curve over a time window
#'
#'@param  t2PROGRESSION  time to progression/death or censoring
#'@param  STATUS_PROGRESSION   binary indicator for progression status: 1 for progression/death; 0 for censoring
#'@param  t2RESPONSE  time to response or censoring
#'@param  STATUS_RESPONSE  binary indicator for response status: 1 for response; 0 for censoring
#'@param  time        user-selected time points at which the PBIR is estimated; the default value is "NULL" and the PBIR will be estimated at all observed time points
#'@param  alpha       coverage level of the point-wise confidence interval for PBIR curve; the default value is 0.95
#'@return a data matrix containing "time", "PBIR estimates", "standard errors of PBIR estimates",  "confidence intervals of the PBIR"
#'@examples
#'
#'library(survival)
#'n=100
#'set.seed(10)
#'
#'# Generate the data
#'
#'trt=rbinom(n, 1, 0.5)
#'error=rnorm(n)
#'tr=exp(rnorm(n)+error-trt*0.5+0.5)
#'tp=exp(rnorm(n)+error+trt*0.25)
#'tr[tp<tr]=Inf
#'tc=runif(n, 3, 8.5)

#'t2response=pmin(tr, tc)
#'delta_response=1*(tr<tc)
#'t2progression=pmin(tp, tc)
#'delta_progression=1*(tp<tc)

#'
#'# Estimate the PBIR in two groups
#'
#'fit1=PBIR1(t2PROGRESSION=t2progression[trt==1],
#'            STATUS_PROGRESSION=delta_progression[trt==1],
#'            t2RESPONSE=t2response[trt==1],
#'            STATUS_RESPONSE=delta_response[trt==1])
#'
#'fit0=PBIR1(t2PROGRESSION=t2progression[trt==0],
#'            STATUS_PROGRESSION=delta_progression[trt==0],
#'            t2RESPONSE=t2response[trt==0],
#'            STATUS_RESPONSE=delta_response[trt==0])
#'
#'
#'# Plot the estimated PBIR by group
#'
#'tt1=c(0, fit1$time)
#'PBIR1=c(0, fit1$PBIR)
#'B1=length(tt1)
#'tt1=rep(tt1, rep(2, B1))[-1]
#'PBIR1=rep(PBIR1, rep(2, B1))[-(2*B1)]

#'tt0=c(0, fit0$time)
#'PBIR0=c(0, fit0$PBIR)
#'B0=length(tt0)
#'tt0=rep(tt0, rep(2, B0))[-1]
#'PBIR0=rep(PBIR0, rep(2, B0))[-(2*B0)]

#'plot(range(c(fit1$time, fit0$time)), range(c(fit1$PBIR, fit0$PBIR)),
#'      xlab="time",  ylab="PBIR",
#'      main="black: group 0; red: group 1", type="n")
#'lines(tt0, PBIR0, col=1)
#'lines(tt1, PBIR1, col=2)

#'@importFrom("survival", survfit)
#'@references Huang, B., Tian, L., Talukder, E., Rothenberg, M., Kim, DY., and Wei, LJ. (2018) Evaluating Treatment Effect Based on Duration of Response for a Comparative Oncology Study. JAMA Oncol, doi: 10.1001/jamaoncol.2018.0275
#'@references Huang, B., Tian, L., McCaw, Z., Luo, Talukder, E., X., Rothenberg, M., Xie, W., Choueiri, T., Kim, DY., & Wei, LJ. (2020). Analysis of Response Data for Assessing Treatment Effects in Comparative Clinical Studies. Ann Intern Med, doi: 10.7326/M20-0104.
#'@export



PBIR1=function(t2PROGRESSION, STATUS_PROGRESSION, t2RESPONSE, STATUS_RESPONSE, time=NULL, alpha=0.95){

  t2RESPONSE[STATUS_RESPONSE==0]=Inf

  y1=pmin(t2PROGRESSION, t2RESPONSE)
  d1=1*(STATUS_RESPONSE+STATUS_PROGRESSION>0)

  y2=t2PROGRESSION
  d2=STATUS_PROGRESSION

  fit1=survfit(Surv(y1, d1)~1)
  n1=length(fit1$time)

  fit2=survfit(Surv(y2, d2)~1)
  n2=length(fit2$time)

  #plot(fit2)
  #lines(fit1, col=2)

  tau.grd=sort(unique(c(y1, y2)))


  if(min(d1[y1==max(y1)])==0)  taumax=min(max(y1), max(y2))
  if(min(d1[y1==max(y1)])==1)  taumax=max(y2)


  tau.grd=c(tau.grd[tau.grd<taumax], taumax)
  m.tau.grd=length(tau.grd)
  n=length(t2PROGRESSION)

  surv1.tot=surv2.tot=rep(NA, m.tau.grd)
  atrisk1.tot=atrisk2.tot=rep(NA, m.tau.grd)
  hazard1.tot=hazard2.tot=rep(NA, m.tau.grd)

  for(i in 1:m.tau.grd)
  {t0=tau.grd[i]

   if(t0>=max(y1))  t0new=max(y1)
   if(t0<=max(y1)) t0new=t0

   id1=max((1:(n1+1))[c(0, fit1$time)<=t0new])
   id2=max((1:(n2+1))[c(0, fit2$time)<=t0])

   surv1.tot[i]=c(1, fit1$surv)[id1]
   surv2.tot[i]=c(1, fit2$surv)[id2]

   atrisk1.tot[i]=sum(y1>=t0new)
   atrisk2.tot[i]=sum(y2>=t0)

   hazard1.tot[i]=sum((y1==t0new)*d1)/sum(y1>=t0new)
   hazard2.tot[i]=sum((y2==t0)*d2)/sum(y2>=t0)
   }
  dsurv.tot=surv2.tot-surv1.tot


  dsd.tot=rep(NA, m.tau.grd)
  for(i in 1:m.tau.grd)
  {t0=tau.grd[i]

   if(t0>max(y1))  t0new=max(y1)
   if(t0<=max(y1)) t0new=t0

   tau1=tau2=rep(0, n)
   for(j in 1:n)
     tau1[j]=(y1[j]<=t0new)*d1[j]/sum(y1>=y1[j])-sum((hazard1.tot/atrisk1.tot)[tau.grd<=min(t0new,y1[j])])
   for(j in 1:n)
     tau2[j]=(y2[j]<=t0)*d2[j]/sum(y2>=y2[j])-sum((hazard2.tot/atrisk2.tot)[tau.grd<=min(t0,y2[j])])
   tau1=surv1.tot[i]*tau1
   tau2=surv2.tot[i]*tau2
   dsd.tot[i]=sqrt(sum((tau1-tau2)^2))
  }

  logd.tot=log(dsurv.tot/(1-dsurv.tot))
  logsd.tot=dsd.tot/(dsurv.tot*(1-dsurv.tot))

  cilogd=cbind(logd.tot-qnorm((1+alpha)/2)*logsd.tot, logd.tot+qnorm((1+alpha)/2)*logsd.tot)
  cid=exp(cilogd)/(1+exp(cilogd))
  cid[dsurv.tot==0,]=0
  cid[dsurv.tot==1,]=1
  ci1=cid[,1]
  ci2=cid[,2]


  if(length(time)>0)
  {if(max(time)>taumax)
     {message("The PBIR is not identifiable at some selected time points,
which are replaced by the maximum time point, where it is identifiable.")
      time=time[time<taumax]
      time=c(time, taumax)
      }
   dsurv.tot=approx(x=tau.grd, y=dsurv.tot, xout=time, method="constant", yleft=0, yright=0)$y
   dsd.tot=approx(x=tau.grd, y=dsd.tot, xout=time, method="constant", yleft=0, yright=0)$y
   ci1=approx(x=tau.grd, y=ci1, xout=time, method="constant", yleft=0, yright=0)$y
   ci2=approx(x=tau.grd, y=ci2, xout=time, method="constant", yleft=0, yright=0)$y
   tau.grd=time
   }

  res.tot=cbind(tau.grd, dsurv.tot, dsd.tot, ci1, ci2)
  colnames(res.tot)=c("time", "PBIR", "std", "ci-low", "ci-up")

  res.tot=data.frame(res.tot)

  return(res.tot)
}
