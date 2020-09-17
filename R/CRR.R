#' Estimate cumulative response rates (CRR) and test their equality between two groups
#'
#'@param  t2PROGRESSION  time to progression/death or censoring
#'@param  STATUS_PROGRESSION   binary indicator for progression status: 1 for progression/death; 0 for censoring
#'@param  t2RESPONSE  time to response or censoring
#'@param  STATUS_RESPONSE  binary indicator for response status: 1 for response; 0 for censoring
#'@param  TRT   binary indicator for treatment assignment: 1 for treatment arm and 0 for control arm
#'@param  time        user-selected time points at which the cumulative response rate is to be estimated; the default value is "NULL" and the cumulative response rate will be estimated at all observed time points
#'@param  alpha       coverage level of the point-wise confidence interval for the cumulative response rate; the default value is 0.95
#'@return A list with following elements
#'\itemize{
#'\item{result0:}   a data matrix containing "time", "CRR estimates (group 0)", "standard error of CRR estimates (group 0)",  "confidence interval of CRR (group 0)"
#'\item{result1:}   a data matrix containing "time", "CRR estimates (group 1)", "standard error of CRR estimates (group 1)",  "confidence interval of CRR (group 1)"
#'\item{pvalue:}    the p-value from two group comparison
#'}
#'@examples
#'
#'library(cmprsk)
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
#'fit=CRR(t2PROGRESSION=t2progression,
#'          STATUS_PROGRESSION=delta_progression,
#'          t2RESPONSE=t2response,
#'          STATUS_RESPONSE=delta_response,
#'          TRT=trt)
#'
#'fit
#'
#'# Plot the estimated PBIR by group
#'
#'tt1=c(0, fit$result1$time)
#'CRR1=c(0, fit$result1$CRR)
#'B1=length(tt1)
#'tt1=rep(tt1, rep(2, B1))[-1]
#'CRR1=rep(CRR1, rep(2, B1))[-(2*B1)]

#'tt0=c(0, fit$result0$time)
#'CRR0=c(0, fit$result0$CRR)
#'B0=length(tt0)
#'tt0=rep(tt0, rep(2, B0))[-1]
#'CRR0=rep(CRR0, rep(2, B0))[-(2*B0)]

#'plot(range(c(fit$result1$time, fit$result0$time)),
#'      range(c(fit$result1$CRR, fit$result0$CRR)),
#'      xlab="time",  ylab="CRR",
#'      main="black: group 0; red: group 1", type="n")
#'lines(tt0, CRR0, col=1)
#'lines(tt1, CRR1, col=2)

#'@importFrom("cmprsk", cuminc)
#'@references Gray, RJ. (1988) A class of K-sample tests for comparing the cumulative incidence of a competing risk, ANNALS OF STATISTICS, 16:1141-1154.
#'@references Aalen, O. (1978) Nonparametric estimation of partial transition probabilities in multiple decrement models, ANNALS OF STATISTICS, 6:534-545.
#'@export



CRR=function(t2PROGRESSION, STATUS_PROGRESSION, t2RESPONSE, STATUS_RESPONSE, TRT, time=NULL, alpha=0.95){

  t2RESPONSE[STATUS_RESPONSE==0]=Inf

  y1=pmin(t2PROGRESSION, t2RESPONSE)
  d1=rep(0, length(y1))
  d1[STATUS_RESPONSE>0 & y1==t2RESPONSE]=1
  d1[STATUS_PROGRESSION>0 & y1==t2PROGRESSION]=2

  fit=cuminc(ftime=y1, fstatus=d1,  group=TRT)

  time0=fit$"0 1"$time
  crr0=fit$"0 1"$est
  var0=fit$"0 1"$var

  time1=fit$"1 1"$time
  crr1=fit$"1 1"$est
  var1=fit$"1 1"$var

  tt0=unique(time0)
  crrnew0=varnew0=tt0
  for(b in 1:length(tt0))
     {crrnew0[b]=max(crr0[time0==tt0[b]])
      varnew0[b]=var0[time0==tt0[b] & crr0==crrnew0[b]][1]
     }
  senew0=sqrt(varnew0)
  log.crr0=log(crrnew0/(1-crrnew0))
  log.se0=senew0/(crrnew0*(1-crrnew0))
  ci.logcrr0=cbind(log.crr0-qnorm((1+alpha)/2)*log.se0, log.crr0+qnorm((1+alpha)/2)*log.se0)
  ci0=exp(ci.logcrr0)/(1+exp(ci.logcrr0))
  ci0[varnew0==0,]=0
  ci01=ci0[,1]
  ci02=ci0[,2]

  tt1=unique(time1)
  crrnew1=varnew1=tt1
  for(b in 1:length(tt1))
    {crrnew1[b]=max(crr1[time1==tt1[b]])
     varnew1[b]=var1[time1==tt1[b] & crr1==crrnew1[b]][1]
     }
  senew1=sqrt(varnew1)
  result1=cbind(tt1, crrnew1, varnew1)
  log.crr1=log(crrnew1/(1-crrnew1))
  log.se1=senew1/(crrnew1*(1-crrnew1))
  ci.logcrr1=cbind(log.crr1-qnorm((1+alpha)/2)*log.se1, log.crr1+qnorm((1+alpha)/2)*log.se1)
  ci1=exp(ci.logcrr1)/(1+exp(ci.logcrr1))
  ci1[varnew0==0,]=0
  ci11=ci1[,1]
  ci12=ci1[,2]

  if(length(time)>0)
     {mtime=max(time)
      if(mtime>max(tt0))
      {message("The CRR in group 0 is not identifiable at some selected time points,
which are replaced by the maximum time point, where it is identifiable.")}
      time0=time[time<max(tt0)]
      if(mtime>max(tt0))
         {time0=c(time0, max(tt0))}
      crrnew0=approx(x=tt0, y=crrnew0, xout=time0, method="constant", yleft=0, yright=0)$y
      senew0=approx(x=tt0, y=senew0, xout=time0, method="constant", yleft=0, yright=0)$y
      ci01=approx(x=tt0, y=ci01, xout=time0, method="constant", yleft=0, yright=0)$y
      ci02=approx(x=tt0, y=ci02, xout=time0, method="constant", yleft=0, yright=0)$y
      tt0=time0
     }

  if(length(time)>0)
    {mtime=max(time)
     if(max(time)>max(tt1))
      {message("The CRR in group 1 is not identifiable at some selected time points,
which are replaced by the maximum time point, where it is identifiable.")}
     time1=time[time<max(tt1)]
     if(mtime>max(tt1))
        {time1=c(time1, max(tt1))}
     crrnew1=approx(x=tt1, y=crrnew1, xout=time1, method="constant", yleft=0, yright=0)$y
     senew1=approx(x=tt1, y=senew1, xout=time1, method="constant", yleft=0, yright=0)$y
     ci11=approx(x=tt1, y=ci11, xout=time1, method="constant", yleft=0, yright=0)$y
     ci12=approx(x=tt1, y=ci12, xout=time1, method="constant", yleft=0, yright=0)$y
     tt1=time1
     }

  result1=cbind(tt1, crrnew1, senew1, ci11, ci12)
  result0=cbind(tt0, crrnew0, senew0, ci01, ci02)


  colnames(result1)=colnames(result0)=c("time", "CRR", "se", "ci-low", "ci-up")

  result0=data.frame(result0)
  result1=data.frame(result1)

  return(list(result0=result0, result1=result1, pvalue=fit$Tests[1,2]))
}
