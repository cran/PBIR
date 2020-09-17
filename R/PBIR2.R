#' Estimate and compare PBIR curves from two groups over a time window
#'
#'@param  t2PROGRESSION  time to progression/death or censoring
#'@param  STATUS_PROGRESSION   binary indicator for progression status: 1 for progression/death; 0 for censoring
#'@param  t2RESPONSE  time to response or censoring
#'@param  STATUS_RESPONSE  binary indicator for response status: 1 for response; 0 for censoring
#'@param  TRT         treatment indicator: 1 for treatment arm; 0 for control arm
#'@param  time        user-selected time points at which PBIRs are to be compared; the default value is "NULL" and PBIRs at all observed time points are compared
#'@param  alpha       coverage level of the point-wise confidence interval for the difference in the PBIR, the default value is 0.95
#'@return a data matrix containing "time", "estimated differences in PBIR (treatment-control)", "standard errors of estimated PBIR differences",  "confidence intervals of the PBIR difference"
#'@examples
#'
#'library(survival)
#'n=100
#'set.seed(10)
#'
#'# Generate the data
#'
#'TRT=trt=rbinom(n, 1, 0.5)
#'
#'error=rnorm(n)
#'tr=exp(rnorm(n)+error-trt*0.5+0.5)
#'tp=exp(rnorm(n)+error+trt*0.25)
#'tr[tp<tr]=Inf
#'tc=runif(n, 3, 8.5)
#'
#'t2response=pmin(tr, tc)
#'delta_response=1*(tr<tc)
#'t2progression=pmin(tp, tc)
#'delta_progression=1*(tp<tc)
#'
#'# Estimate the difference in PBIR
#'# the analysis is truncated at time 8, which is slightly smaller than the largest follow-up time
#'
#'fit=PBIR2(t2PROGRESSION=t2progression,
#'           STATUS_PROGRESSION=delta_progression,
#'           t2RESPONSE=t2response,
#'           STATUS_RESPONSE=delta_response,
#'           TRT=trt)
#'
#'
#'# Plot the estimated differnece in PBIR
#'
#'tt=fit$time
#'diff=fit$diff
#'low=fit$ci.low
#'up=fit$ci.up
#'
#'tt=c(0, tt)
#'diff=c(0, diff)
#'low=c(0, low)
#'up=c(0, up)
#'B=length(tt)
#'
#'tt=rep(tt, rep(2, B))[-1]
#'diff=rep(diff, rep(2, B))[-(2*B)]
#'low=rep(low, rep(2, B))[-(2*B)]
#'up=rep(up, rep(2, B))[-(2*B)]
#'
#'plot(range(c(fit$time, 0)), range(c(low, up)),
#'      xlab="time", ylab="difference in PBIR",
#'      lwd=2, type="n")
#'lines(tt, diff, lwd=2, col=3)
#'lines(tt, low,  col=2)
#'lines(tt, up, col=2)
#'lines(range(fit$time), rep(0, 2), col=4, lty=4)
#'@importFrom("survival", survfit)
#'@references Huang, B., Tian, L., Talukder, E., Rothenberg, M., Kim, DY., and Wei, LJ. (2018) Evaluating Treatment Effect Based on Duration of Response for a Comparative Oncology Study. JAMA Oncol, doi: 10.1001/jamaoncol.2018.0275
#'@references Huang, B., Tian, L., McCaw, Z., Luo, Talukder, E., X., Rothenberg, M., Xie, W., Choueiri, T., Kim, DY., & Wei, LJ. (2020). Analysis of Response Data for Assessing Treatment Effects in Comparative Clinical Studies. Ann Intern Med, doi: 10.7326/M20-0104.
#'@export


PBIR2=function(t2PROGRESSION, STATUS_PROGRESSION, t2RESPONSE, STATUS_RESPONSE, TRT, time=NULL, alpha=0.95){

  res1=PBIR1(t2PROGRESSION[TRT==1], STATUS_PROGRESSION[TRT==1], t2RESPONSE[TRT==1], STATUS_RESPONSE[TRT==1], time=time)
  res0=PBIR1(t2PROGRESSION[TRT==0], STATUS_PROGRESSION[TRT==0], t2RESPONSE[TRT==0], STATUS_RESPONSE[TRT==0], time=time)

  time.max1=max(res1$time)
  time.max0=max(res0$time)

  taumax=min(time.max1, time.max0)

  res1=res1[res1$time<=taumax, ]
  res0=res0[res0$time<=taumax, ]

  n1=length(res1[,1])
  n0=length(res0[,1])

  time0=sort(unique(c(res1$time, res0$time)))
  n.time0=length(time0)

  diff=rep(NA, n.time0)
  sd.diff=rep(NA, n.time0)
  ci1.diff=rep(NA, n.time0)
  ci2.diff=rep(NA, n.time0)

  for(i in 1:n.time0)
  {id1=max((1:(n1+1))[c(0, res1$time)<=time0[i]])
   id0=max((1:(n0+1))[c(0, res0$time)<=time0[i]])
   diff[i]=c(0,res1$PBIR)[id1]-c(0, res0$PBIR)[id0]
   sd.diff[i]=sqrt(c(0,res1$std)[id1]^2+c(0, res0$std)[id0]^2)
  }

  diff.tr=log((1-diff)/(1+diff))
  sd.difftr=2*sd.diff/(1-diff^2)

  ci.tr=cbind(diff.tr-qnorm((1+alpha)/2)*sd.difftr, diff.tr+qnorm((1+alpha)/2)*sd.difftr)
  ci=(1-exp(ci.tr))/(1+exp(ci.tr))

  res=cbind(time0, diff, sd.diff, ci[,2], ci[,1])
  colnames(res)=c("time", "diff in PRIB", "std", "ci-low", "ci-up")

  if(length(time)>0 && taumax<max(time))
    {message("The PBIR difference is only estimated at time points, where it is identifiable.")
    }

  res=data.frame(res)

  return(res)
}
