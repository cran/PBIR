#' Estimate mean duration of response
#'
#'@param  t2PROGRESSION  time to progression/death or censoring
#'@param  STATUS_PROGRESSION   binary indicator for progression/death status: 1 for progression/death; 0 for censoring
#'@param  t2RESPONSE  time to response or censoring
#'@param  STATUS_RESPONSE  binary indicator for response status: 1 for response; 0 for censoring
#'@param  time.max    maximum time point, up to which the mean DOR is to be estimated; the default value corresponds to the maximum time window in which the mean DOR is estimable
#'@details The mean duration of response restricted within a time window is also the area under the PBIR curve over the same time window. The estimated mean duration can be viewed as a global summary of the PBIR curve. One may compare the mean duration of response between two groups, which is also a global comparison between two PBIR curves.
#'@return A list with following elements
#'\itemize{
#'\item{meandor.est:}  the restricted mean DOR estimate
#'\item{meandor.se:}   the standard error of the esimated DOR
#'\item{time.truncation:}     the truncation time point used in DOR.
#'}
#'@examples
#'
#'library(survival)
#'n=100
#'set.seed(10)
#'
#'# Generate the data
#'
#'error=rnorm(n)
#'tr=exp(rnorm(n)+error+0.5)
#'tp=exp(rnorm(n)+error)
#'tr[tp<tr]=Inf
#'tc=runif(n, 3, 8.5)

#'t2response=pmin(tr, tc)
#'delta_response=1*(tr<tc)
#'t2progression=pmin(tp, tc)
#'delta_progression=1*(tp<tc)

#'
#'# Estimate the mean duration of response (point estimator and its standard error)
#'
#'fit=mduration(t2PROGRESSION=t2progression,
#'               STATUS_PROGRESSION=delta_progression,
#'               t2RESPONSE=t2response,
#'               STATUS_RESPONSE=delta_response,
#'               time.max=8)
#'
#'fit
#'@importFrom("survival", survfit)
#'@references Huang, B., Tian, L., Talukder, E., Rothenberg, M., Kim, DY., and Wei, LJ. (2018) Evaluating Treatment Effect Based on Duration of Response for a Comparative Oncology Study. JAMA Oncol, doi: 10.1001/jamaoncol.2018.0275
#'@references Huang, B., Tian, L., McCaw, Z., Luo, Talukder, E., X., Rothenberg, M., Xie, W., Choueiri, T., Kim, DY., & Wei, LJ. (2020). Analysis of Response Data for Assessing Treatment Effects in Comparative Clinical Studies. Ann Intern Med, doi: 10.7326/M20-0104.
#'@export


mduration=function(t2PROGRESSION, STATUS_PROGRESSION, t2RESPONSE, STATUS_RESPONSE, time.max=-1){

  t2RESPONSE[STATUS_RESPONSE==0]=Inf

  y=pmin(t2PROGRESSION, t2RESPONSE)
  delta=1*(STATUS_RESPONSE+STATUS_PROGRESSION>0)
  if(min(delta[y==max(y)])==0)
    taumax=max(pmin(t2PROGRESSION, t2RESPONSE))
  if(min(delta[y==max(y)])==1)
    taumax=max(t2PROGRESSION)


  if(time.max<0)
     time.max=taumax

  if(time.max>taumax)
     {time.max=taumax
      message("time.max is reduced to", taumax)
      }


  y1=pmin(t2PROGRESSION, t2RESPONSE, time.max)
  delta1=1*(STATUS_RESPONSE+STATUS_PROGRESSION>0)
  delta1[y1==time.max]=1

  y2=pmin(t2PROGRESSION, time.max)
  delta2=STATUS_PROGRESSION
  delta2[y2==time.max]=1

  fit1=survfit(Surv(y1, delta1)~1)
  n1=length(fit1$time)
  auc1=sum(c(1, fit1$surv[-n1])*(fit1$time-c(0, fit1$time[-n1])))


  fit2=survfit(Surv(y2, delta2)~1)
  n2=length(fit2$time)
  auc2=sum(c(1, fit2$surv[-n2])*(fit2$time-c(0, fit2$time[-n2])))

  muhat=auc2-auc1

  N=length(y1)
  prisk1=pauc1=rep(0, N)
  for(i in 1:N)
  {k1=sum(fit1$time<=y1[i])
  if(k1>1)
    pauc1[i]=sum(c(1, fit1$surv[1:(k1-1)])*(fit1$time[1:k1]-c(0, fit1$time[1:(k1-1)])))
  if(k1==1)
    pauc1[i]=fit1$time[1]
  prisk1[i]=mean(y1>=y1[i])
  }
  pauc1=auc1-pauc1

  xi1=delta1*pauc1/prisk1

  tau1=rep(0, N)
  for(i in 1:N)
  {tau1[i]=xi1[i]-mean(xi1*(y1[i]>=y1)/prisk1)}

  prisk2=pauc2=rep(0, N)
  for(i in 1:N)
  {k2=sum(fit2$time<=y2[i])
  if(k2>1)
    pauc2[i]=sum(c(1, fit2$surv[1:(k2-1)])*(fit2$time[1:k2]-c(0, fit2$time[1:(k2-1)])))
  if(k2==1)
    pauc2[i]=fit2$time[1]
  prisk2[i]=mean(y2>=y2[i])
  }
  pauc2=auc2-pauc2

  xi2=delta2*pauc2/prisk2

  tau2=rep(0, N)
  for(i in 1:N)
  {tau2[i]=xi2[i]-mean(xi2*(y2[i]>=y2)/prisk2)}

  sigma=sqrt(var(tau1-tau2)/N)

  return(list(meandor.est=muhat, meandor.se=sigma, time.truncation=time.max))
}
