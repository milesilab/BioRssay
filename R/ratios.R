# helpers
#1. function to calculate the CI of the regressions for each strain
CIplot<-function(mods,pmort_min,pmort_max,conf.level){
  summ<-summary(mods)
  zz<-qnorm((1-conf.level)/2,lower.tail=FALSE)
  a<-summ$coefficient[2] # mods slope
  b<-summ$coefficient[1] # mods intercept
  minldose<-((pmort_min-b)/a) # predicted dose for 0% mortality
  maxldose<-((pmort_max-b)/a) # predicted dose for 100% mortality
  datalfit<-seq(minldose-0.2,maxldose+0.2,0.01) # generates a set of doses
  datafit<-data.frame(dose=10^datalfit)
  pred<-predict.glm(mods,newdata=datafit,type="response",se.fit=TRUE) #generates the predicted mortality for the set of doses and SE
  ci<-cbind(pred$fit-zz*pred$se.fit,pred$fit+zz*pred$se.fit) #CI for each dose
  probitci<-cbind(datafit$dose,suppressWarnings(qnorm(ci))) # apply probit transformation to CI
  return(probitci)
}

#2. glm?
LD <- function(mod, conf.level) { ## mod=glm regression, d=dose (untransformed), confidence interval level
  p <- c(50,95) # the LD values of interest
  het = deviance(mod)/df.residual(mod)
  if(het < 1){het = 1} # Heterogeneity cannot be less than 1
  summary <- summary(mod, dispersion=het, cor = F)
  intercept <- summary$coefficients[1]
  interceptSE <- summary$coefficients[3]
  slope <- summary$coefficients[2]
  slopeSE <- summary$coefficients[4]
  z.value <- summary$coefficients[6]

  b0<-intercept # Intercept (alpha)
  b1<-slope # Slope (beta)
  vcov = summary(mod)$cov.unscaled
  var.b0<-vcov[1,1] # Intercept variance
  var.b1<-vcov[2,2] # Slope variance
  cov.b0.b1<-vcov[1,2] # Slope intercept covariance

  # Adjust alpha depending on heterogeneity (Finney, 1971, p. 76)
  alpha<-1-conf.level
  if(het > 1) {talpha <- -qt(alpha/2, df=df.residual(mod))} else {talpha <- -qnorm(alpha/2)}

  ## Calculate g (Finney, 1971, p 78, eq. 4.36)
  # "With almost all good sets of data, g will be substantially smaller than 1.0 and seldom greater than 0.4."
  g <- het * ((talpha^2 * var.b1)/b1^2)

  # Calculate theta.hat for all LD levels based on probits in eta (Robertson et al., 2007, pg. 27; or "m" in Finney, 1971, p. 78)
  eta = family(mod)$linkfun(p/100)  #probit distribution curve
  theta.hat <- (eta - b0)/b1
  # Calculate correction of fiducial limits according to Fieller method (Finney, 1971, p. 78-79. eq. 4.35)
  const1 <- (g/(1-g))*(theta.hat + cov.b0.b1/var.b1) # const1 <- (g/(1-g))*(theta.hat -   cov.b0.b1/var.b1)
  const2a <- var.b0 + 2*cov.b0.b1*theta.hat + var.b1*theta.hat^2 - g*(var.b0 - (cov.b0.b1^2/var.b1))
  const2 <- talpha/((1-g)*b1) * sqrt(het * (const2a))
  #Calculate the confidence intervals LCL=lower, UCL=upper (Finney, 1971, p. 78-79. eq. 4.35)
  LCL <- (theta.hat + const1 - const2)
  UCL <- (theta.hat + const1 + const2)
  #Calculate variance for theta.hat (Robertson et al., 2007, pg. 27)
  var.theta.hat <- (1/(theta.hat^2)) * ( var.b0 + 2*cov.b0.b1*theta.hat + var.b1*theta.hat^2 )

  # Make a vector of the different values of interest
  ECtable <- c(10^theta.hat[1],10^LCL[1],10^UCL[1],var.theta.hat[1],10^theta.hat[2],10^LCL[2],10^UCL[2],var.theta.hat[2],slope,slopeSE,intercept,interceptSE,het,g)
  return(ECtable)
}

#3. test the significance of model pairs of strains
reg.pair<-function(data){
yy<-cbind(data$dead,data$total-data$dead)
mod1<-glm(yy~data$strain/log10(data$dose)-1,family = quasibinomial(link=probit))
mod2<-glm(yy~log10(data$dose),family = quasibinomial(link=probit))
anova(mod1,mod2,test="Chi")
}

#3. get dxt
get.dxt<-function(strains,data,conf.level){dxt<-lapply(strains,function(ss,data,conf.level){
  tmp<-data[data$strain == ss,]
  y<-with(tmp,cbind(dead,total-dead))
  mods<-glm(y~log10(dose),data=tmp,family = quasibinomial(link=probit))
  res.strain <- LD(mods, conf.level)
  dat<-res.strain
  chq<-sum(((mods$fitted.values*tmp$total-tmp$dead)^2)/(mods$fitted.values*tmp$total))
  dat<-c(dat,pchisq(q=chq,df=length(tmp$dead)-2,lower.tail=FALSE))
  return(list(mods,dat))
},data=data,conf.level=conf.level)
return(dxt)
}

#' Calculate regression coefficients and resistance ratios
#'
#' This function calculates coefficients of the linear regression and the LD values for the strains tested using a GLM. Further it also outputs 95 percent confidence intervals of the RR (log-dose), according to Robertson and Preisler 1992 (also check references).
#'
#' @param data a data frame of probit transformed mortality data using the function probit.trans
#' @param conf.level numerical. confidence level for the LD calculation (default 0.95)
#' @param ref.strain character. name of the reference strain if present (see details)
#' @param plot logical. Whether to plot the mortality plot
#' @param plot.conf logical. If plot=TRUE, whether to plot the confidence intervals
#' @param test.validity logical. (If plot=TRUE) Test and plot the validity of the models
#' @param ... parameters to be passed on to graphics for the plot (e.g. col, pch)
#'
#' @importFrom graphics abline axis legend lines mtext
#' @importFrom stats deviance df.residual family pchisq predict.glm qt qnorm
#'
#' @details If a name of the reference strain is provided as it appears in the data, it will be used in the model. If not provided, the function will look for strains with the suffix "ref" in their names and will use it as the reference. If this returns NULL, the function will run the model considering no reference strain.
#'
#' @return returns a data frame of statistical output and if plot=TRUE, plots the mortality on a log transformed scale.
#'
#' @author Pascal Milesi, Piyal Karunarathne
#'
#' @references Finney DJ(1971). Probitanalysis. Cambridge:Cambridge UniversityPress. 350p.
#'
#' HommelG(1988). A stage wise rejective multiple test procedure based on a modified Bonferroni test. Biometrika 75, 383-6.
#'
#' Johnson RM, Dahlgren L, Siegfried BD,EllisMD(2013). Acaricide,fungicide and druginteractions in honeybees (Apis mellifera). PLoSONE8(1): e54092.
#'
#' Robertson, J. L., and H.K. Preisler.1992. Pesticide bioassays with arthropods. CRC, Boca Raton, FL.
#'
#' @examples
#' data(bioassay)
#' transd<-probid.trans(bioassay$assay2)
#' data<-transd$tr.data
#' resist.ratio(data,plot=TRUE)
#'
#' @export
resist.ratio<-function(data,conf.level=0.95,ref.strain=NULL,plot=FALSE,plot.conf=TRUE,test.validity=TRUE,...) {
  data$strain<-as.factor(data$strain)
  strains<-levels(data$strain)
  dxt<-get.dxt(strains,data,conf.level)
  dat<-do.call(rbind,lapply(dxt,function(x){x[[2]]}))
  colnames(dat)<-c("LD50", "LD50min","LD50max","varLD50","LD95", "LD95min","LD95max","varLD95",
                   "Slope", "SlopeSE", "Intercept", "InterceptSE", "h", "g", "Chi(p)")
  rownames(dat)<-strains

  if(is.null(ref.strain)){
    ref <- which(strains == strains[grep("-ref$",as.character(strains))], arr.ind=TRUE)
  } else {
    ref=ref.strain
  }
  if (length(ref)==0) refrow <- which(dat[,"LD50"]==min(dat[,"LD50"]),arr.ind=TRUE) else refrow <-ref

  rr50<-dat[,"LD50"]/dat[refrow,"LD50"]
  rr95<-dat[,"LD95"]/dat[refrow,"LD95"]
  CI50<-1.96*sqrt(dat[,"varLD50"]+dat[refrow,"varLD50"])
  rr50max<-10^(log10(rr50)+CI50);rr50max[refrow]<-0
  rr50min<-10^(log10(rr50)-CI50);rr50min[refrow]<-0
  CI95<-1.96*sqrt(dat[,"varLD95"]+dat[refrow,"varLD95"])
  rr95max<-10^(log10(rr95)+CI50);rr95max[refrow]<-0
  rr95min<-10^(log10(rr95)-CI50);rr95min[refrow]<-0
  if(plot){
    mort.plot(data,strains,plot.conf,test.validity=test.validity,conf.level=conf.level,...)
  }
  dat<-cbind(dat,rr50,rr50max,rr50min,rr95,rr95max,rr95min)
  return(dat)
}



#' Test mortality model significance
#'
#' This function tests whether there are differences between the regression of mortality between the strains (likelihood ratio test); It compares a model with and without different strains (i.e. to compare the significance if the data came from one strain or different strains)
#'
#' @param data a data frame of probit transformed mortality data using the function probit.trans
#' @importFrom stats anova
#' @importFrom utils combn
#'
#' @return a list with model outputs: a chi-square test if there only two strains or bonferroni test significance from a pariwise model comparison if there are more than two strains
#'
#' @author Pascal Milesi, Piyal Karunarathne
#'
#' @examples
#' data(bioassay)
#' transd<-probid.trans(bioassay$assay2)
#' data<-transd$tr.data
#' model.signif(data)
#'
#' @export
model.signif<-function(data){
  data$strain<-as.factor(data$strain)
  strains<-levels(data$strain)
  if (length(strains)>=2) {
    if (length(strains)==2) {
      Test<-reg.pair(data)
    } else {
      Test<-sapply(strains, function(x,data) sapply(strains, function(y,data){if(x!=y){dat<-data[data$strain==x | data$strain==y,];reg.pair(dat)$Pr[2]}},data=data),data=data)
      dv<-sapply(strains, function(x,data) sapply(strains, function(y,data){if(x!=y){dat<-data[data$strain==x | data$strain==y,];reg.pair(dat)$Deviance[2]}},data=data),data=data)
      dff=sapply(strains, function(x,data) sapply(strains, function(y,data){if(x!=y){dat<-data[data$strain==x | data$strain==y,];reg.pair(dat)$Df[2]}},data=data),data=data)
      pval<-unique(sort(unlist(Test)))
      rk<-(length(strains)*(length(strains)-1)/2):1
      rk<-0.05/rk
      toget<-t(combn(rownames(Test),2))
      Test<-data.frame(cbind(toget,unlist(Test[toget]),unlist(dv[toget]),unlist(dff[toget]),ifelse(pval<rk,"sig","non-sig")))
      colnames(Test)<-c("strain1","strain2","model.pval","deviance","df","bonferroni")
    }
  } else {
    message("Only one strain present; check your data")
  }
  return(list(model=Test))
}
# add deviance and degree of freedom to the table






