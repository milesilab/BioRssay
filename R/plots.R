#' Test validity of the probit model
#' @noRd
validity<-function(strains,data){
  ndataf<-do.call(rbind,lapply(strains,function(ss,data){
    tmp<-data[data$strain==ss,]
    dss<-unique(tmp$dose)
    dd<-do.call(rbind,lapply(dss,function(x,tmp){
      tt<-tmp[tmp$dose==x,]
      cbind(sum(tt$dead),sum(tt$total))
    },tmp))
    tx<-cbind(ss,dss,dd)
  },data=data))
  ndataf<-data.frame(strain=ndataf[,1],apply(ndataf[,-1],2,as.numeric))
  colnames(ndataf)<-c("strain","dose","dead","total")
  d.ratio<-ifelse(ndataf$dead/ndataf$total==0,0.006,
                  ifelse(ndataf$dead/ndataf$total==1,1-0.006,ndataf$dead/ndataf$total))
  p.mortality<-sapply(d.ratio,qnorm)
  ndataf<-cbind(ndataf,d.ratio,p.mortality)
  return(ndataf)
}

#' Legend assembly
#' @noRd
lgas<-function(legend.par, ll,strains){llg<-as.list(legend.par)
lpos<-c("bottomleft","bottomright","topleft","topright","top","bottom","center")
ps<-match(lpos,legend.par)
#if(any(lpos==llg[[1]])) {llg$x<-llg[[1]]} else {llg$x <-"bottomleft"}
if(sum(ps,na.rm = TRUE)>0){x<-unlist(llg[na.omit(ps)]);llg<-llg[-na.omit(ps)];llg$x<-x} else {llg$x <-"bottomleft"}
if(any(names(llg)=="y")) llg$y<-llg$y
if(!any(names(llg)=="legend")) llg$legend<-strains
if(is.null(llg$col)) llg$col=rainbow_hcl(length(strains))
#if(is.null(llg$pch)) {if(length(strains)<=6)llg$pch=15:20 else llg$pch=1:20}
if(!is.null(ll$pch)){llg$pch<-ll$pch}
if(is.null(llg$lwd)) llg$lwd=1.5
llg$lwd<-as.numeric(llg$lwd)
if(is.null(llg$cex)) llg$cex=0.8
llg$cex<-as.numeric(llg$cex)
if(any(is.na(ll$pch))) {llg$lty=1;llg$pch=NA}
llg$lty<-as.numeric(llg$lty)
if(is.null(llg$bg)) llg$bg="grey60"
if(is.null(llg$bty)) llg$bty="o"
if(is.null(llg$box.col)) llg$box.col=NA
lnames<-names(formals(legend))

do.call("legend",llg[names(llg)%in%lnames])}

#' Plot dose-mortality response for each strain
#'
#' This function plots the probit-transformed mortalities (probit.trans() function) as a function of the log10 of the dose, the regressions predicted by the resist.ratio() function,  with or without confidence levels, if the dose-mortality responses are linear (option).
#'
#' @param data a data frame of probit transformed mortality data using the function probit.trans()
#' @param strains character. list of test strains to be plotted. If not provided, the function will plot all the strains in the data set.
#' @param plot.conf logical. Whether to plot the confidence intervals for each strain, default TRUE
#' @param conf.level numerical. The confidence interval to be plotted
#' @param LD.value numerical. Level of lethal dose to be tested. default=c(25,50,95)
#' @param test.validity logical. When TRUE (default), if a strain mortality-dose response fails the chi-square test for linearity in the resist.ratio() function, no regression will be plotted, only the observed data.
#' @param legend.par multi-type. Arguments to be passed to the legend as in \code{\link[graphics]{legend}}. default position \code{bottomleft}. If no legend desired use FALSE. Note: if pch, lty, and col are passed to the plot, they don't need to be passed to \code{legend()}
#' @param ... parameters to be passed on to graphics for the plot (e.g. col, pch)
#'
#' @importFrom graphics points layout par plot.default title
#' @importFrom colorspace rainbow_hcl
#'
#' @return A plot of dose-mortality responses for bioassays
#'
#' @author Piyal Karunarathne, Pascal Milesi, Pierrick Labbé
#'
#' @examples
#' data(bioassay)
#' transd<-probit.trans(bioassay$assay2)
#' data<-transd$tr.data
#' strains<-levels(data$strain)
#' mort.plot(data,strains)
#'
#' @export
mort.plot<-function(data,strains=NULL,plot.conf=TRUE,conf.level=0.95,
                    LD.value=c(25,50,95),test.validity=TRUE,legend.par=c("bottomleft"),...){
  #opars<-par(no.readonly = TRUE)
  #on.exit(par(opars))
  data$strain<-as.factor(data$strain)
  if(is.null(strains)){
    strains<-levels(data$strain)
  }
  dmin<-floor(log10(min(data$dose)))
  dmax<-ceiling(log10(max(data$dose)))
  dose_min<- 10^(dmin)
  dose_max<- 10^(dmax)
  pmort_min<- qnorm(0.006)
  pmort_max<- qnorm(1-0.006)
  ll<-list(...)
  if(is.null(ll$col)) ll$col=rainbow_hcl(length(strains))
  if(is.null(ll$pch)) {if(length(strains)<=6)ll$pch=15:20 else ll$pch=1:20}
  if(is.null(ll$conf.level)) ll$conf.level=0.95
  if(is.null(ll$lwd)) ll$lwd=1.5
  if(is.null(ll$cex)) ll$cex=1

  if(is.null(ll$xlim)) {ll$xlim=c(dose_min,dose_max)}
  if(is.null(ll$ylim)) {ll$ylim=c(floor(pmort_min*100)/100,ceiling(pmort_max*100)/100)}
  if(is.null(ll$ylab)) {ll$ylab="mortality"}
  if(is.null(ll$yaxt)) {ll$yaxt="n"}
  if(is.null(ll$xaxt)) {ll$xaxt="n"}
  if(is.null(ll$log)) {ll$log="x"}
  if(is.null(ll$ann)) {ll$ann=FALSE}
  cl<-ll$col
  ph<-ll$pch
  ll<-ll[-which(names(ll)=="pch")]
  ll<-ll[-which(names(ll)=="col")]
  pnames<-c(names(formals(plot.default)),names(par()))

  dxt<-get.dxt(strains,data,ll$conf.level,LD.value=LD.value)

  do.call("plot",c(list(x=data$dose,y=data$probmort),col=list(cl[data$strain]),pch=list(ph[data$strain]),ll[names(ll)%in%pnames]))
  if(!is.null(ll$main)){title(ll$main)}
  ll$col<-cl
  ll$pch<-ph

  abline(v = dose_min, col = "grey95", lwd = 180000)
  points(data$dose,data$probmort,col=ll$col[data$strain],
         pch=ll$pch[data$strain])

  labely<-c(1,5,seq(10,90,10),95,99)
  axis(2, at=qnorm(labely/100),labels=labely,las=2, adj=0)
  axis(4, at=qnorm(labely/100),labels=FALSE)
  mtext("Mortality (%)", side=2, line=3)

  for (i in dmin:dmax) {
    axis(1,at=10^i,labels=substitute(10^k,list(k=i)))
  }
  axis.at <- 10 ^ c(dmin:dmax)
  axis(1, at = 2:9 * rep(axis.at[-1] / 10, each = 8),
       tcl = -0.5, labels = FALSE)
  mtext(expression(Dose (mg.L^-1) ), side=1, line=3)
  abline(h=pmort_min, lty=3)
  abline(h=pmort_max, lty=3)

  if(plot.conf){
    if(test.validity){
      ndataf<-validity(strains,data)
      for(i in 1:length(strains)){
        if(dxt[[i]][[2]][[length(dxt[[i]][[2]])]]>0.05){ # dxt[[i]][[2]][[15]]
          abline(dxt[[i]][[1]], col=ll$col[i],lwd=ll$lwd)
          CIfit<-CIplot(dxt[[i]][[1]],pmort_min,pmort_max,conf.level=conf.level)
          lines(CIfit[,1],CIfit[,2],type="l", lty=3, col=ll$col[i],lwd=ll$lwd)
          lines(CIfit[,1],CIfit[,3],type="l", lty=3, col=ll$col[i],lwd=ll$lwd)
        } else {
          points(ndataf$dose[ndataf$strain==strains[i]],
                 ndataf$p.mortality[ndataf$strain==strains[i]],type="l", col=ll$col[i])
        }
      }
    } else {
      for(i in 1:length(strains)){
        abline(dxt[[i]][[1]], col=ll$col[i],lwd=ll$lwd)
        CIfit<-CIplot(dxt[[i]][[1]],pmort_min,pmort_max,conf.level=conf.level)
        lines(CIfit[,1],CIfit[,2],type="l", lty=3, col=ll$col[i],lwd=ll$lwd)
        lines(CIfit[,1],CIfit[,3],type="l", lty=3, col=ll$col[i],lwd=ll$lwd)
      }
    }
  } else {
    for(i in 1:length(strains)){
      abline(dxt[[i]][[1]], col=ll$col[i],lwd=ll$lwd)
    }
  }

  if(length(legend.par)<2){
    if(!isFALSE(legend.par)){
      lgas(legend.par,ll,strains)
    }
  } else {
    lgas(legend.par,ll,strains)
  }
}

