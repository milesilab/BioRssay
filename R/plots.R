#1. test validity of the model
validity<-function(strains,data){
  ndataf<-do.call(rbind,lapply(strains,function(ss,data){
    tmp<-data[data$strain==ss,]
    dss<-unique(tmp$dose)
    dd<-do.call(rbind,lapply(dss,function(x,tmp){tt<-tmp[tmp$dose==x,];cbind(sum(tt$dead),sum(tt$total))},tmp))
    tx<-cbind(ss,dss,dd)
  },data=data))
  ndataf<-data.frame(strain=ndataf[,1],apply(ndataf[,-1],2,as.numeric));colnames(ndataf)<-c("strain","dose","dead","total")
  d.ratio<-ifelse(ndataf$dead/ndataf$total==0,0.006,ifelse(ndataf$dead/ndataf$total==1,1-0.006,ndataf$dead/ndataf$total))
  p.mortality<-sapply(d.ratio,qnorm)
  ndataf<-cbind(ndataf,d.ratio,p.mortality)
  return(ndataf)
}


#' Plot mortality
#'
#' Plots the mortality of a given data set with and without confidence levels, and adds the validity of the regression.
#'
#' @param data a data frame of probit transformed mortality data using the function probit.trans
#' @param strains character. list of test strains. If not provided, the function will plot all the strains in the data set.
#' @param plot.conf logical. Whether to plot the confidence level on the plot for each strain
#' @param conf.level numerical. The confidence interval to be plotted
#' @param test.validity logical. Test and plot the validity of the models
#' @param ... parameters to be passed on to graphics for the plot (e.g. col, pch)
#'
#' @importFrom graphics points
#' @importFrom colorspace rainbow_hcl
#'
#' @return A plot of mortality
#'
#' @author Piyal Karunarathne, Pascal Milesi
#'
#' @examples
#' data(bioassay)
#' transd<-probid.trans(bioassay$assay2)
#' data<-transd$tr.data
#' strains<-levels(data$strain)
#' mort.plot(data,strains)
#'
#' @export
mort.plot<-function(data,strains=NULL,plot.conf=TRUE,conf.level=0.95,test.validity=TRUE,...){
  if(is.null(strains)){
    strains<-as.character(unique(data$strain))
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
  dxt<-get.dxt(strains,data,ll$conf.level)

  plot(data$dose,data$probmort,log="x",xlim=c(dose_min,dose_max),
       ylim=c(floor(pmort_min*100)/100,ceiling(pmort_max*100)/100),
       ylab="mortalit?",yaxt="n",xaxt="n", ann=FALSE ,col=ll$col[data$strain],
       pch=ll$pch[data$strain])
  abline(v = dose_min, col = "grey95", lwd = 180000)
  points(data$dose,data$probmort,col=ll$col[data$strain],
         pch=ll$pch[data$strain])

  labely<-c(1,5,seq(10,90,10),95,99)
  axis(2, at=qnorm(labely/100),labels=labely,las=2, adj=0)
  axis(4, at=qnorm(labely/100),labels=FALSE)
  mtext("Mortality (%)", side=2, line=3)
  for (i in dmin:dmax) axis(1,at=10^i,labels=substitute(10^k,list(k=i)))
  axis.at <- 10 ^ c(dmin:dmax)
  axis(1, at = 2:9 * rep(axis.at[-1] / 10, each = 8),
       tcl = -0.5, labels = FALSE)
  mtext(expression(Dose (mg.L^-1) ), side=1, line=3)
  legend("bottomright", strains, col = ll$col, pch=ll$pch, cex = 0.8,
         inset=c(0,1), xpd=TRUE, horiz=TRUE, bty="n")
  abline(h=pmort_min, lty=3)
  abline(h=pmort_max, lty=3)

  if(plot.conf){
    if(test.validity){
      ndataf<-validity(strains,data)
      for(i in 1:length(strains)){
        if(dxt[[i]][[2]][[15]]>0.05){
          abline(dxt[[i]][[1]], col=ll$col[i],lwd=ll$lwd)
          CIfit<-CIplot(dxt[[i]][[1]],pmort_min,pmort_max,conf.level=conf.level)
          lines(CIfit[,1],CIfit[,2],type="l", lty=3, col=ll$col[i],lwd=ll$lwd)
          lines(CIfit[,1],CIfit[,3],type="l", lty=3, col=ll$col[i],lwd=ll$lwd)
        } else {
          points(ndataf$dose[ndataf$strain==strains[i]],ndataf$p.mortality[ndataf$strain==strains[i]],type="l", col=ll$col[i])
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
}

