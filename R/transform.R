#' Apply Abbott's correction
#'
#' Apply Abbott's correction to morbidity data
#' @noRd
probit_C <- function(Cx,ii,dataf,x){
  datac<-dataf[dataf$dose==0,]
  data<-dataf[dataf$dose>0,]
  data$mort[ii]<-(data$mort[ii]-Cx)/(1-Cx)
  dataf$dead[ii]<-dataf$mort[ii]*dataf$total[ii]
  y<-cbind(data$dead[ii],data$total[ii]-data$dead[ii])
  moda<-glm(y~log10(data$dose[ii]),family = quasibinomial(link=probit))
  L<--moda$deviance/2
  j<-datac$strain==x
  L<-L+sum(datac$dead[j]*log10(Cx)+(datac$total[j]-datac$dead[j])*log10(1-Cx))
  return(L)
}


#' Probit-transform the data and apply Abbott's correction
#'
#' This function applies probit transformation to the data, after applying
#' Abbott's correction (see reference) when control groups (e.g. unexposed
#' susceptible strain) show non-negligible mortality.
#'
#' @param dataf a data frame of mortality data containing four mandatory
#' columns "strain", "dose", "total", "dead" (not necessarily in that order).
#' @param conf numerical. Threshold for the mortality in the controls above
#' which the correction should be applied (default=0.05)
#'
#' @importFrom stats glm optim qnorm quasibinomial runif
#'
#' @return Returns a list. convrg: with correction values and convergence
#' (NULL if mortality in the controls is below conf.), tr.data: transformed
#' data
#'
#' @author Pascal Milesi, Piyal Karunarathne, Pierrick Labbé
#'
#' @references Abbott, WS (1925). A method of computing the effectiveness of
#' an insecticide. J. Econ. Entomol.;18:265‐267.
#'
#' @examples
#' data(bioassay)
#' transd<-probit.trans(bioassay)
#' @export
probit.trans<-function(dataf,conf=0.05){
  mort<-ifelse(dataf$dead/dataf$total==0,0.006,ifelse(dataf$dead/dataf$total==1,
                                                      1-0.006,dataf$dead/dataf$total))
  dataf<-cbind(dataf,mort)
  if(any(dataf$dose==0)){
    if(any(dataf[dataf$dose==0,"mort"]>conf)){
      st<-unique(as.character(dataf$strain))
      tt<-lapply(st,function(x,dataf){
        data<-dataf[dataf$dose>0,]
        sig<-data$strain==x
        if(!any(data$mort[sig]==0)){
          bottom<-1e-12;top<-min(data$mort[sig])
          pin<-runif(1,min=bottom,max=top)
          opz<-optim(pin,probit_C,ii=sig,dataf=dataf,x=x,control=list(fnscale=-1,trace=1),
                     method="L-BFGS-B",lower=bottom,upper=top)
          val<-c(opz$par,opz$convergence)
        } else {
          val<-c(0,0)
        }
        return(c(x,val))
      },dataf=dataf)
      tt<-data.frame(do.call(rbind,tt))
      colnames(tt)<-c("Strain","ControlMortality","Convergence(OK if 0)")
      data<-dataf[dataf$dose>0,]
      for(i in seq_along(tt[,1])){
        data$mort[data[,"strain"]==tt[i,1]]<-(data$mort[data[,"strain"]==tt[i,1]]-as.numeric(tt[i,2]))/(1-as.numeric(tt[i,2]))
        data$mort[data$mort==0]<-0.000000005
        data$dead[data[,"strain"]==tt[i,1]]<-data$mort[data[,"strain"]==tt[i,1]]*data$total[data[,"strain"]==tt[i,1]]
      }
    } else {
      data<-dataf[dataf$dose>0,]
      tt<-NULL
    }
  } else {
    data<-dataf
    tt<-NULL
  }
  probmort<-sapply(data$mort,qnorm) # apply probit transformation to the data
  data<-cbind(data,probmort)
  data$strain<-as.factor(data$strain)
  outdata<-list(convrg=tt,tr.data=data)
  return(outdata)
}


