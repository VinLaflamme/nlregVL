#' Example data
#' 
#' A dataset used for examples
#' \itemize{
#' \item Duration. Duration in ms.
#' \item Sub. Participant code.
#' \item Nb1. Number of short responses.
#' \item Nb3. Number of long responses.
#' }
#' 
#' @docType data
#' @name exampleData
#' @usage data(exampleData)
#' @format data.table with 105 observation of 4 variables.

NULL

#' Settings for the fitting process
#' 
#' @slot xx quantity on the x axis (e.g. duration). Must be a vector of numerics
#' @slot pp quantity on the y axis (e.g. proportion of long responses). Must be a vector of numerics AND of the same length as xx.
#' @slot getPred Function that returns the predicted value.
#' @slot residFun Function that returns the residuals.
#' @slot fnJacob Function that returns the jacobian of the function to be fitted.
#' @slot fnStart Function that returns the starting parameters.
#' @export

setClass("psyFitset",representation = representation(xx="numeric",pp="numeric",getPred="function",residFun="function",fnJacob="function",fnStart="function"))
.psyFitset.valid<-function(object){
  if(length(object@xx)!=length(object@pp)){
    return("length mismatch")
  }else{
    if(!all(names(formals((object@getPred)))==c("parS","xx"))) return("Bad arguments in getPred")
    if(!all(names(formals((object@residFun)))==c("p","observed","xx"))) return("Bad arguments in residFun")
    if(!all(names(formals((object@fnJacob)))==c("p","observed","xx"))) return("Bad arguments in fnJacob")
    if(!all(names(formals((object@fnStart)))==c("xx"))) return("Bad arguments in fnStart")
  }
  return(T)
}

setValidity("psyFitset",.psyFitset.valid)

#' An S4 class containing the results of a non-linear regression by the levenberg-Marquardt algorithm
#' @slot psyF psyFitset object containing the settings for the regression
#' @slot pars list of the estimated parameters.
#' @slot Chi Chi square statistic.
#' @slot ci matrix containing the confidence intervals for the estimated parameters.
#' @slot R2 R square statistic.
#' @slot info information. See \code{\link[minpack.lm]{nls.lm}}
#' @export

setClass("psyFitted",representation=representation(psyF="psyFitset",pars="list",Chi="numeric",ci="matrix",R2="numeric",info="numeric"))

#' Creates an instance of a psyFitset object
#' 
#' @param type string expliciting which of three available fitting functions will be used. "CG" for cumulative Gaussian, "LN" for lognormal and "PL" for pseudologistic. 
#' @param xx quantity on the x axis (e.g. duration). Must be a vector of numerics
#' @param pp quantity on the y axis (e.g. proportion of long responses). Must be a vector of numerics AND of the same length as xx.
#' @param getPred Function that returns the predicted value.
#' @param residFun Function that returns the residuals.
#' @param fnJacob Function that returns the jacobian of the function to be fitted.
#' @param fnStart Function that returns the starting parameters.
#' @return A psyFitset object for furture fitting.
#' @examples
#' 
#' data("exampleData",package="nlregVL")
#' setkey(exampleData,Sub)
#' setkey(exampleData,Sub)
#' xx<-exampleData[J("TP01")][,Duration]
#' nbb1<-exampleData[J("TP01"),Nb1]
#' nbb3<-exampleData[J("TP01"),Nb3]
#' pp<-nbb3/(nbb1+nbb3)
#' result<-nlregVL::newpsyFitset(xx = xx,pp = pp)
#' 
#' @seealso \code{\link{psyFitset-class}}
#' @export

newpsyFitset<-function(type="CG",xx,pp,getPred,residFun,fnJacob,fnStart){
  if(missing(getPred)||missing(residFun)||missing(fnJacob)||missing(fnStart)){
    if(!(all(class(type)=="character") && length(type)==1)) stop("mauvais param type")
    switch(type,PL={
      getPred <- function(parS, xx) 1/(1+exp(pi*(parS$U-xx)/(sqrt(3)*parS$S)))
      residFun <- function(p, observed, xx) observed - getPred(p,xx)
      fnjacoo<-function(p, observed, xx){
        mumu<- (exp(pi * (p$U - xx)/(sqrt(3) * p$S)) *
                  (pi/(sqrt(3) * p$U))/(1 + exp(pi * (p$S - xx)/(sqrt(3) * p$S)))^2)
        sisi<- -exp(pi * (p$U - xx)/(sqrt(3) * p$S)) * (pi * (p$U - xx) * sqrt(3)/
                                                          (sqrt(3) * p$S)^2)/(1 + exp(pi * (p$U - xx)/(sqrt(3) * p$S)))^2
        return(c(mumu,sisi))
      }
      parStart <- function(xx) list(U=mean(xx),S=mean(xx)/10)
    }, CG={
      getPred <- function(parS, xx) pnorm(xx,mean=parS$U,sd=parS$S)
      residFun <- function(p, observed, xx) observed - getPred(p,xx)
      fnjacoo<-function(p, observed, xx){
        mumu<-exp(-((xx-p$U)/(p$S*sqrt(2)))^2)/(p$S*sqrt(2*pi))
        sisi<-(xx-p$U)*exp(-((xx-p$U)/(p$S*sqrt(2)))^2)/((p$S^2)*sqrt(2*pi))
        return(c(mumu,sisi))
      }
      parStart <- function(xx) list(U=mean(xx),S=mean(xx)/10)
    }, LN={
      getPred <- function(parS, xx) plnorm(xx, meanlog = parS$U, sdlog = parS$S)
      residFun <- function(p, observed, xx) observed - getPred(p,xx)
      fnjacoo<-function(p, observed, xx){
        mumu<- exp(-((log(xx)-p$U)/(p$S*sqrt(2)))^2)/(p$S*sqrt(2*pi))
        sisi<- (log(xx)-p$U)*exp(-((log(xx)-p$U)/(p$S*sqrt(2)))^2)/((p$S^2)*sqrt(2*pi))
        return(c(mumu,sisi))
      }
      parStart <- function(xx) list(U=mean(xx),S=mean(xx)/10)
    },{
      getPred<-NULL
      residFun<-NULL
      fnjacoo<-NULL
      stop("erreur:mauvais parametre typeF")})
  }
  obj<-new(Class = "psyFitset",getPred=getPred,residFun=residFun,fnJacob=fnjacoo,fnStart=parStart,xx=xx,pp=pp)
  return(obj)
}


#' Applies a nonlinear regression on a specific dataset.
#' @inheritParams newpsyFitset
#' @param control a control list for the non linear regression. See \code{\link[minpack.lm]{nls.lm.control}}
#' @return a \code{\link{psyFitted-class}} object containing the fitted data.
#' @examples
#' data("exampleData",package="nlregVL.R")
#' setkey(exampleData,Sub)
#' xx<-exampleData[J("TP01")][,Duration]
#' nbb1<-exampleData[J("TP01"),Nb1]
#' nbb3<-exampleData[J("TP01"),Nb3]
#' pp<-nbb3/(nbb1+nbb3)
#' result<-nlregVL::doPsyFit(xx = xx,pp = pp)
#' @export
#' 

doPsyFit<-function(xx,pp,type="CG",getPred,residFun,fnJacob,fnStart,control){
  obj<-newpsyFitset(type,xx,pp,getPred,residFun,fnJacob,fnStart)
  parStart<-obj@fnStart(obj@xx)
  ci<-matrix(data = c(0,0,0,0),nrow = 2)
  obj2<-new(Class = "psyFitted",psyF=obj,pars=parStart,Chi=0,ci=ci,R2=0,info=666)
  try({
    for ( w in 1:100) {
      nls.out <- nls.lm(par=parStart, fn = obj@residFun,jac=obj@fnJacob, observed = obj@pp, xx = obj@xx,control=control)
      if(sum(abs(unlist(nls.out$par)-unlist(parStart)))<1e-6) break
      parStart <- nls.out$par
    }
    sstot<-(length(xx)-1)*var(pp)
    chic<-deviance(nls.out)/df.residual(nls.out)
    reserr<-deviance(nls.out)
    Rsq<-1-(reserr/sstot)
    
    try(ci<-confint(nls.out))
    obj2<-new(Class = "psyFitted",psyF=obj,pars=nls.out$par,Chi=chic,ci=ci,R2=Rsq,info=nls.out$info)
  }) 
  
  return(obj2)
}

#' Run the non-linear regression on a data table
#' 
#' @param dataF The data table containing the data to be fitted
#' @param cD name (or position) of the column containing the x-axis data. In the case of a temporal bisection task, this corresponds to the name (or position) of the column containing the duration.
#' @param cNbS name (or position) of the column containing the number of "short" responses.
#' @param CNbL name (or position) of the column containing the number of "long" responses.
#' @param cBy name (or position) of the column(s) to group by. If left empty, it will take all other columns.
#' @param rpsyF Return psyF object as a column?
#' @inheritParams doPsyFit
#' @return A data.table with the parameter estimates, their confidence interval, the R square statistic and the Chi square statistic.
#' @examples
#' data("exampleData",package="nlregVL")
#' result<-nlregVL::fitLSDT(dataF = exampleData,cD = "Duration",cNbS = "Nb1",cNbL = "Nb3",cBy = "Sub")
#' 
#' @seealso \code{\link{psyFitted-class}}
#' @export

fitLSDT<-function(dataF,cD,cNbS,cNbL,cBy,type="CG",control,residFun,getPred,fnJacob,fnStart,rpsyF=T){
  if(!any(class(dataF)=="data.table")) stop("mauvais parametre dataF")
  if(is.numeric(x = cD)){
    namecD<-names(dataF)[cD]
  }else{
    namecD<-cD
    cD<-match(x = namecD,table = names(dataF))
  }
  if(any(is.na(namecD))||length(namecD)!=1) stop("mauvais parametre cD")
  if(is.numeric(cNbS)){
    nameNb1<-names(dataF)[cNbS]
  }else{
    nameNb1<-cNbS
    cNbS<-match(x = nameNb1,table = names(dataF))
  }  
  if(any(is.na(nameNb1))||length(nameNb1)!=1) stop("mauvais parametre CNbS")
  if(is.numeric(cNbL)){
    nameNb3<-names(dataF)[cNbL]
  }else{
    nameNb3<-cNbL
    cNbL<-match(x = nameNb3,table = names(dataF))
  }
  if(any(is.na(nameNb3))||length(nameNb3)!=1) stop("mauvais parametre cNbL")
  if(missing(cBy)){
    namesBy<-names(dataF)[-c(cNbS,cNbL,cD)]
  }else{
    if(is.character(cBy)){
      namesBy<-cBy
      if(length(namesBy)==1){
        namesBy<-gsub(pattern = " ",replacement = "",x = namesBy)
        namesBy<-unlist(strsplit(x = namesBy,split = ","))
      }
      if(!(all(namesBy %in% names(dataF))) || any(namesBy %in% c(namecD,nameNb1,nameNb3))){
        namesBy<-names(dataF)[-c(cNbS,cNbL,cD)]
      }
    }else{
      if(is.numeric(cBy)){
        namesBy<-names(dataF)[cBy]
        if(!(all(namesBy %in% names(dataF))) || any(namesBy %in% c(namecD,nameNb1,nameNb3))){
          namesBy<-names(dataF)[-c(cNbS,cNbL,cD)]
        }
      }
    }
  }
  if(missing(control)){
    control<-nls.lm.control(ftol=1e-9,ptol=1e-9,nprint=0,maxiter=100)
  }
  fstr<-paste0("list(",nameNb1,"=sum(",nameNb1,"),",nameNb3,"=sum(",nameNb3,"))")
  dataFm<-dataF[,eval(parse(text=fstr)),by=c(namecD,namesBy)]
  functStr<-"{www<-doPsyFit("
  if(rpsyF){
    if(missing(getPred)||missing(fnJacob)||missing(fnStart)||missing(residFun)){
      functStr<-paste0(functStr,"xx=",namecD,",pp=",nameNb3,"/(",nameNb3,"+",nameNb1,"),type=\"",type,"\",control=control);list(par1=www@pars[[1]],par1inf=www@ci[1,1],par1sup=www@ci[1,2],par2=www@pars[[2]],par2inf=www@ci[2,1],par2sup=www@ci[2,2],Chi=www@Chi,R2=www@R2,info=www@info,psyF=list(www))}")
    }else{
      functStr<-paste0(functStr,"xx=",namecD,",pp=",nameNb3,"/(",nameNb3,"+",nameNb1,"),type=\"",type,"\",control=control,residFun=residFun,getPred=getPred,fnJacob=fnJacob,fnStart=fnStart);list(par1=www@pars[[1]],par1inf=www@ci[1,1],par1sup=www@ci[1,2],par2=www@pars[[2]],par2inf=www@ci[2,1],par2sup=www@ci[2,2],Chi=www@Chi,R2=www@R2,info=www@info,psyF=list(www))}")
    }
  }else{
    if(missing(getPred)||missing(fnJacob)||missing(fnStart)||missing(residFun)){
      functStr<-paste0(functStr,"xx=",namecD,",pp=",nameNb3,"/(",nameNb3,"+",nameNb1,"),type=\"",type,"\",control=control);list(par1=www@pars[[1]],par1inf=www@ci[1,1],par1sup=www@ci[1,2],par2=www@pars[[2]],par2inf=www@ci[2,1],par2sup=www@ci[2,2],Chi=www@Chi,R2=www@R2,info=www@info}")
    }else{
      functStr<-paste0(functStr,"xx=",namecD,",pp=",nameNb3,"/(",nameNb3,"+",nameNb1,"),type=\"",type,"\",control=control,residFun=residFun,getPred=getPred,fnJacob=fnJacob,fnStart=fnStart);list(par1=www@pars[[1]],par1inf=www@ci[1,1],par1sup=www@ci[1,2],par2=www@pars[[2]],par2inf=www@ci[2,1],par2sup=www@ci[2,2],Chi=www@Chi,R2=www@R2,info=www@info}")
    }
  }
  functExp<-parse(text=functStr)
  dataFR<-dataFm[,eval(functExp),by=namesBy]  
  return(dataFR)
}

graFromFitted<-function(psyFittedObj,title="Blank",nFitPoints=1000,Xlab="Duration (seconds)",Ylab="Proportion of Long Responses",Annotate=T){
  dataP<-psyFittedObj@psyF@xx
  propM<-psyFittedObj@psyF@pp
  dataPdt<-data.table(dataP,propM)
  mIL<-mean((c(dataP,0)-c(0,dataP))[2:length(dataP)])
  
  fxx<-seq(from = dataP[1]-(mIL/2),to = dataP[length(dataP)]+(mIL/2),length.out = nFitPoints)
  fpp<-psyFittedObj@psyF@getPred(psyFittedObj@pars,fxx)
  fxxdt<-data.table(fxx,fpp)
  gra<-ggplot(data = fxxdt,aes(x=fxx,y=fpp))+geom_line()+
    geom_point(data = dataPdt,mapping = aes(x=dataP,y=propM))
  if(Annotate) gra<-gra+annotate("text",x=dataP[1]+(mIL/2),y=.8,label=paste0("R2 = ",round(psyFittedObj@R2,digits = 2),"\nPar 1 =",round(psyFittedObj@pars[[1]],2),"\nPar 2 = ", round(psyFittedObj@pars[[2]],2)),colour="red",size=3)
  gra<-gra+
    scale_x_continuous(breaks=dataP,name=Xlab)+
    scale_y_continuous(breaks=0:5/5,name=Ylab)+
    labs(title=title)+
    coord_cartesian(ylim=c(-.04,1.04))+
    theme_classic(base_size = 10)
  return(list(gra))
}

#' Construct graphs from data tables containing psyFitted objects
#' 
#' @param dataT The data.table containing the psyFitted objects
#' @param byC Name of the columns to contain the factors. There should be no duplicate rows.Can be a character vector or can be a single string with each name separated by commas.
#' @param psyFCol Name of the column containing the psyFitted objects
#' @param Xlab X axis label
#' @param Ylab y axis label
#' @param Annotate If TRUE, graphs will be annotated.
#' @return A data table with the byC columns and a column containing ggplots
#' @examples
#' data("exampleData",package="nlregVL")
#' result<-nlregVL::fitLSDT(dataF = exampleData,cD = "Duration",cNbS = "Nb1",cNbL = "Nb3",cBy = "Sub")
#' result2<-graFromFittedDT(dataT = result,byC = "Sub",psyFCol = "psyF")
#' 
#' @export

graFromFittedDT<-function(dataT,byC,psyFCol,Xlab="Duration (seconds)",Ylab="Proportion of Long Responses",Annotate=T){
  if(length(byC)==1){
    byCS<-byC
    byCL<-gsub(pattern = " ",replacement = "",x = byC)
    byCL<-unlist(strsplit(x = byCL,split = ","))
  }else{
    byCL<-byC
    byCS<-paste0(byC,collapse=",")
  }
  if(!(length(unique(c(byCL,psyFCol)))==length(c(byCL,psyFCol)))) stop("Bad parameters")
  if(!(all(c(byCL,psyFCol) %in% names(dataT)))) stop("Bad parameters")
  DTtest<-dataT[,list(T),by=byCS]
  if(!(length(DTtest[[1]])==length(dataT[[1]]))) stop("Bad parameters")
  unEv<-parse(text = paste0("list(graph=graFromFitted(psyFittedObj=",psyFCol,"[[1]],title=paste0(",paste0("\"",byCL," = \",",byCL,",\" \"",collapse=", "),"),Xlab=Xlab,Ylab=Ylab,Annotate=Annotate))"))
  ee<-dataT[,eval(unEv),by=byCS]
  return(ee)
}



loglikLog<-function(par,Dur,NbS,NbL){
  mu<-par[1]
  sigma<-par[2] 
  PD<-plogis(q = Dur,location = mu,scale = sigma)
  out<-sum(lchoose(n = NbS+NbL,k = NbL)+NbL*log(x = PD)+NbS*log(x = (1-PD)))
  return(out)
}
loglikGau<-function(par,Dur,NbS,NbL){
  mu<-par[1]
  sigma<-par[2]
  PD<-pnorm(q = Dur,mean = mu,sd = sigma)
  out<-sum(lchoose(n = NbS+NbL,k = NbL)+NbL*log(x = PD)+NbS*log(x = (1-PD)))
  return(out)
}

loglikLog1<-function(par,Dur,NbS,NbL){
  mu<-par[1]
  sigma<-par[2]
  PD<-plogis(q = Dur,location = mu,scale = sigma)
  dPDdmu<-(-1)*exp(-((Dur-mu)/sigma))/(sigma*((1+exp(-((Dur-mu)/sigma)))^2))
  dPDdsigma<-(mu-Dur)*exp(-((Dur-mu)/sigma))/((sigma*(1+exp(-((Dur-mu)/sigma))))^2)
  dloglik<-(NbL-NbL*PD-NbS*PD)/(PD*(1-PD))
  dloglikdmu<-sum(dloglik*dPDdmu)
  dloglikdsigma<-sum(dloglik*dPDdsigma)
  return(c(dloglikdmu,dloglikdsigma))
}

loglikGau1<-function(par,Dur,NbS,NbL){
  mu<-par[1]
  sigma<-par[2]
  PD<-pnorm(q = Dur,mean = mu,sd = sigma)
  dPDdmu<-(-1)*exp(-((Dur-mu)/(sigma*sqrt(2)))^2)/(sigma*sqrt(2*pi))
  dPDdsigma<-(mu-Dur)*exp(-((Dur-mu)/(sigma*sqrt(2)))^2)/((sigma^2)*sqrt(2*pi))
  dloglik<-(NbL-NbL*PD-NbS*PD)/(PD*(1-PD))
  dloglikdmu<-sum(dloglik*dPDdmu)
  dloglikdsigma<-sum(dloglik*dPDdsigma)
  return(c(dloglikdmu,dloglikdsigma))
}

fnd2PDdmu2<-deriv(~(-1)*exp(-((Dur-mu)/sigma))/(sigma*((1+exp(-((Dur-mu)/sigma)))^2)),"mu",function.arg=c("mu","sigma","Dur"))
fnd2PDdsigma2<-deriv(~(mu-Dur)*exp(-((Dur-mu)/sigma))/((sigma*(1+exp(-((Dur-mu)/sigma))))^2),"sigma",function.arg=c("mu","sigma","Dur"))

loglikLog2<-function(par,Dur,NbS,NbL){
  mu<-par[1]
  sigma<-par[2]
  PD<-plogis(q = Dur,location = mu,scale = sigma)
  dPDdmu<-(-1)*exp(-((Dur-mu)/sigma))/(sigma*((1+exp(-((Dur-mu)/sigma)))^2))
  dPDdsigma<-(mu-Dur)*exp(-((Dur-mu)/sigma))/((sigma*(1+exp(-((Dur-mu)/sigma))))^2) 
  d2PDdmu2<-attr(fnd2PDdmu2(mu = mu,sigma = sigma,Dur = Dur),which = "gradient")[,1]
  d2PDdsigma2<-attr(fnd2PDdsigma2(mu = mu,sigma = sigma,Dur = Dur),which = "gradient")[,1]
  d2llikpart1<-((NbS+NbL)*(PD*(1-PD))+((NbL-(NbS+NbL)*PD)*(1-2*PD)))*(-1)/((PD*(1-PD))^2)
  d2llikpart2<-(NbL-NbL*PD-NbS*PD)/(PD*(1-PD))
  d2llikdmu2<-sum(d2llikpart1*(dPDdmu^2)+d2llikpart2*d2PDdmu2)
  d2llikdsigma2<-sum(d2llikpart1*(dPDdsigma^2)+d2llikpart2*d2PDdsigma2)
  return(c(d2llikdmu2,d2llikdsigma2))
}

loglikGau2<-function(par,Dur,NbS,NbL){
  mu<-par[1]
  sigma<-par[2]
  PD<-pnorm(q = Dur,mean = mu,sd = sigma)
  dPDdmu<-(-1)*exp(-((Dur-mu)/(sigma*sqrt(2)))^2)/(sigma*sqrt(2*pi))
  dPDdsigma<-(mu-Dur)*exp(-((Dur-mu)/(sigma*sqrt(2)))^2)/((sigma^2)*sqrt(2*pi))
  d2PDdmu2<-(mu-Dur)*exp(-((Dur-mu)/(sigma*sqrt(2)))^2)/((sigma^3)*sqrt(2*pi))
  part1d2dsig<-2+((mu-Dur)/(sigma^2))
  d2PDdsigma2<-part1d2dsig*(Dur-mu)*exp(-((Dur-mu)/(sigma*sqrt(2)))^2)/((sigma^3)*sqrt(2*pi))
  d2llikpart1<-((NbS+NbL)*(PD*(1-PD))+((NbL-(NbS+NbL)*PD)*(1-2*PD)))*(-1)/((PD*(1-PD))^2)
  d2llikpart2<-(NbL-NbL*PD-NbS*PD)/(PD*(1-PD))
  d2llikdmu2<-sum(d2llikpart1*(dPDdmu^2)+d2llikpart2*d2PDdmu2)
  d2llikdsigma2<-sum(d2llikpart1*(dPDdsigma^2)+d2llikpart2*d2PDdsigma2)
  return(c(d2llikdmu2,d2llikdsigma2))
}


fitslik<-function(NbS,NbL,Dur,parStartfn,getPred,loglikfn,loglikgr,loglikgr2,maxN=100,printw=F,method="BFGS",parstart.allow=F,returnFitObj=F,ctrlL=list(fnscale=-1,reltol=.Machine$double.eps)){
  parStart<-parStartfn(Dur)
  ctrlL$fnscale<- -1
  if(parstart.allow){
    if(is.nan(loglikfn(parStart,Dur,NbS,NbL))){
      nn<-NbS+NbL
      nnw<-nn/sum(nn)
      uu<-weighted.mean(x = Dur,w = nnw)
      parStart<-c(uu,uu/5)
    }
  }
  while(is.nan(loglikfn(parStart,Dur,NbS,NbL))){
    uu<-runif(n = 1,min = min(Dur),max=max(Dur))
    vv<-runif(n=1,min(Dur)/5,max=max(Dur)/5)
    parStart<-c(uu,vv)
  }
  gotToMaxN<-F
  bbb<-list(par=rep(NA_real_,times=length(parStart)),value=NA_real_,R2=NA_real_,info=NA_integer_)
  for( w in 1:maxN) {
    aaa<-try(optim(par=parStart,fn=loglikfn,gr=loglikgr,NbS=NbS,NbL=NbL,Dur=Dur,method = method,control=ctrlL))
    if(class(aaa)=="try-error"){
      aaa<-bbb
      break
    }
    if(abs(sum(parStart)-sum(aaa$par))<1e-8) {
      if(printw) cat(paste0(as.character(w),"\n"))
      gotToMaxN<-T
      break
    }
    parStart<-aaa$par
    bbb<-aaa
  }
  if((!gotToMaxN) & printw) cat(paste0(as.character(maxN),"\n"))
  if(!(is.null(getPred))){
    yy<-NbL/(NbL+NbS)
    ypred<-getPred(aaa$par,xx=Dur)
    R2<-1-(sum((yy-ypred)^2)/sum((yy-mean(yy))^2))
  }else{
    R2<-NA_real_
  }
  if(!(is.null(loglikgr2))){
    ll2<-loglikgr2(aaa$par,Dur=Dur,NbS=NbS,NbL=NbL)
  }else{
    ll2<-rep(x = NA_real_,times=length(aaa$par))
  }
  outList<-list()
  for(iii in 1:length(aaa$par)){
    outList[[paste0("par",iii,"_value")]]<-aaa$par[iii]
    outList[[paste0("par",iii,"_se")]]<-1/sqrt(-1*ll2[iii])
  }
  if(returnFitObj){
    outList<-c(outList,list(loglik=aaa$value,R2=R2,info=aaa$convergence,OBJ=aaa))
  }else{
    outList<-c(outList,list(loglik=aaa$value,R2=R2,info=aaa$convergence))
  }
  return(outList)
}

#' Maximum likelihood fits for psychometric functions
#' 
#' @inheritParams fitLSDT
#' @param loglikfn Only used if type is not specified. A function that computes the loglikelihood. It must take the following arguments: par (a vector of numerical values), Dur (a vector of durations), NbS (vector of the number of short respnses), NbL (vector of the number of short respnses). It must return a numerical vector of size 1.
#' @param loglikgr Only used if type is not specified. A function that computes the gradient of the loglikelihood. It must take the following arguments: par (a vector of numerical values), Dur (a vector of durations), NbS (vector of the number of short respnses), NbL (vector of the number of short respnses). It must return a numerical vector of the same size as par.
#' @param loglikgr2 Only used if type is not specified. A function that computes the second order gradient of the loglikelihood. It must take the following arguments: par (a vector of numerical values), Dur (a vector of durations), NbS (vector of the number of short respnses), NbL (vector of the number of short respnses). It must return a numerical vector of the same size as par.
#' @param fnStart A function that calculates the starting values for the fitting process. It must take a single argument (Dur).
#' @param method See \code{\link[stats]{optim}}
#' @param maxN Maximum number of iterations.
#' @param parstart.allow If the starting parameters are bad, allows the algorithm to find new ones.
#' @param allow.print Allows the printing of the number of iterations of the fits.
#' @param ctrlL control list for optim. See \code{\link[stats]{optim}}
#' @return A data.table containing the parameter estimates, the standard error of those estimates, the log likelihood as well as an information code (see \code{\link[stats]{optim}})
#' @export


fitslikDT<-function(dataF,cD,cNbS,cNbL,cBy,type="CG",getPred,loglikfn,loglikgr,loglikgr2,fnStart,method="BFGS",maxN=100,parstart.allow=F,allow.print=F,ctrlL=list(fnscale=-1,reltol=.Machine$double.eps)){
  if(!any(class(dataF)=="data.table")) stop("mauvais parametre dataF")
  if(is.numeric(x = cD)){
    namecD<-names(dataF)[cD]
  }else{
    namecD<-cD
    cD<-match(x = namecD,table = names(dataF))
  }
  if(any(is.na(namecD))||length(namecD)!=1) stop("mauvais parametre cD")
  if(is.numeric(cNbS)){
    nameNb1<-names(dataF)[cNbS]
  }else{
    nameNb1<-cNbS
    cNbS<-match(x = nameNb1,table = names(dataF))
  }  
  if(any(is.na(nameNb1))||length(nameNb1)!=1) stop("mauvais parametre CNbS")
  if(is.numeric(cNbL)){
    nameNb3<-names(dataF)[cNbL]
  }else{
    nameNb3<-cNbL
    cNbL<-match(x = nameNb3,table = names(dataF))
  }
  if(any(is.na(nameNb3))||length(nameNb3)!=1) stop("mauvais parametre cNbL")
  if(missing(cBy)){
    namesBy<-names(dataF)[-c(cNbS,cNbL,cD)]
  }else{
    if(is.character(cBy)){
      namesBy<-cBy
      if(length(namesBy)==1){
        namesBy<-gsub(pattern = " ",replacement = "",x = namesBy)
        namesBy<-unlist(strsplit(x = namesBy,split = ","))
      }
      if(!(all(namesBy %in% names(dataF))) || any(namesBy %in% c(namecD,nameNb1,nameNb3))){
        namesBy<-names(dataF)[-c(cNbS,cNbL,cD)]
      }
    }else{
      if(is.numeric(cBy)){
        namesBy<-names(dataF)[cBy]
        if(!(all(namesBy %in% names(dataF))) || any(namesBy %in% c(namecD,nameNb1,nameNb3))){
          namesBy<-names(dataF)[-c(cNbS,cNbL,cD)]
        }
      }
    }
  }
  switch(type,LOGI={
    loglikfn <- loglikLog
    loglikgr<-loglikLog1
    loglikgr2<-loglikLog2
    getPred<-function(Par,xx) plogis(xx,location = Par[1],scale = Par[2])
    parStartfn <- function(Dur) c(mean(Dur),mean(Dur)/5)
  }, CG={
    loglikfn <- loglikGau
    loglikgr<-loglikGau1
    loglikgr2<-loglikGau2
    getPred <- function(Par, xx) pnorm(xx,mean=Par[1],sd=Par[2])
    parStartfn <- function(Dur) c(mean(Dur),mean(Dur)/5)
  }, {
    if(missing(loglikfn)||missing(loglikgr)||missing(loglikgr2)) stop("Custom likelihood functions not specified and Type argument missing")
    loglikfn<-loglikfn
    loglikgr<-loglikgr
    loglikgr2<-loglikgr2
    if(missing(fnStart)){
      parStartfn <- function(Dur) c(mean(Dur),mean(Dur)/5)
    }else{
      parStartfn <- fnStart
    }
    if(missing(getPred)){
      getPred<-NULL
    }
  })
  setkeyv(dataF,cols = c(namesBy,namecD))[,{
    data.table::setattr(x=list(sum(get(nameNb1)),sum(get(nameNb3))),name = "names",value = c(nameNb1,nameNb3))
  },keyby=c(namesBy,namecD)][,{
    fitslik(NbS = get(nameNb1),NbL = get(nameNb3),Dur = get(namecD),
            parStartfn = parStartfn,getPred = getPred,loglikfn = loglikfn,
            loglikgr=loglikgr,loglikgr2=loglikgr2,method=method,maxN=maxN,
            parstart.allow=parstart.allow,printw=allow.print,returnFitObj=F,ctrlL=ctrlL)
  },keyby=namesBy]
}


#' Add a geom, a theme, etc. to ggplots in a data.table
#' 
#' @param dataT The data.table containing the ggplots
#' @param gg2Str String containing the unevaluated object to add to the ggplots. Can contain dataT column names.
#' @param graphCol name of the column containing the ggplot objects
#' @param byC Column to do the by.
#' 
#' @return A data.table whose graphCol has been updated.
#' @export

add2ggDT<-function(dataT,gg2Str,graphCol,byC){
  if(length(byC)==1){
    byCS<-byC
    byCL<-gsub(pattern = " ",replacement = "",x = byC)
    byCL<-unlist(strsplit(x = byCL,split = ","))
  }else{
    byCL<-byC
    byCS<-paste0(byC,collapse=",")
  }
  if(!(length(unique(c(byCL,graphCol)))==length(c(byCL,graphCol)))) stop("Bad parameters")
  if(!(all(c(byCL,graphCol) %in% names(dataT)))) stop("Bad parameters")
  DTtest<-dataT[,list(T),by=byCS]
  if(!(length(DTtest[[1]])==length(dataT[[1]]))) stop("Bad parameters")
  strstr<-paste0("list(graph=list(",graphCol,"[[1]]+",gg2Str,"))")
  dataT<-dataT[,eval(parse(text=strstr)),by=byCS]
  return(dataT)
}

StatNLFit <- ggproto("StatNLFit", Stat, 
                     required_aes = c("x", "nS","nL","y"),
                     
                     setup_params = function(data, params) {
                       if((!is.null(params$mMin))&(!is.null(params$mMax))) return(params)
                       mMin <- min(data$x) - (max(data$x)-min(data$x))/5
                       mMax <- max(data$x) + (max(data$x)-min(data$x))/5
                       
                       list(
                         n = params$n,
                         mMin = mMin,
                         mMax = mMax,
                         type = params$type,
                         na.rm=params$na.rm
                       )
                     },
                     compute_group = function(data, scales, n = 100, mMin, mMax, type="CG") {
                       xx<-seq(mMin, mMax, length = n)
                       dataT<-data.table(data)
                       parS<- nlregVL::fitslikDT(dataF = dataT[,list(x,nS,nL)],cD = "x",cNbS = "nS",cNbL = "nL",type = type)
                       Par<-c(parS[,par1_value][1],parS[,par2_value][1])
                       switch(type,
                              LOGI=getpPred<-function(Par,xx) plogis(xx,location = Par[1],scale = Par[2]),
                              CG=getPred <- function(Par, xx) pnorm(xx,mean=Par[1],sd=Par[2]),
                              getPred <- function(Par, xx) pnorm(xx,mean=Par[1],sd=Par[2]))
                       yy<-getPred(Par=Par,xx=xx)
                       data.frame(x=xx,y=yy)
                     }
)

#' Draws the fitted psychometric function.
#' 
#' @inheritParams ggplot2::stat_identity
#' @param formula The modelling formula passed to \code{lm}. Should only 
#'   involve \code{y} and \code{x}
#' @param n Number of points used for interpolation.
#' @param type Type of fit
#' @param mMin Minimum of the plotted range
#' @param mMax Maximum of the plotted range
#' @export
 
 
stat_nlreg <- function(mapping = NULL, data = NULL, geom = "line",
                       position = "identity", na.rm = FALSE, show.legend = NA, 
                       inherit.aes = TRUE, n = 100, type="CG",mMin=NULL,mMax=NULL, 
                       ...) {
  layer(
    stat = StatNLFit, data = data, mapping = mapping, geom = geom, 
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(n = n, mMin=mMin,mMax=mMax, type=type, na.rm=na.rm, ...)
  )
}

#' Calculate degrees of visual angle
#' 
#' @param stimX Stimulus width in pixels
#' @param stimY Stimulus height in pixels
#' @param dist Distance from screen to user
#' @param distUnit Measurement unit for dist
#' @param scrX screen width
#' @param scrY Screen height
#' @param scrUnit Measurement unit for scrX and scrY
#' @param resX Screen resolution (horizontal)
#' @param resY Screen resolution (vertical)
#' @return A list with the visual angles in radian and in degrees
#' @export

DVA_VL<-function(stimX,stimY,dist,distUnit="cm",scrX,scrY,scrUnit="mm",resX,resY){
  switch(scrUnit,
         mm={
           scrXmm<-scrX
           scrYmm<-scrY
         },cm={
           scrXmm<-scrX*10
           scrYmm<-scrY*10
         },inch={
           scrXmm<-scrX*10*2.54
           scrYmm<-scrY*10*2.54
         },{
           scrXmm<-scrX
           scrYmm<-scrY})
  switch(distUnit,
         mm={
           distmm<-dist
         },cm={
           distmm<-dist*10
         },m={
           distmm<-dist*1000
         },inch={
           distmm<-dist*10*2.54
         },feet={
           distmm<-dist*10*2.54*12
         },{distmm<-dist*10})
  stimXmm<-stimX*scrXmm/resX
  vaX<-2*atan(x = stimXmm/(2*distmm))
  vaXdeg<-vaX*360/(2*pi)
  stimYmm<-stimY*scrYmm/resY
  vaY<-2*atan(x = stimYmm/(2*distmm))
  vaYdeg<-vaY*360/(2*pi)
  return(list(stimXmm=stimXmm,stimYmm=stimYmm,vaX=vaX,vaXdeg=vaXdeg,vaY=vaY,vaYdeg=vaYdeg))
}


#' @import data.table
#' @import minpack.lm
#' @import ggplot2

NULL