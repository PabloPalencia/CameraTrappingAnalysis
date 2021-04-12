### REM FUNCTIONS (not included in packages. Developed by M. Rowcliffe & St. Andrews University)
## 12/04/2021

TRD <- function(P, T, param, strata=NULL, areas=NULL){
  if(length(P)!=length(T)) stop("P and T have unequal lengths")
  if(!("g" %in% names(param))) param <- c(param,g=1)
  if(!("p" %in% names(param))) param <- c(param,p=1)
  if(is.null(strata))
  {	res <- pi * param$g * sum(P) / (sum(T) * param$DR * param$r * (2+param$theta))
  names(res) <- NULL
  } else
  {	if(is.null(areas)) stop("areas are missing")
    if(length(strata)!=length(P)) stop("strata vector is a different length to P/T")
    if(sum(names(areas) %in% levels(strata)) != length(names(areas)) |
       sum(levels(strata) %in% names(areas)) != length(levels(strata)))
      stop("strata levels do not match areas names")
    nstrata <- length(areas)
    areas <- unlist(areas[order(names(areas))])
    locdens <- sapply(1:nstrata, STRD, P, T, param, strata)
    res <- sum(locdens * areas) / sum(areas)
  }
  res
} # Modified to work with DR instead of  speed and activity
TRDsample <- function(i, P, T, param, strata=NULL, areas=NULL){
  if(is.null(strata)){
    x <- sample(1:length(T), replace=TRUE)
  } else{
    nstrata <- length(areas)
    x <- NULL
    for (i in 1:nstrata) x <- c(x, sample(which(strata==levels(strata)[i]), replace=TRUE))
    strata <- sort(strata)
  }
  TRD(P[x], T[x], param, strata, areas)
}
bootTRD <- function(P, T, param, paramSE, strata=NULL, areas=NULL, its=1000){
  BSdens <- sapply(1:its, TRDsample, P, T, param, strata, areas)
  BSse <- sd(BSdens)
  Dens <- TRD(P,T,param,strata,areas)
  prms <- length(param)
  addn <- rep(0,prms)
  addn[which(names(param)=="theta")] <- 2
  Es <- c(Dens,addn+unlist(param))
  SEs <- c(BSse,unlist(paramSE))
  SE <- Dens * sqrt(sum((SEs/Es)^2))
  cbind(Density=Dens, SE=SE)
  
}

EDRtransform <- function(dsobject, alpha=0.05) {
  if(class(dsobject) != "dsmodel") stop("First argument must be a dsmodel object")
  if(!dsobject$ddf$meta.data$point) stop("EDR can only be computed for point transect data")
  summary.ds.model <- summary(dsobject)
  p_a <- summary.ds.model$ds$average.p
  se.p_a <- summary.ds.model$ds$average.p.se
  cv.p_a <- se.p_a / p_a
  w <- summary.ds.model$ds$width
  edr <- sqrt(p_a * w^2)
  se.edr <- cv.p_a/2 * edr
  degfree <- summary.ds.model$ds$n - length(summary.ds.model$ddf$par)
  t.crit <- qt(1 - alpha/2, degfree)
  se.log.edr <- sqrt(log(1 + (cv.p_a/2)^2))
  c.mult <- exp(t.crit * se.log.edr)
  ci.edr <- c(edr / c.mult, edr * c.mult)
  return(list(EDR=edr, se.EDR=se.edr, ci.EDR=ci.edr))
}
