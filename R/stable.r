library.dynam("stable.so")

###########################################################################
# Tools for using the Wilkinson-Rogers notation
wr <- function(formula){
	mt <- terms(formula)
	data <- sys.frame(sys.parent())
	mf <- model.frame(mt, data, na.action=na.fail)
	x <- model.matrix(mt, mf)
	y <- model.response(mf, "numeric")
	z <- list(y,x)
	names(z) <- c("response","design")
	z}

###########################################################################
# Computation of the mode ytilde as a function of (loc,disp,skew,tail)

stable.mode <- function(loc,disp,skew,tail){
 if (tail <1) warning("stable.mode is only reliable for tail in (1,2)") 
 coef1 <- 1.7665114+1.8417675*tail-2.2954390*tail^2+0.4666749*tail^3
 coef2 <- -0.003142967+632.4715*exp(-7.106035*tail)*tail
 return(list(ytilde=loc+disp*coef1*exp(-coef2*abs(skew))*skew,
              c1=coef1,c2=coef2))
}

###########################################################################
# Density of a stable distribution
#
# This function produces the stable density computed at y.
# Note that it is obtained by integrating a function from 0 to Infinity.
# This integral is approximated by the finite integral from 0 to UP
# using the Simpson's method with npt points or Romberg's integration
dstable <- function(y, loc=0, disp=1/sqrt(2), skew=0, tail=2,
		    npt=501, up=10, eps=1.0e-6, integration="Romberg"){
  type <- ifelse(integration=="Simpson",1,2)
  ly <- length(y)
  z0 <- rep(0,ly)
  skew <- skew+z0
  tail <- tail+z0
  yy <- (y-loc)/disp
  z <- .C("stable",
	  as.integer(ly),
	  yy,
	  skew,
	  tail, 
	  as.integer(npt),
	  up,
	  eps,
	  as.integer(type),
	  err=integer(1),
	  ffy=double(ly))
  z$ffy/disp}

# R version of dstable (written in C)
dstable2 <- function(y, loc, disp, skew, tail, npt=501, up=10,
		     integration="Romberg"){
  ffy <- 0*y
  eta <- skew*(1-abs(1-tail))*pi/2
  seta <- sin(eta)
  ceta <- cos(eta)
  yy <- (y-loc)/disp
  if(integration=="Simpson"){
    npt <- npt-npt%%2
    h <- up/npt
    for (i in 0:npt){
      s <- (npt-i)*h
      sa <- s^tail
      ffy <- ffy+(4-2*(i%%2==0)-(i==1||i==npt))*
	cos(-yy*s+sa*seta)*exp(-sa*ceta)}
    ffy <- ffy*h/3/pi}
  else {
    for (i in 1:length(y)) {
      g <- function(s) {
	sa <- s^tail[i]
	cos(-yy[i]*s+sa*seta[i])*exp(-sa*ceta[i])}
      ffy[i] <- int(g,0,"infty")}
    ffy <- ffy/pi}
  ffy/disp}

###########################################################################
# c.d.f of a stable distribution
#
pstable <- function(y, loc=0,disp=1/sqrt(2),skew=0,tail=2,eps=1.0e-6){
  yy <- (y-loc)/disp
  ly <- length(yy)
  z0 <- rep(0,ly)
  skew <- skew+z0
  tail <- tail+z0 
  eta <- skew*(1-abs(1-tail))*pi/2
#        if ((yy==0)&(eta==0)) return(0.5+z0)
  z <- .C("pstable",
	  as.integer(ly),
	  yy,
	  skew,
	  tail, 
	  eps,
	  err=integer(1),
	  ffy=double(ly))
  z$ffy}

###########################################################################
# Quantile function of a stable random variable
#
uniroot2 <- function (f, interval, lower = min(interval), upper = max(interval), 
		      tol = .Machine$double.eps^0.25, ...) 
{
  while ((f(interval[1], ...) * f(interval[2], ...) >= 0)) {
    interval <- 1.5*interval}
  val <- .Internal(zeroin(function(arg) f(arg, ...), lower, 
			  upper, tol))
  list(root = val, f.root = f(val, ...))
}

qstable <- function(q, loc=0, disp=1/sqrt(2), skew=0, tail=2, eps=1.0e-6){
  nn <- length(q);zero.nn <- rep(0,nn)
  loc <- zero.nn+loc ; disp <- zero.nn+disp
  skew <- zero.nn+skew ; tail <- zero.nn+tail  
  tamp <- rep(0,nn)
  for (i in (1:nn)){
    h <- function(x) {
      ttamp <- pstable(x, loc[i],disp[i],skew[i],tail[i],eps=eps)
      return(ttamp-q[i])}
    z0 <- uniroot2(h,c(-20,20))
    tamp[i] <- z0$root}
  return(tamp)}

###########################################################################
# Generation of stable random deviates
#
rstable <- function(n=1,loc=0,disp=1/sqrt(2),skew=0,tail=2,eps=1.0e-6){
 return(qstable(runif(n),loc=loc,disp=disp,skew=skew,tail=tail,eps=eps))}

###########################################################################
# Stable hazard
#
hstable <- function(y, loc=0,disp=1/sqrt(2),skew=0,tail=2,eps=1.0e-6){
 return(dstable(y,loc=loc,disp=disp,skew=skew,tail=tail,eps=eps)/
         (1-pstable(y,loc=loc,disp=disp,skew=skew,tail=tail,eps=eps))
       )}

###########################################################################
# Link and inverse link functions for use in stableglm

loc.g <- function(x) x # link function for loc
loc.h <- function(x) x # inverse link function for disp
disp.g <- function(x) log(x) # link function for disp
disp.h <- function(x) exp(x) # inverse link function for disp
skew.g <- function(x) log((1+x)/(1-x)) # link function for skew 
skew.h <- function(x) 2/(1+exp(-x))-1  # inverse link function for skew
tail.g <- function(x) log((x-1)/(2-x)) # link function for tail in (1,2)
tail.h <- function(x) 1/(1+exp(-x))+1  # inverse link function for tail
#tail.g <- function(x) log(x/(2-x)) # link function for tail in (0,2)
#tail.h <- function(x) 2/(1+exp(-x))  # inverse link function for tail

###########################################################################
# Regression models for the four parameters in the stable distribution 
# Note that the returned optimized parameters are on the normalized
# scale, i.e. that they have to be transformed back to the right scale.

stableglm <- function(y, Delta=1, loc=0, disp=1, skew=0, tail=1.5, 
		     oloc=T,odisp=T, oskew=T, otail=T, noopt=F,
		     iloc=NULL, idisp=NULL,iskew=NULL, itail=NULL,
		     wt=1, exact=F, censor=F,
		     integration="Romberg", eps=1e-6, up=10, npoint=501,
                     hessian=T,msg=1, ndigit=-1, steptol=1e-12,
                     gradtol=0.0001, typf=1, typsize=1, stepmax=1, iterlim=100,
                     output=F, llik.output=F){
  call <- sys.call()
  
  if(censor)
    if(length(dim(y))!=2||dim(y)[2]!=2){
      print("Two column matrix required for response")
      print("Times and censor indicator")
      return(22)}
    else n <- dim(y)[1]
  else {
    if(!is.vector(y)){
      print("y must be a vector")
      return(22)}
    n <- length(y)}
  if(censor){
    y[,2] <- as.integer(y[,2])
    if(y[,2]!=-1&y[,2]!=0&y[,2]!=1){
      print("Censor indicator must be -1s, 0s, and 1s")
      return(23)}
    cc <- ifelse(y[,2]==1,1,0)  # observed
    rc <- ifelse(y[,2]==0,1,ifelse(y[,2]==-1,-1,0))  # right censored
    lc <- ifelse(y[,2]==-1,0,1) # left censored
    if(Delta<=0&y[,2]==1){
      print("All deltas for uncensored data must be positive")
      return(25)}}
#	else {
#		Delta <- ifelse(Delta<=0,0.000001,Delta)
#		Delta <- ifelse(y[,1]-Delta/2<=0,Delta-0.00001,Delta)}}
  else {
    if(min(Delta)<=0){
      print("All deltas for must be positive")
      return(25)}}

# cc  is 1 if y[,2]=1  i.e. if observed
#        0 otherwise   i.e. if censoring
#
# rc is  1 if y[,2]=0  i.e. if right censoring
#       -1 if y[,2]=-1 i.e. if left censoring
#        0 if y[,2]=1  i.e. if observed
#
# lc is  0 if y[,2]=-1  i.e. if left censoring
#        1 otherwise    i.e. if right censoring or observed

  if (censor) y <- y[,1]
  
# If don't want to optimize at all, just set NOOPT=T
  if (noopt){
    oloc <- F
    odisp <- F
    oskew <- F
    otail <- F}
  
  if (length(wt)==1)wt <- rep(wt,n)

# LOC
# From the coming lines, we see that if one does not want to optimize over
# loc (ie. OLOC=F), then there are 2 alternatives: 
# (1) set loc equal to a value identical for all units: use LOC=<scalar>
# (2) set loc equal to values varying through units: use LOC=<language>
#      and ILOC=<corresponding vector of initial conditions>
# Note that we work with the log.g link (identity link by default).
# This means that language description of the loc model and
# initial conditions are understood on that scale!!
  if (is.language(loc)){
    mf <- model.frame(terms(loc),sys.frame(sys.parent()),na.action=na.fail)
    wrloc <- model.matrix(terms(loc), mf)
    nploc <- dim(wrloc)[2]
    if (oloc) {
      if(is.numeric(iloc) & length(iloc)!=nploc){
	cat("iloc must be of size ",nploc,"\n")
	return(1)}
      if (!is.numeric(iloc)) iloc <- rep(0,nploc)
      fnloc <- function(p) loc.h(wrloc %*% p)}
    else {
      if (!is.numeric(iloc)){
	print("Missing initial conditions for loc")
	return(2)}
    else if (length(iloc)!=nploc){
      cat("iloc must be of size ",nploc,"\n")
      return(3)}
    else fnloc <- function(p) loc.h(wrloc %*% iloc)}}
  else if (oloc){
    fnloc <- function(p) loc.h(rep(p[1],n))
    if (!is.numeric(iloc)){
      iloc <- loc[1]
      nploc <- 1}}
  else {
                                        # IMPORTANT
    if (length(loc)==n)fnloc <- function(p) loc.h(loc)
                                        # IMPORTANT
    else fnloc <- function(p) rep(loc.h(loc[1]),n)
    nploc <- 1}

# DISP
# Note that we work with the disp.g link (log link by default).
# This means that language description of the disp model and
# initial conditions are understood on that scale!!
# Non language description of disp are also understood on the disp.g scale 
# and specified using disp=<value>, yielding disp.h(<value>) for the 
# parameter disp.

if(is.language(disp)){
  mf <- model.frame(terms(disp),sys.frame(sys.parent()),na.action=na.fail)
  wrdisp <- model.matrix(terms(disp), mf)	
  npdisp <- dim(wrdisp)[2]
  if (odisp) {
    if(is.numeric(idisp) & length(idisp)!=npdisp){
      cat("idisp must be of size ",npdisp,"\n")
      return(4)}
    else if (!is.numeric(idisp)) idisp <- rep(0,npdisp)
    fndisp <- function(p) disp.h(wrdisp %*% p)}
  else {
    if (!is.numeric(idisp)){
      print("Missing initial conditions for disp")
      return(5)}
    else if(length(idisp)!=npdisp){
      cat("idisp must be of size ",npdisp,"\n")
      return(6)}
    else fndisp <- function(p) disp.h(wrdisp %*% idisp)}}
else if (odisp){
  fndisp <- function(p) disp.h(rep(p[1],n))
  if (!is.numeric(idisp)){
    idisp <- disp[1]
    npdisp <- 1}
}
else {
					# IMPORTANT
  if (length(disp)==n)fndisp <- function(p) disp.h(disp)
					# IMPORTANT
  else fndisp <- function(p) rep(disp.h(disp[1]),n)
  npdisp <- 1}

# SKEW
# Note that, y default, we work with the skew.g([.])=log{(1+[.])/(1-[.])} link.
# This means that language description of the skew model and
# the possible initial conditions are understood on that scale!!
# Non language description of skew are also understood on the skew.g scale 
# and specified using skew=<value>, yielding skew.h(<value>) for the 
# parameter skew.

  if(is.language(skew)){ 
    mf <- model.frame(terms(skew),sys.frame(sys.parent()),na.action=na.fail)
    wrskew <- model.matrix(terms(skew), mf)
    npskew <- dim(wrskew)[2]
    if (oskew){
      if(is.numeric(iskew) & length(iskew)!=npskew){
	cat("iskew must be of size ",npskew,"\n")
	return(8)}
      else if (!is.numeric(iskew)) iskew <- rep(0,npskew)
      fnskew <- function(p) skew.h(wrskew %*% p)}
    else {
      if (!is.numeric(iskew)) {
	print("Missing initial conditions for skew")
	return(9)}
      else if (length(iskew)!=npskew){
	cat("iskew must be of size ",npskew,"\n")
	return(10)}
      else fnskew <- function(p) skew.h(wrskew %*% iskew)}}
  else if (oskew){
    fnskew <- function(p) skew.h(rep(p[1],n))
    if (!is.numeric(iskew)){
      iskew <- skew[1]
      npskew <- 1}
  }
  else {
    if (length(skew)==n) fnskew <- function(p) skew.h(skew)
    else fnskew <- function(p) rep(skew.h(skew[1]),n) # IMPORTANT
    npskew <- 1}

# TAIL
# Note that we work with the tail.g([.])=log{([.]-1)/(2-[.])} link.
# This means that language description of the tail model and
# the possible initial conditions are understood on that scale!!
# Non language description of tail are also understood on the tail.g scale 
# and specified using tail=<value>, yielding tail.h(<value>) for the 
# parameter tail.

  if(is.language(tail)){ 
    mf <- model.frame(terms(tail),sys.frame(sys.parent()),na.action=na.fail)
    wrtail <- model.matrix(terms(tail), mf)
    nptail <- dim(wrtail)[2]
    if (otail){
      if (is.numeric(itail) & length(itail)!=nptail){
	cat("itail must be of size ",nptail,"\n")
	return(12)}
      else if (!is.numeric(itail)) itail <- rep(0,nptail)
      fntail <- function(p) tail.h(wrtail %*% p)}
    else {
		if (!is.numeric(itail)){
		  print("Missing initial conditions for tail")
		  return(13)}
		else if (length(itail)!=nptail){
		  cat("itail must be of size ",nptail,"\n")
		  return(14)}
		else fntail <- function(p) tail.h(wrtail%*%itail)}}
  else if (otail){
    fntail <- function(p) tail.h(rep(p[1],n))
    if (!is.numeric(itail)){
      itail <- tail[1]
      nptail <- 1}
  }
  else {
    if (length(tail)==n)fntail <- function(p) tail.h(tail)
    else fntail <- function(p) rep(tail.h(tail[1]),n) # IMPORTANT
    nptail <- 1}

# Computation of -log-likelihood
  llikstable <- function(p){  # ,up=up,integration=integration) {
    i1 <- 1
    if (oloc) {
      i2 <- i1+nploc
      loc <- fnloc(p[i1:(i2-1)])
      i1 <- i2}
    else loc <- fnloc()
    if (odisp) {
      i2 <- i1+npdisp
      disp <- fndisp(p[i1:(i2-1)])
      i1 <- i2}
    else disp <- fndisp()
    if (oskew) {
      i2 <- i1+npskew
      skew <- fnskew(p[i1:(i2-1)])
      i1 <- i2}
    else skew <- fnskew()
    if (otail) {
      i2 <- i1+nptail
      tail <- fntail(p[i1:(i2-1)])}
    else tail <- fntail()
    
    if (!censor){
      if (exact) {
	tamp <- pstable(y=y+Delta/2,loc=loc,disp=disp,
			skew=skew,tail=tail)-
			  pstable(y=y-Delta/2,loc=loc,disp=disp, 
				  skew=skew,tail=tail)
	llikcomp <- -(log(tamp))*wt}
      else {tamp <- dstable(y=y,loc=loc,disp=disp,
			   skew=skew,tail=tail,
			   npt=npoint, up=up, integration=integration)
	    llikcomp <- -(log(tamp)+log(Delta))*wt}}
    else {
      if (exact) {
	p1 <- pstable(y=y+Delta/2,loc=loc,disp=disp,
		      skew=skew,tail=tail)
	ps <- pstable(y=y-Delta/2,loc=loc,disp=disp,
		      skew=skew,tail=tail)
	llikcomp <- -wt*(cc*log(p1-ps)+log(lc-rc*ps))}
      else {tamp <- dstable(y=y,loc=loc,disp=disp,
			   skew=skew,tail=tail,
			   npt=npoint,up=up, integration=integration)
	    llikcomp <- -wt*(cc*(log(tamp)+log(Delta))
			    +log(lc-rc*pstable(y=y-Delta/2,loc=loc, 
					       disp=disp,skew=skew,
					       tail=tail)))}}
    
    llik <- sum(llikcomp)
    if (output){
      if (length(p)==0) cat("-LogLik: ",sum(llikcomp),"\n")
      else cat("-LogLik: ",sum(llikcomp)," ",p,"\n")}
    z <- list(
	      llik=llik,
	      llikcomp=llikcomp,
	      loc=loc,
	      disp=disp,
	      skew=skew,
	      tail=tail)
    return(z)}

# This function returns -llik (and NOT the deviance) to get the 
# correct s.e's with the hessian returned by nlm (optimization routine).
optstable <- function(p){
  tamp <- llikstable(p)$llik
  if (llik.output) cat("-LogLik: ",tamp," (",p,")","\n")
  return(ifelse(is.na(tamp),1e20,tamp))}

# initial conditions
  p0 <- c()
  if (oloc) {
    p0 <- c(iloc)
    names(p0) <- c(rep("iloc",length(iloc)))}
  if (odisp) {
    tamp <- names(p0) 
    p0 <- c(p0,idisp)
    names(p0) <- c(tamp,rep("idisp",length(idisp)))}
  if (oskew) {
    tamp <- names(p0) 
    p0 <- c(p0,iskew)
    names(p0) <- c(tamp,rep("iskew",length(iskew)))}
  if (otail) {
    tamp <- names(p0) 
    p0 <- c(p0,itail)
    names(p0) <- c(tamp,rep("itail",length(itail)))}
  if (output) {
    cat("No. of parameters: ",nploc,"",npdisp,"",npskew,"",nptail,"\n")
    if(oloc||odisp||oskew||otail){
      cat("Vector of initial conditions on IR^p:","\n")
      print(p0)}}

# -Log-likelihood at p0
  llik0 <- llikstable(p=p0)
  np0 <- length(p0)

  if (np0>0){
    p.opt <- nlm(optstable, p=p0, hessian=hessian, fscale=typf,
		 typsize=rep(1,length(p0)), print.level=msg, ndigit=ndigit,
		 gradtol=gradtol, steptol=steptol, iterlim=iterlim,
		 stepmax=stepmax)
    z <- llikstable(p.opt$estimate)}
  else z <- llikstable(p0)
  
  ytilde.tamp <- stable.mode(z$loc,z$disp,z$skew,z$tail)$ytilde  
  # corresponding mode 
  tamp <- dstable(y=ytilde.tamp, loc=z$loc, disp=z$disp, skew=z$skew,
		 tail=z$tail, npt=npoint, up=up, integration=integration)
  llik.ytilde <- -(log(tamp)+log(Delta))*wt

  ytilde.tamp <- stable.mode(z$loc,z$disp,z$skew,z$tail)$ytilde  

  np <- nploc+npdisp+npskew+nptail
  llik <- z$llik
  aic <- llik+np
  aic.opt <- llik+np0
  nobs <- sum(as.numeric(wt))
  df <- nobs-np0
  loc <- as.vector(z$loc)
  disp <- as.vector(z$disp)
  skew <- as.vector(z$skew)
  tail <- as.vector(z$tail)
  ytilde <- stable.mode(loc,disp,skew,tail)$ytilde  # corresponding mode 
  residuals <- as.vector((y-loc)/disp)

  if (np0>0){
    cov <- diag(np0)
    if(hessian){if(np0==1)cov <- 1/p.opt$hessian
    else cov <- solve(p.opt$hessian)}
    se <- sqrt(diag(cov))
    z1 <- list(
	       call=call,
	       y=y,
	       loc=loc,
	       disp=disp,
	       skew=skew,
	       tail=tail,
	       ytilde=ytilde, 
	       start=p0,
	       weights=wt,
	       llik=llik,
	       llikcomp=z$llikcomp,
	       residuals=residuals,
	       aic=aic,
	       aic.opt=aic.opt,
	       npar=np,
	       npar.opt=np0,
	       n=nobs,
	       df=df,
	       coefficients=p.opt$estimate,
	       se=se,
	       cov=cov,
	       corr=cov/(se%o%se),
	       gradient=p.opt$gradient,
	       error=p.opt$error,
	       code=p.opt$code,
	       integration=integration)}
  else z1 <- list(
		  call=call,
		  y=y,
		  loc=loc,
		  disp=disp,
		  skew=skew,
		  tail=tail,
		  ytilde=ytilde, 
		  start=p0,
		  weights=wt,
		  llik=llik,
		  llikcomp=z$llikcomp,
		  residuals=residuals,
		  aic=aic,
		  aic.opt=aic.opt,
		  npar=np,
		  npar.opt=np0,
		  n=nobs,
		  df=df,
		  integration=integration)

  class(z1) <- "stable"
  return(z1)}

stablelm <- function(...) stableglm(...)

###########################################################################
# Stable class functions
#
residuals.stable <- function(z) z$residuals
fitted.values.stable <- function(z) z$loc
coefficients.stable <- function(z) z$coefficients
weights.stable <- function(z) z$weights
df.residual.stable <- function(z) z$df
llik.stable <- function(z) z$llik
llik.comp.stable <- function(z) z$llikcomp
deviance.stable <- function(z) 2*z$llik
deviance.comp.stable <- function(z) 2*z$llikcomp
aic.opt.stable <- function(z) z$aic.opt
aic.stable <- function(z) z$aic
mode.stable <- function(z) z$ytilde

print.stable <- function(z) {
  np <- z$npar.opt
  cat("\nCall:\n",deparse(z$call),"\n\n",sep="")
  cat("-Log likelihood             ",z$llik,"\n")
  cat("No. of obs                  ",z$n,"\n")
  cat("No. of estimated parameters ",np,"\n")
  cat("No. of parameters           ",z$npar,"\n")
  cat("Degrees of freedom          ",z$df,"\n")
  cat("Optimized AIC               ",z$aic.opt,"\n")
  cat("AIC                         ",z$aic,"\n\n")
  cat("CAREFUL: AIC = -Log-likelihood + No.parameters","\n\n")

  if (np>0){
    cat("Coefficients:\n")
    coef.table <- cbind(z$coefficients, z$se)
    dimnames(coef.table) <- list(seq(1,length(z$coefficients)),
				 c("estimate", "se"))
    print.default(coef.table, digits=4, print.gap=2)}
  if (np>1){
    cat("\nCorrelations:\n")
    dimnames(z$corr) <- list(seq(1,np),seq(1,np))
    print.default(z$corr, digits=4)}
  invisible(z)}
