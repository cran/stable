\name{stableglm}
\title{Stable Generalized Linear Regression Model}
\usage{
stableglm(y, Delta=1,
         loc=0, disp=1, skew=0, tail=1.5, 
         iloc=NULL, idisp=NULL,iskew=NULL, itail=NULL,
         oloc=T,odisp=T, oskew=T, otail=T, noopt=F,
         wt=1, exact=F, censor=F,
         integration="Romberg", eps=1e-6, up=10, npoint=501,
         hessian=T, msg=0, ndigit=-1, steptol=1e-12,
         gradtol=0.0001, typf=1, typsize=1, stepmax=1, iterlim=50,
         output=F, llik.output=T)

mode.stable(obj.stable)
fitted.values.stable(obj.stable)
coefficients.stable(obj.stable)
weights.stable(obj.stable)
df.residual.stable(obj.stable)
llik.stable(obj.stable)
llik.comp.stable(obj.stable)
deviance.stable(obj.stable)
deviance.comp.stable(obj.stable)
aic.stable(obj.stable)
aic.opt.stable(obj.stable)
}
\alias{stableglm}
\alias{mode.stable}
\alias{fitted.values.stable}
\alias{coefficients.stable}
\alias{weights.stable}
\alias{df.residual.stable}
\alias{llik.stable}
\alias{llik.comp.stable}
\alias{deviance.stable}
\alias{deviance.comp.stable}
\alias{aic.stable}
\alias{aic.opt.stable}

\arguments{
\item{y}{The response vector.
For censored data, two columns with the second being the censoring
indicator (1: uncensored, 0: right censored, -1: left censored.)}
\item{Delta}{Scalar or vector giving the unit of measurement
for each response value, set to unity by default. For example,
if a response is measured to two decimals, \code{Delta=0.01}.
If the response is transformed, this must be multiplied by
the Jacobian.  For example, with a log transformation,
\code{Delta=1/y}. (The \code{Delta} values for the censored response are
ignored.) The transformation cannot contain unknown
parameters.}
\item{loc, oloc and iloc}{Describe the linear model
considered for the location parameter of the fitted stable distribution,
after transformation by the link function \code{loc.g} (set to the
identity by default. The inverse link function is denoted by \code{loc.h}.
Note that these functions cannot contain unknown parameters).

Two specifications are possible:

(1) \code{loc} is a language expression beginning with ~, describing
the linear part of the regression function for the location parameter
(after transformation by \code{loc.g}, the link function).

\code{iloc} is a vector of initial conditions for the parameters 
in the linear part.

\code{oloc} is a boolean indicating if
an optimization of the likelihood has to be carried out on these parameters.
If \code{oloc} is set to TRUE, a default zero value is considered for the
starting values \code{iloc}. But if no optimization is desired on the location
parameters, i.e. when the likelihood has to be evaluated or optimized
at a fixed location, then \code{iloc} has to be explicitely specified.

(2) \code{loc} is a numeric expression (i.e. a scalar or a vector of the same 
size as the data vector \code{y}, or \code{y[,1]} when censoring is considered).

If \code{oloc} is set to TRUE, i.e. when an optimization of the likelihood has
to be carried out on the location parameter, then the location parameter
(after transformation by the link function loc.g) is set to an unknown
parameter with initial value equal to \code{iloc[1]} or \code{loc[1]}
when \code{iloc} is not specified.

But when \code{oloc} is set to FALSE, i.e. when the likelihood has to
be evaluated  or optimized at a fixed location, then the transformed
location is assumed to be equal to \code{loc} when it is of the same
length as the data vector \code{y} (or \code{y[,1]} when censoring is
considered), and to \code{loc[1]} otherwise.

Specification (1) is especially useful in ANOVA-like
situations where the location is assumed to change with the levels
of some factor variable.
}

\item{disp, odisp and idisp}{describe the linear model
considered for the dispersion parameter of the fitted stable
distribution, after transformation by the link function \code{disp.g}
(set to the \code{log} function by default).
The inverse link function is denoted by \code{disp.h}.
Again these functions cannot contain unknown parameters.
The same rules as above apply when specifying the generalized linear
model for the dispersion parameter.}

\item{skew, oskew and iskew}{describe the linear model considered for
the skewness parameter of the fitted stable distribution, after
transformation by the link function \code{skew.g} (set to 
\code{log{(1 + [.])/(1 - [.])}} by default). The inverse link function
is denoted by \code{skew.h}. Again these functions cannot contain unknown
parameters.
The same rules as above apply when specifying the generalized linear
model for the skewness parameter.}

\item{tail, otail and itail}{describe the linear model considered for
the tail parameter of the fitted stable distribution, after
transformation by the link function \code{tail.g} (set to 
\code{log{([.] - 1)/(2 - [.])}} by default. The inverse link function
is denoted by \code{tail.h}. Again these functions cannot contain
unknown parameters). The same rules as above apply when specifying
the generalized linear model for the tail parameter.}
\item{noopt}{When set to TRUE, it forces \code{oloc},
\code{odisp}, \code{oskew} and \code{otail} to FALSE, whatever the user choice
for these last three arguments. It is especially useful when looking 
for appropriate initial values for the linear model parameters, before
undertaking the optimization of the likelihood.}
\item{wt}{Weight vector.}
\item{exact}{If TRUE, fits the exact likelihood function for
continuous data by integration over intervals of observation,
i.e. interval censoring.}
\item{censor}{TRUE if data are censored, in which case \code{y} must be
a two column matrix.}
\item{integration, eps, up and npoint}{\code{integration} indicates
which algorithm must be used to evaluate the stable density when
the likelihood is computed with \code{exact} set to FALSE.
See the man page on \code{dstable} for extra information.}
\item{llik.output}{is TRUE when the likelihood has to be displayed 
at each iteration of the optimization.}
\item{others}{Arguments controlling the optimization procedure \code{nlm}.}
\item{obj.stable}{Object of class stable, as produced by \code{stableglm}}
}


\description{
\code{stableglm} fits user specified generalized linear regression
models based on the stable distribution to (uncensored, right and/or
left censored) data. This allows the location, the dispersion, the skewness
and the tails of the fitted stable distribution to vary with explanatory
variables. 
\code{mode.stable} gives the mode of the fitted stable distribution for each
component of \code{y}. \cr
\code{fitted.values.stable} gives the fitted values, i.e. the value of the 
location parameter \code{loc} of the fitted stable distribution for each
component of \code{y}. \cr
\code{coefficients.stable} provides the MLEs of the parameters in the 
corresponding model. \cr
\code{weights.stable} gives the weight vector specified by \code{wt}. \cr
\code{df.residual.stable} gives the residual degree of freedom. \cr
\code{llik.stable(obj.stable)}
\code{llik.comp.stable(obj.stable)}
\code{deviance.stable(obj.stable)}
\code{deviance.comp.stable(obj.stable)}
\code{aic.stable(obj.stable)}
\code{aic.opt.stable(obj.stable)}



\value{
A list of class \code{sable} is returned.
The printed output includes the -log-likelihood, the corresponding AIC,
the maximum likelihood estimates, standard errors, and correlations.
It also include  all the relevant information calculated, including error codes.
}

\reference{
Lambert, P. and Lindsey, J.K. (1998) Generalized regression models
for heavy-tailed processes based on non-symmetric stable
distributions: a likelihood approach. (To appear in Applied Statistics).
}

\author{
Philippe Lambert (University of Liege, Belgium, \email{plambert@ulg.ac.be})
and Jim Lindsey.
}

\seealso{
\code{\link{lm}}, \code{\link{glm}}, \code{\link{stable}} and 
\code{\link{stable.mode}}.
}
\examples{
library(stable)
## Share return over a 50 day period (see reference above)
# share
y <- c(296,296,300,302,300,304,303,299,293,294,294,293,295,287,288,297,305,
307,307,304,303,304,304,309,309,309,307,306,304,300,296,301,298,295,
295,293,292,297,294,293,306,303,301,303,308,305,302,301,297,299)  

#return
ret <- (y[2:50]-y[1:49])/y[1:49]
hist(ret, breaks=seq(-0.035,0.045,0.01))

day <- seq(0,0.48,by=0.01) # time measured in days/100
x <- seq(1,length(ret))-1

# Classic stationary normal model tail=2
print(z1  <- stableglm(y = ret, Delta = 1/y[1:49],
     loc = ret ~ 1, disp= ret ~ 1, skew = ret ~ 1, tail = tail.g(1.9999999),
     iloc = 0,idisp = 0, iskew = 0,
     otail = F, oskew = F))

# Normal model (tail=2) with dispersion=disp.h(b0+b1*day)
print(z2 <- stableglm(y = ret,Delta = 1/y[1:49],
      loc = ~ day,disp = ret ~ 1, skew = ret ~ 1,tail = tail.g(1.999999),
      iloc = c(0,0),idisp = 0,iskew = 0,
      oskew = F,otail =F))

# Stable model with loc(ation)=loc.h(b0+b1*day)
print(z3 <- stableglm(y = ret, Delta = 1/y[1:49],
      loc = ~ day, disp = ret ~ 1, skew = ret ~ 1, tail = ret ~ 1,
      iloc = c(0.001,0),idisp = 0,iskew = 0, itail = 0))

# Stable model with disp(ersion)=disp.h(b0+b1*day)
print(z4 <- stableglm(y = ret, Delta = 1/y[1:49],
      loc = ret ~ 1, disp = ret ~ day, skew = ret ~ 1, tail = ret ~ 1,
      iloc = 0,idisp = c(-4.5,0),iskew = -2, itail = 1))

# Stable model with skew(ness)=skew.h(b0+b1*day)
# Evaluation at fixed parameter values (since noopt is set to TRUE)
print(z5 <- stableglm(y = ret, Delta = 1/y[1:49],
      loc = ret ~ 1, disp = ret ~ 1, skew = ~ day, tail = ret ~ 1,
      iloc = 5.557e-04,idisp = -4.957,iskew = c(2.811,-2.158),itail = -5.261e-1,
      noopt=T))

# Stable model with tail=tail.h(b0+b1*day)
print(z6 <- stableglm(y = ret, Delta = 1/y[1:49],
      loc = ret ~ 1, disp = ret ~ 1, skew = ret ~ 1, tail = ret ~ day,
      iloc = 0.001,idisp = -5,iskew = -3, itail = c(2,-7),hessian=F))
}
