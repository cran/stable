\name{stable}
\title{The Stable Distribution}
\usage{
dstable(x, loc=0, disp=1/sqrt(2), skew=0, tail=2, eps=1e-6)
pstable(q, loc=0, disp=1/sqrt(2), skew=0, tail=2, eps=1e-6)
qstable(p, loc=0, disp=1/sqrt(2), skew=0, tail=2, eps=1e-6)
hstable(x, loc=0, disp=1/sqrt(2), skew=0, tail=2, eps=1e-6)
rstable(n, loc=0, disp=1/sqrt(2), skew=0, tail=2, eps=1e-6)
}
\alias(stable)
\alias{dstable}
\alias{pstable}
\alias{qstable}
\alias{hstable}
\alias{rstable}
\arguments{
\item{y,q}{vector of quantiles.}
\item{p}{vector of probabilites.}
\item{n}{number of observations.}
\item{loc}{vector of (real) location parameters.}
\item{disp}{vector of (positive) dispersion parameters.}
\item{skew}{vector of skewness parameters (in [-1,1]).}
\item{tail}{vector of parameters (in [0,2]) related to the tail thickness.}
\item{eps}{scalar giving the required precision in computation.}
}
\description{
These functions provide information about the stable distribution
with the location, the dispersion, the skewness and the tail thickness
respectively modelled by the parameters \code{loc}, \code{disp},
\code{skew} and  \code{tail}.

\code{dstable}, \code{pstable}, \code{qstable} and \code{hstable}
compute the density, the distribution, the quantile and the hazard functions
of a stable variate. \code{rstable} generates random deviates with 
the prescribed stable distribution.

\code{loc} is a location parameter in the same way as the mean
in the normal distribution: it can take any real value. 

\code{disp} is a dispersion parameter in the same way as the standard 
deviation in the normal distribution: it can take any positive value. 

\code{skew} is a skewness parameter: it can take any value in \eqn{(-1,1)}.
The distribution is right-skewed, symmetric and left-skewed
when \code{skew} is negative, null or positive respectively. 

\code{tail} is a tail parameter (often named the characteristic exponent):
it can take any value in \eqn{(0,2)} (with \code{tail=1} and \code{tail=2}
yielding the Cauchy and the normal distributions respectively
when symmetry holds). 

If \code{loc} or \code{disp} or \code{skew} or \code{tail} are not
specified they assume the default values of \eqn{0}, \eqn{1/sqrt(2)},
\eqn{0} and \eqn{2} respectively. This corresponds to a normal 
variate with mean\eqn{=0} and variance\eqn{=1=2 disp^2}.

The stable characteristic function is given by
\deqn{greekphi(t) = i loca t - disp {|t|}^{tail}  [1+i skew sign(t) greekomega(t,tail)]}{phi(t) = i loc t - disp |t|^tail [1+i skew sign(t) omega(t,tail)]}
where
\deqn{greekomega(t,tail) = \frac{2}{\pi} LOG(ABS(t))}{omega(t,tail) = (2/pi) log|t|}
when \code{tail=1}, and
\deqn{greekomega(t,tail) = tan(\frac{\pi tail}{2})}{omega(t,tail) = tan(pi alpha / 2)}
otherwise.

The characteristic function can be inverted using Fourier's transform
to obtain the corresponding stable density. This inversion requires the
numerical evaluation of an integral from \eqn{0} to \eqn{\infty}{infinity}.
Two algorithms
are proposed to do this. The default is the Romberg's method 
(\code{integration}="Romberg") which is used to evaluate the integral
with an error bounded by \code{eps}.
The alternative method is Simpson's integration 
(\code{integration}="Simpson"): it approximates the
integral from \eqn{0} to \eqn{\infty}{infinity} by an integral 
from \eqn{0} to \code{up}
with \code{npt} points subdividing \eqn{(O, up)}. 
These three extra arguments -- namely \code{integration}, \code{up} and
\code{npt} -- are only available when using \code{dstable}.
The other functions are all based on Romberg's algorithm.
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
\code{\link{stableglm}} to fit generalized linear models for the stable 
distribution parameters. \cr
\code{\link{stable.mode}} to compute the mode of a stable
distribution.
}

\examples{
par(mfrow=c(2,2))
x <- seq(-5,5,by=0.1)

# Influence of loc (location)
plot(x,dstable(x,loc=-2,disp=1/sqrt(2),skew=-0.8,tail=1.5),
  type="l",ylab="",title("Varying LOCation"))
lines(x,dstable(x,loc=0,disp=1/sqrt(2),skew=-0.8,tail=1.5))
lines(x,dstable(x,loc=2,disp=1/sqrt(2),skew=-0.8,tail=1.5))

# Influence of disp (dispersion)
plot(x,dstable(x,loc=0,disp=0.5,skew=0,tail=1.5),
  type="l",ylab="",title("Varying DISPersion"))
lines(x,dstable(x,loc=0,disp=1/sqrt(2),skew=0,tail=1.5))
lines(x,dstable(x,loc=0,disp=0.9,skew=0,tail=1.5))

# Influence of skew (skewness)
plot(x,dstable(x,loc=0,disp=1/sqrt(2),skew=-0.8,tail=1.5),
  type="l",ylab="",title("Varying SKEWness"))
lines(x,dstable(x,loc=0,disp=1/sqrt(2),skew=0,tail=1.5))
lines(x,dstable(x,loc=0,disp=1/sqrt(2),skew=0.8,tail=1.5))

# Influence of tail (tail)
plot(x,dstable(x,loc=0,disp=1/sqrt(2),skew=0,tail=0.8),
  type="l",ylab="",title("Varying TAIL thickness"))
lines(x,dstable(x,loc=0,disp=1/sqrt(2),skew=0,tail=1.5))
lines(x,dstable(x,loc=0,disp=1/sqrt(2),skew=0,tail=2))
}
