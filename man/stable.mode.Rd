\name{stable.mode}
\title{Mode of a stable distribution}
\usage{
stable.mode(loc,disp,skew,tail)
}
\alias(stable.mode)

\arguments{
\item{loc}{vector of (real) location parameters.}
\item{disp}{vector of (positive) dispersion parameters.}
\item{skew}{vector of skewness parameters (in [-1,1]).}
\item{tail}{vector of parameters (in [1,2]) related to the tail thickness.}
}

\description{
This function gives a reliable approximation to the mode of a stable
distribution with location, dispersion, skewness and tail thickness 
specified by the parameters \code{loc}, \code{disp}, \code{skew} and \code{tail}.
Note that \code{tail} must be in (1,2).

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

The simplest empirical formula found to give a satisfactory approximation
to the mode for values of \code{tail} in \eqn{(1,2)} is
\deqn{
   loc+disp*a*skew*exp(-b*abs(skew))
}
with
\deqn{
 a = 1.7665114+1.8417675*tail-2.2954390*tail^2+0.4666749*tail^3 \cr
 b = -0.003142967+632.4715*tail*exp(-7.106035*tail)
}
}

\value{
A list of size 3 giving the mode, \eqn{a} and \eqn{b}.
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
\code{\link{stable}} for more details on the stable distribution. \cr
\code{\link{stableglm}} to fit generalized linear models for the stable 
distribution parameters.
}

\examples{
x <- seq(-5,5,by=0.1)
plot(x,dstable(x,loc=0,disp=1,skew=-1,tail=1.5),type="l",ylab="f(x)")
xhat <- stable.mode(loc=0,disp=1,skew=-1,tail=1.5)$ytilde
fxhat <- dstable(xhat,loc=0,disp=1,skew=-1,tail=1.5)
lines(c(xhat,xhat),c(0,fxhat),lty="dotted")
}
