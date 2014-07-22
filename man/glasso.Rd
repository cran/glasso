\name{glasso}
\alias{glasso}
\title{Graphical lasso}
\description{
Estimates a sparse inverse covariance matrix using a lasso (L1) penalty
}
\usage{
glasso(s, rho, zero=NULL, thr=1.0e-4, maxit=1e4,  approx=FALSE, 
penalize.diagonal=TRUE, start=c("cold","warm"), 
w.init=NULL,wi.init=NULL, trace=FALSE)
}
\arguments{
\item{s}{Covariance matrix:p by p matrix (symmetric)}
\item{rho}{(Non-negative) regularization parameter for lasso. rho=0 means no regularization. Can be a scalar (usual) or a symmetric p by p matrix, or a vector of length p. In the latter case, the penalty matrix has jkth element sqrt(rho[j]*rho[k]).}
\item{zero}{(Optional) indices of entries of inverse covariance to be constrained to be zero. The input should be a matrix with two columns, each row indicating
the indices of elements to be constrained to be zero. The solution must be symmetric, so you need only  specify one of (j,k) and (k,j). An entry in the zero matrix
overrides any entry in the rho matrix for a given element.}
\item{thr}{Threshold for convergence. Default value is 1e-4.  Iterations stop when average absolute parameter change is less than thr * ave(abs(offdiag(s)))}
\item{maxit}{Maximum number of iterations of outer loop. Default 10,000}
\item{approx}{Approximation flag: if true, computes Meinhausen-Buhlmann(2006)
approximation}
\item{penalize.diagonal}{Should diagonal of inverse covariance be penalized?
Dafault TRUE.}
\item{start}{Type of start. Cold start is default. Using Warm start, can provide starting values for w and wi}
\item{w.init}{Optional starting values for estimated covariance matrix (p by p).
Only needed when start="warm" is specified}
\item{wi.init}{Optional starting values for estimated inverse covariance matrix (p by p)
Only needed when start="warm" is specified}
\item{trace}{Flag for printing out information as iterations proceed.
Default FALSE}
}

\details{
Estimates a sparse inverse covariance matrix using a lasso (L1) penalty,
using the approach of Friedman, Hastie and Tibshirani (2007).
The Meinhausen-Buhlmann (2006) approximation is also implemented. 
The algorithm can also be used to estimate a graph with missing edges,
by specifying which edges to omit in the zero argument, and setting rho=0.
Or both fixed zeroes for some elements and regularization on the other elements
can be specified.

This version 1.7 uses a block diagonal screening rule to speed up
computations considerably. Details are given in the paper "New insights
and fast computations for the graphical lasso" by Daniela Witten, Jerry
Friedman, and Noah Simon, to appear in "Journal of Computational and
Graphical Statistics". The idea is as follows: it is possible to quickly
check whether the solution to the graphical lasso problem will be block
diagonal, for a given value of the tuning parameter. If so, then one can
simply apply the graphical lasso algorithm to each block separately,
leading to massive speed improvements.
}
\value{
A list with components
\item{w}{Estimated covariance matrix}
\item{wi}{Estimated inverse covariance matrix}
\item{loglik}{Value of maximized log-likelihodo+penalty}
\item{errflag}{Memory allocation error flag: 0 means no error; !=0 means
 memory allocation error - no output returned}
\item{approx}{Value of input argument approx}
\item{del}{Change in parameter value at convergence}
\item{niter}{Number of iterations of outer loop used by algorithm}
}


 \references{
Jerome Friedman, Trevor Hastie and Robert Tibshirani (2007).
Sparse inverse covariance estimation with the lasso.
Biostatistics 2007. http://www-stat.stanford.edu/~tibs/ftp/graph.pdf

Meinshausen, N. and  Buhlmann, P.(2006)
 High dimensional graphs
and variable selection with the lasso.
Annals of Statistics,34, p1436-1462.

Daniela Witten, Jerome Friedman, and Noah Simon (2011). New insights and
faster computations for the graphical lasso. To appear in Journal of
Computational and Graphical Statistics.
}
\examples{


set.seed(100)

x<-matrix(rnorm(50*20),ncol=20)
s<- var(x)
a<-glasso(s, rho=.01)
aa<-glasso(s,rho=.02, w.init=a$w, wi.init=a$wi)

# example with structural zeros and no regularization,
# from Whittaker's Graphical models book  page xxx.

s=c(10,1,5,4,10,2,6,10,3,10)
S=matrix(0,nrow=4,ncol=4)
S[row(S)>=col(S)]=s
S=(S+t(S))
diag(S)<-10
zero<-matrix(c(1,3,2,4),ncol=2,byrow=TRUE)
a<-glasso(S,0,zero=zero)
}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}

