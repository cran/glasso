\name{glassopath}
\alias{glassopath}
\title{Compute the Graphical lasso  along a path}
\description{
Estimates a sparse inverse covariance matrix using a lasso (L1) penalty, along a path
of values for the regularization parameter
}
\usage{
glassopath(s, rholist=NULL, thr=1.0e-4, maxit=1e4,  approx=FALSE, 
penalize.diagonal=TRUE, w.init=NULL,wi.init=NULL, trace=1)
}
\arguments{
\item{s}{Covariance matrix:p by p matrix (symmetric)}
\item{rholist}{Vector of non-negative regularization parameters for the lasso.
Should be  increasing from smallest to largest; actual path is computed from largest
 to smallest value of rho).
 If NULL, 10 values in a   (hopefully reasonable) range are used. Note that
the same parameter rholist[j] is used for all entries of the inverse covariance matrix;
different penalties for different entries are not allowed. }
\item{thr}{Threshold for convergence. Default value is 1e-4.  Iterations stop when average absolute parameter change is less than thr * ave(abs(offdiag(s)))}
\item{maxit}{Maximum number of iterations of outer loop. Default 10,000}
\item{approx}{Approximation flag: if true, computes Meinhausen-Buhlmann(2006)
approximation}
\item{penalize.diagonal}{Should diagonal of inverse covariance be penalized?
Dafault TRUE.}
\item{w.init}{Optional starting values for estimated covariance matrix (p by p).
Only needed when start="warm" is specified}
\item{wi.init}{Optional starting values for estimated inverse covariance matrix (p by p)
Only needed when start="warm" is specified}
\item{trace}{Flag for printing out information as iterations proceed. trace=0 means no printing; trace=1 means outer level printing;  trace=2 means  full printing
Default FALSE}
}

\details{
Estimates a sparse inverse covariance matrix using a lasso (L1) penalty, along a path of
regularization paramaters,
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
\item{w}{Estimated covariance matrices, an array  of dimension (nrow(s),ncol(n), length(rholist))}
\item{wi}{Estimated inverse covariance matrix, an array  of dimension (nrow(s),ncol(n), length(rholist))}
\item{approx}{Value of input argument approx}
\item{rholist}{Values of regularization parameter used}
\item{errflag}{values of error flag (0 means no memory allocation error)}
}


 \references{
Jerome Friedman, Trevor Hastie and Robert Tibshirani (2007).
Sparse inverse covariance estimation with the lasso.
Biostatistics 2007. http://www-stat.stanford.edu/~tibs/ftp/graph.pdf

Meinshausen, N. and  Buhlmann, P.(2006)
 High dimensional graphs
and variable selection with the lasso.
Annals of Statistics,34, p1436-1462.

Daniela Witten, Jerome Friedman, Noah Simon (2011).
New insights and fast computation for the graphical lasso. To appear in
Journal of Computational and Graphical Statistics.
}
\examples{


set.seed(100)

x<-matrix(rnorm(50*20),ncol=20)
s<- var(x)
a<-glassopath(s)

}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}

