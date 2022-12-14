% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/allfuns.r
\name{fun.4plm}
\alias{fun.4plm}
\title{The definition of four parameter log-logistic model function}
\usage{
fun.4plm(alp, x, logarg)
}
\arguments{
\item{alp}{the parameter of 4PLM, 4-dim}

\item{x}{the independent variable}

\item{logarg}{to choose the parameterzation form of 4PLM}
}
\value{
a 4PLM value of specified point
}
\description{
The definition of four parameter log-logistic model(4PLM) function. The log-logistic models are the most used dose-response models.
}
\examples{
\dontrun{
fx<-fun.4plm(c(2,2,5,1),seq(0.1, 20, len = 10))
} 
}
