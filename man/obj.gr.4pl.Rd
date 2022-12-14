% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/allfuns.r
\name{obj.gr.4pl}
\alias{obj.gr.4pl}
\title{4PLM's gradient}
\usage{
obj.gr.4pl(p, x)
}
\arguments{
\item{p}{the parameter of 4PLM, 4-dim}

\item{x}{the independent variable}
}
\description{
Calculate 4PLM's gradient of specified point
}
\examples{
\dontrun{
4plm.grad<-obj.gr.4pl(c(2,2,5,1),seq(0.1, 20, len = 10))
}
}
