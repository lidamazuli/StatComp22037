##############################################################################
#####################################################

#' @title A Gibbs sampler using R
#' @description A Gibbs sampler using R
#' @param N the number of samples
#' @param burn burn-in length, the length of deserted samples from the first
#' @param mu1 the mean of first marginal value
#' @param mu2 the mean of second marginal value
#' @param sigma1 the variance of first marginal value
#' @param sigma2 the variance of second marginal value
#' @param rho correlation
#' @return a random sample of size \code{n}
#' @importFrom stats rnorm rgamma ts fft
#' @examples
#' \dontrun{
#' rnR <- Gibbsinr(5000,0,0,0,1,1,.9)
#' par(mfrow=c(2,1));
#' plot(rnR[,1],type='l')
#' plot(rnR[,2],type='l')
#' }
#' @export
Gibbsinr<-function(N,burn,mu1,mu2,sigma1,sigma2,rho){
  #initialize constants and parameters
  #N: length of chain
  #burn: burn-in length
  #rho: correlation
  X <- matrix(0, N, 2)    #the chain, a bivariate sample
  s1 <- sqrt(1-rho^2)*sigma1
  s2 <- sqrt(1-rho^2)*sigma2

  ###### generate the chain #####

  X[1, ] <- c(mu1, mu2)            #initialize

  for (i in 2:N) {
    x2 <- X[i-1, 2]
    m1 <- mu1 + rho * (x2 - mu2) * sigma1/sigma2
    X[i, 1] <- rnorm(1, m1, s1)
    x1 <- X[i, 1]
    m2 <- mu2 + rho * (x1 - mu1) * sigma2/sigma1
    X[i, 2] <- rnorm(1, m2, s2)
  }

  b <- burn + 1
  x <- X[b:N, ]
  return(x)
}
#' @title The definition of two parameter log-logistic model function
#' @description The definition of two parameter log-logistic model(2PLM) function. The log-logistic models are the most used dose-response models. 
#' @param alp the parameter of 2PLM, 2-dim
#' @param x the independent variable
#' @param logarg to choose the parameterzation form of 2PLM
#' @return a 2PLM value of specified point
#' @examples
#' \dontrun{
#' fx<-fun.2plm(c(2,2),seq(0.1, 20, len = 10),1)
#' } 
#' 
#' @export
fun.2plm<-function(alp,x, logarg)
{
  b=alp[1]; e=alp[2]
  if(logarg==1){y = 1/(1+exp(b*(log(x)-log(e))))}  #--> LL.2()
  if(logarg!=1){y = 1/(1+exp(b*(log(x)-e)))}       #--> LL2.3()
  return(y)
}


# calculate gradient and Hessian
der.2plm<-deriv(~(1/(1+exp(b*(log(x)-log(e))))),
                c("b","e"),
                function(b,e,x) NULL,hessian = TRUE)

der.2plm2<-deriv(~(1/(1+exp(b*(log(x)-e)))),
                 c("b","e"),
                 function(b,e,x) NULL,hessian = TRUE)

#' @title 2PLM's gradient
#' @description Calculate 2PLM's gradient of specified point
#' @param p the parameter of 2PLM, 2-dim
#' @param x the independent variable
#' @return a gradient of 2PLM at the specified point
#' @examples
#' \dontrun{
#' 2plm.grad<-obj.gr.2pl(c(2,2),seq(0.1, 20, len = 10))
#' } 
#' @export
obj.gr.2pl<-function(p, x){
  b<-p[1]; e<-p[2];
  attr(der.2plm(b,e,x), "gradient")
  #attr(obj.der.2plm2(b,e,x), "gradient")
}

#' @title 2PLM's Hessian matrix
#' @description Calculate 2PLM's Hessian matrix of specified point
#' @param p the parameter of 2PLM, 2-dim
#' @param x the independent variable
#' @examples
#' \dontrun{
#' 2plm.hess<-obj.hess.2pl(c(2,2),seq(0.1, 20, len = 10))
#' } 
#' @export
obj.hess.2pl<-function(p, x){
  b<-p[1]; e<-p[2];
  attr(der.2plm(b,e,x), "hessian")
}

# end 2PLM
#####################################################


#' @title The definition of three parameter log-logistic model function
#' @description The definition of three parameter log-logistic model(3PLM) function. The log-logistic models are the most used dose-response models. 
#' @param alp the parameter of 3PLM, 3-dim
#' @param x the independent variable
#' @param logarg to choose the parameterzation form of 3PLM
#' @return a 3PLM value of specified point
#' @examples
#' \dontrun{
#' fx<-fun.3plm(c(2,2,5),seq(0.1, 20, len = 10))
#' } 
#' @export

fun.3plm<-function(alp,x, logarg)
{
  b=alp[1]; d=alp[2]; e=alp[3]
  if(logarg==1){y = d/(1+exp(b*(log(x)-log(e))))}       #--> LL.3()
  if(logarg!=1){y = d/(1+exp(b*(log(x)-e)))}            #--> LL2.3()
  return(y)
}

# calculate gradient and Hessian
der.3plm<-deriv(~(d/(1+exp(b*(log(x)-log(e))))),
                c("b","d","e"),
                function(b,d,e,x) NULL,hessian = TRUE)

der.3plm2<-deriv(~(d/(1+exp(b*(log(x)-e)))),
                 c("b","d","e"),
                 function(b,d,e,x) NULL,hessian = TRUE)

#' @title 3PLM's gradient
#' @description Calculate 3PLM's gradient of specified point
#' @param p the parameter of 3PLM, 3-dim
#' @param x the independent variable
#' @examples
#' \dontrun{
#' 3plm.grad<-obj.gr.3pl(c(2,2,5),seq(0.1, 20, len = 10))
#' } 
#' @export
obj.gr.3pl<-function(p, x){
  b<-p[1];d<-p[2]; e<-p[3];
  attr(der.3plm(b,d,e,x), "gradient")
  #attr(obj.der.3plm2(b,d,e,x), "gradient")
}

#' @title 3PLM's Hessian matrix
#' @description Calculate 3PLM's Hessian matrix of specified point
#' @param p the parameter of 3PLM, 3-dim
#' @param x the independent variable
#' @examples
#' \dontrun{
#' 3plm.hess<-obj.hess.2pl(c(2,2,5),seq(0.1, 20, len = 10))
#' }
#' @export
obj.hess.3pl<-function(p, x){
  b<-p[1];d<-p[2]; e<-p[3];
  attr(der.3plm(b,d,e,x), "hessian")
}

# end 3PLM
#####################################################


#' @title The definition of four parameter log-logistic model function
#' @description The definition of four parameter log-logistic model(4PLM) function. The log-logistic models are the most used dose-response models. 
#' @param alp the parameter of 4PLM, 4-dim
#' @param x the independent variable
#' @param logarg to choose the parameterzation form of 4PLM
#' @return a 4PLM value of specified point
#' @examples
#' \dontrun{
#' fx<-fun.4plm(c(2,2,5,1),seq(0.1, 20, len = 10))
#' } 
#' @export

fun.4plm<-function(alp,x, logarg)
{
  b=alp[1]; c=alp[2]; d=alp[3]; e=alp[4]
  if(logarg==1){y = c+(d-c)/(1+exp(b*(log(x)-log(e)))) } #--> LL.4()
  if(logarg!=1){y = c+(d-c)/(1+exp(b*(log(x)-e)))}      #--> LL2.4()
  return(y)
}

# calculate gradient and Hessian
der.4plm<-deriv(~(c+(d-c)/(1+exp(b*(log(x)-log(e))))),
                c("b","c","d","e"),
                function(b,c,d,e,x) NULL,hessian = TRUE)

der.4plm2<-deriv(~(c+(d-c)/(1+exp(b*(log(x)-e)))),
                 c("b","c","d","e"),
                 function(b,c,d,e,x) NULL,hessian = TRUE)

#' @title 4PLM's gradient
#' @description Calculate 4PLM's gradient of specified point
#' @param p the parameter of 4PLM, 4-dim
#' @param x the independent variable
#' @examples
#' \dontrun{
#' 4plm.grad<-obj.gr.4pl(c(2,2,5,1),seq(0.1, 20, len = 10))
#' }
#' @export
obj.gr.4pl<-function(p, x){
  b<-p[1];c<-p[2]; d<-p[3]; e<-p[4]
  attr(der.4plm(b,c,d,e,x), "gradient")
  #attr(obj.der.4plm2(b,c,d,e,x), "gradient")
}


#' @title 4PLM's Hessian matrix
#' @description Calculate 4PLM's Hessian matrix of specified point
#' @param p the parameter of 4PLM, 4-dim
#' @param x the independent variable
#' @examples
#' \dontrun{
#' 4plm.hess<-obj.hess.4pl(c(2,2,5,1),seq(0.1, 20, len = 10))
#' }
#' @export
obj.hess.4pl<-function(p, x){
  b<-p[1];c<-p[2]; d<-p[3]; e<-p[4]
  attr(der.4plm(b,c,d,e,x), "hessian")
}

# end 4PLM
#####################################################

#' @title The definition of five parameter log-logistic model function
#' @description The definition of five parameter log-logistic model(5PLM) function. The log-logistic models are the most used dose-response models. 
#' @param alp the parameter of 5PLM, 5-dim
#' @param x the independent variable
#' @param logarg to choose the parameterzation form of 5PLM
#' @return a 5PLM value of specified point
#' @examples
#' \dontrun{
#' fx<-fun.5plm(c(2,2,5,1,3),seq(0.1, 20, len = 10))
#' } 
#' @export

fun.5plm<-function(alp,x, logarg)
{
  b=alp[1]; c=alp[2]; d=alp[3]; e=alp[4]; f=alp[5]
  if(logarg==1){y = c+(d-c)/((1+exp(b*(log(x)-log(e))))^f)}   #--> LL.5()
  if(logarg!=1){y = c+(d-c)/((1+exp(b*(log(x)-e)))^f)}        #--> LL2.5()
  return(y)
}

# calculate gradient and Hessian
der.5plm<-deriv(~(c+(d-c)/((1+exp(b*(log(x)-log(e))))^f)),
                c("b","c","d","e","f"),
                function(b,c,d,e,f,x) NULL,hessian = TRUE)

der.5plm2<-deriv(~(c+(d-c)/((1+exp(b*(log(x)-e)))^f)),
                 c("b","c","d","e","f"),
                 function(b,c,d,e,f,x) NULL,hessian = TRUE)

#' @title 5PLM's gradient
#' @description Calculate 5PLM's gradient of specified point
#' @param p the parameter of 5PLM, 5-dim
#' @param x the independent variable
#' @examples
#' \dontrun{
#' 5plm.grad<-obj.gr.5pl(c(2,2,5,1,3),seq(0.1, 20, len = 10))
#' }
#' @export
obj.gr.5pl<-function(p, x){
  b<-p[1];c<-p[2]; d<-p[3]; e<-p[4]; f<-p[5]
  attr(der.5plm(b,c,d,e,f,x), "gradient")
  #attr(obj.der.5plm2(b,c,d,e,f,x), "gradient")
}

#' @title 5PLM's Hessian matrix
#' @description Calculate 5PLM's Hessian matrix of specified point
#' @param p the parameter of 5PLM, 5-dim
#' @param x the independent variable
#' @examples
#' \dontrun{
#' 5plm.hess<-obj.hess.5pl(c(2,2,5,1,3),seq(0.1, 20, len = 10))
#' }
#' @export
obj.hess.5pl<-function(p, x){
  b<-p[1];c<-p[2]; d<-p[3]; e<-p[4]; f<-p[5]
  attr(der.5plm(b,c,d,e,f,x), "hessian")
}

# end 5PLM
#####################################################




