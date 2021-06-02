##' Fit a von Bertalanffy growth model to length-at-age fisheries data,
##' correcting for gear selectivity, length-stratified sampling, and ageing
##' errors.
##'
##' @title von Bertalanffy model for length selected data.
##' @param Age Vector of integer ages (zone counts)
##' @param Length Vector of measured lengths
##' @param Aerr Ageing error matrix
##' @param Lbreaks Internal length bin boundaries
##' @param fraction Sampling fraction for each bin
##' @param selectivity Selectivity function
##' @param Aoffset Ageing offset
##' @param Amin Minimum age
##' @param start Initial estimates of model parameters
##' @param control List of control parameters for nlminb
##' @param verbose Enable tracing information
##' @seealso \code{\link{lsgmFit}}, \code{\link{lsGompertz}}, \code{\link{lsLogistic}},
##' \code{\link{lsSchnuteRichards}}, \code{\link{nlminb}}, \code{\link[TMB]{MakeADFun}}
##' @return An object of class \code{lsVB} with elements
##' \item{\code{coefficients}}{a named vector of estimated parameters}
##' \item{\code{cov}}{the covariance of the estimated parameters}
##' \item{\code{logL}}{log likelihood}
##' \item{\code{aic}}{Akaike information criteria}
##' \item{\code{bic}}{Bayesian information criteria}
##' \item{\code{opt}}{return value from \code{nlminb}}
##' \item{\code{obj}}{TMB object}
##' \item{\code{obj}}{selectivity function}
##' \item{\code{model}}{model name}
##' \item{\code{call}}{the matched call}
##' \item{\code{call}}{function to evaluate the growth curve}
##' \item{\code{call}}{function to evaluate the gradient of the growth curve}
##' @useDynLib lsgm
##' @importFrom TMB MakeADFun sdreport
##' @importFrom stats setNames
##' @export
lsVB <- function(Age,Length,Aerr,Lbreaks,fraction=1,
                 selectivity=function(L) 1,
                 Aoffset=0,Amin=0,
                 start=list(Linf=1,K=1,t0=0,CV=1),
                 control=list(),verbose=FALSE) {

  cl <- match.call()
  ## VB growth model
  growth <- function(age,par) {
    pnames <- c("Linf","K","t0","CV")
    with(as.list(setNames(par,pnames[seq_len(par)])),
         Linf*(1-exp(-K*(age-t0))))
  }
  ## Gradient of the growth model wrt model parameters
  gradient <- function(age,par) {
    pnames <- c("Linf","K","t0","CV")
    with(as.list(setNames(par,pnames[seq_len(par)])),{
         .expr2 <- age-t0
         .expr4 <- exp(-K*.expr2)
         .expr5 <- 1-.expr4
         .value <- Linf*.expr5
         cbind(Linf=.expr5,
               K=Linf*(.expr4*.expr2),
               t0=-(Linf*(.expr4*K)),
               CV=0)})
  }

  growth <- function(age,par) par[1]*(1-exp(-par[2]*(age-par[3])))
  ## Gradient of the growth curve wrt parameters
  gradient <- function(age,par)
    cbind((1-exp(-par[2]*(age-par[3]))),
          par[1]*(age-par[3])*exp(-par[2]*(age-par[3])),
          -par[1]*par[2]*exp(-par[2]*(age-par[3])),
          0)

  ## von Bertalanffy parameters for TMB
  tmbpars <- list(logp1=log(as.numeric(start[c("Linf","K")])),
                  p2 = start$t0,
                  logCV=log(start$CV))

  ## Fit the model with TMB
  fit <- lsgmFit(Age,Length,Aerr,Lbreaks,fraction,selectivity,Aoffset,Amin,model=1L,tmbpars,control,verbose)

  ## Set VB specific return parameters
  names(fit$coefficients) <- c("Linf","K","t0","CV")
  dimnames(fit$cov) <- list(names(fit$coefficients),names(fit$coefficients))
  fit$call <- cl

  ## Growth curve and its gradient
  formals(growth)$par <- fit$coefficients[-length(fit$coefficients)]
  fit$growth <- growth
  formals(gradient)$par <- fit$coefficients[-length(fit$coefficients)]
  fit$gradient <- gradient
  fit$model <- "Length selected von Bertalanffy"

  class(fit) <- c("lsVB","lsgm")
  fit
}



##' Fit a Gompertz growth model to length-at-age fisheries data,
##' correcting for gear selectivity, length-stratified sampling, and
##' ageing errors.
##'
##'
##' @title Gompertz model for length-sampled data.
##' @param Age Vector of measured ages
##' @param Length Vector of measured lengths
##' @param Aerr Ageing error matrix
##' @param Lbreaks Internal length bin boundaries
##' @param fraction Sampling fraction for each bin
##' @param selectivity Selectivity function
##' @param Aoffset Ageing offset
##' @param Amin Minimum age
##' @param start Initial estimates of model parameters
##' @param control List of control parameters for nlminb
##' @param verbose Enable tracing information
##' @seealso \code{\link{lsgmFit}}, \code{\link{lsVB}}, \code{\link{lsLogistic}},
##' \code{\link{lsSchnuteRichards}}, \code{\link{nlminb}}, \code{\link[TMB]{MakeADFun}}
##' @return An object of class \code{lsVB} with elements
##' \item{\code{coefficients}}{a named vector of estimated parameters}
##' \item{\code{cov}}{the covariance of the estimated parameters}
##' \item{\code{logL}}{log likelihood}
##' \item{\code{aic}}{Akaike information criteria}
##' \item{\code{bic}}{Bayesian information criteria}
##' \item{\code{opt}}{return value from \code{nlminb}}
##' \item{\code{obj}}{TMB object}
##' \item{\code{obj}}{selectivity function}
##' \item{\code{model}}{model name}
##' \item{\code{call}}{the matched call}
##' \item{\code{call}}{function to evaluate the growth curve}
##' \item{\code{call}}{function to evaluate the gradient of the growth curve}
##' @useDynLib lsgm
##' @importFrom TMB MakeADFun sdreport
##' @importFrom stats setNames
##' @export
lsGompertz <- function(Age,Length,Aerr,Lbreaks,fraction=1,
                       selectivity=function(L) 1,
                       Aoffset=0,Amin=0,
                       start=list(Linf=1,K=1,t0=0,CV=1),
                       control=list(),verbose=FALSE) {

  cl <- match.call()

  ## Gompertz growth model
  growth <- function(age,par) {
    pnames <- c("Linf","K","t0","CV")
    with(as.list(setNames(par,pnames[seq_len(par)])),
         Linf*exp(-exp(-K*(age-t0))/K))
  }
  ## Gradient of the growth model wrt model parameters
  gradient <- function(age,par) {
    pnames <- c("Linf","K","t0","CV")
    with(as.list(setNames(par,pnames[seq_len(par)])),{
         .expr2 <- age-t0
         .expr4 <- exp(-K*.expr2)
         .expr7 <- exp(-.expr4/K)
         .value <- Linf*.expr7
         cbind(Linf=.expr7,
               K=Linf*(.expr7*(.expr4*.expr2/K+.expr4/K^2)),
               t0=-(Linf*(.expr7*(.expr4*K/K))),
               CV=0)})
  }

  ## Gompertz parameters for TMB
  tmbpars <- list(logp1=log(as.numeric(start[c("Linf","K")])),
                  p2 = start$t0,
                  logCV=log(start$CV))

  ## Fit the model with TMB
  fit <- lsgmFit(Age,Length,Aerr,Lbreaks,fraction,selectivity,Aoffset,Amin,model=2L,tmbpars,control,verbose)

  ## Set Gompertz specific return parameters
  names(fit$coefficients) <- c("Linf","K","t0","CV")
  dimnames(fit$cov) <- list(names(fit$coefficients),names(fit$coefficients))
  fit$call <- cl

  ## Growth curve and its gradient
  formals(growth)$par <- fit$coefficients[-length(fit$coefficients)]
  fit$growth <- growth
  formals(gradient)$par <- fit$coefficients[-length(fit$coefficients)]
  fit$gradient <- gradient
  fit$model <- "Length selected Gompertz"

  class(fit) <- c("lsGompertz","lsgm")
  fit
}

##' Fit a logistic growth model to length-at-age fisheries data, correcting for
##' gear selectivity, length-stratified sampling, and ageing errors.
##'
##'
##' @title Logistic model for length-sampled data.
##' @param Age Vector of measured ages
##' @param Length Vector of measured lengths
##' @param Aerr Ageing error matrix
##' @param Lbreaks Internal length bin boundaries
##' @param fraction Sampling fraction for each bin
##' @param selectivity Selectivity function
##' @param Aoffset Ageing offset
##' @param Amin Minimum age
##' @param start Initial estimates of model parameters
##' @param control List of control parameters for nlminb
##' @param verbose Enable tracing information
##' @seealso \code{\link{lsgmFit}}, \code{\link{lsVB}}, \code{\link{lsGompertz}},
##' \code{\link{lsSchnuteRichards}}, \code{\link{nlminb}}, \code{\link[TMB]{MakeADFun}}
##' @return An object of class \code{lsVB} with elements
##' \item{\code{coefficients}}{a named vector of estimated parameters}
##' \item{\code{cov}}{the covariance of the estimated parameters}
##' \item{\code{logL}}{log likelihood}
##' \item{\code{aic}}{Akaike information criteria}
##' \item{\code{bic}}{Bayesian information criteria}
##' \item{\code{opt}}{return value from \code{nlminb}}
##' \item{\code{obj}}{TMB object}
##' \item{\code{obj}}{selectivity function}
##' \item{\code{model}}{model name}
##' \item{\code{call}}{the matched call}
##' \item{\code{call}}{function to evaluate the growth curve}
##' \item{\code{call}}{function to evaluate the gradient of the growth curve}
##' @useDynLib lsgm
##' @importFrom TMB MakeADFun sdreport
##' @importFrom stats setNames
##' @export
lsLogistic <- function(Age,Length,Aerr,Lbreaks,fraction=1,
                       selectivity=function(L) 1,
                       Aoffset=0,Amin=0,
                       start=list(Linf=1,K=1,t0=1,CV=1),
                       control=list(),verbose=FALSE) {

  cl <- match.call()

  ## logistic growth model.
  growth <- function(age,par) {
    pnames <- c("Linf","K","t0","CV")
    with(as.list(setNames(par,pnames[seq_len(par)])),
         Linf/(1+exp(-K*(age-t0))))
  }
  ## Gradient of the growth curve wrt model parameters
  gradient <- function(age,par) {
    pnames <- c("Linf","K","t0","CV")
    with(as.list(setNames(par,pnames[seq_len(par)])),{
         .expr2 <- age-t0
         .expr4 <- exp(-K*.expr2)
         .expr5 <- 1+.expr4
         .expr10 <- .expr5^2
         .value <- Linf/.expr5
         cbind(Linf=1/.expr5,
               K=Linf*(.expr4*.expr2)/.expr10,
               t0=-(Linf*(.expr4*K)/.expr10),
               CV=0)})
  }

  ## logistic parameters for TMB
  tmbpars <- list(logp1=log(as.numeric(start[c("Linf","K")])),
                  p2 = start$t0,
                  logCV=log(start$CV))

  ## Fit the model with TMB
  fit <- lsgmFit(Age,Length,Aerr,Lbreaks,fraction,selectivity,Aoffset,Amin,model=3L,tmbpars,control,verbose)

  ## Set logistic specific return parameters
  names(fit$coefficients) <- c("Linf","K","t0","CV")
  dimnames(fit$cov) <- list(names(fit$coefficients),names(fit$coefficients))
  fit$call <- cl

  ## Growth curve and its gradient
  formals(growth)$par <- fit$coefficients[-length(fit$coefficients)]
  fit$growth <- growth
  formals(gradient)$par <- fit$coefficients[-length(fit$coefficients)]
  fit$gradient <- gradient
  fit$model <- "Length selected logistic"

  class(fit) <- c("lslogistic","lsgm")
  fit
}


##' Fit a Schnute-Richards growth model to length-at-age fisheries
##' data, correcting for gear selectivity, length-stratified sampling,
##' and ageing errors.
##'
##' @title Schnute-Richards model for length-sampled data.
##' @param Age Vector of measured ages
##' @param Length Vector of measured lengths
##' @param Aerr Ageing error matrix
##' @param Lbreaks Internal length bin boundaries
##' @param fraction Sampling fraction for each bin
##' @param selectivity Selectivity function
##' @param Aoffset Ageing offset
##' @param Amin Minimum age
##' @param start Initial estimates of model parameters
##' @param control List of control parameters for nlminb
##' @param verbose Enable tracing information
##' @seealso \code{\link{lsgmFit}}, \code{\link{lsVB}}, \code{\link{lsGompertz}},
##' \code{\link{lsSchnuteRichards}}, \code{\link{nlminb}}, \code{\link[TMB]{MakeADFun}}
##' @return An object of class \code{lsVB} with elements
##' \item{\code{coefficients}}{a named vector of estimated parameters}
##' \item{\code{cov}}{the covariance of the estimated parameters}
##' \item{\code{logL}}{log likelihood}
##' \item{\code{aic}}{Akaike information criteria}
##' \item{\code{bic}}{Bayesian information criteria}
##' \item{\code{opt}}{return value from \code{nlminb}}
##' \item{\code{obj}}{TMB object}
##' \item{\code{obj}}{selectivity function}
##' \item{\code{model}}{model name}
##' \item{\code{call}}{the matched call}
##' \item{\code{call}}{function to evaluate the growth curve}
##' \item{\code{call}}{function to evaluate the gradient of the growth curve}
##' @useDynLib lsgm
##' @importFrom TMB MakeADFun sdreport
##' @importFrom stats setNames
##' @export
lsSchnuteRichards <- function(Age,Length,Aerr,Lbreaks,fraction=1,
                       selectivity=function(L) 1,
                       Aoffset=0,Amin=0,
                       start=list(Linf=1,K=1,alpha=1,nu=1,gamma=1,CV=1),
                       control=list(),verbose=FALSE) {

  cl <- match.call()

  ## Schnute-Richards growth model
  growth <- function(age,par) {
    pnames <- c("Linf","K","alpha","nu","gamma","CV")
    with(as.list(setNames(par,pnames[seq_len(par)])),
         Linf*(1+alpha*exp(-K*(age^nu)))^(-1/gamma))
  }
  ## Gradient of the growth curve wrt model parameters
  gradient <- function(age,par) {
    pnames <- c("Linf","K","alpha","nu","gamma","CV")
    with(as.list(setNames(par,pnames[seq_len(par)])),{
         .expr2 <- age^nu
         .expr4 <- exp(-K*.expr2)
         .expr6 <- 1+alpha*.expr4
         .expr8 <- -1/gamma
         .expr9 <- .expr6^.expr8
         .expr12 <- .expr6^(.expr8 - 1)
         .value <- Linf*.expr9
         cbind(Linf=.expr9,
               K=-(Linf*(.expr12*(.expr8*(alpha*(.expr4*.expr2))))),
               alpha=Linf*(.expr12*(.expr8*.expr4)),
               nu=-(Linf*(.expr12*(.expr8*(alpha*(.expr4*(K*(.expr2*log(age)))))))),
               gamma=Linf*(.expr9*(log(.expr6)*(1/gamma^2))),
               CV=0)})
  }

  ## Schnute-Richards parameters for TMB
  tmbpars <- list(logp1=log(as.numeric(start[c("Linf","K","alpha","nu","gamma")])),
                  p2 = c(),
                  logCV=log(start$CV))

  ## Fit the model with TMB
  fit <- lsgmFit(Age,Length,Aerr,Lbreaks,fraction,selectivity,Aoffset,Amin,model=4L,tmbpars,control,verbose)

  ## Set Schnute Richards specific return parameters
  names(fit$coefficients) <- c("Linf","K","alpha","nu","gamma","CV")
  dimnames(fit$cov) <- list(names(fit$coefficients),names(fit$coefficients))
  fit$call <- cl

  ## Growth curve and its gradient
  formals(growth)$par <- fit$coefficients[-length(fit$coefficients)]
  fit$growth <- growth
  formals(gradient)$par <- fit$coefficients[-length(fit$coefficients)]
  fit$gradient <- gradient
  fit$model <- "Length selected Schnute-Richards"

  class(fit) <- c("lsSchnuteRichards","lsgm")
  fit
}




##' Fit a growth model to length-at-age fisheries data,
##' correcting for gear selectivity, length-stratified sampling, and ageing
##' errors.
##'
##' @title Growth model for length selected data.
##' @param Age Vector of measured ages
##' @param Length Vector of measured lengths
##' @param Aerr Ageing error matrix
##' @param Lbreaks Internal length bin boundaries
##' @param fraction Sampling fraction for each bin
##' @param selectivity Selectivity function
##' @param Aoffset Ageing offset
##' @param Amin Minimum age
##' @param model Model indicator
##' @param tmbpars Model parameter list
##' @param control List of control parameters for nlminb
##' @param verbose Enable tracing information
##' @seealso \code{\link{nlminb}}, \code{\link[TMB]{MakeADFun}}
##' @return An object of class \code{lsgm}, with elements
##' \item{\code{coefficients}}{a named vector of estimated parameters}
##' \item{\code{cov}}{the covariance of the estimated parameters}
##' \item{\code{logL}}{log likelihood}
##' \item{\code{aic}}{Akaike information criteria}
##' \item{\code{bic}}{Bayesian information criteria}
##' \item{\code{opt}}{return value from \code{nlminb}}
##' \item{\code{obj}}{TMB object}
##' \item{\code{selectivity}}{selectivity function}
##' \item{\code{model}}{model name}
##' @useDynLib lsgm
##' @importFrom TMB MakeADFun sdreport
##' @importFrom stats nlminb
##' @export
lsgmFit <- function(Age,Length,Aerr,Lbreaks,fraction=1,
                    selectivity=function(L) 1,
                    Aoffset=0,Amin=0,
                    model,tmbpars,
                    control=list(),verbose=FALSE) {

  ## Bin midpoints
  midpts <- c(Lbreaks[1],
              0.5*(Lbreaks[-1]+Lbreaks[-length(Lbreaks)]),
              Lbreaks[length(Lbreaks)])
  ## Selectivity and sampling fraction for each bin
  Sbin <- selectivity(midpts)
  Fbin <- rep(fraction,length.out=length(Lbreaks)+1)

  tmbdata <- list(A=Age,
                  L=Length,
                  ## Selectivity of the observations
                  S=selectivity(Length),
                  ## Length bin of each observation
                  bin=.bincode(Length,sort(unique(c(0,Lbreaks,Inf))),right=FALSE),
                  ## Ageing offset for each individual
                  Aoff=rep(Aoffset,length.out=length(Age)),
                  ## Length bin internal breaks
                  Lbrk=Lbreaks,
                  ## Selectivity for each bin
                  Sbin=Sbin,
                  ## Sampling fraction for each bin
                  Fbin=Fbin,
                  ## Ageing error matrix
                  Aerr=Aerr,
                  ## Minimum age
                  Amin=Amin,
                  ## Model type
                  model=model)

  ## Create TMB object
  obj <- MakeADFun(tmbdata,tmbpars,DLL="lsgm",silent=!verbose)
  ## Minimize objective function
  opt <- suppressWarnings(nlminb(obj$par, obj$fn, obj$gr,control=control))
  rep <- sdreport(obj)

  ## Extract estimates and full covariance
  coef <- rep$value
  cov <- rep$cov

  structure(
    list(coefficients=rep$value,
         cov=rep$cov,
         logL=-opt$objective,
         aic=2*length(coef)+2*opt$objective,
         bic=log(length(Length))*length(coef)+2*opt$objective,
         opt=opt,
         obj=obj,
         selectivity=selectivity,
         model="lsgm"),
    class=c("lsgm"))
}


##' Diagnostic quantities for an \code{lsgm} fit
##'
##' detials
##' @title Diagnostics for lsgm objects
##' @param fit an fitted \code{lsgm} object
##' @importFrom stats dnorm optimize
##' @return yes
##' @export
diagnostics <- function(fit) {

  ## Import objects from fit
  coef <- fit$coefficients
  CV <- coef[length(coef)]
  growth <- fit$growth
  selectivity <- fit$selectivity

  ## Import data
  Age <- fit$obj$env$data$A
  Length <- fit$obj$env$data$L
  Aoffset <- fit$obj$env$data$Aoff
  Lbreaks <- c(0,fit$obj$env$data$Lbrk,Inf)
  Sbin <- fit$obj$env$data$Sbin
  Fbin <- fit$obj$env$data$Fbin
  Aerr <- fit$obj$env$data$Aerr
  Amin <- fit$obj$env$data$Amin

  ## Length selectivity adjusted likelihood
  adjustedLik <- function(L,mu,cv) {
    b <- .bincode(L,Lbreaks,right=FALSE)
    N <- colSums(diff(outer(Lbreaks,as.vector(mu),function(L,mu) pnorm(L,mu,cv*mu)))*(Sbin*Fbin))
    dim(N) <- dim(mu)
    dnorm(L,mu,cv*mu)*selectivity(L)*Fbin[b]/N
  }

  fitted <- double(length(Age))
  deviance <- double(length(Age))
  meanLen <- double(length(Age))
  meanAge <- double(length(Age))
  xs <- c(0,fit$obj$env$data$Lbrk,2*max(fit$obj$env$data$Lbrk))
  xs <- sort(c(xs,(xs[-1]+xs[-length(xs)])/2))
  for(i in seq_along(Age)) {
    ## Compute means and weights for the non-zero weights
    w <- Aerr[Age[i]-Amin+1,]
    sub <- which(w>0)
    w <- w[sub]
    ages <- sub+(Amin-1)
    mu <- growth(ages+Aoffset[i],coef)
    ## Estimate fitted values as the maximum likelihood response for the observed age.  The
    ## likelihood is discontinuous at the bin boundaries, so we evaluate on a grid to find the
    ## approximate minimum
    nlik <- function(L) -sum(w*adjustedLik(L,mu,CV))
    j <- which.min(sapply(xs,nlik))
    opt <- optimize(nlik,xs[c(max(1,j-1),min(length(xs),j+1))],tol=1.0E-6)
    fitted[i] <- opt$minimum

    ps <- w*adjustedLik(Length[i],mu,CV)
    meanAge[i] <- sum(ps*ages)/sum(ps)
    meanLen[i] <- sum(ps*mu)/sum(ps)

    ## Fit a separate mean for each observation to calculate the deviance
    opt <- optimize(function(mu) -adjustedLik(Length[i],mu,CV),range(Length),tol=1.0E-6)
    deviance[i] <- 2*log(-opt$objective)-2*log(sum(ps))
  }
  data.frame(age=Age,length=Length,ageOffset=Aoffset,
             fitted=fitted,deviance=deviance,meanAge=meanAge,meanLength=meanLen)
}


##' Predicted values from a length selected growth model fit
##'
##' The computes predictions from a fitted length selected growth
##' model. The supplied ages are assumed to contain no error (but an
##' offset can be supplied).  The standard error of the mean is
##' computed by the delta method, and the standard deviation for the
##' future observations are computed conditional on the estimated CV.
##'
##' @title Predict method for \code{lsgm} fits
##' @param object an \code{lsgm} object
##' @param Age a vector of (true) ages to predict
##' @param Aoffset a vector of age offsets
##' @param ... currently ignored
##' @return a dataframe with columns
##' \item{\code{age}}{the age}
##' \item{\code{length}}{the predicted length}
##' \item{\code{se}}{the standard error of the predicted length}
##' \item{\code{sd}}{the standard deviation of a future observation}
##' @export
predict.lsgm <- function(object,Age,Aoffset=0,...) {

  age <- Age+Aoffset
  coef <- object$coefficients
  ## Fitted means and their gradient wrt to the parameters
  mu <- object$growth(age,coef)
  gr <- object$gradient(age,coef)
  ## Compute approximate variances by the delta method
  v <- rowSums((gr%*%object$cov)*gr)
  CV <- coef[length(coef)]
  data.frame(age=age,length=mu,se=sqrt(v),sd=CV*mu)
}



##' @export
print.lsgm <- function(x, digits = max(3L, getOption("digits")-3L), ...) {
  cat(x$model,"growth model\n\nCall:\n")
  print(x$call)
  invisible(x)
}



##' Summary method for \code{lsgm} objects
##'
##' @title Summarizing length selected growth model fits.
##' @param object an lsgm object
##' @param ... currently ignored
##' @return An object of class \code{summary.lsgm}
##' @importFrom stats pnorm
##' @export
summary.lsgm <- function(object,...) {
  coef <- object$coefficients
  se <- sqrt(diag(object$cov))
  zval <- abs(coef)/se
  coefficients <- cbind(Estimate = coef,
                        `Std. Error` = se,
                        `z value` = zval,
                        `Pr(>|z|)` = 2*pnorm(zval,lower.tail = FALSE))
  structure(
    list(model=object$model,
         logLik=data.frame(AIC=object$aic,BIC=object$bic,logLik=object$logL),
         coefficients=coefficients),
    class=c("summary.lsgm"))
}




##' @importFrom stats printCoefmat
##' @export
print.summary.lsgm <- function(x,digits=3,signif.stars = getOption("show.signif.stars"),...) {
  cat(x$model,"growth model\n\n")
  print(x$logLik, row.names = FALSE)
  cat("\nCoefficients:\n")
  printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars)
  invisible(x)
}


##' Calculate the log likelihood profile with respect to a single
##' parameter for a length selected growth model fit.
##'
##' This function is a wrapper for \code{\link[TMB]{tmbprofile}}.
##'
##' @title Summarizing length selected growth model fits.
##' @param fitted an lsgm object
##' @param which name or index of a parameter to profile
##' @param ... additional parameters to pass to \code{\link[TMB]{tmbprofile}}.
##' @return data.frame with parameter and function values.
##' @seealso \code{\link[TMB]{plot.tmbprofile}}, \code{\link[TMB]{confint.tmbprofile}}.
##' @importFrom TMB tmbprofile
##' @export
profile.lsgm <- function(fitted,which,...) {
  nms <- names(fitted$coefficients)
  if(is.character(which))
    which <- which(which==nms)
  pr <- tmbprofile(fitted$obj,which,...)
  names(pr)[1] <- nms[which]
  if(names(fitted$obj$par)[which] %in% c("logp1","logCV")) pr[[1]] <- exp(pr[[1]])
  pr
}
