"qqt" <-
function(p,df) sign(p-0.5)*sqrt(qf(1-2*pmin(p,1-p),1,df))

"mmake.link" <-
function (link,nu = 1, a = 0, b = 0) 
{
    if (is.character(link) && length(grep("^power", link) > 0)) 
        return(eval(parse(text = link)))
    else if (!is.character(link) && !is.na(lambda <- as.numeric(link))) {
        linkfun <- function(mu) mu^lambda
        linkinv <- function(eta) pmax(.Machine$double.eps, eta^(1/lambda))
        mu.eta <- function(eta) pmax(.Machine$double.eps, (1/lambda) * 
            eta^(1/lambda - 1))
        valideta <- function(eta) all(eta > 0)
    }
    else switch(link, logit = {
        linkfun <- function(mu) log(mu/(1 - mu))
        linkinv <- function(eta) {
            thresh <- -log(.Machine$double.eps)
            eta <- pmin(thresh, pmax(eta, -thresh))
            exp(eta)/(1 + exp(eta))
        }
        mu.eta <- function(eta) {
            thresh <- -log(.Machine$double.eps)
            res <- rep.int(.Machine$double.eps, length(eta))
            res[abs(eta) < thresh] <- (exp(eta)/(1 + exp(eta))^2)[abs(eta) < 
                thresh]
            res
        }
        valideta <- function(eta) TRUE
    }, cauchit = {
        linkfun <- function(mu) qcauchy(mu)
        linkinv <- function(eta) {
            thresh <- -qcauchy(.Machine$double.eps)
            eta <- pmin(thresh, pmax(eta, -thresh))
            pcauchy(eta)
        }
        mu.eta <- function(eta) pmax(dcauchy(eta), .Machine$double.eps)
        valideta <- function(eta) TRUE
    }, Gosset = {
        linkfun <- function(mu) qqt(mu,nu)
        linkinv <- function(eta) {
            thresh <- -qqt(.Machine$double.eps,nu)
            eta <- pmin(thresh, pmax(eta, -thresh))
            pt(eta,nu)
        }
        mu.eta <- function(eta) pmax(dt(eta,nu), .Machine$double.eps)
        valideta <- function(eta) TRUE
}, Pregibon = {
        linkfun <- function(mu) -qPregibon(1 - mu,a = a, b = b)
        linkinv <- function(eta) {
            tlo <- qPregibon(.Machine$double.eps^.5, a = a, b = b)
            thi <- qPregibon(1 - .Machine$double.eps^.5, a = a, b = b)
            eta <- -pmin(thi, pmax(-eta, tlo))
            1 - pPregibon(-eta, a = a, b = b)
        }
        mu.eta <- function(eta)
                pmax(dPregibon(-eta, a = a, b = b), .Machine$double.eps^.5)
        valideta <- function(eta) TRUE
    }, probit = {
        linkfun <- function(mu) qnorm(mu)
        linkinv <- function(eta) {
            thresh <- -qnorm(.Machine$double.eps)
            eta <- pmin(thresh, pmax(eta, -thresh))
            pnorm(eta)
        }
        mu.eta <- function(eta) pmax(dnorm(eta), .Machine$double.eps)
        valideta <- function(eta) TRUE
    }, cloglog = {
        linkfun <- function(mu) log(-log(1 - mu))
        linkinv <- function(eta) pmax(.Machine$double.eps, pmin(1 - 
            .Machine$double.eps, -expm1(-exp(eta))))
        mu.eta <- function(eta) {
            eta <- pmin(eta, 700)
            pmax(.Machine$double.eps, exp(eta) * exp(-exp(eta)))
        }
        valideta <- function(eta) TRUE
    }, identity = {
        linkfun <- function(mu) mu
        linkinv <- function(eta) eta
        mu.eta <- function(eta) rep.int(1, length(eta))
        valideta <- function(eta) TRUE
    }, log = {
        linkfun <- function(mu) log(mu)
        linkinv <- function(eta) pmax(.Machine$double.eps, exp(eta))
        mu.eta <- function(eta) pmax(.Machine$double.eps, exp(eta))
        valideta <- function(eta) TRUE
    }, sqrt = {
        linkfun <- function(mu) mu^0.5
        linkinv <- function(eta) eta^2
        mu.eta <- function(eta) 2 * eta
        valideta <- function(eta) all(eta > 0)
    }, "1/mu^2" = {
        linkfun <- function(mu) 1/mu^2
        linkinv <- function(eta) 1/eta^0.5
        mu.eta <- function(eta) -1/(2 * eta^1.5)
        valideta <- function(eta) all(eta > 0)
    }, inverse = {
        linkfun <- function(mu) 1/mu
        linkinv <- function(eta) 1/eta
        mu.eta <- function(eta) -1/(eta^2)
        valideta <- function(eta) all(eta != 0)
    }, stop(paste(link, "link not recognised")))
    list(linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta, 
        valideta = valideta)
}

"mbinomial" <-
function (link = "logit", ...) 
{
    linktemp <- substitute(link)
    if (!is.character(linktemp)) {
        linktemp <- deparse(linktemp)
        if (linktemp == "link") 
            linktemp <- eval(link)
    }
    if (any(linktemp == c("logit", "probit", "cloglog", "log", "Pregibon", "Gosset", "cauchit"))) 
        stats <- mmake.link(linktemp, ...)
    else stop(paste(linktemp, "link not available for binomial", 
        "family, available links are \"logit\", ", "\"probit\", \"cloglog\" and \"log\""))
    variance <- function(mu) mu * (1 - mu)
    validmu <- function(mu) all(mu > 0) && all(mu < 1)
    dev.resids <- function(y, mu, wt) 2 * wt * (y * log(ifelse(y == 
        0, 1, y/mu)) + (1 - y) * log(ifelse(y == 1, 1, (1 - y)/(1 - 
        mu))))
    aic <- function(y, n, mu, wt, dev) {
        m <- if (any(n > 1)) 
            n
        else wt
        -2 * sum(ifelse(m > 0, (wt/m), 0) * dbinom(round(m * 
            y), round(m), mu, log = TRUE))
    }
    initialize <- expression({
        if (NCOL(y) == 1) {
            if (is.factor(y)) y <- y != levels(y)[1]
            n <- rep.int(1, nobs)
            if (any(y < 0 | y > 1)) stop("y values must be 0 <= y <= 1")
            mustart <- (weights * y + 0.5)/(weights + 1)
            m <- weights * y
            if (any(abs(m - round(m)) > 0.001)) warning("non-integer #successes in a binomial glm!")
        } else if (NCOL(y) == 2) {
            if (any(abs(y - round(y)) > 0.001)) warning("non-integer counts in a binomial glm!")
            n <- y[, 1] + y[, 2]
            y <- ifelse(n == 0, 0, y[, 1]/n)
            weights <- weights * n
            mustart <- (n * y + 0.5)/(n + 1)
        } else stop(paste("For the binomial family, y must be", 
            "a vector of 0 and 1's or a 2 column", "matrix where col 1 is no. successes", 
            "and col 2 is no. failures"))
    })
    structure(list(family = "binomial", link = linktemp, linkfun = stats$linkfun, 
        linkinv = stats$linkinv, variance = variance, dev.resids = dev.resids, 
        aic = aic, mu.eta = stats$mu.eta, initialize = initialize, 
        validmu = validmu, valideta = stats$valideta), class = "family")
}

# pqd functions for Pregibon link from gld package of Robert King
# Note that we are fixing the scale so that we have logistic scale for all a,b 

qPregibon <- function(x,a = 0,b = 0){
        scale <- (qgl(3/4,c(0,1,a-b,a+b))- qgl(1/4,c(0,1,a-b,a+b)))/2.197224
	#scale = 1
        #qgl(x,c(0,1,a-b,a+b)) -  qgl(.5,c(0,1,a-b,a+b))
        qgl(x,c(0,1,a-b,a+b))/scale
        }
pPregibon <- function(x,a = 0,b = 0,tol=1e-12){
        scale <- (qgl(3/4,c(0,1,a-b,a+b))- qgl(1/4,c(0,1,a-b,a+b)))/2.197224
	#scale = 1
        #pgl(x + qgl(.5,c(0,1,a-b,a+b)),c(0,1,a-b,a+b),inverse.eps=tol)
        pgl(x*scale, c(0,1,a-b,a+b),inverse.eps=tol)
        }
dPregibon <- function(x,a = 0,b = 0,tol=1e-12){
        scale <- (qgl(3/4,c(0,1,a-b,a+b))- qgl(1/4,c(0,1,a-b,a+b)))/2.197224
	#scale = 1
        #dgl(x + qgl(.5,c(0,1,a-b,a+b)), c(0,1,a-b,a+b),inverse.eps=tol)
        dgl(x*scale, c(0,1,a-b,a+b),inverse.eps=tol)*scale
        }

rPregibon <- function(n,a = 0,b = 0){
	qPregibon(runif(n),a=a,b=b)
	}

pglm <- function(formula,link="Pregibon", method = "optim", Warm = TRUE, plevel = 0){
  switch(link,
    Gosset = {
      nuInt <- c(.15,30)
      g <- function(nu, x, y)
        glm(formula, family=mbinomial(link="Gosset",nu))$aic
      h <- function(nu, x, y, gopt)
        g(nu, x, y) - gopt -  3.84
      f <- optimize(g, nuInt, x = x, y = y)
      nuhat <- f$minimum
      gopt <- f$objective
      f <- glm(formula,family=mbinomial(link="Gosset",nuhat))
      if(h(nuInt[1],x,y,gopt) > 0 )
        nulo <- uniroot(h,c(nuInt[1],nuhat),x = x, y = y, gopt = gopt)$root
      else
        nulo <-  nuInt[1]
      if(h(nuInt[2],x,y,gopt) > 0 )
        nuhi <- uniroot(h,c(nuhat,nuInt[2]),x = x, y = y, gopt = gopt)$root
      else
        nuhi <- nuInt[2]
      return(list(f=f, nuhat=nuhat, nulo=nulo, nuhi=nuhi)) }, 
    Pregibon = {
      g <- function(theta, x, y)
        glm(y ~ x, family=mbinomial(link="Pregibon",a=theta[1],b=theta[2]))$aic
      if(Warm){ #Pregibon one-step starting values
        f <- glm(formula,family=binomial(link="logit"))
        m <- model.frame(f)
        u <- f$fitted.values
        m$ga <-  - .5*(log(u)^2 - log(1-u)^2)
        m$gd <-  - .5*(log(u)^2 + log(1-u)^2)
        fu <- update(f, . ~ . + ga + gd, data = m)
        #plot(ellipse(f,which=c("ga","gd"))
	GOLtest <- anova(f,fu)$Deviance[2]
        theta <- coef(fu)[c("ga","gd")] 
	}
      else
	theta <- rep(0,2)
      if(method == "nlm"){
        f <- nlm(g, theta, typsize = c(.5,2), 
		print.level = plevel, ndigit = 5, x = x, y = y)
        coef <- f$estimate 
        gopt <- f$minimum 
	} 
      else{
	#NB:  box limits changed from +\- 2 to +/- .5 May 18, 2006
        f <- optim(theta,g,method = "L-BFGS-B", lower = c(-.5,-.5), 
		upper = c(.5,.5),x=x,y=y)
	mess <- f$message
        coef <- f$par
        gopt <- f$value
	}
      gnull <- glm(y ~ x, family=binomial(link="logit"))$aic
      LRtest <- gopt - gnull 
      if(Warm)
      	return(list(theta = theta, GOLtest = GOLtest, coef = coef, LRtest = LRtest))
      else
      	return(list(coef = coef, LRtest = LRtest)) 
    }, stop(paste(link, "link not recognized")))
}
