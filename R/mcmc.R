#' @title mcmc function
#' @description MCMC routine for the strucatural equation model
#' 
#' @references {Maity, A. K., Lee, S. C., Mallick, B. K., & Sarkar, T. R. (2020). 
#' Bayesian structural equation modeling in multiple omics data integration 
#' with application to circadian genes. Bioinformatics, 36(13), 3951-3958.}
#' 
#' @param ct survival response, a \eqn{n*2} matrix with first column as response and second column as right censored indicator,
#' 1 is event time and 0 is right censored.
#' @param X Matrix of covariates, dimension \eqn{n*p}.
#' @param u1 Matix of predictors from the first platform, dimension \eqn{n*q_1}
#' @param u2 Matix of predictors from the first platform, dimension \eqn{n*q_2}
#' @param nburnin number of burnin samples
#' @param nmc number of markov chain samples
#' @param nthin thinning parameter. Default is 1 (no thinning)
#' 
#' @return \item{pMean.beta.t}{}
#' \item{pMean.beta.t}{} 
#' \item{pMean.alpha.t}{}
#' \item{pMean.alpha.t}{} 
#' \item{pMean.phi.t}{} 
#' \item{pMean.phi.t}{}
#' \item{pMean.alpha.u1}{}
#' \item{pMean.alpha.u2}{}
#' \item{pMean.alpha.u2}{}
#' \item{pMean.phi.u1}{}
#' \item{pMean.eta1}{}
#' \item{pMean.eta2}{}
#' \item{pMean.sigma.t.square}{}
#' \item{pMean.sigma.u1.square}{}
#' \item{pMean.sigma.u2.square}{}
#' \item{alpha.t.samples}{}
#' \item{phi.t.samples}{}
#' \item{beta1.t.samples}{}
#' \item{beta2.t.samples}{}
#' \item{beta.t.samples}{}
#' \item{alpha.u1.samples}{}
#' \item{alpha.u2.samples}{}
#' \item{phi.u1.samples}{}
#' \item{phi.u2.samples}{}
#' \item{eta1.samples}{}
#' \item{eta2.samples}{}
#' \item{sigma.t.square.samples}{}
#' \item{sigma.u1.square.samples}{}
#' \item{sigma.u2.square.samples}{}
#' \item{pMean.logt.hat}{}
#' \item{DIC}{}
#' \item{WAIC}{}
#'
#' @importFrom stats dnorm median pnorm rnorm runif
#' 
#' @export
#'
#'
#' @examples 
#' 
#' require(MASS)  
#' # for random number from multivariate normal distribution
#' 
#' n <- 100  # number of individuals
#' p <- 5    # number of variables
#' q1 <- 20    # dimension of the response
#' q2 <- 20    # dimension of the response
#' 
#' ngrid   <- 1000
#' nburnin <- 100
#' nmc     <- 200
#' nthin   <- 5
#' niter   <- nburnin + nmc
#' effsamp <- (niter - nburnin)/nthin
#' 
#' alpha.tt  <- runif(n = 1, min = -1, max = 1)  # intercept term
#' alpha.u1t <- runif(n = 1, min = -1, max = 1)  # intercept term
#' alpha.u2t <- runif(n = 1, min = -1, max = 1)  # intercept term
#' beta.tt   <- runif(n = p, min = -1, max = 1)  # regression parameter 
#' gamma1.t  <- runif(n = q1, min = -1, max = 1)
#' gamma2.t  <- runif(n = q2, min = -1, max = 1)
#' phi.tt    <- 1 
#' phi.u1t   <- 1
#' phi.u2t   <- 1
#' 
#' sigma.tt    <- 1
#' sigma.u1t   <- 1
#' sigma.u2t   <- 1
#' sigma.etat1 <- 1
#' sigma.etat2 <- 1
#' 
#' x <- mvrnorm(n = n, mu = rep(0, p), Sigma = diag(p))
#' 
#' 
#' 
#' eta2 <- rnorm(n = 1, mean = 0, sd = sigma.etat2)
#' eta1 <- rnorm(n = 1, mean = eta2, sd = sigma.etat1)
#' logt <- rnorm(n = n, mean = alpha.tt + x %*% beta.tt + eta1 * phi.tt, 
#' sd = sigma.tt)
#' u1   <- matrix(rnorm(n = n * q1, mean = alpha.u1t + eta1 * phi.u1t, 
#' sd = sigma.u1t), nrow = n, ncol = q1)
#' u2   <- matrix(rnorm(n = n * q2, mean = alpha.u2t + eta2 * phi.u2t, 
#' sd = sigma.u2t), nrow = n, ncol = q2)
#' logt <- rnorm(n = n, mean = alpha.tt + x %*% beta.tt + u1 %*% gamma1.t + 
#' u2 %*% gamma2.t, sd = sigma.tt)
#'   
#' # Survival time generation
#' T <- exp(logt)   # AFT model
#' C <- rgamma(n, shape = 1, rate = 1)  # 50% censor
#' time <- pmin(T, C)  # observed time is min of censored and true
#' status = time == T   # set to 1 if event is observed
#' 1 - sum(status)/length(T)   # censoring rate
#' censor.rate <- 1 - sum(status)/length(T)    # censoring rate
#' censor.rate
#' summary(C)
#' summary(T)
#' ct <- as.matrix(cbind(time = time, status = status))  # censored time
#' logt.grid <- seq(from = min(logt) - 1, to = max(logt) + 1, length.out = ngrid)
#'   
#' index1 <- which(ct[, 2] == 1)  # which are NOT censored
#' ct1    <- ct[index1, ]
#'   
#' posterior.fit.sem <- mcmc(ct, u1, u2, x, nburnin = nburnin, 
#' nmc = nmc, nthin = nthin)
#'   
#' pMean.beta.t   <- posterior.fit.sem$pMean.beta.t
#' pMean.alpha.t  <- posterior.fit.sem$pMean.alpha.t
#' pMean.phi.t    <- posterior.fit.sem$pMean.phi.t
#' pMean.alpha.u1 <- posterior.fit.sem$pMean.alpha.u1
#' pMean.alpha.u2 <- posterior.fit.sem$pMean.alpha.u2
#' pMean.phi.u1   <- posterior.fit.sem$pMean.phi.u1
#' pMean.phi.u2   <- posterior.fit.sem$pMean.phi.u2
#' pMean.eta1     <- posterior.fit.sem$pMean.eta1
#' pMean.eta2     <- posterior.fit.sem$pMean.eta2
#' pMean.logt.hat <- posterior.fit.sem$posterior.fit.sem
#' 
#' pMean.sigma.t.square  <- posterior.fit.sem$pMean.sigma.t.square
#' pMean.sigma.u1.square <- posterior.fit.sem$pMean.sigma.u1.square
#' pMean.sigma.u2.square <- posterior.fit.sem$pMean.sigma.u2.square
#' pMean.logt.hat        <- posterior.fit.sem$pMean.logt.hat
#' 
#' DIC.sem  <- posterior.fit.sem$DIC
#' WAIC.sem <- posterior.fit.sem$WAIC
#' mse.sem  <- mean(pMean.logt.hat[index1] - log(ct1[, 1]))^2
#' 
#' alpha.t.samples         <- posterior.fit.sem$alpha.t.samples
#' beta1.t.samples         <- posterior.fit.sem$beta1.t.samples
#' beta2.t.samples         <- posterior.fit.sem$beta2.t.samples
#' beta.t.samples          <- posterior.fit.sem$beta.t.samples
#' phi.t.samples           <- posterior.fit.sem$phi.t.samples          
#' alpha.u1.samples        <- posterior.fit.sem$alpha.u1.samples
#' alpha.u2.samples        <- posterior.fit.sem$alpha.u2.samples
#' phi.u1.samples          <- posterior.fit.sem$phi.u1.samples
#' phi.u2.samples          <- posterior.fit.sem$phi.u2.samples
#' sigma.t.square.samples  <- posterior.fit.sem$sigma.t.square.samples
#' sigma.u1.square.samples <- posterior.fit.sem$sigma.u1.square.samples
#' sigma.u2.square.samples <- posterior.fit.sem$sigma.u2.square.samples
#' eta1.samples            <- posterior.fit.sem$eta1.samples
#' eta2.samples            <- posterior.fit.sem$eta2.samples
#'   
#' inv.cpo <- matrix(0, nrow = effsamp, ncol = n)  
#' # this will store inverse cpo values
#' log.cpo <- rep(0, n)                        # this will store log cpo  
#' for(iter in 1:effsamp)  # Post burn in
#' {
#'   inv.cpo[iter, ] <- 1/(dnorm(ct[, 1], mean = alpha.t.samples[iter] + 
#'   x %*% beta.t.samples[, iter] + 
#'                                 + eta1.samples[iter] * phi.t.samples[iter], 
#'                               sd = sqrt(sigma.t.square.samples[iter]))^ct[, 2] * 
#'                           pnorm(ct[, 1], mean = alpha.t.samples[iter] + 
#'                           x %*% beta.t.samples[, iter] + 
#'                                   + eta1.samples[iter] * phi.t.samples[iter], 
#'                                 sd = sqrt(sigma.t.square.samples[iter]), 
#'                                 lower.tail = FALSE)^(1 - ct[, 2]))
#' }                   # End of iter loop
#' for (i in 1:n){
#'   log.cpo[i]   <- -log(mean(inv.cpo[, i]))    
#'   # You average invcpo[i] over the iterations,
#'   # then take 1/average and then take log.
#'   # Hence the negative sign in the log
#' }
#' lpml.sem <- sum(log.cpo)
#'   
#'   
#'   
#'   
#' DIC.sem
#' WAIC.sem
#' mse.sem

mcmc <- function(ct, u1, u2, X, nburnin = 1000, 
                 nmc = 2000, nthin = 1)
{
  
  n  <- nrow(X)  # sample size
  p  <- ncol(X)  # dimension
  q1 <- ncol(u1)
  q2 <- ncol(u2)
  
  
  # Survival response
  time         <- ct[, 1]
  status       <- ct[, 2]
  censored.id  <- which(status == 0)
  n.censored   <- length(censored.id)  # number of censored observations
  X.censored   <- X[censored.id, ]
  X.observed   <- X[-censored.id, ]
  logt <- logtime <- log(time)   # for coding convenience, since the whole code is written with y
  logt.censored <- logt[censored.id]
  logt.observed <- logt[-censored.id]
  
  
  alpha.t  <- runif(n = 1, min = -1, max = 1)  # intercept term
  alpha.u1 <- runif(n = q1, min = -1, max = 1)  # intercept term
  alpha.u2 <- runif(n = q2, min = -1, max = 1)  # intercept term
  beta.t   <- runif(n = p, min = -1, max = 1)  # regression parameter 
  phi.t    <- 1
  phi.u1   <- rep(1, q1)
  phi.u2   <- rep(1, q2)
  
  sigma.t.square    <- 1
  sigma.u1.square   <- rep(1, q1)
  sigma.u2.square   <- rep(1, q2)
  sigma.eta1.square <- 1
  sigma.eta2.square <- 1
  
  alpha.sigmau <- 0.5
  beta.sigmau  <- 0.5
  
  eta2 <- rnorm(n = 1, mean = 0, sd = 1)
  eta1 <- rnorm(n = 1, mean = eta2, sd = 1)
  
  nburnin <- nburnin
  nmc     <- nmc
  niter   <- nburnin + nmc
  effsamp <- (niter - nburnin)/nthin
  
  # output
  beta.tout   <- matrix(0, nrow = p, ncol = effsamp)
  alpha.tout  <- rep(0, effsamp)
  phi.tout    <- rep(0, effsamp)
  alpha.u1out <- matrix(0, nrow = q1, ncol = effsamp)
  alpha.u2out <- matrix(0, nrow = q2, ncol = effsamp)
  phi.u1out   <- matrix(0, nrow = q1, ncol = effsamp)
  phi.u2out   <- matrix(0, nrow = q2, ncol = effsamp)
  eta1.out    <- rep(0, effsamp)
  eta2.out    <- rep(0, effsamp)
  
  sigma.t.square.out  <- rep(0, effsamp)
  sigma.u1.square.out <- rep(1, effsamp)
  sigma.u2.square.out <- rep(1, effsamp)
  
  logt.out          <- matrix(0, nrow = n, ncol = effsamp)
  logt.hat.out      <- matrix(0, nrow = n, ncol = effsamp)
  loglikelihood.out <- rep(0, effsamp)
  likelihood.out    <- matrix(0, nrow = n, ncol = effsamp)
  
  for(iter in 1:niter)  # MCMC
  {
    
    ## Update survival latent variable ##
    mean.impute       <- alpha.t + as.matrix(X.censored) %*% as.matrix(beta.t) + eta1 * phi.t
    sd.impute         <- sqrt(sigma.t.square)
    time.censored     <- msm::rtnorm(n.censored, mean = mean.impute, sd = sd.impute, lower = logt.censored)
    logt[censored.id] <- time.censored  # truncated at log(time) for censored data
    
    
    # Sample $ \beta.t $
    A            <- crossprod(x = X) + chol2inv(chol(100 * diag(p)))
    Ainv         <- chol2inv(chol(A))
    Sigma.beta.t <- sigma.t.square * Ainv
    mean.beta.t  <- as.vector(Ainv %*% t(X) %*% (logt - alpha.t - eta1 * phi.t))
    beta.t       <- as.vector(MASS::mvrnorm(n = 1, mu = mean.beta.t, Sigma = Sigma.beta.t))
    
    # Sample $ \alpha_t $
    A             <- crossprod(x = rep(1, n)) + chol2inv(chol(diag(1)))
    Ainv          <- chol2inv(chol(A))
    Sigma.alpha.t <- sigma.t.square * Ainv
    mean.alpha.t  <- as.vector(Ainv %*% t(rep(1, n)) %*% (logt - X %*% beta.t - eta1 * phi.t))
    alpha.t       <- as.vector(MASS::mvrnorm(n = 1, mu = mean.alpha.t, Sigma = Sigma.alpha.t))
    
    # Sample $ \phi_t $
    A           <- crossprod(x = rep(eta1, n)) + chol2inv(chol(diag(1)))
    Ainv        <- chol2inv(chol(A))
    Sigma.phi.t <- Ainv
    mean.phi.t  <- as.vector(Ainv %*% t(rep(eta1, n)) %*% (logt - alpha.t - X %*% beta.t))
    phi.t       <- as.vector(MASS::mvrnorm(n = 1, mu = mean.phi.t, Sigma = Sigma.phi.t))
    
    # Sample $ \alpha_u1 $
    A             <- crossprod(x = rep(1, n)) + chol2inv(chol(diag(1)))
    Ainv          <- chol2inv(chol(A))
    Sigma.alpha.u <- sigma.u1.square * Ainv
    mean.alpha.u  <- as.vector(Ainv %*% t(rep(1, n)) %*% (u1 - eta1 * phi.u1))
    alpha.u1      <- as.vector(rnorm(n = q1, mean = mean.alpha.u, sd = sqrt(Sigma.alpha.u)))
    
    # Sample $ \alpha_u2 $
    A             <- crossprod(x = rep(1, n)) + chol2inv(chol(diag(1)))
    Ainv          <- chol2inv(chol(A))
    Sigma.alpha.u <- as.vector(sigma.u2.square) * Ainv
    mean.alpha.u  <- as.vector(Ainv %*% t(rep(1, n)) %*% (u2 - eta2 * phi.u2))
    alpha.u2      <- as.vector(rnorm(n = q2, mean = mean.alpha.u, sd = sqrt(Sigma.alpha.u)))
    
    # Sample $ \phi_u1 $
    A           <- crossprod(x = rep(eta1, n)) + chol2inv(chol(diag(1)))
    Ainv        <- chol2inv(chol(A))
    Sigma.phi.u <- Ainv
    mean.phi.u  <- as.vector(Ainv %*% t(rep(eta1, n)) %*% (u1 - alpha.u1))
    phi.u1      <- as.vector(rnorm(n = q1, mean = mean.phi.u, sd = sqrt(Sigma.phi.u)))
    
    # Sample $ \phi_u2 $
    A           <- crossprod(x = rep(eta2, n)) + chol2inv(chol(diag(1)))
    Ainv        <- chol2inv(chol(A))
    Sigma.phi.u <- Ainv
    mean.phi.u  <- as.vector(Ainv %*% t(rep(eta2, n)) %*% (u2 - alpha.u2))
    phi.u2      <- as.vector(rnorm(n = q2, mean = mean.phi.u, sd = sqrt(Sigma.phi.u)))
    
    
    # Sample $ \eta1 $
    sum.phi.eta1 <- 0
    for(k in 1:q1)
    {
      sum.phi.eta1 <- sum.phi.eta1 + sum((phi.u1[k] * u1[, k])/sigma.u1.square[k])
    }
    Sigma.eta <- 1/(1/sigma.eta1.square + sum(phi.u1^2/sigma.u1.square))
    mean.eta  <- Sigma.eta * (eta2/sigma.eta1.square + sum.phi.eta1 +
                                phi.t * sum(logt)/sigma.t.square - sum((phi.u1 * alpha.u1)/sigma.u1.square) -
                                phi.t * alpha.t/sigma.t.square -
                                (t(beta.t) %*% t(X) %*% rep(phi.t, n))/sigma.t.square)
    eta1      <- as.vector(MASS::mvrnorm(n = 1, mu = mean.eta, Sigma = Sigma.eta))
    
    # Sample $ \eta2 $
    sum.phi.eta2 <- 0
    for(l in 1:q2)
    {
      sum.phi.eta2 <- sum.phi.eta2 + sum((phi.u2[l] * u2[, l])/sigma.u2.square[l])
    }
    Sigma.eta <- 1/(1/sigma.eta1.square + 1/sigma.eta2.square + sum(phi.u2^2/sigma.u2.square))
    mean.eta  <- Sigma.eta * (eta1/sigma.eta1.square + sum.phi.eta2 +
                                sum((phi.u2 * alpha.u2)/sigma.u2.square))
    eta2      <- as.vector(MASS::mvrnorm(n = 1, mu = mean.eta, Sigma = Sigma.eta))
    
    
    
    
    # Sample $ \sigma_t^2
    E1 <- crossprod(x = logt - alpha.t - X %*% beta.t  - eta1 * phi.t)
    E2 <- crossprod(x = beta.t)
    E3 <- crossprod(x = alpha.t)
    # E4 <- crossprod(x = phi.t)
    sigma.t.square <- 1/stats::rgamma(1, shape = (n + p + 1)/2, scale = 2/(E1 + E2 + E3))
    # if(sigma.t.square.propesed < 10)
    # {
    #   sigma.t.square <- sigma.t.square.propesed
    # }
    
    # Sample $ \sigma_u1^2
    E1 <- rep(NA, q1)
    for(k in 1:q1)
    {
      E1[k] <- crossprod(x = u1[k] - alpha.u1 - eta1 * phi.u1)
    }
    
    E3 <- crossprod(x = alpha.u1)
    E4 <- crossprod(x = phi.u1)
    sigma.u1.square <- 1/stats::rgamma(q1, alpha.sigmau + (n + 1)/2, 
                                       scale = 1/(beta.sigmau + (as.vector(E1) + as.vector(E3) + as.vector(E4))/2))
    
    # Sample $ \sigma_u2^2
    E1 <- rep(NA, q2)
    for(k in 1:q2)
    {
      E1[k] <- crossprod(x = u2[k] - alpha.u2 - eta2 * phi.u2)
    }
    
    E3 <- crossprod(x = alpha.u2)
    E4 <- crossprod(x = phi.u2)
    sigma.u2.square <- 1/stats::rgamma(q2, alpha.sigmau + (n + 1)/2, 
                                       scale = 1/(beta.sigmau + (as.vector(E1) + as.vector(E3) + as.vector(E4))/2))
    
    
    logt.hat <- alpha.t + X %*% beta.t + eta1 * phi.t
    
    loglikelihood <- sum(status * dnorm(logt, mean = alpha.t + X %*% beta.t + eta1 * phi.t, 
                                        sd = sqrt(sigma.t.square), log = TRUE) + 
                           (1 - status) * pnorm(logt, mean = alpha.t + X %*% beta.t + eta1 * phi.t, 
                                                sd = sqrt(sigma.t.square), 
                                                lower.tail = FALSE, log.p = TRUE)) 
    likelihood    <- exp(loglikelihood)
    
    
    # stores the MCMC samples after discarding burnin samples
    if (iter > nburnin)
    {
      beta.tout[, (iter - nburnin)/nthin]   <- beta.t
      alpha.tout[(iter - nburnin)/nthin]    <- alpha.t
      phi.tout[(iter - nburnin)/nthin]      <- phi.t
      alpha.u1out[, (iter - nburnin)/nthin] <- alpha.u1
      alpha.u2out[, (iter - nburnin)/nthin] <- alpha.u2
      phi.u1out[, (iter - nburnin)/nthin]   <- phi.u1
      phi.u2out[, (iter - nburnin)/nthin]   <- phi.u2
      eta1.out[(iter - nburnin)/nthin]      <- eta1
      eta2.out[(iter - nburnin)/nthin]      <- eta2
      
      sigma.t.square.out[(iter - nburnin)/nthin]  <- sigma.t.square
      # sigma.u1.square.out[(iter - nburnin)/nthin] <- sigma.u1.square
      # sigma.u2.square.out[(iter - nburnin)/nthin] <- sigma.u2.square
      
      logt.out[, (iter - nburnin)/nthin]        <- logt
      logt.hat.out[, (iter - nburnin)/nthin]    <- logt.hat
      loglikelihood.out[(iter - nburnin)/nthin] <- loglikelihood
      likelihood.out[, (iter - nburnin)/nthin]  <- likelihood
    }
  }
  
  
  
  # Posterior Mean
  pMean.beta.t   <- apply(beta.tout, 1, mean)
  pMean.alpha.t  <- mean(alpha.tout)
  pMean.phi.t    <- mean(phi.tout)
  pMean.alpha.u1 <- apply(alpha.u1out, 1, mean)
  pMean.alpha.u2 <- apply(alpha.u2out, 1, mean)
  pMean.phi.u1   <- apply(phi.u1out, 1, mean)
  pMean.phi.u2   <- apply(phi.u2out, 1, mean)
  pMean.eta1     <- mean(eta1.out)
  pMean.eta2     <- mean(eta2.out)
  
  pMean.sigma.t.square  <- median(sigma.t.square.out)
  # pMean.sigma.u1.square <- mean(sigma.u1.square.out)
  # pMean.sigma.u2.square <- mean(sigma.u2.square.out)
  
  pMean.logt     <- apply(logt.out, 1, mean)
  pMean.logt.hat <- apply(logt.hat.out, 1, mean)
  pLoglikelihood <- mean(loglikelihood.out)
  plikelihood    <- apply(likelihood.out, 1, mean)
  
  
  loglikelihood.posterior <- sum(status * dnorm(pMean.logt, mean = pMean.alpha.t + X %*% pMean.beta.t + 
                                                  pMean.eta1 * pMean.phi.t, 
                                                sd = sqrt(pMean.sigma.t.square), log = TRUE) + 
                                   (1 - status) * pnorm(pMean.logt, mean = pMean.alpha.t + X %*% pMean.beta.t + 
                                                          pMean.eta1 * pMean.phi.t, 
                                                        sd = sqrt(pMean.sigma.t.square), 
                                                        lower.tail = FALSE, log.p = TRUE))
  
  DIC  <- -4 * pLoglikelihood + 2 * loglikelihood.posterior
  lppd <- sum(log(plikelihood))
  WAIC <- -2 * (lppd - 2 * (loglikelihood.posterior - pLoglikelihood))
  # WAIC <- -2 * (lppd - 2 * (lppd - pLoglikelihood))
  
  
  result = list(pMean.beta.t = pMean.beta.t, pMean.alpha.t = pMean.alpha.t, pMean.phi.t = pMean.phi.t,
                pMean.alpha.u1 = pMean.alpha.u1, pMean.alpha.u2 = pMean.alpha.u2, 
                pMean.phi.u1 = pMean.phi.u1, pMean.phi.u2 = pMean.phi.u2, 
                pMean.eta1 = pMean.eta1, pMean.eta2 = pMean.eta2, 
                pMean.sigma.t.square = pMean.sigma.t.square, 
                # pMean.sigma.u1.square = pMean.sigma.u1.square, pMean.sigma.u2.square = pMean.sigma.u2.square,
                alpha.t.samples = alpha.tout, phi.t.samples = phi.tout, 
                beta1.t.samples = beta.tout[1, ], #beta2.t.samples = beta.tout[2, ],
                beta.t.samples = beta.tout,
                alpha.u1.samples = alpha.u1out, alpha.u2.samples = alpha.u2out,
                phi.u1.samples = phi.u1out, phi.u2.samples = phi.u2out, 
                eta1.samples = eta1.out, eta2.samples = eta2.out,
                sigma.t.square.samples = sigma.t.square.out, 
                sigma.u1.square.samples = sigma.u1.square.out, sigma.u2.square.samples = sigma.u2.square.out,
                pMean.logt.hat = pMean.logt.hat, DIC = DIC, WAIC = WAIC)
  
  return(result)
}