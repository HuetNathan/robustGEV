### R code for influence functions plots ###

library(mev)

### Score functions of the GEV ###

# for the location parameter

S_mu <- Vectorize(function(x, mu, sigma, xi) {
  
  z <- (x - mu) / sigma
  term <- 1 + xi * z
  
  if (abs(xi) >= 0.01){
    result <- (1 / sigma) * (term)^(-1) * (xi + 1 - term^(-1 / xi))
  }
  
  if (abs(xi) < 0.01){
    result <- (1 / sigma) * (1-exp(-z))
  }
  
  return(result)
}, vectorize.args = c("x", "mu", "sigma", "xi"))

# for the scale parameter

S_sigma <- Vectorize(function(x, mu, sigma, xi) {
  
  z <- (x - mu) / sigma
  term <- 1 + xi * z
  
  if (abs(xi) >= 0.01){
    result <- -1 / sigma + (x - mu) / sigma^2 * (term)^(-1) * ((xi + 1) - term^(-1 / xi))
  }
  
  if (abs(xi) < 0.01){
    result <- (1 / sigma) *(-1 + z*(1-exp(-z)))
  }
  
  return(result)
}, vectorize.args = c("x", "mu", "sigma", "xi"))

# for the shape parameter

S_xi <- Vectorize(function(x, mu, sigma, xi) {
  
  z <- (x - mu) / sigma
  term <- 1 + xi * z
  
  if (abs(xi) >= 0.01){
    
    part1 <- (1/xi^2)*log(term)*(1-term^(-1/xi))
    
    part2 <- (1/xi)*z*term^(-1)*(term^(-1/xi)-1)
    
    part3 <- -z*term^(-1)
    
    result <- part1+part2+part3
    
  }
  
  if (abs(xi) < 0.01){
    
    result <- (1/2)*z^2*(1-exp(-z))-z
    
  }
  return(result)
}, vectorize.args = c("x", "mu", "sigma", "xi"))

# Global score

S <- Vectorize(function(x, mu, sigma, xi) {
  c(S_mu(x, mu, sigma, xi),S_sigma(x, mu, sigma, xi),S_xi(x, mu, sigma, xi))
}, vectorize.args = c("x", "mu", "sigma", "xi"))


### Influence functions of the GEV ###

# for the location parameter

function_denominator_mu <- function(x, mu, sigma, xi, alpha){
  return(S_mu(x, mu, sigma, xi)^2*dgev(x = x,loc = mu, scale = sigma,shape = xi)^(1+alpha))
}

function_numerator_mu <- function(x, mu, sigma, xi, alpha){
  return(S_mu(x, mu, sigma, xi)*dgev(x = x,loc = mu, scale = sigma,shape = xi)^(1+alpha))
}

IF_mu <- function(x, mu, sigma, xi, alpha){
  
  if (xi >= 0.01){
    cste_denominator <- integrate(function_denominator_mu,
                                  mu-sigma/xi+0.001, mu+100,
                                  mu = mu,sigma=sigma,xi=xi,alpha=alpha)$value
    
    cste_numerator <- integrate(function_numerator_mu,
                                mu-sigma/xi, mu+100,
                                mu = mu,sigma=sigma,xi=xi,alpha=alpha)$value
    return((1/cste_denominator)*(S_mu(x, mu, sigma, xi)*dgev(x = x,loc = mu, scale = sigma,shape = xi)^alpha - cste_numerator))
  }
  
  if (abs(xi) < 0.01){
    cste_denominator <- integrate(function_denominator_mu,
                                  mu-100, mu+100,
                                  mu = mu,sigma=sigma,xi=0,alpha=alpha)$value
    
    cste_numerator <- integrate(function_numerator_mu,
                                mu-100, mu+100,
                                mu = mu,sigma=sigma,xi=0,alpha=alpha)$value
    return((1/cste_denominator)*(S_mu(x, mu, sigma, 0)*dgev(x = x,loc = mu, scale = sigma,shape = 0)^alpha - cste_numerator))
  }
  
  if (xi <= -0.01){
    cste_denominator <- integrate(function_denominator_mu,
                                  mu-100, mu-sigma/xi-0.001,
                                  mu = mu,sigma=sigma,xi=xi,alpha=alpha)$value
    
    cste_numerator <- integrate(function_numerator_mu,
                                -Inf,mu-sigma/xi,
                                mu = mu,sigma=sigma,xi=xi,alpha=alpha)$value
    
    return((1/cste_denominator)*(S_mu(x, mu, sigma, xi)*dgev(x = x,loc = mu, scale = sigma,shape = xi)^alpha - cste_numerator))
    
  }
}

# for the scale parameter

function_denominator_sigma <- function(x, mu, sigma, xi, alpha){
  return(S_sigma(x, mu, sigma, xi)^2*dgev(x = x,loc = mu, scale = sigma,shape = xi)^(1+alpha))
}

function_numerator_sigma <- function(x, mu, sigma, xi, alpha){
  return(S_sigma(x, mu, sigma, xi)*dgev(x = x,loc = mu, scale = sigma,shape = xi)^(1+alpha))
}

IF_sigma <- function(x, mu, sigma, xi, alpha){
  
  if (xi >= 0.01){
    cste_denominator <- integrate(function_denominator_sigma,
                                  mu-sigma/xi+0.001, mu+100,
                                  mu = mu,sigma=sigma,xi=xi,alpha=alpha)$value
    
    cste_numerator <- integrate(function_numerator_sigma,
                                mu-sigma/xi+0.001, mu+100,
                                mu = mu,sigma=sigma,xi=xi,alpha=alpha)$value
    return((1/cste_denominator)*(S_sigma(x, mu, sigma, xi)*dgev(x = x,loc = mu, scale = sigma,shape = xi)^alpha - cste_numerator))
  }
  
  if (abs(xi) < 0.01){
    cste_denominator <- integrate(function_denominator_sigma,
                                  -Inf, Inf,
                                  mu = mu,sigma=sigma,xi=0,alpha=alpha)$value
    
    cste_numerator <- integrate(function_numerator_sigma,
                                mu-100, mu+100,
                                mu = mu,sigma=sigma,xi=0,alpha=alpha)$value
    return((1/cste_denominator)*(S_sigma(x, mu, sigma, 0)*dgev(x = x,loc = mu, scale = sigma,shape = 0)^alpha - cste_numerator))
  }
  
  if (xi <= -0.01){
    cste_denominator <- integrate(function_denominator_sigma,
                                  mu-100, mu-sigma/xi-0.001,
                                  mu = mu,sigma=sigma,xi=xi,alpha=alpha)$value
    
    cste_numerator <- integrate(function_numerator_sigma,
                                mu-100,mu-sigma/xi-0.001,
                                mu = mu,sigma=sigma,xi=xi,alpha=alpha)$value
    
    return((1/cste_denominator)*(S_sigma(x, mu, sigma, xi)*dgev(x = x,loc = mu, scale = sigma,shape = xi)^alpha - cste_numerator))
    
  }
}

# for the shape parameter

function_denominator_xi <- function(x, mu, sigma, xi, alpha){
  return(S_xi(x, mu, sigma, xi)^2*mev::dgev(x = x,loc = mu, scale = sigma,shape = xi)^(1+alpha))
}

function_numerator_xi <- function(x, mu, sigma, xi, alpha){
  return(S_xi(x, mu, sigma, xi)*mev::dgev(x = x,loc = mu, scale = sigma,shape = xi)^(1+alpha))
}

IF_xi <- function(x, mu, sigma, xi, alpha){

  if (xi >= 0.01){
    cste_denominator <- integrate(function_denominator_xi,
                                  mu-sigma/xi+0.001, mu+100,
                                  mu = mu,sigma=sigma,xi=xi,alpha=alpha)$value
    
    cste_numerator <- integrate(function_numerator_xi,
                                mu-sigma/xi+0.001, mu+100,
                                mu = mu,sigma=sigma,xi=xi,alpha=alpha)$value
    return((1/cste_denominator)*(S_xi(x, mu, sigma, xi)*mev::dgev(x = x,loc = mu, scale = sigma,shape = xi)^alpha - cste_numerator))
  }
  
  if (abs(xi) < 0.01){
    cste_denominator <- integrate(function_denominator_xi,
                                  mu-100, mu+100,
                                  mu = mu,sigma=sigma,xi=0,alpha=alpha)$value
    
    cste_numerator <- integrate(function_numerator_xi,
                                mu-100, mu+100,
                                mu = mu,sigma=sigma,xi=0,alpha=alpha)$value
    return((1/cste_denominator)*(S_xi(x, mu, sigma, 0)*dgev(x = x,loc = mu, scale = sigma,shape = 0)^alpha - cste_numerator))
  }
  
  if (xi <= -0.01){
    cste_denominator <- integrate(function_denominator_xi,
                                  mu-100, mu-sigma/xi-0.001,
                                  mu = mu,sigma=sigma,xi=xi,alpha=alpha)$value
    
    cste_numerator <- integrate(function_numerator_xi,
                                mu-100,mu-sigma/xi-0.001,
                                mu = mu,sigma=sigma,xi=xi,alpha=alpha)$value
    
    return((1/cste_denominator)*(S_xi(x, mu, sigma, xi)*mev::dgev(x = x,loc = mu, scale = sigma,shape = xi)^alpha - cste_numerator))
  }
}
