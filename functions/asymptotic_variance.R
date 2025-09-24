library(latex2exp)
library(ggplot2)
library(mev)
library(matlib)
library(MASS)

options(scipen = 999)


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


### Quantities involved in the asymptotic variance ###

### J_\alpha

function_J_alpha <- function(x, mu, sigma, xi, alpha){
  s1 <- S_mu(x, mu, sigma, xi)
  s2 <- S_sigma(x, mu, sigma, xi)
  s3 <- S_xi(x, mu, sigma, xi)
  d  <- dgev(x = x, loc = mu, scale = sigma, shape = xi)^(1 + alpha)
  
  M <- c(s1, s2, s3) %*% t(c(s1, s2, s3)) * d
  return(M)
}

J_alpha <- function(mu, sigma, xi, alpha){
  
  J <- matrix(nrow = 3, ncol = 3)
  
  if (xi >= 0.01){
    lower <- mu - sigma/xi+0.001
    upper <- mu+100
  } else if (abs(xi) < 0.01){
    lower <- mu-100
    upper <- mu+100
  } else {
    lower <- mu-100
    upper <- mu + sigma/abs(xi)-0.001
  }
  
  for (i in 1:3){
    for (j in 1:3){
      integrand_ij <- function(x) {
        sapply(x, function(xx) {
          M <- function_J_alpha(xx, mu, sigma, xi, alpha)
          M[i, j]
        })
      }
      res <- tryCatch(
        integrate(integrand_ij, lower, upper, stop.on.error = T),
        error = function(e) list(value = NA, abs.error = NA)
      )
      J[i, j] <- res$value
    }
  }
  
  return(J)
}


### U_\alpha

function_U_alpha <- function(x, mu, sigma, xi, alpha){
  s1 <- S_mu(x, mu, sigma, xi)
  s2 <- S_sigma(x, mu, sigma, xi)
  s3 <- S_xi(x, mu, sigma, xi)
  d  <- dgev(x = x, loc = mu, scale = sigma, shape = xi)^(1 + alpha)
  
  # for scalar x returns length-3 vector
  M <- c(s1, s2, s3) * d
  return(M)
}

U_alpha <- function(mu, sigma, xi, alpha){
  
  U <- numeric(3)
  
  if (xi >= 0.01){
    lower <- mu - sigma/xi+0.001
    upper <- mu+100
  } else if (abs(xi) < 0.01){
    lower <- mu-100
    upper <- mu+100
  } else {
    lower <- mu-100
    upper <- mu - sigma/xi-0.001
  }
  
  for (i in 1:3){
    integrand_i <- function(x) {
      sapply(x, function(xx) {
        M <- function_U_alpha(xx, mu, sigma, xi, alpha)
        M[i]  
      })
    }
    res <- tryCatch(
      integrate(integrand_i, lower, upper, stop.on.error = FALSE),
      error = function(e) list(value = NA, abs.error = NA)
    )
    U[i] <- res$value
  }
  
  return(U)
}


### K_\alpha

K_alpha <- function(mu, sigma, xi, alpha){
  return(J_alpha(mu, sigma, xi, 2*alpha)-U_alpha(mu, sigma, xi, alpha) %*% t(U_alpha(mu, sigma, xi, alpha)))
}

### Asymptotic variance MDPDE ###

asymptotic_variance_mdpde <- function(mu, sigma, xi, alpha){
  
  J <- J_alpha(mu, sigma, xi, alpha)
  
  inv_J_alpha <- tryCatch(
    ginv(J),
    error=function(e) {
      matrix(NA, nrow=nrow(J), ncol=ncol(J))
    }
  )
  
  return(inv_J_alpha %*% K_alpha(mu, sigma, xi, alpha) %*% inv_J_alpha)
}

asymptotic_variance_mdpde_componentwise <- function(mu, sigma, xi, alpha,coord=c(i,j)){
  
  asymptotic_variance_mdpde <- asymptotic_variance_mdpde(mu, sigma, xi, alpha)
  
  return(asymptotic_variance_mdpde[coord[1],coord[2]])
}


