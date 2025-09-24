
R code for influence functions plots

```{r}

library(latex2exp)
library(ggplot2)
library(mev)

```

Score functions of the GEV

```{r}

# for the location parameter

S_mu <- function(x, mu, sigma, xi) {
  
  z <- (x - mu) / sigma
  term <- 1 + xi * z
  
  result <- (1 / sigma) * (term)^(-1) * (xi + 1 - term^(-1 / xi))
  
  return(result)
}

# for the scale parameter

S_sigma <- function(x, mu, sigma, xi) {

  z <- (x - mu) / sigma
  term <- 1 + xi * z
  
  result <- -1 / sigma + 
    (x - mu) / sigma^2 * 
    (term)^(-1) * 
    ((xi + 1) - term^(-1 / xi))
  
  return(result)
}

# for the shape parameter

S_xi <- function(x, mu, sigma, xi) {

  z <- (x - mu) / sigma
  term <- 1 + xi * z
  
  log_term <- log(term)
  pow_term <- term^(-1 / xi)
  
  part1 <- (1 / xi^2) * log_term * (1 - pow_term)
  part2 <- (1 / xi) * (term^(-1)) * ((xi + 1) - pow_term)
  
  result <- part1 + part2
  
  return(result)
}
```

Influence functions of the GEV

```{r}
# for the location parameter

function_denominator_mu <- function(x, mu, sigma, xi, alpha){
  return(S_mu(x, mu, sigma, xi)^2*dgev(x = x,loc = mu, scale = sigma,shape = xi)^(1+alpha))
}

function_numerator_mu <- function(x, mu, sigma, xi, alpha){
  return(S_mu(x, mu, sigma, xi)*dgev(x = x,loc = mu, scale = sigma,shape = xi)^(1+alpha))
}

IF_mu <- function(x, mu, sigma, xi, alpha){
  # no xi = 0
  
  if (xi>0){
    cste_denominator <- integrate(function_denominator_mu,
                                  mu-sigma/xi, Inf,
                                  mu = mu,sigma=sigma,xi=xi,alpha=alpha)$value
    
    cste_numerator <- integrate(function_numerator_mu,
                                mu-sigma/xi, Inf,
                                mu = mu,sigma=sigma,xi=xi,alpha=alpha)$value
    return((1/cste_denominator)*(S_mu(x, mu, sigma, xi)*dgev(x = x,loc = mu, scale = sigma,shape = xi)^alpha - cste_numerator))
  }
  if (xi<0){
    cste_denominator <- integrate(function_denominator_mu,
                                  -Inf, mu-sigma/xi,
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
  # no xi = 0
  
  if (xi>0){
    cste_denominator <- integrate(function_denominator_sigma,
                                  mu-sigma/xi, Inf,
                                  mu = mu,sigma=sigma,xi=xi,alpha=alpha)$value
    
    cste_numerator <- integrate(function_numerator_sigma,
                                mu-sigma/xi, Inf,
                                mu = mu,sigma=sigma,xi=xi,alpha=alpha)$value
    return((1/cste_denominator)*(S_sigma(x, mu, sigma, xi)*dgev(x = x,loc = mu, scale = sigma,shape = xi)^alpha - cste_numerator))
  }
  if (xi<0){
    cste_denominator <- integrate(function_denominator_sigma,
                                  -Inf, mu-sigma/xi,
                                  mu = mu,sigma=sigma,xi=xi,alpha=alpha)$value
    
    cste_numerator <- integrate(function_numerator_sigma,
                                -Inf,mu-sigma/xi,
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
  # no xi = 0
  
  if (xi>0){
    cste_denominator <- integrate(function_denominator_xi,
                                  mu-sigma/xi, Inf,
                                  mu = mu,sigma=sigma,xi=xi,alpha=alpha)$value
    
    cste_numerator <- integrate(function_numerator_xi,
                                  mu-sigma/xi, Inf,
                                  mu = mu,sigma=sigma,xi=xi,alpha=alpha)$value
    return((1/cste_denominator)*(S_xi(x, mu, sigma, xi)*mev::dgev(x = x,loc = mu, scale = sigma,shape = xi)^alpha - cste_numerator))
  }
  if (xi<0){
    cste_denominator <- integrate(function_denominator_xi,
                                  -Inf, mu-sigma/xi,
                                  mu = mu,sigma=sigma,xi=xi,alpha=alpha)$value
    
    cste_numerator <- integrate(function_numerator_xi,
                                -Inf,mu-sigma/xi,
                                mu = mu,sigma=sigma,xi=xi,alpha=alpha)$value
    
    return((1/cste_denominator)*(S_xi(x, mu, sigma, xi)*mev::dgev(x = x,loc = mu, scale = sigma,shape = xi)^alpha - cste_numerator))
    
  }
}

```

Plots influence functions for positive shape parameter (Figure 2)

```{r}

# Parameters

mu <- 0
sigma <- 1
xi <-0.3
alphas <- c(0.0, 0.1, 0.3, 0.5, 1.0)
```

```{r}

# plot influence function location parameter (Figure 2, left)

p <- c(0.001,0.1,0.5,0.8,0.9)
list_quantile <- mev::qgev(p, loc = mu, scale = sigma, shape = xi, lower.tail = TRUE, log.p = FALSE)

x_mu <- seq(mu - sigma/xi + 1.8, 4, by = 0.1)  

plot_data <- data.frame()

for (alpha in alphas) {
  y_vals <- IF_mu(x_mu, mu, sigma, xi, alpha)
  df <- data.frame(x = x_mu, y = y_vals, alpha = as.factor(alpha))
  plot_data <- rbind(plot_data, df)
}

plot_if <- ggplot(plot_data, aes(x = x, y = y, color = alpha)) +
  geom_line(size = 1.2) +
  theme_minimal()+
    scale_colour_viridis_d(option = "plasma",
                           name = expression(~~~~alpha))+
  labs(
    x = TeX("$q$"),
    y = TeX("$\\textit{IF}_{\\alpha,\\mu}(q,\\mu_0,\\sigma_0,\\xi_0)$")
  ) +
  scale_x_continuous(
    breaks = list_quantile,
    labels = p
  ) +
  theme(legend.position = "none",
        legend.title=element_text(size=26,face="bold"),
        legend.text =element_text(size=24))+
  theme(
    axis.title = element_text(size = 24),   
    axis.text = element_text(size = 18)
  )
plot_if
```

```{r}
ggsave(plot = plot_if,
       filename = paste0("figures/IF_mu_positive_xi.png"),
       width = 1760 / 180,
       height = 1060 / 180,
       dpi = 180,
       units = "in")

```

```{r}

# plot influence function scale parameter (Figure 2, middle)

p <- c(0.001,0.2,0.5,0.8,0.9)
list_quantile <- mev::qgev(p, loc = mu, scale = sigma, shape = xi, lower.tail = TRUE, log.p = FALSE)


x_sig <- seq(mu - sigma/xi + 1.6, 10, by = 0.1)  

plot_data <- data.frame()

for (alpha in alphas) {
  y_vals <- IF_sigma(x_sig, mu, sigma, xi, alpha)
  df <- data.frame(x = x_sig, y = y_vals, alpha = as.factor(alpha))
  plot_data <- rbind(plot_data, df)
}

plot_if <- ggplot(plot_data, aes(x = x, y = y, color = alpha)) +
  geom_line(size = 1.2) +
  theme_minimal() +
    scale_colour_viridis_d(option = "plasma",
                           name = expression(~~~~alpha))+
  labs(
    x = TeX("$q$"),
    y = TeX("$\\textit{IF}_{\\alpha,\\sigma}(q,\\mu_0,\\sigma_0,\\xi_0)$")
  )+
  scale_x_continuous(
    breaks = list_quantile,
    labels = p
  ) +
  theme(legend.position ="none",
        legend.title=element_text(size=26,face="bold"),
        legend.text =element_text(size=24))+
  theme(
    axis.title = element_text(size = 24),   
    axis.text = element_text(size = 18)
  )
plot_if
```

```{r}
ggsave(plot = plot_if,
       filename = paste0("figures/IF_sig_positive_xi.png"),
       width = 1760 / 180,
       height = 1060 / 180,
       dpi = 180,
       units = "in")

```

```{r}

# plot influence function shape parameter (Figure 2, right)

p <- c(0.001,0.1,0.5,0.8,0.9,0.95)
list_quantile <- mev::qgev(p, loc = mu, scale = sigma, shape = xi, lower.tail = TRUE, log.p = FALSE)

x_xi <- seq(mu - sigma/xi + 1.4, 7, by = 0.1)  

plot_data <- data.frame()

for (alpha in alphas) {
  y_vals <- IF_xi(x_xi, mu, sigma, xi, alpha)
  df <- data.frame(x = x_xi, y = y_vals, alpha = as.factor(alpha))
  plot_data <- rbind(plot_data, df)
}

plot_if <- ggplot(plot_data, aes(x = x, y = y, color = alpha)) +
  geom_line(size = 1.2) +
  theme_minimal() +
  scale_colour_viridis_d(option = "plasma",
                         name = expression(~~~~alpha)) +
  labs(
    x = TeX("$q$"),
    y = TeX("$\\textit{IF}_{\\alpha,\\xi}(q,\\mu_0,\\sigma_0,\\xi_0)$")
  ) +
  scale_x_continuous(
    breaks = list_quantile,
    labels = p
  ) +
  #scale_y_log10() +  # <-- log scale for y-axis
  theme(
    legend.position = c(.7, .7),
    legend.title = element_text(size = 26, face = "bold"),
    legend.text = element_text(size = 24),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 18)
  )
plot_if
```

```{r}
ggsave(plot = plot_if,
       filename = paste0("figures/IF_xi_positive_xi.png"),
       width = 1760 / 180,
       height = 1060 / 180,
       dpi = 180,
       units = "in")

```

Plots influence functions for negative shape parameter (Figure 1)

```{r}

# Parameters

mu <- 0
sigma <- 1
xi <- -0.3
alphas <- c(0.0, 0.1, 0.3, 0.5, 1.0)
```

```{r}

# plot influence function location parameter (Figure 1, left)

p <- c(0.001,0.1,0.5,0.8,0.9)
list_quantile <- mev::qgev(p, loc = mu, scale = sigma, shape = xi, lower.tail = TRUE, log.p = FALSE)

x_mu <- seq(-3, mu - sigma/xi - 0.4, by = 0.1) 

plot_data <- data.frame()

for (alpha in alphas) {
  y_vals <- IF_mu(x_mu, mu, sigma, xi, alpha)
  df <- data.frame(x = x_mu, y = y_vals, alpha = as.factor(alpha))
  plot_data <- rbind(plot_data, df)
}

plot_if <- ggplot(plot_data, aes(x = x, y = y, color = alpha)) +
  geom_line(size = 1.2) +
  theme_minimal()+
    scale_colour_viridis_d(option = "plasma",
                           name = expression(~~~~alpha))+
  labs(
    x = TeX("$q$"),
    y = TeX("$\\textit{IF}_{\\alpha,\\mu}(q,\\mu_0,\\sigma_0,\\xi_0)$")
  ) +
  scale_x_continuous(
    breaks = list_quantile,
    labels = p
  ) +
  theme(legend.position = "none",
        legend.title=element_text(size=26,face="bold"),
        legend.text =element_text(size=24))+
  theme(
    axis.title = element_text(size = 24),   
    axis.text = element_text(size = 18)
  )
plot_if
```

```{r}
ggsave(plot = plot_if,
       filename = paste0("figures/IF_mu_negative_xi.png"),
       width = 1760 / 180,
       height = 1060 / 180,
       dpi = 180,
       units = "in")

```

```{r}

# plot influence function sclae parameter (Figure 1, middle)

p <- c(0.001,0.1,0.5,0.7,0.9)
list_quantile <- mev::qgev(p, loc = mu, scale = sigma, shape = xi, lower.tail = TRUE, log.p = FALSE)


x_sig <- seq(-3, mu - sigma/xi - 0.4, by = 0.1)  

plot_data <- data.frame()

for (alpha in alphas) {
  y_vals <- IF_sigma(x_sig, mu, sigma, xi, alpha)
  df <- data.frame(x = x_sig, y = y_vals, alpha = as.factor(alpha))
  plot_data <- rbind(plot_data, df)
}

plot_if <- ggplot(plot_data, aes(x = x, y = y, color = alpha)) +
  geom_line(size = 1.2) +
  theme_minimal() +
    scale_colour_viridis_d(option = "plasma",
                           name = expression(~~~~alpha))+
  labs(
    x = TeX("$q$"),
    y = TeX("$\\textit{IF}_{\\alpha,\\sigma}(q,\\mu_0,\\sigma_0,\\xi_0)$")
  )+
  scale_x_continuous(
    breaks = list_quantile,
    labels = p
  ) +
  theme(legend.position ="none",
        legend.title=element_text(size=26,face="bold"),
        legend.text =element_text(size=24))+
  theme(
    axis.title = element_text(size = 24),   
    axis.text = element_text(size = 18)
  )
plot_if
```

```{r}
ggsave(plot = plot_if,
       filename = paste0("figures/IF_sig_negative_xi.png"),
       width = 1760 / 180,
       height = 1060 / 180,
       dpi = 180,
       units = "in")

```


```{r}

# plot influence function location parameter (Figure 1, right)

p <- c(0.001,0.1,0.5,0.8,0.9)
list_quantile <- mev::qgev(p, loc = mu, scale = sigma, shape = xi, lower.tail = TRUE, log.p = FALSE)

x_xi <- seq(-3, mu - sigma/xi - 0.5, by = 0.1)  

plot_data <- data.frame()

for (alpha in alphas) {
  y_vals <- IF_xi(x_xi, mu, sigma, xi, alpha)
  df <- data.frame(x = x_xi, y = y_vals, alpha = as.factor(alpha))
  plot_data <- rbind(plot_data, df)
}

plot_if <- ggplot(plot_data, aes(x = x, y = y, color = alpha)) +
  geom_line(size = 1.2) +
  theme_minimal() +
  scale_colour_viridis_d(option = "plasma",
                         name = expression(~~~~alpha)) +
  labs(
    x = TeX("$q$"),
    y = TeX("$\\textit{IF}_{\\alpha,\\xi}(q,\\mu_0,\\sigma_0,\\xi_0)$")
  ) +
  scale_x_continuous(
    breaks = list_quantile,
    labels = p
  ) +
  theme(
    legend.position = c(.6, .4),
    legend.title = element_text(size = 26, face = "bold"),
    legend.text = element_text(size = 24),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 18)
  )
plot_if
```

```{r}
ggsave(plot = plot_if,
       filename = paste0("figures/IF_xi_negative_xi.png"),
       width = 1760 / 180,
       height = 1060 / 180,
       dpi = 180,
       units = "in")

```


