### Simulate reproductive success data for analysis using conventional and Bayesian aster models
require(ggplot2) 
require(dplyr)
source("R/rpoisson0.R")
source("R/invlogit.R")

### Set parameters ----
### Number of nests to simulate
  n_plots <- 12
  n_nests <- 40

### Parasitism parameters
  alpha.p <- 0.1      # Mean parastism rate
  sd.p <- 0            # Plot-level variation in parasitism rate
  beta1 <- -1.2          # Slope parameter controlling effect of clutch completion on parasitism rate
  beta2 <- -3       # Slope parameter controlling effect of distance to edge on parasitism rate
  beta3 <- 3           # Slope parameter controlling effect of cowbird density on parasitism rate

### BHCO hatching success parameters
  alpha.ch <- 0.8      # Mean BHCO hatching success probability
  sd.ch <- 0           # Plot-level variation in BHCO hatching success
  beta8 <- 0           # Slope parameter controlling effect of nest height on BHCO hatching success
  beta9 <- 0           # Slope parameter controlling effect of vegetation structure on BHCO hatching success

### WOTH hatching success parameters
  alpha.wh <- 0.85     # Mean WOTH hatching success probability
  sd.wh <- 0           # Plot-level variation in WOTH hatching success
  beta4 <- -4          # Slope parameter controlling effect of BHCO hatching success on WOTH hatching success
  beta10 <- 0           # Slope parameter controlling effect of nest height on WOTH hatching success
  beta11 <- 0          # Slope parameter controlling effect of vegetation structure on WOTH hatching success

### BHCO offspring parameters
  mu.cf <- 0.8         # Mean number of BHCO offspring
  sd.cf <- 0           # Plot-level variation in number of BHCO offspring
  beta5 <- -2.5          # Slope parameter controlling effect of WOTH hatching success on number of BHCO offspring

### WOTH offspring parameters
  mu.wf <- 1.5           # Mean number of WOTH offspring
  sd.wf <- 0           # Plot-level variation in number of WOTH offspring
  beta6 <- -0.3          # Slope parameter controlling effect of clutch completion date on number of WOTH offspring
  beta7 <- -0.4          # Slope parameter controlling effect of number of BHCO offspring on number of WOTH offspring

  
### Simulation function ----
nest_success <- function(nPlots = n_plots, nNests = n_nests, 
                           para.mean = log(alpha.p/(1 - alpha.p)), para.sd = sd.p,
                           slope.c.dens = beta3, slope.dist = beta2, slope.para.cc = beta1, 
                           bhco.hatch.mean = log(alpha.ch/(1 - alpha.ch)), bhco.hatch.sd = sd.ch, 
                           bhco.ht.slope = beta8, bhco.veg.slope = beta9, 
                           woth.hatch.mean = log(alpha.wh/(1 - alpha.wh)), woth.hatch.sd = sd.wh, 
                           bhco.hatch.slope = beta4, woth.ht.slope = beta10, woth.veg.slope = beta11, 
                           bhco.fledge.mean = log(mu.cf), bhco.fledge.sd = sd.cf, 
                           woth.hatch.slope = beta5, 
                           woth.fledge.mean = log(mu.wf), woth.fledge.sd = sd.wf, 
                           slope.cc.woth = beta6, bhco.fledge.slope = beta7) 
  {

### Input variables ----
    N <- nPlots * nNests # Number of nests
    plot <- gl(n = nPlots, k = nNests)

### Randomly generate nest and plot covariates ----
    clutch.comp <- scale(runif(n = N, min = 0, max = 30))[,1]         # Clutch completion dates (0 - 30 days)
    ht <- scale(runif(n = N, min = 2, max = 6))[,1]                 # Nest height (2 - 6 meters)
    veg <- scale(runif(n = N, min = 0, max = 1))[,1]                  # Vegetation index (0 - 1)
    dist <- scale(runif(n = N, min = 0, max = 100))[,1]               # Distance to patch edge (0 - 100m)
    c.dens <- scale(rep(runif(n = nPlots, min = 0, max = 20), each = nNests))[,1] # Plot-level cowbird density (0 - 20 Cowbirds/unit)

### Generate parasitism data ----
    bhco.intercept.effects <- rnorm(n = nPlots, mean = para.mean, sd = para.sd)
    parasitism.effects <- c(bhco.intercept.effects, slope.c.dens, slope.dist, slope.para.cc)
    
    Xmat.para <- model.matrix(~ plot + c.dens + dist + clutch.comp - 1)
    lin.pred.para <- Xmat.para[,] %*% parasitism.effects
    exp.p.para <- invlogit(lin.pred.para)
    
    C.para <- rbinom(n = N, size = 1, prob = exp.p.para)
  
    df <- data.frame(woth.eggs = woth.eggs, bhco.eggs = bhco.eggs, C.para = C.para, P.para = exp.p.para, dist = dist, c.dens = c.dens, plot = plot, clutch.comp = clutch.comp)

### BHCO hatching success ----
  bhco.intercept.effects <- rnorm(n = nPlots, mean = bhco.hatch.mean, sd = bhco.hatch.sd)
  bhco.hatch.effects <- c(bhco.intercept.effects, bhco.ht.slope, bhco.veg.slope)
  
  Xmat.bhco.hatch <- model.matrix(~ plot + ht + veg - 1)
  
  lin.pred.bhco.hatch <- Xmat.bhco.hatch[,] %*% bhco.hatch.effects
  exp.p.bhco.hatch <- invlogit(lin.pred.bhco.hatch) * C.para
  
  C.bhco.hatch <- rbinom(n = N, size = 1, prob = exp.p.bhco.hatch) 
  
  df2 <- df %>% mutate(C.bhco.hatch = C.bhco.hatch, P.bhco.hatch = exp.p.bhco.hatch[,1], ht = ht, veg = veg)
  
### WOTH hatching success ----
  woth.intercept.effects <- rnorm(n = nPlots, mean = woth.hatch.mean, sd = woth.hatch.sd)
  
  woth.hatch.effects <- c(woth.intercept.effects, bhco.hatch.slope, woth.ht.slope, woth.veg.slope)
  
  Xmat.woth.hatch <- model.matrix(~ plot + C.bhco.hatch + ht + veg - 1)
  
  lin.pred.woth.hatch <- Xmat.woth.hatch[,] %*% woth.hatch.effects
  exp.p.woth.hatch <- invlogit(lin.pred.woth.hatch)
  
  C.woth.hatch <- rbinom(n = N, size = 1, prob = exp.p.woth.hatch)

  df3 <- df2 %>% mutate(C.woth.hatch = C.woth.hatch, P.woth.hatch = exp.p.woth.hatch[,1])
  
### BHCO fledging success ----
  bhco.intercept.effects <- rnorm(n = nPlots, mean = bhco.fledge.mean, sd = bhco.fledge.sd)
  bhco.fledge.effects <- c(bhco.intercept.effects, woth.hatch.slope)
  
  Xmat.bhco.fledge <- model.matrix(~ plot + C.woth.hatch - 1)
  
  lin.pred.bhco.fledge <- Xmat.bhco.fledge[,] %*% bhco.fledge.effects
  exp.n.bhco.fledge <- exp(lin.pred.bhco.fledge) 
  
  # Simulate zero-truncated Poisson distribution for expected number of BHCO fledged
  C.bhco.fledge <- rpoisson0(n = N, lambda = exp.n.bhco.fledge) * C.bhco.hatch

  df4 <- df3 %>% mutate(C.bhco.fledge = C.bhco.fledge, P.bhco.fledge = exp.n.bhco.fledge[,1])
  
### WOTH fledging success ----
  woth.intercept.effects <- rnorm(n = nPlots, mean = woth.fledge.mean, sd = woth.fledge.sd)
  woth.fledge.effects <- c(woth.intercept.effects, bhco.fledge.slope, slope.cc.woth)
  
  Xmat.woth.fledge <- model.matrix(~ plot + C.bhco.fledge + clutch.comp - 1)
  
  lin.pred.woth.fledge <- Xmat.woth.fledge[,] %*% woth.fledge.effects
  exp.n.woth.fledge <- exp(lin.pred.woth.fledge) 
  
  # Simulate zero-truncated Poisson distribution for expected number of BHCO fledged
  C.woth.fledge <- rpoisson0(n = N, lambda = exp.n.woth.fledge) * C.woth.hatch

  
  df5 <- df4 %>% mutate(C.woth.fledge = C.woth.fledge, P.woth.fledge = exp.n.woth.fledge[,1])

  return(df5)
}

### Simulate and visualize data ----
  data_sim <- nest_success()

  data_sim %>%
    #filter(C.para == 1) %>%
    ggplot(., aes(x = C.bhco.hatch, y = C.woth.fledge)) + geom_point(alpha = 0.5) +
      #stat_smooth(method = "glm", method.args = list(family = "binomial"))
      stat_smooth(method = "glm", method.args = list(family = "poisson")) +
    theme_bw()

