
### Load data on bill length for the adult population in both habitats
Jays <- read.csv("Island_Scrub_Jay_Phenotypes.csv", header=TRUE)

## Match bill length phenotypes to individuals in known breeding pairs in Pine
Pairs <- read.csv("Mated_Pairs.csv", header=TRUE)
Pine <- Jays[Jays$Habitat=="Pine",]
Pairs$malenares <- with(Pairs,Pine$Ave_Nares_mm[match(Pairs$M_ID, Pine$ID)])
Pairs$femalenares <- with(Pairs,Pine$Ave_Nares_mm[match(Pairs$F_ID, Pine$ID)])
Pine_Pairs<-na.omit(Pairs)
mPine<- subset(Pine, Sex=="M")
fPine <- subset(Pine, Sex=="F")
# Remove paired pine individuals from pine singles dataframe
ID_f <- intersect(fPine$ID, Pine_Pairs$F_ID)
ID_m <- intersect(mPine$ID, Pine_Pairs$M_ID)
fPine <- fPine[!fPine$ID %in% ID_f, ]
mPine <- mPine[!mPine$ID %in% ID_m, ]

## Match bill length phenotypes to individuals in known breeding pairs in Oak
Pairs <- read.csv("Mated_Pairs.csv", header=TRUE)
Oak <- Jays[Jays$Habitat=="Oak",]
Pairs$malenares <- with(Pairs,Oak$Ave_Nares_mm[match(Pairs$M_ID, Oak$ID)])
Pairs$femalenares <- with(Pairs,Oak$Ave_Nares_mm[match(Pairs$F_ID, Oak$ID)])
Oak_Pairs<-na.omit(Pairs)
mOak <- subset(Oak, Sex=="M")
fOak <- subset(Oak, Sex=="F")
# Remove paired oak individuals from oak singles dataframe
ID_f <- intersect(fOak$ID, Oak_Pairs$F_ID)
ID_m <- intersect(mOak$ID, Oak_Pairs$M_ID)
fOak <- fOak[!fOak$ID %in% ID_f, ]
mOak <- mOak[!mOak$ID %in% ID_m, ]



##### Prepare the list of phenotypes by subpopulation 
xs = fPine$Ave_Nares_mm
ys = mPine$Ave_Nares_mm
xp = Pine_Pairs$femalenares
yp = Pine_Pairs$malenares
# 
xs = fOak$Ave_Nares_mm
ys = mOak$Ave_Nares_mm
xp = Oak_Pairs$femalenares
yp = Oak_Pairs$malenares


Jaydata <- vector(mode = "list", length = 2)
Jaydata[[1]] <- list(xs = fPine$Ave_Nares_mm,
                     ys = mPine$Ave_Nares_mm,
                     xp = Pine_Pairs$femalenares,
                     yp = Pine_Pairs$malenares)

Jaydata[[2]] <- list(xs = fOak$Ave_Nares_mm,
                     ys = mOak$Ave_Nares_mm,
                     xp = Oak_Pairs$femalenares,
                     yp = Oak_Pairs$malenares)
# Function to calculate the proportionality constant from Clancey et al (2022)
mgf <- function(t, lambda) {
  exp(lambda * t / (1 - 2*t))/sqrt(1 - 2*t)
}
# Function to calculate the likelihood function
nloglik <- function(theta, dlst) {
  K <- length(dlst)
  alpha <- theta[1]
  delta <- theta[2]
  mx <- theta[3:(K + 2)]
  my <- theta[(K + 3):(2*K + 2)]
  sx <- theta[(2*K + 3):(3*K + 2)]
  sy <- theta[(3*K + 3):(4*K + 2)]
  
  loglik <- 0
  for (k in 1:K) {
    xp <- dlst[[k]]$xp
    yp <- dlst[[k]]$yp
    xs <- dlst[[k]]$xs
    ys <- dlst[[k]]$ys
    mxy <- my[k] - mx[k] - delta
    sxy <- sx[k]^2 + sy[k]^2
    loglik <- loglik + sum(dnorm(xp, mx[k], sx[k], log = TRUE)) + sum(dnorm(yp, my[k], sy[k], log = TRUE)) -
      sum(alpha * (yp - xp - delta)^2) - length(xp) * log(mgf(- alpha * sxy, mxy^2 / sxy))
    if (!is.null(xs)) {
      loglik <- loglik + sum(dnorm(xs, mx[k], sx[k], log = TRUE))
    }
    if (!is.null(ys)) {
      loglik <- loglik + sum(dnorm(ys, my[k], sy[k], log = TRUE))
    }
  }
  return(-loglik)
}
# Function to calculate the maximize likelihood function
alphamle <- function(dlst, start = NA) {
  
  K <- length(dlst)
  
  if (is.na(start)) {
    start <- rep(NA, 2 + 4*K)
    start[1] <- 0.0001
    start[2] <- 0.0001
    for (k in 1:K) {
      start[k + 2] <- mean(c(dlst[[k]]$xp, dlst[[k]]$xs))
      start[K + k + 2] <- mean(c(dlst[[k]]$yp, dlst[[k]]$ys))
      start[2*K + k + 2] <- sd(c(dlst[[k]]$xp, dlst[[k]]$xs))
      start[3*K + k + 2] <- sd(c(dlst[[k]]$yp, dlst[[k]]$ys))
    }
  }
  m <- optim(start, nloglik, method = "L-BFGS-B", lower = c(0, -Inf, rep(-Inf, K*2), rep(1e-7, K*2)), hessian = TRUE, dlst = dlst)
  out <- cbind(m$par, sqrt(diag(solve(m$hessian))))
  rownames(out) <- c("alpha", "delta", paste("mx", 1:K, sep = ""), paste("my", 1:K, sep = ""), paste("sx", 1:K, sep = ""), paste("sy", 1:K, sep = ""))
  colnames(out) <- c("estimate","se")
  return(out)
}

# Read the results output
mle <- alphamle(Jaydata)

