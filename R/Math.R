# This code checks math done with 4D array and matrix

# Faster computation for interaction surface
# Main effect + interaction effect 
# = exposure * effect + exposure * surface * exposure

# ------------------------------------------------
# n = 1, two exposures, lag = 12
# ------------------------------------------------
e1 <- matrix(rpois(12, 1), nrow = 1, ncol = 12)
e2 <- matrix(rpois(12, 1), nrow = 1, ncol = 12)

eff1 <- matrix(rnorm(12), nrow = 12, ncol = 1)
eff2 <- matrix(rnorm(12), nrow = 12, ncol = 1)

surf <- matrix(rnorm(12^2), nrow = 12, ncol = 12)

# Raw
mtx1 <- e1 %*% eff1 + e2 %*% surf %*% t(e1)

# matrix x matrix
mtx2 <- e1 %*% eff1 + sum(outer(c(e2), c(e1)) * surf)


# check equal
all(mtx1 == mtx2)




# ------------------------------------------------
# n = 7, three exposures, lag = 12
# ------------------------------------------------
e1 <- matrix(rpois(7*12, 1), nrow = 7, ncol = 12)
e2 <- matrix(rpois(7*12, 1), nrow = 7, ncol = 12)
e3 <- matrix(rpois(7*12, 1), nrow = 7, ncol = 12)

eff1 <- matrix(rnorm(12), nrow = 12, ncol = 1)
eff2 <- matrix(rnorm(12), nrow = 12, ncol = 1)
eff3 <- matrix(rnorm(12), nrow = 12, ncol = 1)

surf1 <- matrix(rnorm(12^2), nrow = 12, ncol = 12)
surf2 <- matrix(rnorm(12^2), nrow = 12, ncol = 12)
surf3 <- matrix(rnorm(12^2), nrow = 12, ncol = 12)

# Raw
mtx1 <- e1 %*% eff1 + e2 %*% eff2 + e3 %*% eff3 +
          + diag(e2 %*% surf1 %*% t(e1)) +
          + diag(e3 %*% surf2 %*% t(e2)) + 
          + diag(e3 %*% surf3 %*% t(e1))


# for loop one-by-one brute force
mtx2 <- matrix(0, nrow = 7, ncol = 1)
for(i in 1:7){
  mtx2[i] <- 
    e1[i,] %*% eff1 + t(e2[i,]) %*% surf1 %*% e1[i,] + 
    e2[i,] %*% eff2 + t(e3[i,]) %*% surf2 %*% e2[i,] + 
    e3[i,] %*% eff3 + t(e3[i,]) %*% surf3 %*% e1[i,]
}

# check equal
all(round(mtx1, 3) == round(mtx2, 3))
