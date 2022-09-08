#based on Sherman-Morrison formula

outer_prod <- function(vec1, vec2) {
  outer(vec1, vec2, function(x,y) x*y)
}

inverse_perturbation <- function(M, i, j, lambda) {
  m_ij <- M[i,j]
  m_ii <- M[i,i]
  m_jj <- M[j,j]
  M_i <- M[i,]
  M_j <- M[j,]
  
  M_ij <- outer_prod(M_i, M_j)
  M_ji <- outer_prod(M_j, M_i)
  M_ii <- outer_prod(M_i, M_i)
  M_jj <- outer_prod(M_j, M_j)
  
  
  inv_new <- M - ((lambda * (1 + lambda*m_ij)) / ((1 + lambda * m_ij)^2 - lambda^2* m_ii * m_jj)) * (M_ij + M_ji -
                                                              ((lambda * m_jj) / (1+lambda*m_ij)) * M_ii -
                                                              ((lambda * m_ii) / (1+lambda*m_ij)) * M_jj)
  inv_new
}

#test for correctness
for (i in (1:10000)) {
  m <- matrix(runif(20*20, min = -1, max = 1), nrow = 20, ncol = 20)
  A <- m + t(m)
  A.inv <- solve(A)
  
  delta <- matrix(0, nrow = 20, ncol = 20)
  delta[1,4] = 2
  delta[4,1] = 2
  B <- A + delta
  
  B.inv <- solve(B)
  B.inv2 <- inverse_perturbation(A.inv, 1, 4, 2)
  if (sum(abs(B.inv - B.inv2)) > 10^(-3)) {
    cat("Perturbation method gives wrong result", sum(abs(B.inv - B.inv2)))
    break
  }
}