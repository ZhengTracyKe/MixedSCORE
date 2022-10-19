Estimation <- function(R, vertices, K, eig.values, eig.vectors){
    n = nrow(R)
    memberships = matrix(0, nrow = n, ncol = K)
    
    for(i in 1:n){
        out = limSolve::lsei(A = t(vertices), 
                             B = R[i,], 
                             E = rep(1,K), F = 1,
                             G = diag(K), H = rep(0,K),
                             type = 2)
        memberships[i, ] = out$X
    }
    
    # truncate and normalize
    memberships[memberships > 1] = 1
    memberships[memberships < 0] = 0
    for(i in 1:n){
        memberships[i,] = memberships[i,]/sum(memberships[i,])
    }
    
    # recover the original memberships
    tildeV = cbind(rep(1,K), vertices)
    b1.inv = sqrt(diag(tildeV%*%diag(eig.values[1:K])%*%t(tildeV))) 
    
    
    # get degrees
    degrees = abs(eig.vectors[,1] * rowSums(memberships %*% diag(b1.inv)) )
    
    # estimate P
    P = diag(1/b1.inv) %*% tildeV %*% diag(eig.values[1:K]) %*% t(tildeV) %*% diag(1/b1.inv)
    
    
    # get tilted memberships
    memberships = memberships%*%diag( b1.inv) 
    

    # normalize again
    for(i in 1:n){
        memberships[i,] = memberships[i,]/sum(memberships[i,])
    }

    return(list(memberships = memberships, degrees = degrees, P=P))
}

